#include "../../include/matting/CGMattingSolver.h"

CGMattingSolver::CGMattingSolver(const ImageMatrix& image, const ImageMatrix& trimap, const size_t kernel_radius)
	: m_image(image), m_trimap(trimap), m_kernel_radius(kernel_radius)
{
	if (m_image.width() != m_trimap.width() ||
	    m_image.height() != m_trimap.height())
	{
		throw std::invalid_argument("image and trimap are of different sizes");
	}

	m_image.normalizeIfNot();
	m_trimap.normalizeIfNot();

	m_constraint = 65536.0;
	m_reg_param = 0.01;

	initConstraints();
	initImageSAT();
	initImageCovariance();
	m_def_constraint = m_constraint * (m_constraint_mat * m_trimap.colormap());
}

void CGMattingSolver::initConstraints()
{
	size_t img_area = m_image.area();
	m_constraint_mat = Eigen::SparseMatrix<double>(img_area, img_area);
	for (size_t i = 0; i < img_area; ++i)
	{
		Point coord = m_image.coords2D(i);
		double val = m_trimap.getAt(coord.first, coord.second, 0);
		m_constraint_mat.coeffRef(i, i) = (val > 1 - KThreshold || val < KThreshold) ? 1 : 0;
	}
}

void CGMattingSolver::initImageSAT()
{
	size_t channels = m_image.numComponents();
	m_image_sat.clear();
	for (size_t channel = 0; channel < channels; ++channel)
	{
		auto channel_data = m_image.colormap().col(channel).data();
		MatrixXdRow image_channel = Eigen::Map<MatrixXdRow>(channel_data, m_image.height(), m_image.width());
		m_image_sat.emplace_back(image_channel);
	}
}

void CGMattingSolver::initImageCovariance()
{
	m_win_cov = std::vector<std::vector<Eigen::MatrixXd>>(m_image.height(),
														  std::vector<Eigen::MatrixXd>(m_image.width()));
	size_t rows = m_image.height(), cols = m_image.width();
	for (size_t row = 0; row < rows; ++row)
	{
		for (size_t col = 0; col < cols; ++col)
		{
			m_win_cov[row][col] = calcCov(row, col);
		}
	}
}

Eigen::MatrixXd CGMattingSolver::calcCov(const size_t row, const size_t col)
{
	size_t channels = m_image.numComponents();
	size_t rows = m_image.height(), cols = m_image.width();
	Vector mean(channels);
	ImageWindow cur_window({rows, cols}, {row, col}, m_kernel_radius);
	size_t win_area = cur_window.getArea(), index = 0;
	Eigen::MatrixXd observations(win_area, channels);

	Point left_top, right_bottom;
	left_top = cur_window.getOriginBound();
	right_bottom = cur_window.getEndBoundIncl();

	for (size_t channel = 0; channel < channels; ++channel)
	{
		mean[channel] = m_image_sat[channel].windowSum(left_top, right_bottom) / win_area;
	}
	for (auto &coord: cur_window)
	{
		size_t flat_coord = m_image.flatCoords(coord);
		observations.row(index) = m_image.colormap().row(flat_coord);
		++index;
	}
	auto deviation = observations - Vector::Ones(win_area) * mean.transpose();
	Eigen::Matrix3d covariance = deviation.transpose() * deviation / (win_area - 1);
	return covariance;
}

CGMattingSolver::Vector CGMattingSolver::alphaMatting(size_t iterations, double precision)
{
	m_cg_solution = Vector::Zero(m_image.area());
	m_residual = m_def_constraint;
	m_conjugate = m_def_constraint;
	for (size_t iter = 0; iter < iterations; ++iter)
	{
		std::cout << "Iteration: #" << iter << std::endl;
		CGIterate();
		if (m_residual.norm() < precision) break;
	}
	return m_cg_solution;
}

void CGMattingSolver::CGIterate()
{
	Vector residual_prev = m_residual;
	MatrixXdRow conjugate_mat = Eigen::Map<MatrixXdRow>(m_conjugate.data(),
												 m_image.height(),
												 m_image.width());
	SATMatrix conjugate_sat(conjugate_mat);

	std::vector<SATMatrix> Ip_cwise_sat;
	size_t rows = m_image.height(), cols = m_image.width();
	for (size_t channel = 0; channel < m_image.numComponents(); ++channel)
	{
		MatrixXdRow image_channel = Eigen::Map<MatrixXdRow>(m_image.colormap().col(channel).data(), rows, cols);
		image_channel = image_channel.cwiseProduct(conjugate_mat);
		Ip_cwise_sat.emplace_back(image_channel);
	}
	std::vector<SATMatrix> slope_sat = calcSlopeAndBiasSAT(Ip_cwise_sat, conjugate_sat);
	SATMatrix bias_sat = slope_sat.back();
	slope_sat.pop_back();
	Vector Lp = laplacianProduct(slope_sat, bias_sat);

	Vector coeff_product = Lp + m_constraint * (m_constraint_mat * m_conjugate);
	double alpha = m_residual.dot(m_residual) / handleZeroDiv(coeff_product.dot(m_conjugate));
	m_cg_solution += alpha * m_conjugate;
	m_residual -= alpha * coeff_product;
	double beta = m_residual.dot(m_residual) / handleZeroDiv(residual_prev.dot(residual_prev));
	m_conjugate = m_residual + beta * m_conjugate;
}

std::vector<SATMatrix> CGMattingSolver::calcSlopeAndBiasSAT(const std::vector<SATMatrix>& Ip_sat,
													 		const SATMatrix& conjugate_sat)
{
	size_t rows = m_image.height(), cols = m_image.width(), channels = m_image.numComponents();
	std::vector<SATMatrix> sat_slope;
	Eigen::MatrixXd bias_mat(rows, cols);
	std::vector<Eigen::MatrixXd> slope_channels(channels, Eigen::MatrixXd(rows, cols));
	Vector mean(channels), Ip(channels), cur_slope(channels);
	Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(channels, channels);
	for (size_t row = 0; row < rows; ++row)
	{
		for (size_t col = 0; col < cols; ++col)
		{
			Point shape, center;
			center = { row, col };
			shape = { rows, cols };
			ImageWindow cur_window(shape, center, m_kernel_radius);
			size_t win_area = cur_window.getArea();
			double conjugate_mean;

			Point left_top, right_bottom;
			left_top = cur_window.getOriginBound();
			right_bottom = cur_window.getEndBoundIncl();

			conjugate_mean = conjugate_sat.windowSum(left_top, right_bottom) / win_area;
			for (size_t channel = 0; channel < channels; ++channel)
			{
				Ip[channel] = Ip_sat[channel].windowSum(left_top, right_bottom);
				mean[channel] = m_image_sat[channel].windowSum(left_top, right_bottom) / win_area;
				cur_slope[channel] = (Ip[channel] / win_area) - conjugate_mean * mean[channel];
			}

			Eigen::MatrixXd delta = m_win_cov[row][col] + identity * m_reg_param / win_area;
			Eigen::MatrixXd delta_inv = delta.inverse();
			cur_slope = delta_inv * cur_slope;
			bias_mat(row, col) = conjugate_mean - cur_slope.dot(mean);
			for (size_t channel = 0; channel < channels; ++channel)
				slope_channels[channel](row, col) = cur_slope[channel];
		}
	}
	for (size_t channel = 0; channel < channels; ++channel)
		sat_slope.emplace_back(slope_channels[channel]);
	sat_slope.emplace_back(bias_mat);
	return sat_slope;
}

CGMattingSolver::Vector
CGMattingSolver::laplacianProduct(const std::vector<SATMatrix> &slope_sat, const SATMatrix &bias_sat)
{
	size_t rows = m_image.height(), cols = m_image.width();
	Vector lap_product(rows * cols);
	for (int i = 0; i < lap_product.size(); ++i)
	{
		Vector slope_sum(m_image.numComponents());
		double b_star_sum;

		Point shape, center;
		center = m_image.coords2D(i);
		shape = { rows, cols };
		ImageWindow cur_window(shape, center, m_kernel_radius);

		Point left_top, right_bottom;
		left_top = cur_window.getOriginBound();
		right_bottom = cur_window.getEndBoundIncl();
		for (size_t channel = 0; channel < m_image.numComponents(); ++channel)
		{
			slope_sum[channel] = slope_sat[channel].windowSum(left_top, right_bottom);
		}
		b_star_sum = bias_sat.windowSum(left_top, right_bottom);
		lap_product[i] = cur_window.getArea() * m_conjugate[i] - (slope_sum.dot(m_image.colormap().row(i)) + b_star_sum);
	}
	return lap_product;
}
