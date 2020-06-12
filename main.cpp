#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <jpeglib.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCore>

#include "include/cfmatting/ImageMatrix.h"
#include "include/cfmatting/ImageWindow.h"
#include "include/cfmatting/SATMatrix.h"

// TODO introduce class with these functions as
//  methods to avoid argument copying and provide access constraints
Eigen::VectorXd CGLaplacianSolver(Eigen::MatrixXd &image,
								  const Eigen::VectorXd &trimap,
								  size_t rows, size_t cols, double eps=1.e-5,
								  double constraint_coeff=100.0);

ImageMatrix read_jpeg(FILE *jpeg_file);

Eigen::VectorXd laplacian_product(const Eigen::MatrixXd &image,
								  size_t rows, size_t cols,
		  					  	  const Eigen::VectorXd &conjugate,
								  const SATMatrix &red, const SATMatrix &green, const SATMatrix &blue,
		  					  	  const size_t &window_radius=1);

double laplacian_product_entry(const Eigen::MatrixXd &image, ImageWindow &window,
							   const Eigen::VectorXd &conjugate);

Eigen::Vector3d laplacian_slope(const Eigen::MatrixXd &image, ImageWindow &window,
					   const Eigen::VectorXd &conjugate, double eps);

double laplacian_bias(const Eigen::MatrixXd &image, ImageWindow &window,
					  const Eigen::VectorXd &conjugate,
					  const Eigen::Vector3d &slope);

Point deploy_flat_coords(size_t flat_coords, const Point &img_shape);


int main(int argc, char* argv[])
{
    std::string file_path = argv[1];
    FILE *in_file;
    if ((in_file = fopen(file_path.c_str(), "rb")) == nullptr)
    {
        std::cout << "Cannot open file " << file_path << std::endl;
        return -1;
    }
    ImageMatrix im = read_jpeg(in_file);
    fclose(in_file);

	Eigen::MatrixXd mat(im.width() * im.height(), 3), trimap(im.width() * im.height(), 1);
	for (size_t row = 0; row < im.height(); ++row)
	{
		for (size_t col = 0; col < im.width(); ++col)
		{
			mat(im.width() * row + col, 0) = im.getAt(row, col, 0);
			mat(im.width() * row + col, 1) = im.getAt(row, col, 1);
			mat(im.width() * row + col, 2) = im.getAt(row, col, 2);
		}
	}

//	Eigen::MatrixXd mat(9,3), trimap(9, 1);
//	mat << 31,48,132,
//	       17,32,112,
//	       194,59,36,
//	       0,15,95,
//	       194,59,36,
//	       255,36,0,
//	       0,19,120,
//	       0,15,95,
//	       0,19,120;
//	trimap << 0,100,100,0,255,255,0,0,0;


	file_path = argv[2];
	if ((in_file = fopen(file_path.c_str(), "rb")) == nullptr)
	{
		std::cout << "Cannot open file " << file_path << std::endl;
		return -1;
	}
	ImageMatrix jpeg_trimap = read_jpeg(in_file);
	fclose(in_file);
	for (size_t row = 0; row < im.height(); ++row)
	{
		for (size_t col = 0; col < im.width(); ++col)
		{
			trimap(im.width() * row + col, 0) = jpeg_trimap.getAt(row, col, 0);
		}
	}

	mat /= 255;
	trimap /= 255;

	Eigen::VectorXd alpha = CGLaplacianSolver(mat, trimap, im.height(), im.width());
	alpha = (alpha - Eigen::VectorXd::Ones(alpha.rows()) * alpha.minCoeff());
	alpha = alpha / alpha.maxCoeff() * 255;


	std::ofstream out_fs("output.txt", std::fstream::out);

	for (size_t row = 0; row < im.height(); ++row)
	{
		for (size_t col = 0; col < im.width(); ++col)
		{
			out_fs << alpha[im.width() * row + col] << " ";
		}
		out_fs << std::endl;
	}

	out_fs.close();


	return 0;
}

// TODO add exception handlers
ImageMatrix read_jpeg(FILE *jpeg_file)
{
    const size_t buffer_sample_n = 1;
    size_t row_sample_len;
    JSAMPARRAY dcomp_buffer = new JSAMPROW[buffer_sample_n], image;
    struct jpeg_decompress_struct dcomp{};
    struct jpeg_error_mgr jerr{};
    dcomp.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&dcomp);

    jpeg_stdio_src(&dcomp, jpeg_file);
    jpeg_read_header(&dcomp, true);
    row_sample_len = dcomp.image_width * dcomp.num_components;
    image = new JSAMPROW[dcomp.image_height];
    for (size_t jsample_index = 0;
                    jsample_index < dcomp.image_height;
                        ++jsample_index)
    {
        image[jsample_index] = new JSAMPLE[row_sample_len];
    }
    jpeg_start_decompress(&dcomp);
    while (dcomp.output_scanline < dcomp.image_height)
    {
        dcomp_buffer[0] = image[dcomp.output_scanline];
        jpeg_read_scanlines(&dcomp, dcomp_buffer, buffer_sample_n);
    }
    ImageMatrix image_mat(dcomp.image_width, dcomp.image_height, dcomp.num_components, image);
    jpeg_finish_decompress(&dcomp);
    jpeg_destroy_decompress(&dcomp);
    delete[] dcomp_buffer;
    return image_mat;
}

Eigen::VectorXd CGLaplacianSolver(Eigen::MatrixXd &image,
								  const Eigen::VectorXd &trimap,
								  const size_t rows, const size_t cols,
								  double eps,
								  double constraint_coeff)
{
	// TODO global variables
	const double zero_eps = 0.00001;
	const size_t radius = 1;
	using vec = Eigen::VectorXd;

	std::vector<SATMatrix> sat_image, sat_Ip_cwise, sat_slope;
	std::vector<std::vector<Eigen::Matrix3d>> window_covariances(rows, std::vector<Eigen::Matrix3d>(cols));
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> image_channel;

	for (size_t channel = 0; channel < 3; ++channel)
	{
		image_channel = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(image.col(channel).data(), rows, cols);
		sat_image.emplace_back(image_channel);
	}

	for (size_t row = 0; row < rows; ++row)
	{
		for (size_t col = 0; col < cols; ++col)
		{
			Eigen::Vector3d mean;
			Point shape, center;
			center << row, col;
			shape << rows, cols;
			ImageWindow cur_window(shape, center, radius);
			size_t win_area = cur_window.getArea(), index = 0;
			Eigen::MatrixXd observations(win_area, 3);

			std::pair<size_t, size_t> left_top, right_bottom;
			left_top.first = cur_window.getOriginBound()[0];
			left_top.second = cur_window.getOriginBound()[1];
			right_bottom.first = cur_window.getEndBound()[0] - 1;
			right_bottom.second = cur_window.getEndBound()[1] - 1;

			for (size_t channel = 0; channel < 3; ++channel)
			{
				mean[channel] = sat_image[channel].windowSum(left_top, right_bottom) / win_area;
			}
			for (auto &coord: cur_window)
			{
				observations.row(index) = image.row(cur_window.flatCoords(coord));
				++index;
			}
			auto x = observations - Eigen::VectorXd::Ones(win_area) * mean.transpose();
			Eigen::Matrix3d covariance = x.transpose() * x / (win_area - 1);
			window_covariances[row][col] = covariance;
		}
	}

	size_t aaa = rows * cols;
	Eigen::SparseMatrix<double> constraints_diag(aaa, aaa);
	for (size_t i = 0; i < rows * cols; ++i)
	{
		constraints_diag.coeffRef(i, i) = (trimap[i] > 0.99 || trimap[i] < 0.01) ? 1 : 0;
	}
	Eigen::VectorXd bias = constraint_coeff * (constraints_diag * trimap);







	size_t n_size = bias.size();
	double alpha, beta;
	vec solution(n_size);
	solution.fill(0.);
	vec residual = bias, conjugate = bias;


	// CG START

	for (size_t i = 0; i < n_size; ++i)
	{
		vec residual_prev = residual;

		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> conj_matrix(conjugate.data(), rows, cols);
		SATMatrix sat_conj(conj_matrix);

		// TODO vector type cwise product; less lines
		sat_Ip_cwise.clear();
		for (size_t channel = 0; channel < 3; ++channel)
		{
			image_channel = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(image.col(channel).data(), rows, cols);
			image_channel = image_channel.cwiseProduct(conj_matrix);
			sat_Ip_cwise.emplace_back(image_channel);
		}

		std::vector<Eigen::MatrixXd> slope_channels(3, Eigen::MatrixXd(rows, cols));
		Eigen::MatrixXd biases(rows, cols);
		for (size_t row = 0; row < rows; ++row)
		{
			for (size_t col = 0; col < cols; ++col)
			{
				Eigen::Vector3d mean, Ip, cur_slope;
				Point shape, center;
				center << row, col;
				shape << rows, cols;
				ImageWindow cur_window(shape, center, radius);
				size_t win_area = cur_window.getArea();
				double conjugate_mean;

				std::pair<size_t, size_t> left_top, right_bottom;
				left_top.first = cur_window.getOriginBound()[0];
				left_top.second = cur_window.getOriginBound()[1];
				right_bottom.first = cur_window.getEndBound()[0] - 1;
				right_bottom.second = cur_window.getEndBound()[1] - 1;

				conjugate_mean = sat_conj.windowSum(left_top, right_bottom) / win_area;
				for (size_t channel = 0; channel < 3; ++channel)
				{
					Ip[channel] = sat_Ip_cwise[channel].windowSum(left_top, right_bottom);
					mean[channel] = sat_image[channel].windowSum(left_top, right_bottom) / win_area;
					cur_slope[channel] = Ip[channel] / win_area - conjugate_mean * mean[channel];
				}


				Eigen::Matrix3d delta = window_covariances[row][col] + Eigen::Matrix3d::Identity() * eps / win_area;
				cur_slope = delta.inverse() * cur_slope;
				for (size_t channel = 0; channel < 3; ++channel)
					slope_channels[channel](row, col) = cur_slope[channel];

				biases(row, col) = conjugate_mean - cur_slope.dot(mean);
			}
		}
		sat_slope.clear();
		for (size_t channel = 0; channel < 3; ++channel)
		{
			sat_slope.emplace_back(slope_channels[channel]);
		}
		SATMatrix sat_bias(biases);

		Eigen::VectorXd lap_product(rows * cols);
		Point shape;
		shape << rows, cols;

		for (size_t j = 0; j < rows * cols; ++j)
		{
			Point center;
			center = deploy_flat_coords(j, shape);
			Eigen::Vector3d slope_sum;
			double b_star_sum;
			ImageWindow cur_window(shape, center, radius);

			std::pair<size_t, size_t> left_top, right_bottom;
			left_top.first = cur_window.getOriginBound()[0];
			left_top.second = cur_window.getOriginBound()[1];
			right_bottom.first = cur_window.getEndBound()[0] - 1;
			right_bottom.second = cur_window.getEndBound()[1] - 1;

			for (size_t channel = 0; channel < 3; ++channel)
			{
				slope_sum[channel] = sat_slope[channel].windowSum(left_top, right_bottom);
			}
			b_star_sum = sat_bias.windowSum(left_top, right_bottom);
			lap_product[j] = cur_window.getArea() * conjugate[j] - (slope_sum.dot(image.row(j)) + b_star_sum);
		}

		// TODO constraint_coeff * constraints_diag as variable
//		vec conj_solution = laplacian_product(image, rows, cols, conjugate, sat_image[0], sat_image[1], sat_image[2], radius) +
//								(constraint_coeff * (constraints_diag * conjugate));

		vec conj_solution = lap_product +
							(constraint_coeff * (constraints_diag * conjugate));

		alpha = residual.dot(residual) / conj_solution.dot(conjugate);
		solution += alpha * conjugate;
		residual -= alpha * conj_solution;

		std::cout << i << std::endl << std::endl;
//		if (i == 20)
//			break;

		if (residual.norm() < eps)
			break;
		// TODO zero division
		beta = residual.dot(residual) / residual_prev.dot(residual_prev);
		conjugate = residual + beta * conjugate;
	}
	return solution;
}

Eigen::VectorXd laplacian_product(const Eigen::MatrixXd &image,
								  const size_t rows, const size_t cols,
								  const Eigen::VectorXd &conjugate,
								  const SATMatrix &red, const SATMatrix &green, const SATMatrix &blue,
								  const size_t &window_radius)
{
	size_t flat_img_size = rows * cols;
	Eigen::VectorXd product(flat_img_size);
	Point shape, center;
	shape << rows, cols;
	ImageWindow window;
	for (size_t i = 0; i < flat_img_size; ++i)
	{
		center = deploy_flat_coords(i, shape);
		window = ImageWindow(shape, center, window_radius);
		product(i) = laplacian_product_entry(image, window, conjugate);
	}
	return product;
}

// TODO const iterator for window and const argument passing
double laplacian_product_entry(const Eigen::MatrixXd &image, ImageWindow &window,
							   const Eigen::VectorXd &conjugate)
{
	size_t flat_coord;
	Eigen::Vector3d slope;
	ImageWindow window_kernel;
	double bias, entry = 0;
	for (const auto &coord: window)
	{
		Eigen::Vector3d img_i;
		if (!window.inBounds(coord))
		{
			continue;
		}
		flat_coord = window.flatCoords(coord);
		img_i = image.row(flat_coord);
		window_kernel = ImageWindow(window.getImageShape(), coord, window.getRadius());
		// TODO adaptive epsilon substitution
		slope = laplacian_slope(image, window_kernel, conjugate, 0.000001);
		bias = laplacian_bias(image, window_kernel, conjugate, slope);
		entry += conjugate(window.flatCoords(window.getCenter())) - (slope.dot(img_i) + bias);
	}
	return entry;
}

// TODO refactor code
Eigen::Vector3d laplacian_slope(const Eigen::MatrixXd &image, ImageWindow &window,
							   const Eigen::VectorXd &conjugate, const double eps)
{
//	size_t window_area = (2 * window.getRadius() + 1) * (2 * window.getRadius() + 1);
	size_t window_area = window.getArea();
	Eigen::Vector3d comulative_image_conj = Eigen::Vector3d::Zero(3);
	Eigen::Vector3d win_mean = Eigen::Vector3d::Zero(3);
	Eigen::MatrixXd observations = Eigen::MatrixXd::Zero(0, 3);
	double conjugate_mean = 0.0;
	size_t flat_coord, index = 0;
	for (const auto &coord: window)
	{
		double conjugate_i;
		Eigen::Vector3d img_i;
//		if (!window.inBounds(coord))
//		{
//			continue;
//		}
		flat_coord = window.flatCoords(coord);
		conjugate_i = conjugate[flat_coord];
		img_i = image.row(flat_coord);

		comulative_image_conj += img_i * conjugate_i;
		win_mean += img_i;
		conjugate_mean += conjugate_i;

//		observations.conservativeResize(window_area, Eigen::NoChange);
		observations.row(index) = img_i;
		++index;
	}
	auto x = observations - Eigen::VectorXd::Ones(window_area) * win_mean.transpose() / static_cast<double>(window_area);
	Eigen::Matrix3d covariance = x.transpose() * x / static_cast<double>(window_area - 1);
	Eigen::Matrix3d delta = covariance + Eigen::Matrix3d::Identity() * eps / static_cast<double>(window_area);
	return  delta.inverse() * (comulative_image_conj - win_mean * conjugate_mean) / static_cast<double>(window_area);
}

double laplacian_bias(const Eigen::MatrixXd &image, ImageWindow &window,
					  const Eigen::VectorXd &conjugate,
					  const Eigen::Vector3d &slope)
{
	size_t window_area = 0;
	Eigen::Vector3d win_mean = Eigen::Vector3d::Zero(3);
	double conjugate_mean = 0.0;
	size_t flat_coord;
	for (const auto &coord: window)
	{
		double conjugate_i;
		Eigen::Vector3d img_i;
		if (!window.inBounds(coord))
		{
			continue;
		}
		++window_area;
		flat_coord = window.flatCoords(coord);
		conjugate_i = conjugate[flat_coord];

		img_i = image.row(flat_coord);
		conjugate_mean += conjugate_i;
		win_mean += img_i;
	}
	conjugate_mean /= window_area;
	win_mean /= window_area;
	return conjugate_mean - slope.dot(win_mean);
}

Point deploy_flat_coords(const size_t flat_coords, const Point &img_shape)
{
	Point coords;
	coords << flat_coords / img_shape[1], flat_coords % img_shape[1];
	return coords;
}