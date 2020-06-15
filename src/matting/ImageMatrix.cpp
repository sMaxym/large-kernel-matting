#include "../../include/matting/ImageMatrix.h"

ImageMatrix::ImageMatrix(const size_t& width,
                         const size_t& height,
                         const size_t& num_components)
        : m_width(width), m_height(height), m_comps_n(num_components)
{
	if (width == 0 || height == 0)
	{
		throw std::invalid_argument("zero size");
	}
    m_colormap = Eigen::MatrixXd(m_height * m_width, m_comps_n);
    updateNormalized();
}

ImageMatrix::ImageMatrix(const Eigen::MatrixXd &matrix)
	: ImageMatrix(matrix.cols(), matrix.rows(), 1)
{
	for (size_t row = 0; row < m_height; ++row)
	{
		for (size_t col = 0; col < m_width; ++col)
		{
			setAt(row, col, 0, matrix.coeffRef(row, col));
		}
	}
	updateNormalized();
}

ImageMatrix::~ImageMatrix() = default;

ImageMatrix::ImageMatrix(const ImageMatrix &other)
{
    m_width = other.m_width;
    m_height = other.m_height;
    m_comps_n = other.m_comps_n;
    m_colormap = other.m_colormap;
    m_is_normalized = other.m_is_normalized;
}

ImageMatrix &ImageMatrix::operator= (ImageMatrix other)
{
	swap(*this, other);
	return *this;
}

void ImageMatrix::swap(ImageMatrix &first, ImageMatrix &second)
{
	std::swap(first.m_width, second.m_width);
	std::swap(first.m_height, second.m_height);
	std::swap(first.m_comps_n, second.m_comps_n);
	std::swap(first.m_colormap, second.m_colormap);
	std::swap(first.m_is_normalized, second.m_is_normalized);
}

ImageMatrix::Point ImageMatrix::coords2D(const size_t &flat_coords) const
{
	size_t row, col;
	row = flat_coords / m_width;
	col = flat_coords % m_width;
	assertCoords(row, col, 0);
	return ImageMatrix::Point(row, col);
}

void ImageMatrix::setAt(const size_t &row, const size_t &col, const size_t &component, const double val)
{
	assertCoords(row, col, component);
	size_t flat = flatCoords(row, col);
	m_colormap(flat, component) = val;
}

double ImageMatrix::getAt(const size_t &row, const size_t &col, const size_t &component) const
{
	assertCoords(row, col, component);
	size_t flat = flatCoords(row, col);
	return m_colormap(flat, component);
}

void ImageMatrix::toImageFormat()
{
	auto min_matrix = Eigen::MatrixXd::Ones(m_width * m_height, m_comps_n) * m_colormap.minCoeff();
	m_colormap = m_colormap - min_matrix;
	m_colormap  = m_colormap  / m_colormap.maxCoeff() * KMax_brightness;
}

void ImageMatrix::updateNormalized()
{
	m_is_normalized = (1. >= m_colormap.array() &&
					   m_colormap.array() >= 0.).all();
}

bool ImageMatrix::isNormalized()
{
	updateNormalized();
	return m_is_normalized;
}

void ImageMatrix::assertCoords(size_t row, size_t col, size_t comp) const
{
	if (row >= m_height || col >= m_width)
	{
		throw std::invalid_argument("coords not in range of image");
	}
	if (comp >= m_comps_n)
	{
		throw std::invalid_argument("component does not exist");
	}
}

ImageMatrix ImageMatrix::grayscale() const
{
	ImageMatrix gray(m_width, m_height, 1);
	double red, green, blue, brightness;
	for (size_t row = 0; row < m_height; ++row)
	{
		for (size_t col = 0; col < m_width; ++col)
		{
			red = getAt(row, col, 0);
			green = getAt(row, col, 1);
			blue = getAt(row, col, 2);
			brightness = 0.3 * red + 0.59 * green + 0.11 * blue;
			gray.setAt(row, col, 0, brightness);
		}
	}
	return gray;
}

void ImageMatrix::expandColorspace()
{
	m_colormap.conservativeResize(m_height * m_width, ++m_comps_n);
	m_colormap.col(m_comps_n - 1) = m_colormap.col(0);
}

ImageMatrix::ImageMatrix(const Eigen::VectorXd &matrix, const size_t &width, const size_t &height)
	: ImageMatrix(width, height, 1)
{
	m_colormap = matrix;
}
