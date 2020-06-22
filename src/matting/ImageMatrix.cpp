#include "../../include/matting/ImageMatrix.h"

ImageMatrix::ImageMatrix(const size_t width,
                         const size_t height,
                         const size_t num_components)
        : m_width(width), m_height(height), m_comps_n(num_components)
{
	if (width == 0 || height == 0)
	{
		throw std::invalid_argument("zero size");
	}
    m_colormap = Eigen::MatrixXd(m_height * m_width, m_comps_n);
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
}

ImageMatrix::~ImageMatrix() = default;

ImageMatrix::ImageMatrix(const ImageMatrix &other)
{
    m_width = other.m_width;
    m_height = other.m_height;
    m_comps_n = other.m_comps_n;
    m_colormap = other.m_colormap;
}

ImageMatrix& ImageMatrix::operator=(const ImageMatrix& other)
{
	ImageMatrix other_copy = other;
	swap(*this, other_copy);
	return *this;
}

void ImageMatrix::swap(ImageMatrix &first, ImageMatrix &second)
{
	std::swap(first.m_width, second.m_width);
	std::swap(first.m_height, second.m_height);
	std::swap(first.m_comps_n, second.m_comps_n);
	std::swap(first.m_colormap, second.m_colormap);
}

ImageMatrix::Point ImageMatrix::coords2D(const size_t flat_coords) const
{
	size_t row, col;
	row = flat_coords / m_width;
	col = flat_coords % m_width;
	return ImageMatrix::Point(row, col);
}

void ImageMatrix::setAt(const size_t row, const size_t col, const size_t component, const double val)
{
	size_t flat = flatCoords(row, col);
	m_colormap(flat, component) = val;
}

double ImageMatrix::getAt(const size_t row, const size_t col, const size_t component) const
{
	size_t flat = flatCoords(row, col);
	return m_colormap(flat, component);
}

void ImageMatrix::toImageFormat()
{
	auto min_matrix = Eigen::MatrixXd::Ones(m_width * m_height, m_comps_n) * m_colormap.minCoeff();
	m_colormap = m_colormap - min_matrix;
	m_colormap  = m_colormap  / m_colormap.maxCoeff() * KMax_brightness;
}

bool ImageMatrix::isNormalized() const
{
	return (1. >= m_colormap.array() &&
			m_colormap.array() >= 0.).all();
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

ImageMatrix::ImageMatrix(const Eigen::VectorXd &matrix, const size_t width, const size_t height)
	: ImageMatrix(width, height, 1)
{
	m_colormap = matrix;
}

void ImageMatrix::normalize()
{
	double colors_min = m_colormap.minCoeff();
	m_colormap -= colors_min * Eigen::MatrixXd::Ones(m_height, m_width);
	double colors_max = m_colormap.maxCoeff();
	m_colormap /= colors_max;
}
