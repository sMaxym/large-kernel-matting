#include "../../include/cfmatting/ImageMatrix.h"

ImageMatrix::ImageMatrix(const size_t& width,
                         const size_t& height,
                         const size_t& num_components)
        : m_width(width), m_height(height), m_comps_n(num_components)
{
    m_colormap = Eigen::MatrixXd(m_height * m_width, m_comps_n);
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

ImageMatrix::ImageMatrix(ImageMatrix &&other) noexcept
	: ImageMatrix()
{
	swap(*this, other);
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
	return ImageMatrix::Point(row, col);
}

void ImageMatrix::setAt(const size_t &row, const size_t &col, const size_t &component, const double val)
{
	size_t flat = flatCoords(row, col);
	m_colormap(flat, component) = val;
}

double ImageMatrix::getAt(const size_t &row, const size_t &col, const size_t &component) const
{
	size_t flat = flatCoords(row, col);
	return m_colormap(flat, component);
}

void ImageMatrix::toImageFormat()
{
	auto min_matrix = Eigen::MatrixXd::Ones(m_width * m_height, m_comps_n) * m_colormap.minCoeff();
	m_colormap = m_colormap - min_matrix;
	m_colormap  = m_colormap  / m_colormap.maxCoeff() * KMax_brightness;
	// TODO round
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
