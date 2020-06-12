#include "../../include/cfmatting/ImageMatrix.h"

ImageMatrix::ImageMatrix(const size_t& width,
                         const size_t& height,
                         const size_t& num_components)
        : m_width(width), m_height(height), m_comps_n(num_components)
{
    size_t row_length = m_width * m_comps_n;
    m_colormap = new JSAMPROW[m_height];
    for (size_t row = 0; row < m_height; ++row)
    {
        m_colormap[row] = new JSAMPLE[row_length];
    }
}

ImageMatrix::ImageMatrix(const size_t &width,
                         const size_t &height,
                         const size_t &num_components,
                         JSAMPARRAY colormap)
        : m_width(width), m_height(height), m_comps_n(num_components),
            m_colormap(colormap) { }

ImageMatrix::~ImageMatrix()
{
    m_free_colormap();
}

ImageMatrix::ImageMatrix(const ImageMatrix &other)
{
    m_width = other.m_width;
    m_height = other.m_height;
    m_comps_n = other.m_comps_n;
    m_colormap = new JSAMPROW[m_height];
    for (size_t row = 0; row < m_height; ++row)
    {
        m_colormap[row] = new JSAMPLE[m_width];
        for (size_t col = 0; col < m_width; ++col)
        {
            m_colormap[row][col] = other.m_colormap[row][col];
        }
    }
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

JSAMPLE ImageMatrix::getAt(const size_t &row,
                           const size_t &col,
                           const size_t &component) const
{
    return m_colormap[row][col * m_comps_n + component];
}

void ImageMatrix::m_free_colormap()
{
    for (size_t row = 0; row < m_height; ++row)
    {
        delete [] m_colormap[row];
    }
    delete [] m_colormap;
}

void ImageMatrix::swap(ImageMatrix &first, ImageMatrix &second)
{
	std::swap(first.m_width, second.m_width);
	std::swap(first.m_height, second.m_height);
	std::swap(first.m_comps_n, second.m_comps_n);
	std::swap(first.m_colormap, second.m_colormap);
}

