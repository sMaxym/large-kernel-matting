#ifndef CFIMAGEMATTING_IMAGEMATRIX_H
#define CFIMAGEMATTING_IMAGEMATRIX_H

#include <stdio.h>
#include <iostream>
#include <jpeglib.h>

class ImageMatrix
{
public:
    ImageMatrix() = default;
    ImageMatrix(const size_t& width,
                const size_t& height,
                const size_t& num_components);
    ImageMatrix(const size_t &width,
                const size_t &height,
                const size_t &num_components,
                JSAMPARRAY colormap);
    ~ImageMatrix(void);

    ImageMatrix(const ImageMatrix& other);
	ImageMatrix &operator= (ImageMatrix other);

    ImageMatrix(ImageMatrix&& other) noexcept;
//    ImageMatrix& operator= (ImageMatrix&& other) noexcept;

    void swap(ImageMatrix &first, ImageMatrix &second);

    inline size_t width(void) { return m_width; }
    inline size_t height(void) { return m_height; }
    inline size_t num_components(void) { return m_comps_n; }
    JSAMPLE getAt(const size_t &row, const size_t &col, const size_t &component) const;

private:
    size_t m_width{}, m_height{}, m_comps_n{};
    JSAMPARRAY m_colormap{};

    void m_free_colormap(void);
};


#endif
