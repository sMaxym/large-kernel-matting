#ifndef CFIMAGEMATTING_IMAGEMATRIX_H
#define CFIMAGEMATTING_IMAGEMATRIX_H

#include <iostream>
#include <Eigen/Core>

class ImageMatrix
{
	typedef std::pair<int, int> Point;
public:
    ImageMatrix() = default;
    ImageMatrix(const size_t& width,
                const size_t& height,
                const size_t& num_components);
    ~ImageMatrix(void);

    ImageMatrix(const ImageMatrix& other);
	ImageMatrix &operator= (ImageMatrix other);

    ImageMatrix(ImageMatrix&& other) noexcept;
//    ImageMatrix& operator= (ImageMatrix&& other) noexcept;

    void swap(ImageMatrix &first, ImageMatrix &second);

    inline size_t width(void) const { return m_width; }
    inline size_t height(void) const { return m_height; }
    inline size_t numComponents(void) const { return m_comps_n; }
    inline size_t area(void) const { return m_width * m_height; }
    inline Eigen::MatrixXd& colormap(void) { return m_colormap; }
    inline size_t flatCoords(const size_t &row, const size_t &col) const { return row * m_width + col; }
	inline size_t flatCoords(const Point &point) const { return flatCoords(point.first, point.second); }
    Point coords2D(const size_t &flat_coords) const;

	void setAt(const size_t &row, const size_t &col, const size_t &component, const double val);
    double getAt(const size_t &row, const size_t &col, const size_t &component) const;

    void toImageFormat(void);
    inline void normalize(void) { m_colormap /= KMax_brightness; };
	bool isNormalized(void);

private:
	const size_t KMax_brightness = 255.;
    size_t m_width{}, m_height{}, m_comps_n{};
    Eigen::MatrixXd m_colormap;
    bool m_is_normalized;

    void updateNormalized(void);
};


#endif
