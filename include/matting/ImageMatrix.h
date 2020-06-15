#ifndef CFIMAGEMATTING_IMAGEMATRIX_H
#define CFIMAGEMATTING_IMAGEMATRIX_H

#include <exception>
#include <Eigen/Core>

class ImageMatrix
{
	typedef std::pair<int, int> Point;
public:
    ImageMatrix() = default;
    ImageMatrix(const size_t& width,
                const size_t& height,
                const size_t& num_components);
    ImageMatrix(const Eigen::MatrixXd& matrix);
	ImageMatrix(const Eigen::VectorXd& matrix,
				const size_t& width,
				const size_t& height);
    ~ImageMatrix();

    ImageMatrix(const ImageMatrix& other);
	ImageMatrix &operator= (ImageMatrix other);

    void swap(ImageMatrix &first, ImageMatrix &second);

    inline size_t width() const { return m_width; }
    inline size_t height() const { return m_height; }
    inline size_t numComponents() const { return m_comps_n; }
    inline size_t area() const { return m_width * m_height; }
    inline Eigen::MatrixXd& colormap() { return m_colormap; }
    inline size_t flatCoords(const size_t &row, const size_t &col) const { return row * m_width + col; }
	inline size_t flatCoords(const Point &point) const { return flatCoords(point.first, point.second); }
    Point coords2D(const size_t &flat_coords) const;

	void setAt(const size_t &row, const size_t &col, const size_t &component, const double val);
    double getAt(const size_t &row, const size_t &col, const size_t &component) const;

    ImageMatrix grayscale() const;
	void expandColorspace();

    void toImageFormat();
    inline void normalize() { m_colormap /= KMax_brightness; };
    inline void normalizeIfNot() { if (!m_is_normalized) normalize(); }
	bool isNormalized();

private:
	const size_t KMax_brightness = 255.;
    size_t m_width{}, m_height{}, m_comps_n{};
    Eigen::MatrixXd m_colormap;
    bool m_is_normalized{};

    void updateNormalized();
    void assertCoords(const size_t row, const size_t col, const size_t comp) const;
};

#endif
