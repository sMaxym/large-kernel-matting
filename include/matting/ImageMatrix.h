#ifndef CFIMAGEMATTING_IMAGEMATRIX_H
#define CFIMAGEMATTING_IMAGEMATRIX_H

#include <exception>
#include <Eigen/Core>

class ImageMatrix
{
	typedef std::pair<int, int> Point;
public:
    ImageMatrix() = default;
    ImageMatrix(const size_t width,
                const size_t height,
                const size_t num_components);
    ImageMatrix(const Eigen::MatrixXd& matrix);
	ImageMatrix(const Eigen::VectorXd& matrix,
				const size_t width,
				const size_t height);
    ~ImageMatrix();

    ImageMatrix(const ImageMatrix& other);
	ImageMatrix& operator=(const ImageMatrix& other);

    void swap(ImageMatrix &first, ImageMatrix &second);

    size_t width() const { return m_width; }
    size_t height() const { return m_height; }
    size_t numComponents() const { return m_comps_n; }
    size_t area() const { return m_width * m_height; }
    Eigen::MatrixXd& colormap() { return m_colormap; }
    size_t flatCoords(const size_t row, const size_t col) const { return row * m_width + col; }
	size_t flatCoords(const Point &point) const { return flatCoords(point.first, point.second); }
    Point coords2D(const size_t flat_coords) const;

	void setAt(const size_t row, const size_t col, const size_t component, const double val);
    double getAt(const size_t row, const size_t col, const size_t component) const;

    ImageMatrix grayscale() const;
	void expandColorspace();

    void toImageFormat();
    void normalize() { m_colormap /= KMax_brightness; };
    void normalizeIfNot() { if (!isNormalized()) normalize(); }
	bool isNormalized() const;

private:
	const double KMax_brightness = 255.;
    size_t m_width{}, m_height{}, m_comps_n{};
    Eigen::MatrixXd m_colormap;
};

#endif
