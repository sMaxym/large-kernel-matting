#ifndef CFIMAGEMATTING_IMAGEWINDOW_H
#define CFIMAGEMATTING_IMAGEWINDOW_H

#include <Eigen/Core>
#include <utility>

typedef Eigen::Vector2i Point;

// TODO add copy/move constructors and operators
// TODO const constraints for methods
class ImageWindow {
public:
	ImageWindow() = default;
	ImageWindow(const Point &image_shape, const Point &center, const size_t radius);

	class iterator
	{
	public:
		iterator(Point ptr, Point upper_b, Point lower_b)
			: m_position(ptr), m_upper_b(upper_b), m_lower_b(lower_b) { }
		~iterator();
		// TODO postfix iterator increment
		iterator operator++();
		inline Point &operator*() { return m_position; }
		inline Point *operator->() { return &m_position; }
		bool operator==(const iterator& rhs) { return m_position == rhs.m_position; }
		bool operator!=(const iterator& rhs) { return !(*this == rhs); }
	private:
		Point m_position, m_upper_b, m_lower_b;
	};

	inline iterator begin() { return iterator(m_origin_bound, m_origin_bound, m_end_bound); }
	inline iterator end() {
		Point end_position;
		end_position << m_end_bound[0], m_origin_bound[1];
		return iterator(end_position, m_origin_bound, m_end_bound);
	}

	bool inBounds(const Point &coord);
	size_t getArea(void);

	inline size_t flatCoords(const Point &coord) { return m_im_shape[1] * coord[0] + coord[1]; }
	inline Point getImageShape(void) { return m_im_shape; }
	inline Point getCenter(void) { return m_center; }
	inline size_t getRadius(void) { return m_radius; }
	inline Point getOriginBound(void) { return m_origin_bound; }
	inline Point getEndBound(void) { return m_end_bound; }

private:
	// m_end_bound not included
	Point m_im_shape, m_center, m_origin_bound, m_end_bound;
	size_t m_radius;
};


#endif
