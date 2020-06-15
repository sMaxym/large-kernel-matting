#ifndef CFIMAGEMATTING_IMAGEWINDOW_H
#define CFIMAGEMATTING_IMAGEWINDOW_H

#include <Eigen/Core>
#include <utility>
#include <algorithm>
#include <cmath>

// TODO add copy/move constructors and operators
// TODO const constraints for methods
class ImageWindow {
	typedef std::pair<int, int> Point;
public:
	ImageWindow() = default;
	ImageWindow(const Point &image_shape, const Point &center, int radius);

	ImageWindow(const ImageWindow& other);
	ImageWindow& operator=(const ImageWindow& other);

	void swap(ImageWindow& other);

	class iterator
	{
	public:
		iterator(Point ptr, Point upper_b, Point lower_b)
			: m_position(std::move(ptr)), m_upper_b(std::move(upper_b)), m_lower_b(std::move(lower_b)) { }
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
		return iterator({m_end_bound.first, m_origin_bound.second},
				     m_origin_bound,
				     m_end_bound);
	}

	bool inBounds(const Point &coord);
	size_t getArea();

	inline Point getImageShape() { return m_im_shape; }
	inline Point getCenter() { return m_center; }
	inline int getRadius() { return m_radius; }
	inline Point getOriginBound() { return m_origin_bound; }
	inline Point getEndBound() { return m_end_bound; }
	inline Point getEndBoundIncl() { return {m_end_bound.first - 1, m_end_bound.second - 1}; }

private:
	Point m_im_shape, m_center, m_origin_bound, m_end_bound; // end bound not included
	int m_radius;
};

#endif
