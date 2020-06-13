#include "../../include/cfmatting/ImageWindow.h"

ImageWindow::ImageWindow(const Point &image_shape, const Point &center, int radius)
{
	m_im_shape = image_shape;
	m_center = center;
	m_radius = radius;

	m_origin_bound = { std::max<int>(0, m_center.first - m_radius),
					 std::max<int>(0, m_center.second - m_radius) };
	m_end_bound = { std::min<int>(image_shape.first, m_center.first + m_radius + 1),
					std::min<int>(image_shape.second, m_center.second + m_radius + 1)};
}

bool ImageWindow::inBounds(const Point &coord)
{
	// TODO DRY
	if (coord.first >= m_im_shape.first || coord.first < 0)
		return false;
	return !(coord.second >= m_im_shape.second || coord.second < 0);
}

size_t ImageWindow::getArea()
{
	Point window_shape = {m_end_bound.first - m_origin_bound.first,
						  m_end_bound.second - m_origin_bound.second};
	return window_shape.first * window_shape.second;
}

ImageWindow::iterator::~iterator() = default;

ImageWindow::iterator ImageWindow::iterator::operator++()
{
	int cur_col = m_position.second;
	if (cur_col == m_lower_b.second - 1)
	{
		m_position.second = m_upper_b.second;
		++m_position.first;
	} else
	{
		++m_position.second;
	}
	return *this;
}
