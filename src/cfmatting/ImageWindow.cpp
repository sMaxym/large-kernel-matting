#include "../../include/cfmatting/ImageWindow.h"

ImageWindow::ImageWindow(const Point &image_shape, const Point &center, size_t radius)
{
	m_im_shape = image_shape;
	m_center = center;
	m_radius = radius;

	m_origin_bound = Point::Zero().cwiseMax(m_center - m_radius * Point::Ones());
	m_end_bound = image_shape.cwiseMin(m_center + (m_radius + 1) * Point::Ones());
}

bool ImageWindow::inBounds(const Point &coord)
{
	// TODO DRY
	if (coord[0] >= m_im_shape[0] || coord[0] < 0)
		return false;
	return !(coord[1] >= m_im_shape[1] || coord[1] < 0);
}

size_t ImageWindow::getArea()
{
	Point window_shape = m_end_bound - m_origin_bound;
	return window_shape[0] * window_shape[1];
}

ImageWindow::iterator::~iterator() = default;

ImageWindow::iterator ImageWindow::iterator::operator++()
{
	Point next_position;
	int cur_col = m_position[1];
	if (cur_col == m_lower_b[1] - 1)
	{
		m_position[1] = m_upper_b[1];
		++m_position[0];
	} else
	{
		++m_position[1];
	}
	return *this;
}
