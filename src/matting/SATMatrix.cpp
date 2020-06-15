#include "../../include/matting/SATMatrix.h"

SATMatrix::SATMatrix(const Eigen::MatrixXd &matrix) {
	size_t n_rows = matrix.rows(), n_cols = matrix.cols();
	m_sat.conservativeResize(n_rows, n_cols);
	m_sat(0, 0) = matrix(0, 0);
	for (size_t row = 1; row < n_rows; ++row)
		m_sat(row, 0) = m_sat(row - 1, 0) + matrix(row, 0);
	for (size_t col = 1; col < n_cols; ++col)
		m_sat(0, col) = m_sat(0, col - 1) + matrix(0, col);
	for (size_t row = 1; row < n_rows; ++row)
	{
		for (size_t col = 1; col < n_cols; ++col)
		{
			m_sat(row, col) = m_sat(row - 1, col) + m_sat(row, col - 1) -
							  m_sat(row - 1, col - 1) + matrix(row, col);
		}
	}
}

SATMatrix::SATMatrix(const SATMatrix &other) : m_sat(other.m_sat) { }

SATMatrix& SATMatrix::operator=(const SATMatrix& other)
{
	SATMatrix new_sat(other);
	swap(new_sat);
	return *this;
}

void SATMatrix::swap(SATMatrix& other)
{
	std::swap(m_sat, other.m_sat);
}

double SATMatrix::windowSum(const std::pair<size_t, size_t> &origin,
						    const std::pair<size_t, size_t> &end) const {
	size_t end_row = end.first, end_col = end.second,
		   orig_row = origin.first, orig_col = origin.second;
	double sum = m_sat(end_row, end_col);
	sum -= orig_row > 0 ? m_sat(orig_row - 1, end_col) : 0;
	sum -= orig_col > 0 ? m_sat(end_row, orig_col - 1) : 0;
	sum += (orig_col > 0 && orig_row > 0) ? m_sat(orig_row - 1, orig_col - 1) : 0;
	return sum;
}
