#ifndef CFIMAGEMATTING_SATMATRIX_H
#define CFIMAGEMATTING_SATMATRIX_H

#include <Eigen/Core>

class SATMatrix {
public:
	SATMatrix() = default;
	SATMatrix(const Eigen::MatrixXd &matrix);

	SATMatrix(const SATMatrix& other);
	SATMatrix& operator=(const SATMatrix& other);

	void swap(SATMatrix& other);

	double windowSum(const std::pair<size_t, size_t> &origin,
				     const std::pair<size_t, size_t> &end) const;

private:
	Eigen::MatrixXd m_sat;
};



#endif
