#ifndef CFIMAGEMATTING_SATMATRIX_H
#define CFIMAGEMATTING_SATMATRIX_H

#include <Eigen/Core>
#include <iostream>

// TODO template Eigen matrix
class SATMatrix {
public:
	SATMatrix() = default;
	SATMatrix(const Eigen::MatrixXd &matrix);

	double windowSum(const std::pair<size_t, size_t> &origin,
				     const std::pair<size_t, size_t> &end);

	inline void print(void) { std::cout << m_sat << std::endl; }

private:
	// TODO T** type m_sat
	Eigen::MatrixXd m_sat;
};



#endif //CFIMAGEMATTING_SATMATRIX_H
