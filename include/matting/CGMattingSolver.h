#ifndef CFIMAGEMATTING_CGMATTINGSOLVER_H
#define CFIMAGEMATTING_CGMATTINGSOLVER_H

#include <iostream>
#include <utility>
#include <exception>
#include <vector>
#include <Eigen/LU>
#include <Eigen/SparseCore>

#include "ImageMatrix.h"
#include "ImageWindow.h"
#include "SATMatrix.h"

class CGMattingSolver {
	typedef Eigen::VectorXd Vector;
	typedef std::pair<int, int> Point;
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdRow;
public:
	CGMattingSolver(const ImageMatrix& image, const ImageMatrix& trimap, const size_t kernel_radius=10);

	CGMattingSolver(const CGMattingSolver& other) = delete;
	CGMattingSolver& operator=(const CGMattingSolver& other) = delete;

	Vector alphaMatting(const size_t iterations=20, const double precision=1.e-6);

	void setConstraintCoeff(const double coeff) { m_constraint = coeff; }
	void setRegParameter(const double eps) { m_reg_param = eps; }


private:
	const double KThreshold = 0.01;

	ImageMatrix m_image, m_trimap;
	std::vector<SATMatrix> m_image_sat;
	Eigen::SparseMatrix<double> m_constraint_mat;
	std::vector<std::vector<Eigen::MatrixXd>> m_win_cov;
	Vector m_cg_solution, m_conjugate, m_residual, m_def_constraint;
	double m_constraint = 65536.0, m_reg_param = 1.;
	size_t m_kernel_radius;


	void initConstraints();
	void initImageCovariance();
	Eigen::MatrixXd calcCov(const size_t row, const size_t col);
	void initImageSAT();
	void CGIterate();
	std::vector<SATMatrix> calcSlopeAndBiasSAT(const std::vector<SATMatrix>& Ip_sat,
												const SATMatrix& conjugate_sat);
	Vector laplacianProduct(const std::vector<SATMatrix>& slope_sat, const SATMatrix& bias_sat);

	double handleZeroDiv(double value) { return value == 0 ? value + KThreshold : value; }
};

#endif
