#ifndef CFIMAGEMATTING_CGMATTINGSOLVER_H
#define CFIMAGEMATTING_CGMATTINGSOLVER_H

#include <utility>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/LU>


#include "ImageMatrix.h"
#include "ImageWindow.h"
#include "SATMatrix.h"

class CGMattingSolver {
	typedef Eigen::VectorXd Vector;
	typedef std::pair<int, int> Point;
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdRow;
public:
	CGMattingSolver(ImageMatrix  image, ImageMatrix  trimap, const size_t kernel_radius=10);

	Vector alphaMatting(const size_t iterations=20, const double precision=1.e-6);

	inline void setConstraintCoeff(const double coeff) { m_constraint = coeff; }
	inline void setRegParameter(const double eps) { m_reg_param = eps; }


private:
	ImageMatrix m_image, m_trimap;
	std::vector<SATMatrix> m_image_sat;
	Eigen::SparseMatrix<double> m_constraint_mat;
	std::vector<std::vector<Eigen::MatrixXd>> m_win_cov;
	Vector m_cg_solution, m_conjugate, m_residual, m_def_constraint;
	double m_constraint = 65536.0, m_reg_param = 1.;
	size_t m_kernel_radius;


	void calcImageSAT(void);
	void CGIterate(void);
	void calcImageCovariance(void);
	std::vector<SATMatrix> calcSlopeAndBiasSAT(const std::vector<SATMatrix>& Ip_sat,
												const SATMatrix& conjugate_sat);
	Vector laplacianProduct(const std::vector<SATMatrix>& slope_sat, const SATMatrix& bias_sat);
};

#endif
