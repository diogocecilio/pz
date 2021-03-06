/**
 * @file
 * @brief Contains TPZCopySolve class which solves clones of the given matrix.
 */
/* Generated by Together */

#ifndef TPZCOPYSOLVE_H
#define TPZCOPYSOLVE_H
#include "pzsolve.h"

template <class TVar>
class TPZFMatrix;

/**
 * @ingroup solver
 * @brief To solve clones of the given matrix. \ref solver "Solver"
 */
template<class TVar>
class TPZCopySolve: public TPZMatrixSolver<TVar>
{
public:
	
	TPZCopySolve(TPZMatrix<TVar> *mat) :
    TPZMatrixSolver<TVar>(mat)
	{
	}
	TPZCopySolve(const TPZCopySolve &other) :
    TPZMatrixSolver<TVar>(other)
	{
	}
	/**
	 * @brief Solves the system of linear equations stored in current matrix \n
	 * As this class implements only a copy operation, it just copies u to F;
	 * @param F contains Force vector
	 * @param result contains the solution
	 * @param residual [out] residual computed
	 */
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual)
	{
		result = F;
	}
	
	/** @brief Clones the current object returning a pointer of type TPZSolver */
	TPZSolver<TVar> *Clone() const
	{
		return new TPZCopySolve(*this);
	}
};

#endif //TPZCOPYSOLVE_H
