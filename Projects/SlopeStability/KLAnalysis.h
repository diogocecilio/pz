/**
 * @file
 * @brief Contains KLAnalysis class which implements the sequence of actions to perform a finite element analysis.
 */

#ifndef KLANALYSIS_H
#define KLANALYSIS_H
#include "KLStrMatrix.h"
#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>
class TPZCompMesh;

template<class TVar> 
class TPZFMatrix;
#include <Eigen/Dense>
 
using namespace std;
using namespace Eigen;
class TPZAnalysis;

/**
 * @ingroup analysis
 * @brief Implements the sequence of actions to perform a finite element analysis. \ref analysis "Analysis"
 */
/** This class will renumerate the nodes upon construction
 */
class KLAnalysis: TPZAnalysis {
		
	public :

	KLAnalysis(TPZCompMesh *mesh);
	
	KLAnalysis();
	virtual ~KLAnalysis(void);

	TPZFMatrix<STATE> &Rhs() { return fRhs;}

	TPZFMatrix<STATE> &Solution() { return fSolution;}
	

	TPZCompMesh *Mesh()const { return fCompMesh;}

	
	void SetSolver(TPZMatrixSolver<STATE> &solver);
	
	void AssembleC();
	
	void AssembleB();
	
	void Solve();
	
	void FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat )
    {

		int rows = eigenmat.rows();
		int cols = eigenmat.cols();
		pzmat.Resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				pzmat(irow,icol)=eigenmat(irow,icol);
			}
		}
	
    }

	void ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat )
    {
		TPZFMatrix<REAL> intpz(pzmat);
		int rows = pzmat.Rows();
		int cols = pzmat.Cols();
		eigenmat.resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				eigenmat(irow,icol)=intpz(irow,icol);
			}
		}
	
    }
    
    REAL IntegrateSqrEigenVecSolution();
	
	TPZFMatrix<REAL> GetItEigenVec(int it)
	{
		return fEigenVectors[it];
	}
	
	REAL GetItEigenVal(int it)
	{
		return fEigenValues(it,0);
	}

	
// 	virtual void DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile)=0;
// 
// 	virtual void PostProcess(int resolution)=0;
	
private:
	KLStrMatrix * fStrMatrix; 
	TPZVec<TPZFMatrix<REAL>> fEigenVectors;
	TPZFMatrix<REAL> fEigenValues;

};


#endif
