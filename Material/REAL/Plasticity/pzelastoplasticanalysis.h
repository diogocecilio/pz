/**
 * @file
 * @brief Contains the declaration of TPZElastoPlasticAnalysis class
 */

#ifndef ELASTOPLASTICANALYSIS_H
#define ELASTOPLASTICANALYSIS_H

#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "TPZGeoElement.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzpostprocanalysis.h"
#include <iostream>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/LU>
#include <Eigen/Core>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
using namespace Eigen;
class TPZElastoPlasticAnalysis : public TPZNonLinearAnalysis {

public:
	/** @brief Constructor */
	TPZElastoPlasticAnalysis(TPZCompMesh *mesh,std::ostream &out=std::cout);
	/** @brief Default constructor */
	TPZElastoPlasticAnalysis();
	/** @brief Default destructor */
	virtual ~TPZElastoPlasticAnalysis();
	
// 	/**
// 	 * @brief It process a Newton's method to solve the non-linear problem.
// 	 * In this implementation, the line search is temporarily disabled.
// 	 */
// 	virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv,bool &ConvOrDiverg);
// 
//     	virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv);
//     
//     /// Iterative process using the linear elastic material as tangent matrix
//     virtual void IterativeProcess(std::ostream &out, TPZAutoPointer<TPZMatrix<STATE> > linearmatrix, REAL tol, int numiter, bool linesearch);
    
	//virtual REAL LineSearch(const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter);
    
    //Improved:A verification is made in order to check convergence.
    virtual REAL LineSearch(const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter, bool &converging );
	

	/**
	 * @brief The code below manages the update of a certain boundary condition (BCId)
	 * to assume values progressing from begin to end values within nsteps, such
	 * that the last step length is a 'lastStepRatio' fraction of the global
	 * load in a geometric progression.
	 * @param out [in] output device to write the output to
	 * @param tol [in] tolerance desired in each loading step
	 * @param numiter [in] maximum iteration steps in each loading step
	 * @param BCId [in] boundary condition Id
	 * @param nsteps [in] number of steps
	 * @param PGRatio [in] ratio for the PG progression
	 * @param val1Begin [in] beginning value of BC.Val1
	 * @param val1End [in] final value of BC.Val1
	 * @param val2Begin [in] beginning value of BC.Val2
	 * @param val2End [in] final value of BC.Val2
	 * @param ppAnalysis [in] TPZPostProcAnalysis object to write the output to
	 * @param res [in] output postprocessing resolution
	 */

	
	/**
	 * @brief Informs the materials to update the plastic memory, assembles the
	 * rhs in order to update the memory and unsets the materials
	 * to update the memory.
	 * @param ResetOutputDisplacements [in] Informs whether to add or reset the deltaSolution to the cumulative solution.
	 */
	virtual REAL AcceptSolution(const int ResetOutputDisplacements = 0);
    
    /** @brief Load the solution into the computable grid, transfering it to the multi physics meshes */
    virtual void LoadSolution();
	
	virtual void LoadSolution(TPZFMatrix<STATE> & loadsol);
	
	virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv);
	void IterativeProcess(std::ostream &out,REAL tol,int numiter);
	
	bool IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv,int &iters);
	
	
	//Implements the cylindrical arc length method ref- Souza Neto 2009
	virtual void IterativeProcessArcLength(std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv);

    REAL LineSearch(const TPZFMatrix<STATE> &Wn, TPZFMatrix<STATE> DeltaW, TPZFMatrix<STATE> &NextW, REAL tol, int niter);
	
    TPZFMatrix<REAL> &CumulativeSolution()
    {
        return fCumSol;
    }
	
	void SetPrecond(TPZMatrixSolver<REAL> &precond);
	
	void SetBiCGStab(int numiter, REAL tol);
	
	void SetBiCGStab_Jacobi(int numiter, REAL tol);
	
	void SetLU();
	
	void TransferSolution(TPZPostProcAnalysis & ppanalysis);
	
    void AddNoPenetration(int matid, int direction)
    {
        fMaterialIds.insert(std::pair<int,int>(matid, direction));
    }
    
    /// build the fEquationstoZero datastructure based on the fMaterialIds data structure
    void IdentifyEquationsToZero();
    
    /// return the vector of active equation indices
    void GetActiveEquations(TPZVec<long> &activeEquations);
    
protected:	
	
	/**
	 * @brief Forces the materials with memory to update the internal
	 * plastic memory during the subsequent assemble/contribute calls
	 * when update is set to 1. When set to 0, the conventional
	 * contribute routine takes place
	 * @param update [in] 1 - updates the material memory during the next assemble calls; 0 - conventional assemble
	 */
	void SetUpdateMem(int update);
	
	/**
	 * @brief Updates block diagonal preconditioning matrix.
	 */
	void UpdatePrecond();
	
	
public:
	
	void CheckConv(std::ostream &out, REAL range);
	
	virtual void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase);
	
	virtual int NumCases();
	
	virtual void Residual(TPZFMatrix<REAL> &residual, int icase);
	
	static void SetAllCreateFunctionsWithMem(TPZCompMesh *cmesh);
	
    /// Structure defining a multiphysics configuration
    bool IsMultiPhysicsConfiguration()
    {
        return fMultiPhysics != NULL;
    }
    
    /// Initialize the multiphysics data structure
    void SetMultiPhysics(TPZCompMesh *mphysics, TPZVec<TPZCompMesh *> &meshvec)
    {
        fMultiPhysics = mphysics;
        fMeshVec = meshvec;
    }
    
    /// Reset the multiphysics data structure
    void ResetMultiPhysics()
    {
        fMultiPhysics = 0;
        fMeshVec.Resize(0);
    }
    
    
 
inline void SolveEigenSparse ( int type, TPZAutoPointer<TPZMatrix<REAL> > A, TPZFMatrix<REAL> b, TPZFMatrix<REAL>& x )
{

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    int sz=A->Rows();


    tripletList.reserve ( sz*10000 );
    // tripletList.reserve(80000);

    x.Resize ( sz, 1);
    SparseMatrix<double> AA ( sz, sz );
    VectorXd bbb ( sz );
    for ( int i = 0; i < sz; i++ ) {
        for ( int j = 0; j < sz; j++ ) {
            if ( fabs ( A->Get(i,j) ) >1.e-12) {
                tripletList.push_back ( T ( i,j,A->Get(i,j) ) );
            }
        }
        bbb ( i ) = b(i,0);
    }

    AA.setFromTriplets ( tripletList.begin(), tripletList.end() );

    AA.makeCompressed();

    VectorXd xx;
    if ( type==0 ) {
        SimplicialLLT< SparseMatrix<double> > solver;
		solver.analyzePattern ( AA );
        solver.factorize ( AA );
		xx = solver.solve ( bbb );
       // xx = solver.compute ( AA ).solve ( bbb );
    } else if ( type==1 ) {
        SparseLU< SparseMatrix<double> > solver2;
        solver2.analyzePattern ( AA );
        solver2.factorize ( AA );
        xx = solver2.solve ( bbb );
    } else if ( type==2 ) {
        ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
		cg.analyzePattern ( AA );
        cg.compute ( AA );
        xx = cg.solve ( bbb );
        //std::cout << "#iterations:     " << cg.iterations() << std::endl;
        //std::cout << "estimated error: " << cg.error()      << std::endl;
    }
    else if ( type==3 ) {
        SimplicialLDLT<SparseMatrix<double> > solver;
		solver.analyzePattern ( AA );
        solver.factorize ( AA );
		xx = solver.solve ( bbb );
    }
    else if(type ==4)
	{
		PartialPivLU<SparseMatrix<double> > solver3;
		//solver3.analyzePattern ( AA );
        //solver3.factorize ( AA );
		xx = solver3.solve ( bbb );
	}
	else if(type==5)
 	{
  		BiCGSTAB<SparseMatrix<double> > solver4;
		solver4.analyzePattern ( AA );
  		solver4.compute(AA);
  		xx = solver4.solve(bbb);
	}
    for ( int i=0; i<sz; i++ ) {
        x(i,0)=xx ( i );
    }
}
	

void SolveEigen ( TPZAutoPointer<TPZMatrix<REAL> > A, TPZFMatrix<REAL> b, TPZFMatrix<REAL>& x )
{

    x.Resize ( A->Rows(), 1 );
    MatrixXd AA ( A->Rows(), A->Rows() );
    VectorXd bbb ( A->Rows());
    for ( int i = 0; i < A->Rows(); i++ ) {
        for ( int j = 0; j < A->Rows(); j++ ) {
            AA ( i, j ) = A->Get(i,j);
        }
        bbb ( i ) = b(i,0);
    }
    LLT<MatrixXd> llt;
    //LU<MatrixXd> lu;
    llt.compute ( AA );
    VectorXd xxx=llt.solve ( bbb );
    for ( int i = 0; i < A->Rows(); i++ ) {
        x(i,0) = xxx ( i );
    }
}
	
void FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat )
    {

        int rows = eigenmat.rows();
        int cols = eigenmat.cols();
        pzmat.Resize ( rows,cols );
        for ( int irow=0; irow<rows; irow++ ) {
            for ( int icol=0; icol<cols; icol++ ) {
                pzmat ( irow,icol ) =eigenmat ( irow,icol );
            }
        }

    }

    void ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat )
    {
        TPZFMatrix<REAL> intpz ( pzmat );
        int rows = pzmat.Rows();
        int cols = pzmat.Cols();
        eigenmat.resize ( rows,cols );
        for ( int irow=0; irow<rows; irow++ ) {
            for ( int icol=0; icol<cols; icol++ ) {
                eigenmat ( irow,icol ) =intpz ( irow,icol );
            }
        }

    }
protected:
	
	/* @brief Cumulative solution vector*/
	TPZFMatrix<REAL> fCumSol;
	
	TPZMatrixSolver<REAL> * fPrecond;
    
    /// Equations with zero dirichlet boundary condition
    std::set<long> fEquationstoZero;
    
    /// Materials with no penetration boundary conditions
    // the second value of the map indicates x (0) or y (1) restraint
    std::multimap<int,int> fMaterialIds;
    
    /// The multiphysics mesh
    TPZCompMesh *fMultiPhysics;
    
    /// The elasticity mesh and vertical deformation mesh
    TPZManVector<TPZCompMesh *,2> fMeshVec;
	
	/**
	 * TPZCompElWithMem<TBASE> creation function setup
	 * These functions should be defined as static members of the TPZCompElWithMem<TBASE> class
	 * but were defined here because the TPZCompElWithMem is a template class and would
	 * require a dummy template argumet in order to be called. That woudn't be elegant.
	 */
	
	static TPZCompEl * CreateCubeElWithMem(  TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	static TPZCompEl * CreateLinearElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	static TPZCompEl * CreatePointElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	static TPZCompEl * CreatePrismElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	static TPZCompEl * CreatePyramElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	static TPZCompEl * CreateQuadElWithMem(  TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	static TPZCompEl * CreateTetraElWithMem( TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	static TPZCompEl * CreateTriangElWithMem(TPZGeoEl *gel, TPZCompMesh &mesh, long &index);
	
};

#endif

