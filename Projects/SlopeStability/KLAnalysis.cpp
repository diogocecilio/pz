#include "KLAnalysis.h"

#include "pzstrmatrix.h"

#include "tpznodesetcompute.h"
#include "tpzsparseblockdiagonal.h"
#include "pzseqsolver.h"
#include "pzbdstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"


#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <sstream>
#include <set>
#include <random>
#include <limits>
//#include "pzfmatrix.h"
#include "pzmatrix.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzsolve.h"
#include "fadType.h"
#include "fad.h"

#include <exception>
#include "pzlog.h"
#include <complex>
#include <pzinterpolationspace.h>
#include <tpzintpoints.h>
using namespace std;


KLAnalysis::KLAnalysis ( TPZCompMesh *mesh ) : TPZAnalysis ( mesh )
{

  KLStrMatrix *mat = new KLStrMatrix ( mesh );
  mat->SetMesh ( mesh );
  fStrMatrix = mat;
  int numeq = fCompMesh->NEquations();
  fSolution.Redim ( numeq,1 );
  fSolution.Zero();
  LoadSolution();
}


KLAnalysis::~KLAnalysis ( void )
{

}


void KLAnalysis::Solve()
{

  TPZFMatrix<REAL> invB,invBC,B,C;
  MatrixXd eigenInvBC,eigenInvB,eigenC,eigenB;
  std::cout << "start assemble C  " << std::endl; 
  fStrMatrix->AssembleC ( C );
	std::cout << "end assembling C. Starting Assembling B  " << std::endl; 
  fStrMatrix->AssembleB ( B );
std::cout << "end assembling B.  " << std::endl;
  ToEigen ( B,eigenB );
  ToEigen ( C,eigenC );

  eigenInvB = eigenB.inverse();
  //std::cout << eigenInvB <<std::endl;
  eigenInvBC = eigenInvB*eigenC;

  MatrixXd val,vec;

  ComplexEigenSolver<MatrixXd> ces;
  std::cout << "Computing Eigenvalues  " << std::endl;
  ces.compute ( eigenInvBC );
  std::cout << "End Computing Eigenvalues  " << std::endl;


  int ncols = ces.eigenvectors().cols();
  int M = fExpansionOrder;
  //int M = ncols;
  int nrows = ces.eigenvectors().rows();
  val.resize ( M,1 );
  vec.resize ( nrows,M );
//   for ( int irow=0; irow<nrows; irow++ )
//   {
//     val ( irow,0 ) =ces.eigenvalues() [nrows-irow-1].real();
//     for ( int icol=0; icol<M; icol++ )
//       {
//         vec ( irow,icol ) =ces.eigenvectors() ( irow,ncols-icol-1 ).real();
//       }
//   }
	for ( int icol=0; icol< M; icol++ )
    {
	  val ( icol,0 ) =ces.eigenvalues() [nrows-icol-1].real();
      for(int irow=0;irow<nrows;irow++)
      {
        vec ( irow,icol ) =ces.eigenvectors() ( irow,ncols-icol-1 ).real();
      }
    }

  //std::cout <<"VT * D"   << "\n" << (vec* val.asDiagonal()* vec.inverse() ) << endl;
  bool print =true;
  if ( print==true )
    {
      std::cout << "Eigenvalues "<<endl;
      std::cout << val << endl;
      std::cout << "Eigenvectors "<<endl;
      std::cout << vec << endl;
     // std::cout << " A "<<endl;
     // std::cout << eigenInvBC <<std::endl;
     //std::cout << " A- V * D * V^(-1) = " << "\n" << (vec * val.asDiagonal() * vec.inverse()) << endl;
    }

    TPZFMatrix<REAL> vecpz(nrows,M);
    fEigenVectors.Resize(M);
	VectorXd intphisqr(M);
    for(int i=0;i<M;i++)fEigenVectors[i].Resize(nrows,1);
    for ( int icol=0; icol< M; icol++ )
    {
      for(int irow=0;irow<nrows;irow++)
      {
        fEigenVectors[icol](irow,0)=vec.col(icol)(irow);
      }

      //fEigenVectors[icol].Print(std::cout);
      LoadSolution(fEigenVectors[icol]);
      REAL integral = IntegrateSqrEigenVecSolution();
	  intphisqr(icol)=integral*integral;
      fEigenVectors[icol]*=1./integral;
    }
    
    for ( int icol=0; icol< M; icol++ )
    {
      for(int irow=0;irow<nrows;irow++)
      {
        vecpz(irow,icol)=fEigenVectors[icol](irow,0);
      }
    }
    LoadSolution(vecpz);
	
	REAL err = 0.;
	REAL err2=0.;
    for ( int i = 0; i < val.size(); i++ ) {
        err += fabs ( val(i)*intphisqr(i) ) ;
		err2 += fabs ( val(i) ) ;
    }
    
    REAL totalarea = ComputeTotalArea();
	std::cout << "total area = "<<totalarea << std::endl;
    //std::cout << " err / (fsig * fsig) = " <<err / (fsig * fsig) << std::endl;
    std::cout << "mean error 1 = " <<1. - 1./totalarea *   err2 << std::endl;
	std::cout << "mean error 2 = " <<1. - 1./totalarea *   err << std::endl;
	
	FromEigen(val,fEigenValues);
// 	
// 	vecpz.Print(std::cout);
	
}
REAL KLAnalysis::IntegrateSqrEigenVecSolution()
{


  int varid=6;//EVECSQR
  int nvar=1;
  int nels = fCompMesh->NElements(), iel;
  REAL sum=0.;
  TPZVec<REAL> vecsol ( nvar );


  for ( iel = 0; iel < nels; iel++ )
    {
      TPZCompEl * cel = fCompMesh->Element ( iel );

      vecsol = cel->IntegrateSolution ( varid );
      sum+=vecsol[0];
    }

  return sqrt ( sum );

}
void KLAnalysis::SetSolver ( TPZMatrixSolver<STATE> &solver )
{
  if ( fSolver ) delete fSolver;
  fSolver = ( TPZMatrixSolver<STATE> * ) solver.Clone();
}
/// Compute the area of the domain
REAL KLAnalysis::ComputeTotalArea()
{
    REAL area = 0.;
    long nelem = fCompMesh->NElements();
    TPZMaterial *pMatWithMem2 = fCompMesh->MaterialVec()[1];
	
    if (!pMatWithMem2) {
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCompMesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if (gel->MaterialId() != 1) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZIntPoints &rule = intel->GetIntegrationRule();
            int np = rule.NPoints();
            for (int ip = 0; ip<np; ip++) {
                TPZManVector<REAL,3> point(2,0.);
                REAL weight;
                rule.Point(ip, point, weight);
                TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
                TPZFNMatrix<9,REAL> axes(2,3);
                REAL detjac;
                gel->Jacobian(point, jac, axes, detjac, jacinv);
                area += weight*fabs(detjac);
            }
        }
    }
    return area;    
}
