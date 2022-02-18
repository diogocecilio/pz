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
  fStrMatrix->AssembleC ( C );
  fStrMatrix->AssembleB ( B );

  ToEigen ( B,eigenB );
  ToEigen ( C,eigenC );

  eigenInvB = eigenB.inverse();
  //std::cout << eigenInvB <<std::endl;
  eigenInvBC = eigenInvB*eigenC;

  MatrixXd val,vec;

  ComplexEigenSolver<MatrixXd> ces;
  ces.compute ( eigenInvBC );
  


  int ncols = ces.eigenvectors().cols();
  int nrows = ces.eigenvectors().rows();
  val.resize ( nrows,1 );
  vec.resize ( nrows,ncols );
  for ( int irow=0; irow<nrows; irow++ )
  {
    val ( irow,0 ) =ces.eigenvalues() [irow].real();
    for ( int icol=0; icol<ncols; icol++ )
      {
        vec ( irow,icol ) =ces.eigenvectors() ( irow,icol ).real();
      }
  }
  
  bool print =false;
  if ( print==true )
    {
      std::cout << "Eigenvalues "<<endl;
      std::cout << val << endl;
      std::cout << "Eigenvectors "<<endl;
      std::cout << vec << endl;
      std::cout << " A "<<endl;
      std::cout << eigenInvBC <<std::endl;
      std::cout << " V * D * V^(-1) = " << endl
                << eigenInvBC-(vec * val.asDiagonal() * vec.inverse()) << endl;
    }

    TPZFMatrix<REAL> vecpz(nrows,ncols);
    fEigenVectors.Resize(ncols);
    for(int i=0;i<ncols;i++)fEigenVectors[i].Resize(nrows,1);
    for ( int icol=0; icol< ncols; icol++ )
    {
      for(int irow=0;irow<nrows;irow++)
      {
        fEigenVectors[icol](irow,0)=vec.col(icol)(irow);
      }

      //fEigenVectors[icol].Print(std::cout);
      LoadSolution(fEigenVectors[icol]);
      REAL integral = IntegrateSqrEigenVecSolution();
      fEigenVectors[icol]*=1./integral;
    }
    
    for ( int icol=0; icol< ncols; icol++ )
    {
      for(int irow=0;irow<nrows;irow++)
      {
        vecpz(irow,icol)=fEigenVectors[icol](irow,0);
      }
    }
    LoadSolution(vecpz);
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
//     int nels = fCompMesh->NElements(), iel;
//     const TPZIntPoints &intrule = this->GetIntegrationRule();
//     for ( iel = 0; iel < nels; iel++ ) {
//
//       TPZCompEl * cel = fCompMesh->Element(iel);
//       const TPZIntPoints &intrule = cel->GetIntegrationRule();
//         int npts = intrule.nrows();
//         for ( Int ipt = 0; ipt < npts; ipt++ ) {
//             xi = intrule[ipt][0];
//             eta = intrule[ipt][1];
//             w = intrule[ipt][2];
//             fshape->shapes ( psis, GradPsi, xi, eta );
//             GetElCoords ( allcoords, iel, elcoords );
//             GradPsi.Mult ( elcoords, Jac );
//             Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
//             for ( Int inode = 0; inode < psis.nrows(); inode++ ) {
//                 sum += psis[inode][0] * pow ( Vec[meshtopology[iel][inode]][0], 2 ) *DetJ*w;
//             }
//         }
//     }
//     solu = sqrt ( sum );
//     cout << "\n Inte = " << solu << endl;
//
//     return solu;



}
void KLAnalysis::SetSolver ( TPZMatrixSolver<STATE> &solver )
{
  if ( fSolver ) delete fSolver;
  fSolver = ( TPZMatrixSolver<STATE> * ) solver.Clone();
}

