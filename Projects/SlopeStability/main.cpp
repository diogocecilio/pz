#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "KLMaterial.h"
#include "KLStrMatrix.h"
#include "KLAnalysis.h"
#include "KLInterpolationSpace.h"
//#include "TPZPlasticityTest.h"
#include <iostream>
#include <cstdlib>
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZTensor.h"
#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "pzfstrmatrix.h"
#include "pzbndmat.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "tpzgeoelrefpattern.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZTensor.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
#include "TPZDruckerPrager.h"
#include "KLRandomField.h"
#include "pzelastoplastic2D.h"
#include <pzmathyperelastic.h>
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"
#include "clock_timer.h"

#include "TPZVTKGeoMesh.h"
using namespace pzshape; // needed for TPZShapeCube and related classes


#include "pzlog.h"
#include "SlopeStabilityAnalysis.h"
#include "pzbfilestream.h"
#include "TPZProjectEllipse.h"
#include "arglib.h"
#include "run_stats_table.h"
#include "pzmganalysis.h"
#include "pzcompel.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pztransfer.h"
#include "tpzgeoelrefpattern.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzsolve.h"
#include "pzgengrid.h"
#include "tpzarc3d.h"

#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

#include "TPZVTKGeoMesh.h"


#include "TPZGeoLinear.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "tpzquadraticline.h"
#include "TPZWavyLine.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"

#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "tpzcompmeshreferred.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include "TPZGeoElement.h"
#include "pzgeoel.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzskylstrmatrix.h"

#include <time.h>
#include <stdio.h>
#include "pzl2projection.h"
#include "tpzgeoelmapped.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>

#include "pzgeotriangle.h"
#include "tpzarc3d.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoelrefpattern.h.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "TPZGeoCube.h"
#include <pzcompel.h>

#include <sstream>
using namespace std;
using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
#include "pzlog.h"
#include "tpzarc3d.h"
#include "pzelasmat.h"
using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;
#ifdef LOG4CXX
static LoggerPtr logger ( Logger::getLogger ( "plasticity.main" ) );
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse ( Logger::getLogger ( "LogEllipse" ) );
#endif

void PrintSol ( TPZCompMesh * plasticCmesh );
TPZGeoMesh  ConfigSlope();
TPZCompMesh* CreateKLCMesh ( TPZGeoMesh * gmesh,int porder,REAL Lx, REAL Ly,REAL Lz, int type, int expansionorder );
TPZCompMesh * KLConfig ( TPZSlopeStabilityAnalysis &slopeA,int porder,int ref );
TPZVec<TPZCompMesh *> CreateFields ( TPZSlopeStabilityAnalysis &slopeA,int porder,int ref );
void SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZVec<TPZCompMesh*> vecmesh,int isol,REAL factor );

TPZGeoMesh * CreteGeoMesh ( TPZSlopeStabilityAnalysis &slopeA,int porder,int ref );
TPZCompMesh * CreteCompMesh ( TPZGeoMesh* gmesh,int porder );
TPZGeoMesh * SimpleGMesh ( int ref );
TPZCompMesh * SimpleCMesh ( TPZGeoMesh * gmesh, int porder );
void SolveSimpleElasticProlem();
void ComputeAdptiveMontecarlo();
void   ConfigSlopeShearRed();
void ShearRed ( TPZSlopeStabilityAnalysis &slopeA, TPZVec<TPZCompMesh *> vecmesh );
void ShearRed ( TPZSlopeStabilityAnalysis &slopeB );
void ShearRed ( );
void Refinine ( TPZSlopeStabilityAnalysis &slope,REAL sqj2_refine,int ref,int order );
TPZSlopeStabilityAnalysis CreateSlopeObj ( REAL porder,REAL bodyforcefac, int ref );
typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;

int main ( int argc, char **argv )
{


//     int porder=1;
//     int ref=4;
//     TPZSlopeStabilityAnalysis slopeA;
// 	std::string output = "slopeA.vtk";
//     slopeA.SetVtkOutPutName ( output );
//
//
//
//
//     TPZVec<TPZCompMesh *> vecmesh;
//     vecmesh=CreateFields ( slopeA, porder, ref );
//
//
//     slopeA.GetCurrentConfig()->CreatePostProcessingMesh();
//
//     slopeA.PostProcess ( 0 );
//
// 	slopeA.LoadingRamp(1);
//
//     slopeA.ExecuteSimulation();


    //SetMaterialParamenters ( &slopeA.GetCurrentConfig()->fCMesh,vecmesh,0,1 );

    //int nsamples = vecmesh[0]->Solution().Cols();

    //vecmesh[0]->Solution().Print ( "SOLUTION PHI" );

    // PrintSol ( vecmesh[0] );

    //ShearRed ( slopeA,vecmesh );

//   SolveSimpleElasticProlem();


    //  ComputeAdptiveMontecarlo();


//     // Mohr Coulomb data
//     REAL mc_cohesion    = 10.0;
//     REAL mc_phi         = ( 30.0*M_PI/180 );
//     REAL mc_psi         = mc_phi;
//
//     TPZElasticResponse ER;
//     REAL nu = 0.49;
//     REAL E = 20000.;
//
//     TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
//     ER.SetUp( E, nu );
// 	LEMC.fER =ER;
//     LEMC.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );
//     int PlaneStrain = 1;
//
//     TPZFMatrix<REAL> Dep;
//
//     TPZTensor<REAL> tensoresp,tensorsig;
//   //      REAL exx=-1.2560424036224354E-003 ,exy=-3.0306513266847217E-004,exz=0.,eyy=7.5770078874204600E-004 ,eyz=0.;
//   //  REAL ezz= 1.6664811388913960E-003  ;
//
//      REAL exx=-2.9780912404287953E-003 ,exy=3.0354177217632269E-002,exz=0.,eyy=4.2960017010575877E-003 ,eyz=0.;
//     REAL ezz= -1.9353008778866254E-017 ;
//
//
//     tensoresp.XX() =exx;
//     tensoresp.XY() =exy;
//     tensoresp.XZ() =exz;
//     tensoresp.YY() =eyy;
//     tensoresp.YZ() =eyz;
//     tensoresp.ZZ() =ezz;
//     LEMC.ApplyStrainComputeDep(tensoresp,tensorsig,Dep);
//   //  tensorsig.Print(std::cout);
// /*    Dep.Print("DEP");
//  */
//
//        		cout << "\nDep " <<endl;
// 		cout << Dep(_XX_,_XX_) << "\t" << Dep(_XX_,_YY_) << "\t" << Dep(_XX_,_XY_) <<"\n";
// 		cout << Dep(_YY_,_XX_) << "\t" << Dep(_YY_,_YY_) << "\t" << Dep(_YY_,_XY_) <<"\n";
// 		cout << Dep(_XY_,_XX_) << "\t" << Dep(_XY_,_YY_) << "\t" << Dep(_XY_,_XY_) <<"\n";
//
//  ConfigSlope();
//ComputeAdptiveMontecarlo();
//   ConfigSlopeShearRed();

    ShearRed();
    cout <<"SUCESS"<<endl;
    return 0;

}

TPZElastoPlasticAnalysis  CreateAnal ( TPZCompMesh *fCMesh )
{
    TPZElastoPlasticAnalysis  analysis ( fCMesh,std::cout );

    TPZSkylineStructMatrix full ( fCMesh );

    TPZStepSolver<REAL> step;

    step.SetDirect ( ELDLt );

    analysis.SetStructuralMatrix ( full );

    analysis.SetSolver ( step );

    return analysis;
}

void Refinine ( TPZSlopeStabilityAnalysis &slope,REAL sqj2_refine,int ref,int order )
{
    cout <<" ----  Refining  ---- "<< endl;
    for ( int i=1; i<=ref; i++ ) {
        std::cout << "ref = " << i<< std::endl;
        slope.PRefineElementAbove ( sqj2_refine, order );
        slope.DivideElementsAbove ( sqj2_refine );
        slope.ExecuteSimulation ( );
        // slope.PostProcess ( 0 );
    }
}

void   ConfigSlopeShearRed()
{

    gRefDBase.InitializeUniformRefPattern ( EOned );
    gRefDBase.InitializeUniformRefPattern ( ETriangle );
    gRefDBase.InitializeUniformRefPattern ( EQuadrilateral );
    std::cout << std::setprecision ( 15 );

    TPZSlopeStabilityAnalysis slope;

    std::string output = "ConfigSlopeShearRed.vtk";

    slope.SetVtkOutPutName ( output );

    int ref=2;
    int porder=2;

    slope.GetCurrentConfig()->CreateGeometricMeshSlope2 ( ref );

    slope.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );

    slope.GetCurrentConfig()->CreatePostProcessingMesh();

    slope.PostProcess ( 0 );

    slope.LoadingRamp ( 1 );

    cout <<" ----  Initial Simulation  ---- "<< endl;

    slope.ExecuteInitialSimulation ( 1 );

    slope.PostProcess ( 0 );

    cout <<" ----  Shear Red  ---- "<< endl;


    ShearRed ( );

    //slope.PostProcess ( 0 );

}


TPZSlopeStabilityAnalysis CreateSlopeObj ( REAL porder,REAL bodyforcefac, int ref )
{
    TPZSlopeStabilityAnalysis slope;

    slope.GetCurrentConfig()->CreateGeometricMeshSlope2 ( ref );

    slope.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );

    slope.GetCurrentConfig()->CreatePostProcessingMesh();

    //slopeB.PostProcess ( 0 );

    slope.LoadingRamp ( bodyforcefac );

    return slope;
}


void ShearRed ( )
{


    REAL FS=1,FSmax=10000.,FSmin=0.,tol=1.e-3;
    REAL c;
    REAL phi;
    int ref=3;
    int porder=2;

    int refinit = 2;

    REAL loadfac=1.;

    int counterout=0;
    do {

        REAL norm = 1000.;
        REAL tol2 = 1.e-5;
        REAL sqj2_refine = 0.005;

        TPZSlopeStabilityAnalysis slopeB;



        slopeB = CreateSlopeObj ( porder,loadfac,refinit );

        TPZCompMesh * cmesh =  &slopeB.GetCurrentConfig()->fCMesh;
        plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
        TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();

        TPZElasticResponse ER = LEMC.fER;
        REAL phi0 = LEMC.fYC.Phi();
        REAL cohesion0 = LEMC.fYC.Cohesion();


        c=cohesion0/FS;
        phi=atan ( tan ( phi0 ) /FS );

        LEMC.fYC.SetUp ( phi, phi, c, ER );
        material->SetPlasticity ( LEMC );


        std::cout << " ------------------------------ | FS = " << FS  << "  | phi = " << phi*180/M_PI << " | c = " <<  c << std::endl;

        norm = slopeB.ExecuteInitialSimulation2 ( 1 );


        if ( norm<tol2 ) {


            for ( int i=1; i<=ref; i++ ) {
				
                std::cout << "nref = " << i<< std::endl;
				
                slopeB.PRefineElementAbove ( sqj2_refine, porder+1 );
				
                slopeB.DivideElementsAbove ( sqj2_refine );
				
                norm = slopeB.ExecuteSimulation2();
				
                if ( norm>tol2 ) {
                    break;
                } else if(i==ref){
					
                    std::string output = "aConfigSlopeShearRed-FS-";
					
                    auto iters =to_string ( FS );
					
                    output+=iters;
					
                    string exts = "-.vtk";
					
                    output+=exts;
					
                    slopeB.SetVtkOutPutName ( output );

                    slopeB.PostProcess ( 0 );
                }
            }

        }

        if ( norm>= tol2 ) {

            FSmax = FS;

            FS = ( FSmin + FSmax ) / 2.;

        } else {

            FSmin = FS;

            FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );

		}

        counterout++;

    }  while ( ( FSmax - FSmin ) / FS > tol );

    std::cout << " ------------------------------ | FS = " << FS  << "  | phi = " << phi*180/M_PI << " | c = " <<  c << std::endl;
}

/*
void ShearRed ( TPZSlopeStabilityAnalysis &slopeBb )
{
    TPZSlopeStabilityAnalysis slopeA ( slopeBb );

    REAL FS=1,FSmax=10000.,FSmin=0.,tol=1.e-3;

    TPZCompMesh * cmesh =  &slopeBb.GetCurrentConfig()->fCMesh;

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();

    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi=phi0;
    REAL c=cohesion0;
    int ref=4;
    int porder=2;

    do {


        std::cout << " ------------------------------ | FS = " << FS  << "  | phi = " << phi*180/M_PI << " | c = " <<  c << std::endl;

        REAL norm = 1000.;
        REAL tol2 = 1.e-5;
        REAL sqj2_refine = 0.001;


        TPZSlopeStabilityAnalysis slopeB;


        slopeB.GetCurrentConfig()->CreateGeometricMeshSlope2 ( 1 );

        slopeB.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );

        slopeB.GetCurrentConfig()->CreatePostProcessingMesh();

        //slopeB.PostProcess ( 0 );

        slopeB.LoadingRamp ( 1 );

        TPZCompMesh * cmesh =  &slopeB.GetCurrentConfig()->fCMesh;
        plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
        TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();

        c=cohesion0/FS;
        phi=atan ( tan ( phi0 ) /FS );

        LEMC.fYC.SetUp ( phi, phi, c, ER );
        material->SetPlasticity ( LEMC );


        norm = slopeB.ExecuteInitialSimulation2 ( 1 );

        //norm = slopeB.ExecuteSimulation2();

        // if(norm<tol2)
        //  {




        if ( norm<tol2 ) {
            std::string output = "aConfigSlopeShearRed-FS-";
            auto iters =to_string ( FS );
            output+=iters;
            string exts = "-.vtk";
            output+=exts;
            slopeB.SetVtkOutPutName ( output );
            //slopeB.PostProcess ( 0 );
            for ( int i=1; i<=ref; i++ ) {
                std::cout << "nref = " << i<< std::endl;
                slopeB.PRefineElementAbove ( sqj2_refine, porder+1 );
                slopeB.DivideElementsAbove ( sqj2_refine );
                norm = slopeB.ExecuteSimulation2();
                if ( norm<tol2 ) {
                    //slopeB.PostProcess ( 0 );
                } else {
                    break;
                }

            }

        }


        //displace=slopeB.GetCurrentConfig()->fSolution;

        // deltau=displace;
        //deltau-=displace0;
        // unorm=Norm(deltau);

        if ( norm>= tol2 ) {

            // displace = displace0;

            FSmax = FS;

            FS = ( FSmin + FSmax ) / 2.;


        } else {
            //displace0 = displace;

            FSmin = FS;

            FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );

            //Refinine(slopeB,sqj2_refine, ref, porder);

            //slopeB.PostProcess ( 0 );
            slopeB.PostProcess ( 0 );

        }

//         TPZCompMesh * cmesh =  &slopeB.GetCurrentConfig()->fCMesh;
//         plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
//         TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
//
//         c=cohesion0/FS;
//         phi=atan ( tan ( phi0 ) /FS );
//
//         LEMC.fYC.SetUp ( phi, phi, c, ER );
//         material->SetPlasticity ( LEMC );
        counterout++;

    }  while ( ( FSmax - FSmin ) / FS > tol );

    std::cout << " ------------------------------ | FS = " << FS  << "  | phi = " << phi*180/M_PI << " | c = " <<  c << std::endl;
}*/
// void ShearRed ( TPZSlopeStabilityAnalysis &slopeB )
// {
//     TPZSlopeStabilityAnalysis slopeA(slopeB);
//     REAL FS=1,FSmax=10000.,FSmin=0.,tol=1.e-3;
// 	int neq = slopeB.GetCurrentConfig()->fCMesh.NEquations();
//     TPZCompMesh * cmesh =  &slopeB.GetCurrentConfig()->fCMesh;
//     TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);
// 	//TPZFNMatrix<10000,REAL> mat(neq,1);
// 	// mat(neq,1);
// 	int counterout = 0, maxcountout = 20;
// 	TPZFMatrix<REAL>  deltau;
// 	REAL unorm =1000.;
//
//   // plasticmat pMatWithMem2 =  cmesh->MaterialVec() [1];
//     plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
//     TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
//
//     TPZElasticResponse ER = LEMC.fER;
//     REAL phi0 = LEMC.fYC.Phi();
//     REAL psi0 = phi0;
//     REAL cohesion0 = LEMC.fYC.Cohesion();
//     REAL phi,psi,c;
//     int ref=3;
//     int order=3;
//
//     do {
//
//
//
//         REAL norm = 1000.;
//         REAL tol2 = 1.e-3;
//         REAL sqj2_refine = 0.001;
//
//         norm = slopeB.ExecuteSimulation2();
//
//
//      //   slopeB.ExecuteSimulation();
//         Refinine(slopeB,sqj2_refine, ref, order);
//
// 		displace=slopeB.GetCurrentConfig()->fSolution;
//
// 		deltau=displace;
// 		deltau-=displace0;
// 		unorm=Norm(deltau);
// 		//std::cout << "| Load step = " << counterout << " | Rhs norm = " << norm << "| delta displacement norm = "<< unorm << std::endl;
//
//         if ( norm>= tol2) {
//
//             displace = displace0;
//
//             FSmax = FS;
//
//             FS = ( FSmin + FSmax ) / 2.;
//
//
//         } else {
//                 displace0 = displace;
//
//                 FSmin = FS;
//
//                 FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
//
//                 slopeB.PostProcess ( 0 );
//
//         }
//
//         slopeB=slopeA;
//
//         std::cout << "FS = " << FS << std::endl;
//         c=cohesion0/FS;
//         phi=atan ( tan ( phi0 ) /FS );
//         psi=phi;
//         LEMC.fYC.SetUp ( phi, psi, c, ER );
//         material->SetPlasticity ( LEMC );
//         counterout++;
//
//     }  while ( ( FSmax - FSmin ) / FS > tol );
// }

/*
void ShearRed ( TPZSlopeStabilityAnalysis &slopeB )
{
    REAL FS=1,FSmax=10000.,FSmin=0.,tol=1.e-3;
	int neq = slopeB.GetCurrentConfig()->fCMesh.NEquations();
    TPZCompMesh * cmesh =  &slopeB.GetCurrentConfig()->fCMesh;
    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);
	//TPZFNMatrix<10000,REAL> mat(neq,1);
	// mat(neq,1);
	int counterout = 0, maxcountout = 20;
	TPZFMatrix<REAL>  deltau;
	REAL unorm =1000.;

  // plasticmat pMatWithMem2 =  cmesh->MaterialVec() [1];
    plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();

    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL psi0 = phi0;
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    do {


		//displace0 = slopeB.GetCurrentConfig()->fCMesh.Solution();
		TPZElastoPlasticAnalysis  anal = CreateAnal(&slopeB.GetCurrentConfig()->fCMesh);

		REAL norm = 1000.;
        REAL tol2 = 1.e-3;
        int NumIter = 30;
        bool linesearch = true;
        bool checkconv = false;



        anal.IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
		norm = Norm(anal.Rhs());
		displace=anal.Solution();
		deltau=displace;
		deltau-=displace0;
		unorm=Norm(deltau);
		//std::cout << "| Load step = " << counterout << " | Rhs norm = " << norm << "| delta displacement norm = "<< unorm << std::endl;

        if ( norm>= tol2) {

            displace = displace0;

            FSmax = FS;

            FS = ( FSmin + FSmax ) / 2.;


        } else {
                displace0 = displace;

                FSmin = FS;

                FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );

        }

        std::cout << "FS = " << FS << std::endl;
        c=cohesion0/FS;
        phi=atan ( tan ( phi0 ) /FS );
        psi=phi;
        LEMC.fYC.SetUp ( phi, psi, c, ER );
        material->SetPlasticity ( LEMC );
        counterout++;

    }  while ( ( FSmax - FSmin ) / FS > tol );
}*/

// void ShearRed ( TPZSlopeStabilityAnalysis &slopeB, TPZVec<TPZCompMesh *> vecmesh )
// {
//     REAL FS=1,FSmax=10000.,FSmin=0.,tol=1.e-3;
// 	int neq = slopeB.GetCurrentConfig()->fCMesh.NEquations();
//     TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);
// 	//TPZFNMatrix<10000,REAL> mat(neq,1);
// 	// mat(neq,1);
// 	int counterout = 0, maxcountout = 20;
// 	TPZFMatrix<REAL>  deltau;
// 	REAL unorm =1000.;
//     do {
//
//
// 		//displace0 = slopeB.GetCurrentConfig()->fCMesh.Solution();
// 		TPZElastoPlasticAnalysis  anal = CreateAnal(&slopeB.GetCurrentConfig()->fCMesh);
//
// 		REAL norm = 1000.;
//         REAL tol2 = 1.e-3;
//         int NumIter = 30;
//         bool linesearch = true;
//         bool checkconv = false;
//
//
//
//         anal.IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
// 		norm = Norm(anal.Rhs());
// 		displace=anal.Solution();
// 		deltau=displace;
// 		deltau-=displace0;
// 		unorm=Norm(deltau);
// 		//std::cout << "| Load step = " << counterout << " | Rhs norm = " << norm << "| delta displacement norm = "<< unorm << std::endl;
//
//         if ( norm>= tol2) {
//
//             displace = displace0;
//
//             FSmax = FS;
//
//             FS = ( FSmin + FSmax ) / 2.;
//
//
//         } else {
//                 displace0 = displace;
//
//                 FSmin = FS;
//
//                 FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
//
//         }
//
//         std::cout << "FS = " << FS << std::endl;
//         SetMaterialParamenters ( &slopeB.GetCurrentConfig()->fCMesh,vecmesh,0,FS );
//         counterout++;
//     }  while ( ( FSmax - FSmin ) / FS > tol );
// }

void PrintSol ( TPZCompMesh * plasticCmesh )
{
    REAL mean=0.;
    int count=1;
    int nels = plasticCmesh->NElements();
    for ( int iel=0; iel<nels; iel++ ) {
        TPZCompEl *celplastic = plasticCmesh->Element ( iel );
        if ( celplastic->Material()->Id() !=1 ) continue;
        TPZInterpolationSpace *intelplastic = dynamic_cast<TPZInterpolationSpace *> ( celplastic );
        const TPZIntPoints &intpoints = intelplastic->GetIntegrationRule();
        TPZManVector<REAL,3> point ( 3,0. );
        TPZMaterialData dataplastic;
        dataplastic.fNeedsSol = true;
        intelplastic->InitMaterialData ( dataplastic );
        REAL weight=0;

        int nint = intpoints.NPoints();
        for ( long ip =0; ip<nint; ip++ ) {
            intpoints.Point ( ip, point, weight );
            //cout << " point = " << point << endl;

            intelplastic->ComputeRequiredData ( dataplastic, point );
            cout << "point = " << point <<endl;
            cout << "dataplastic.x = " << dataplastic.x <<endl;
            cout << "dataplastic.phi " << dataplastic.phi  <<endl;
            cout << "dataplastic.sol " << dataplastic.sol  <<endl;
            cout << "dataplastic.detjac " << dataplastic.detjac <<endl;
            REAL sum=0.;
            for ( int iphi=0; iphi<dataplastic.phi.Rows(); iphi++ ) {
                sum+=dataplastic.phi[iphi];
            }
            cout << "sum = "<< sum << endl ;
            //mean+=dataplastic.sol[0][0];
            count++;
        }

    }

    // mean = mean/count;
    //cout << "MEAN = "<< mean << endl ;
}
void ComputeAdptiveMontecarlo()
{
    //create random fields and put the random material parameters as a solution in two meshes. One for cohesion and other for friction angle.
    int porder=1;
    int ref=2;
    TPZSlopeStabilityAnalysis slopeA;
    TPZVec<TPZCompMesh *> vecmesh;
    vecmesh=CreateFields ( slopeA, porder, ref );




    int nsamples = vecmesh[0]->Solution().Cols();

    vecmesh[0]->Solution().Print ( "SOLUTION PHI" );

    PrintSol ( vecmesh[0] );
    //DebugStop();
    for ( int imc=0; imc<nsamples; imc++ ) {
        std::cout << "--------------Monte Carlo Realization Number: " << imc << std::endl;

        std::string output = "MCSlope";
        TPZSlopeStabilityAnalysis slopeB;
        auto iterstring=std::to_string ( imc );
        output+=iterstring;
        string ext=".vtk";
        output+=ext;
        slopeB.SetVtkOutPutName ( output );

        slopeB.GetCurrentConfig()->CreateGeometricMeshSlope2 ( ref );

        slopeB.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );

        //Insert the two meshes material parameters in the elastoplastic mesh.
        REAL fac=1.;
        SetMaterialParamenters ( &slopeB.GetCurrentConfig()->fCMesh,vecmesh,imc,fac );
        std::cout << "Creating Post Process mesh = " << std::endl;
        slopeB.GetCurrentConfig()->CreatePostProcessingMesh();

        std::cout << "Post Processing solution = " << std::endl;
        slopeB.PostProcess ( 0 );

        //DebugStop();
        std::cout << "Start to solve deterministic = " << std::endl;
        REAL sqj2_refine=0.0025;
        int steps =1;
        REAL delta=1.7/steps;
        REAL factor=delta;
        for ( int i=1; i<=steps; i++ ) {
            std::cout << "fac = " << factor<< std::endl;
            slopeB.LoadingRamp ( factor );
            slopeB.ExecuteInitialSimulation ( 1 );
            factor+=delta;
        }

        const int nsubsteps = 5;
        int nref=2;
        for ( int i=1; i<=nref; i++ ) {
            std::cout << "nref = " << i<< std::endl;
            slopeB.PRefineElementAbove ( sqj2_refine, 2 );
            slopeB.DivideElementsAbove ( sqj2_refine );
            slopeB.ExecuteSimulation ();
            //well.PostProcess ( 0 );
        }
        slopeB.PostProcess ( 0 );
        std::cout << "End deterministic = " << std::endl;
        //slopeB.PopConfiguration();
    }



}




void SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZVec<TPZCompMesh*> vecmesh,int isol,REAL factor )
{

    int nels = plasticCmesh->NElements();
    int nels0 = vecmesh[0]->NElements();
    int nels1 = vecmesh[1]->NElements();
    //vecmesh[1]->Solution().Print ( "SOLUTION PHI" );
    if ( nels!=nels0 || nels!=nels1 || nels0!=nels1 ) {
        //   cout << "nels "<< nels << " | nels0 = "<< nels0 <<" | nels1 = "<< nels1 <<endl;
        //DebugStop();
    }

    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( plasticCmesh->MaterialVec() [1] );
    TPZAdmChunkVector<TPZElastoPlasticMem>  &mem = pMatWithMem2->GetMemory();

    // cout << "mem.NElements() "<< mem.NElements()  <<endl;
    if ( pMatWithMem2 ) {
        pMatWithMem2->SetUpdateMem ( true );
    }

    int globpt=0;
    for ( int iel=0; iel<nels0; iel++ ) {

        TPZCompEl *celplastic = plasticCmesh->Element ( iel );
        TPZInterpolationSpace *intelplastic = dynamic_cast<TPZInterpolationSpace *> ( celplastic );

        TPZCompEl *celrandom1 = vecmesh[0]->Element ( iel );
        TPZInterpolationSpace *intelrandom1 = dynamic_cast<TPZInterpolationSpace *> ( celrandom1 );

        TPZCompEl *celrandom2 = vecmesh[1]->Element ( iel );
        TPZInterpolationSpace *intelrandom2 = dynamic_cast<TPZInterpolationSpace *> ( celrandom2 );

        const TPZIntPoints &intpoints = intelplastic->GetIntegrationRule();
        const TPZIntPoints &intpoints1 = intelrandom1->GetIntegrationRule();
        const TPZIntPoints &intpoints2 = intelrandom2->GetIntegrationRule();

        TPZManVector<REAL,3> point ( 3,0. );
        TPZManVector<REAL,3> point1 ( 3,0. );
        TPZManVector<REAL,3> point2 ( 3,0. );

        TPZMaterialData dataplastic,datarandom1,datarandom2;

        intelplastic->InitMaterialData ( dataplastic );
        intelrandom1->InitMaterialData ( datarandom1 );
        intelrandom2->InitMaterialData ( datarandom2 );

        REAL weight=0;
        int nint = intpoints.NPoints();
        //cout << "intpoints.NPoints() "<< nint  <<endl;
        TPZTensor<REAL> epst,epsp;
        for ( long ip =0; ip<nint; ip++ ) {
            intpoints.Point ( ip, point, weight );
            intpoints1.Point ( ip, point1, weight );
            intpoints2.Point ( ip, point2, weight );
            dataplastic.intLocPtIndex = ip;
            intelplastic->ComputeRequiredData ( dataplastic, point );
            intelrandom1->ComputeRequiredData ( datarandom1, point1 );
            intelrandom2->ComputeRequiredData ( datarandom2, point2 );
            mem[globpt].fPlasticState.fmatprop.Resize ( 2 );

            TPZVec<STATE> Sol1 ;
            TPZVec<STATE> Sol2;
            int varid=0;
            celrandom1->Solution ( point1,varid,Sol1 );
            celrandom2->Solution ( point2,varid,Sol2 );
            // cout << "Sol2 = " << datarandom2.sol[isol][0] <<endl;
            // cout << "Sol2B = " << Sol2[0] <<endl;
            //                    hhatcopy[irow][0] = hhat0[irow][0] / FS;
            //          hhatcopy[irow][1] = atan ( tan ( hhat0[irow][1] ) / FS );
            mem[globpt].fPlasticState.fmatprop[0]=Sol1[0]/factor;
            mem[globpt].fPlasticState.fmatprop[1]=atan ( tan ( Sol2[0] ) /factor );
// 			mem[globpt].fPlasticState.fmatprop[0]=10;
//             mem[globpt].fPlasticState.fmatprop[1]=30*M_PI/180.;
            TPZTensor<REAL> epsTotal,sigma;
            //pMatWithMem2->ApplyStrainComputeSigma(epsTotal,sigma);
            globpt++;

        }

        //pMatWithMem2->SetUpdateMem ( false );
    }
    pMatWithMem2->SetUpdateMem ( false );

    plasticCmesh->Solution().Zero();


}

TPZGeoMesh * CreteGeoMesh ( TPZSlopeStabilityAnalysis &slopeA,int porder,int ref )
{
    slopeA.GetCurrentConfig()->CreateGeometricMeshSlope2 ( ref );
    slopeA.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );
    TPZGeoMesh  * gmesh = ( &slopeA.GetCurrentConfig()->fGMesh );
    std::ofstream file ( "sdda.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,file,true );
    return gmesh;
}

TPZCompMesh * CreteCompMesh ( TPZGeoMesh* gmesh,int porder )
{
    int expansionorder=20;
    REAL Lx=20;
    REAL Ly=2;
    REAL Lz=1.;
    int type=3;
    int id=1;
    int dim = gmesh->Dimension();

    //gmesh->ResetReference();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    KLMaterial * klmat = new KLMaterial ( id,Lx,Ly,Lz,dim,type,expansionorder );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    // cmesh->AdjustBoundaryElements();
    //  cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}


void  InsertMat ( TPZCompMesh * cmesh,int porder )
{
    int expansionorder=20;
    REAL Lx=20;
    REAL Ly=2;
    REAL Lz=1.;
    int type=3;
    int id=1;
    int dim = cmesh->Reference()->Dimension();

    KLMaterial * klmat = new KLMaterial ( id,Lx,Ly,Lz,dim,type,expansionorder );

    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
}

TPZVec<TPZCompMesh *> CreateFields ( TPZSlopeStabilityAnalysis &slopeA,int porder,int ref )
{

    TPZGeoMesh  * gmesh =  CreteGeoMesh ( slopeA,porder,ref );
    TPZCompMesh * cmesh =  CreteCompMesh ( gmesh,porder );
    TPZCompMesh * cmesh2 =  CreteCompMesh ( gmesh,porder );

    //TPZCompMesh * cmesh (&slopeA.GetCurrentConfig()->fCMesh);
    //TPZCompMesh * cmesh2(&slopeA.GetCurrentConfig()->fCMesh);
    cmesh->Print ( std::cout );
    InsertMat ( cmesh, porder );
    cmesh->Print ( std::cout );
    InsertMat ( cmesh2, porder );
    int dim = cmesh->Reference()->Dimension();
    KLAnalysis * klanal = new KLAnalysis ( cmesh );
    KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );
    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );
    klanal->Solve();

    string file ="out-eigenfunctions.vtk";
    klanal->Post ( file,dim,0 );

    TPZVec<TPZCompMesh *> vecmesh ( 2 );
    int samples = 10;
    TPZVec<REAL> mean ( 2 );
    TPZVec<REAL> cov ( 2 );
    mean[0]=10.;
    mean[1]=30.*M_PI/180.;
    cov[0]=0.2;
    cov[1]=0.3;
    REAL corsscorrelation = -0.5;
    KLRandomField * klf = new KLRandomField ( cmesh,klanal,mean,cov,samples,corsscorrelation );
    TPZVec<TPZFMatrix<REAL>> hhat;
    cout << "Creating Log-Normal Random Fields "<<endl;
    hhat = klf->CreateLogNormalRandomField();

    hhat[1].Print ( "PHI" );
    cmesh->LoadSolution ( hhat[0] );
    //TPZCompMesh * cmesh2(cmesh);
    cmesh2->LoadSolution ( hhat[1] );
    vecmesh[0]=cmesh;
    vecmesh[1]=cmesh2;

    file ="out-coes.vtk";
    klanal->Post ( file,dim,0 );
    KLAnalysis * klanal2 = new KLAnalysis ( cmesh2 );
    cmesh2->LoadSolution ( hhat[1] );
    file ="out-phi.vtk";
    klanal2->Post ( file,dim,0 );
    return vecmesh;
}




TPZGeoMesh   ConfigSlope()
{

    gRefDBase.InitializeUniformRefPattern ( EOned );
    gRefDBase.InitializeUniformRefPattern ( ETriangle );
    gRefDBase.InitializeUniformRefPattern ( EQuadrilateral );
    std::cout << std::setprecision ( 15 );

    TPZSlopeStabilityAnalysis well;

    std::string output = "ConfigSlope2.vtk";
    well.SetVtkOutPutName ( output );


    //well.LoadingRamp(factor);
    REAL sqj2_refine=0.005;
    int Startfrom=0;
    const int nsubsteps = 5;
    int porder=1;
    std::cout << "\n ------- 0 -------- "<<std::endl;
    well.GetCurrentConfig()->CreateGeometricMeshSlope2 ( 4 );

    well.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );

    well.GetCurrentConfig()->CreatePostProcessingMesh();

    well.PostProcess ( 0 );

    int steps =6;
    REAL delta=1.73/steps;
    REAL factor=delta;
    for ( int i=1; i<=steps; i++ ) {
        std::cout << "fac = " << factor<< std::endl;
        well.LoadingRamp ( factor );
        well.ExecuteInitialSimulation ( 1 );
        factor+=delta;
    }

    int nref=3;
    for ( int i=1; i<=nref; i++ ) {
        std::cout << "nref = " << i<< std::endl;
        well.PRefineElementAbove ( sqj2_refine, porder+i );
        well.DivideElementsAbove ( sqj2_refine );
        well.ExecuteSimulation ( );
        well.PostProcess ( 0 );
    }
    well.PostProcess ( 0 );
//   std::cout << "last = " << std::endl;
//   well.GetCurrentConfig()->fCMesh.Solution().Zero();
//   well.GetCurrentConfig()->fMCPV.ResetPlasticMem();
//   well.LoadingRamp ( 0.5 );
//   well.ExecuteSimulation ( nsubsteps,0.5 );
//   well.PostProcess ( 0 );
//   well.LoadingRamp ( 0.75 );
//   well.ExecuteSimulation ( nsubsteps,0.75 );
//   well.PostProcess ( 0 );
//   well.LoadingRamp ( 1.0 );
//   well.ExecuteSimulation ( nsubsteps,1.0 );
//   well.PostProcess ( 0 );
//   well.LoadingRamp ( 1.25 );
//   well.ExecuteSimulation ( nsubsteps,1.25 );
//   well.PostProcess ( 0 );
//   well.LoadingRamp ( 1.5 );
//   well.ExecuteSimulation ( nsubsteps,1.5 );
//   well.PostProcess ( 0 );
//   well.LoadingRamp ( 1.77 );
//   well.ExecuteSimulation ( nsubsteps,1.77 );
//   well.PostProcess ( 0 );



    std::ofstream file ( "dsaadfadfd.vtk" );
    cout << "well.GetConfigListSize()" << well.GetConfigListSize() << std::endl;
    TPZGeoMesh *gmesh= &well.GetCurrentConfig()->fGMesh;

    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh, file,true );

    TPZCompMesh c=well.GetCurrentConfig()->fCMesh;
    TPZGeoMesh *g = c.Reference();

    return *g;
}


TPZGeoMesh * Simple3DGMesh ( int ref )
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;

//   std::vector<std::vector<double>> co = {	{0.,0.},{1., 0. },{2., 0. },{2., 2. },{1.,2.},{0.,2.}};
//   std::vector<std::vector<int>> topopology = {{0,1,4,5},{1,2,3,4}};

    gmesh->SetDimension ( 3 );
    gmesh->NodeVec().Resize ( 12 );

    TPZVec<REAL> coord ( 3 );
    coord[0] = -8.;
    coord[1] =  0.;
    coord[2] =  0.;
    gmesh->NodeVec() [0] = TPZGeoNode ( 0, coord, *gmesh );

    coord[0] = 0.;
    coord[1] = 8.;
    coord[2] = 0.;
    gmesh->NodeVec() [1] = TPZGeoNode ( 1, coord, *gmesh );

    coord[0] = 8.;
    coord[1] = 0.;
    coord[2] = 0.;
    gmesh->NodeVec() [2] = TPZGeoNode ( 2, coord, *gmesh );

    coord[0] = 10.;
    coord[1] = 0.;
    coord[2] = 0.;
    gmesh->NodeVec() [3] = TPZGeoNode ( 3, coord, *gmesh );

    coord[0] = 0.;
    coord[1] = 10.;
    coord[2] = 0.;
    gmesh->NodeVec() [4] = TPZGeoNode ( 4, coord, *gmesh );

    coord[0] = -10.;
    coord[1] =  0.;
    coord[2] =  0.;
    gmesh->NodeVec() [5] = TPZGeoNode ( 5, coord, *gmesh );

    ////

    coord[0] = -8.;
    coord[1] =  0.;
    coord[2] =  15.;
    gmesh->NodeVec() [6] = TPZGeoNode ( 6, coord, *gmesh );

    coord[0] = 0.;
    coord[1] = 8.;
    coord[2] = 15.;
    gmesh->NodeVec() [7] = TPZGeoNode ( 7, coord, *gmesh );

    coord[0] = 8.;
    coord[1] = 0.;
    coord[2] = 15.;
    gmesh->NodeVec() [8] = TPZGeoNode ( 8, coord, *gmesh );

    coord[0] = 10.;
    coord[1] = 0.;
    coord[2] = 15.;
    gmesh->NodeVec() [9] = TPZGeoNode ( 9, coord, *gmesh );

    coord[0] = 0.;
    coord[1] = 10.;
    coord[2] = 15.;
    gmesh->NodeVec() [10] = TPZGeoNode ( 10, coord, *gmesh );

    coord[0] = -10.;
    coord[1] =  0.;
    coord[2] =  15.;
    gmesh->NodeVec() [11] = TPZGeoNode ( 11, coord, *gmesh );

    TPZVec<long> TopolQuad ( 4 );
    TopolQuad[0] = 0;
    TopolQuad[1] = 2;
    TopolQuad[2] = 3;
    TopolQuad[3] = 5;

    new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > ( 0,TopolQuad,-1,*gmesh );


    TopolQuad[0] = 6;
    TopolQuad[1] = 8;
    TopolQuad[2] = 9;
    TopolQuad[3] = 11;

    new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > ( 1,TopolQuad,-1,*gmesh );

///laterais
    TopolQuad[0] = 2;
    TopolQuad[1] = 3;
    TopolQuad[2] = 9;
    TopolQuad[3] = 8;

    new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > ( 2,TopolQuad,-1,*gmesh );


    TopolQuad[0] = 0;
    TopolQuad[1] = 5;
    TopolQuad[2] = 11;
    TopolQuad[3] = 6;

    new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > ( 3,TopolQuad,-1,*gmesh );


    TopolQuad[0] = 0;
    TopolQuad[1] = 2;
    TopolQuad[2] = 8;
    TopolQuad[3] = 6;

    new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > ( 4,TopolQuad,-1,*gmesh );


    TopolQuad[0] = 5;
    TopolQuad[1] = 3;
    TopolQuad[2] = 9;
    TopolQuad[3] = 11;

    new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > ( 5,TopolQuad,-1,*gmesh );


    TPZVec<long> TopolArc ( 3 );
    TopolArc[0] = 0;
    TopolArc[1] = 2;
    TopolArc[2] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D> ( 6,TopolArc,-1,*gmesh );

    TopolArc[0] = 5;
    TopolArc[1] = 3;
    TopolArc[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 7,TopolArc,-1,*gmesh );

    TopolArc[0] = 6;
    TopolArc[1] = 8;
    TopolArc[2] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 8,TopolArc,-1,*gmesh );

    TopolArc[0] = 11;
    TopolArc[1] = 9;
    TopolArc[2] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 9,TopolArc,-1,*gmesh );

    TPZVec<long> TopolCube ( 8 );
    TopolCube[0] = 0;
    TopolCube[1] = 2;
    TopolCube[2] = 3;
    TopolCube[3] = 5;

    TopolCube[4] = 6;
    TopolCube[5] = 8;
    TopolCube[6] = 9;
    TopolCube[7] = 11;

    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > ( 10,TopolCube,1,*gmesh );
    gmesh->BuildConnectivity();
// 			//------------------ fazer refinamento direcional ------------------------------------
// 		int  ndirectdivR =1;
// 		set<int> SETmatRefDir11;
// 		for(int j = 0; j < ndirectdivR; j++)
// 		{
// 			long nel = gmesh->NElements();
// 			for (long iref = 0; iref < nel; iref++)
// 			{
// 				TPZVec<TPZGeoEl*> filhos;
// 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
// 				if(!gelP11) continue;
// 				SETmatRefDir11.insert(1);
// 				int matid= gelP11->MaterialId();
//                 TPZRefPatternTools::RefineDirectional ( gelP11, SETmatRefDir11 );
//
// 			}
// 		}

    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    std::ofstream file ( "cubesh222.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh, file,true );
    return gmesh;
}

TPZGeoMesh * SimpleGMesh ( int ref )
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;

    gmesh->SetDimension ( 2 );
//   std::vector<std::vector<double>> co = {	{0.,0.},{1., 0. },{2., 0. },{2., 2. },{1.,2.},{0.,2.}};
//   std::vector<std::vector<int>> topopology = {{0,1,4,5},{1,2,3,4}};

    gmesh->NodeVec().Resize ( 4 );

    TPZVec<REAL> coord ( 2 );
    coord[0] = 0.;
    coord[1] = 0.;
    gmesh->NodeVec() [0] = TPZGeoNode ( 0, coord, *gmesh );


    coord[0] = 1.;
    coord[1] = 0.;
    gmesh->NodeVec() [1] = TPZGeoNode ( 1, coord, *gmesh );

    coord[0] = 1.;
    coord[1] = 1.;
    gmesh->NodeVec() [2] = TPZGeoNode ( 2, coord, *gmesh );

    coord[0] = 0.;
    coord[1] = 1.;
    gmesh->NodeVec() [3] = TPZGeoNode ( 3, coord, *gmesh );


    TPZVec <long> TopoQuad ( 4 );
    TopoQuad[0] = 0;
    TopoQuad[1] = 1;
    TopoQuad[2] = 2;
    TopoQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 0, TopoQuad, 1,*gmesh );

    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 1, TopoLine, -1,*gmesh );//Bottom

    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 2, TopoLine, -2,*gmesh );//Rigth

    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 3, TopoLine, -3,*gmesh );//Top

    TopoLine[0] = 3;
    TopoLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 4, TopoLine, -4,*gmesh );//Left

    TPZVec <long> TopoNode ( 1 );
    TopoLine[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 5, TopoLine, -5,*gmesh );//Bottom

    gmesh->BuildConnectivity();
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }
    return gmesh;
}

TPZCompMesh * SimpleCMesh ( TPZGeoMesh * gmesh, int porder )
{
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    int id=1;
    REAL E=20000;
    REAL nu=0.3;
    REAL fx=0.;
    REAL fy=0.;
    int plainstress=0.;
    TPZElasticityMaterial *elastmat = new TPZElasticityMaterial ( id,  E,  nu,  fx,  fy,  plainstress );
    elastmat->SetId ( 1 );

    REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
    elastmat->SetPreStress ( Sigmaxx,Sigmayx,Sigmayy,Sigmazz );

    //BCS
    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );

    int dirichet =0;
    int directionaldirichet =3;
    int newman =1;
    val2 ( 0,0 ) = 1.;
    val2 ( 1,0 ) = 0.;
    auto * bc_left = elastmat->CreateBC ( elastmat,-4,directionaldirichet, val1, val2 );


    val2 ( 0,0 ) = 0.;
    val2 ( 1,0 ) = 1.;
    auto * bc_node = elastmat->CreateBC ( elastmat,-5,directionaldirichet, val1, val2 );

    val2 ( 0,0 ) = 1000.;
    val2 ( 1,0 ) = 0.;
    auto * loadc_top = elastmat->CreateBC ( elastmat,-2,newman, val1, val2 );

    cmesh->InsertMaterialObject ( elastmat );
    cmesh->InsertMaterialObject ( bc_left );
    cmesh->InsertMaterialObject ( bc_node );
    cmesh->InsertMaterialObject ( loadc_top );


    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( 2 );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    return cmesh;
}
void SolveSimpleElasticProlem( )
{
    int uniformref=0;
    int porder=2;

    //Create Geometry
    TPZGeoMesh * gmesh =  SimpleGMesh ( uniformref );

    //Create Approx Space
    TPZCompMesh * cmesh = SimpleCMesh ( gmesh,porder );

    //Create Analysis
    bool optimizeband =false;
    TPZAnalysis * anal = new TPZAnalysis ( cmesh,optimizeband );

    //Set Matrix Type
    TPZFStructMatrix full ( cmesh );
    full.SetNumThreads ( 0 );
    anal->SetStructuralMatrix ( full );

    //Set Solver Type
    TPZStepSolver<REAL> step;
    step.SetDirect ( ECholesky );
    anal->SetSolver ( step );

    anal->Run();

    auto exactSol = [] ( const TPZVec<REAL> &loc,
                         TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU ) {
        const auto &x=loc[0];
        const auto &y=loc[1];
        u[0]=x*4.55e-5;
        u[1]=-y*1.95e-5;
        gradU ( 0,0 ) = 4.55e-5;
        gradU ( 1,0 ) = -1.95e-5;
        gradU ( 2,0 ) = 0; //optional
    };

    // Post processing
    TPZManVector<std::string> scalarnames ( 3 ), vecnames ( 2 );
    scalarnames[0] = "SigmaX";
    scalarnames[1] = "SigmaY";
    scalarnames[1] = "SigmaZ";
    scalarnames[2] = "Pressure";
    vecnames[0] = "displacement";
    vecnames[1] = "Strain";
    anal->DefineGraphMesh ( 2,scalarnames,vecnames,"SimpleElasticProlem.vtk" );

    anal->PostProcess ( 0 );

    //let us set the exact solution and suggest an integration rule
    //for calculating the error
    anal->SetExact ( exactSol );

    ///Calculating approximation error
    TPZManVector<REAL,3> error;
    std::ofstream anPostProcessFile ( "postprocess.txt" );
    anal->PostProcess ( error,anPostProcessFile );

    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << error[0]<<'\n';
    std::cout << "L1 Norm = " << error[1]<<'\n';
    std::cout << "H1 Seminorm = " << error[2] << "\n\n";

    cmesh->LoadSolution ( anal->Solution() );

    cmesh->Solution().Print ( std::cout );


    PrintSol ( cmesh );

}
