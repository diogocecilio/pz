
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
using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;
#ifdef LOG4CXX
static LoggerPtr logger ( Logger::getLogger ( "plasticity.main" ) );
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse ( Logger::getLogger ( "LogEllipse" ) );
#endif

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb;
// If you have issues with: dyld: Library not loaded: libtbb.dylib
// try setting the LD path. Ex:
//   export DYLD_FALLBACK_LIBRARY_PATH=/Users/borin/Desktop/neopz/tbb40_297oss/lib/
#endif

TPZGeoMesh * MisesPressure(int h, int dir);
RunStatsTable plast_tot ( "-tpz_plast_tot", "Raw data table statistics for the main execution." );
clarg::argInt NumberOfThreads ( "-nt", "Number of threads for WellBoreAnalysis", 8 );
void ScaleVec(TPZVec<REAL> &NodeIni, TPZVec<REAL> &NodeFin, double Norm, TPZVec<REAL> &OriginVec, TPZVec<REAL> &OutputVec);
TPZGeoMesh  ConfigSlope();

TPZCompMesh* CreateKLCMesh ( TPZGeoMesh * gmesh,int porder,REAL Lx, REAL Ly,REAL Lz, int type, int expansionorder );

void KLConfig();
TPZVec<TPZFMatrix<REAL>> CreateFields ( TPZCompMesh * cmesh,KLAnalysis * klanal );
TPZGeoMesh * SimpleGMesh ( int ref );
TPZGeoMesh * Simple3DGMesh ( int ref );
int main ( int argc, char **argv )
{

	
	//TPZGeoMesh * cubegmesh =  Simple3DGMesh (4);
	//TPZGeoMesh * cubegmesh =   MisesPressure(2, 1);


   // std::ofstream file ( "cube.vtk" );
	
   // TPZVTKGeoMesh::PrintGMeshVTK ( cubegmesh, file,true );
	//return 0;
    KLConfig();
    std::cout << " SUCESS " << std::endl;
    //ConfigSlope();
    return 0;

//     TPZGeoMesh  gmesh = ConfigSlope();
//   gmesh.Print();
//
//   TPZCompMesh * cmesh = CreateKLCMesh(&gmesh);
//   KLConfig(cmesh);
//   //ConfigSlope();
//   return 0;

}
void KLConfig()
{

    int porder=2;
    int expansionorder=100;

    TPZSlopeStabilityAnalysis slopeA,slopeB;

//     slopeA.GetCurrentConfig()->CreateGeometricMeshSlope(2);
//     slopeA.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );
//     TPZGeoMesh * gmesh = &slopeA.GetCurrentConfig()->fGMesh;

    int ref=3;
   // TPZGeoMesh * gmesh = SimpleGMesh ( ref );
	
	TPZGeoMesh * gmesh =  Simple3DGMesh (ref);
	//gmesh->Print(std::cout);
// TPZGeoMesh  gmesh = ConfigSlope();
// gmesh.Print();


    REAL Lx=10;
    REAL Ly=0.5;
    REAL Lz=0.;
    int type=4;
//     REAL Lx=1;
//     REAL Ly=1;
//     REAL Lz=1.;
//     int type=3;

    TPZCompMesh * cmesh = CreateKLCMesh ( gmesh,porder,Lx,Ly,Lz,type, expansionorder );
    KLAnalysis * klanal = new KLAnalysis ( cmesh );
    KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );
    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );
    klanal->Solve();


	int dim = gmesh->Dimension();
	ref=4;
	string file ="oute.vtk";
	klanal->Post(file,dim,ref);
//     // Post processing
//     TPZManVector<std::string> scalarnames ( 6 ), vecnames ( 0 );
//     scalarnames[0] = "vec";
//     scalarnames[1] = "vec1";
//     scalarnames[2] = "vec2";
//     scalarnames[3] = "vec3";
//     scalarnames[4] = "vec4";
//     scalarnames[5] = "vec5";
// 
// 	KLMaterial * mat2 = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );
// 	
// 	mat2->Print(std::cout);
//     TPZAnalysis * anal = new TPZAnalysis ( cmesh );
// 	anal->SetCompMesh(cmesh,false);
// 
//     klanal->DefineGraphMesh ( 3,scalarnames,vecnames,"klsol3D.vtk" );
// 
//     klanal->PostProcess ( 0 );

    CreateFields ( cmesh,klanal );

}


TPZGeoMesh   ConfigSlope()
{

    gRefDBase.InitializeUniformRefPattern ( EOned );
    gRefDBase.InitializeUniformRefPattern ( ETriangle );
    gRefDBase.InitializeUniformRefPattern ( EQuadrilateral );
    std::cout << std::setprecision ( 15 );

    TPZSlopeStabilityAnalysis well;

    std::string output = "ConfigSlope.vtk";
    well.SetVtkOutPutName ( output );


    //well.LoadingRamp(factor);
    REAL sqj2_refine=0.005;
    int Startfrom=0;
    const int nsubsteps = 5;
    int porder=2;
    std::cout << "\n ------- 0 -------- "<<std::endl;
    well.GetCurrentConfig()->CreateGeometricMeshSlope ( 2 );

    well.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );

    well.GetCurrentConfig()->CreatePostProcessingMesh();

    well.PostProcess ( 0 );

    int steps =2;
    REAL delta=1.77/steps;
    REAL factor=delta;
    for ( int i=1; i<=steps; i++ ) {
        std::cout << "fac = " << factor<< std::endl;
        well.LoadingRamp ( factor );
        well.ExecuteInitialSimulation ( 5 );
        factor+=delta;
    }

    int nref=2;
    for ( int i=1; i<=nref; i++ ) {
        std::cout << "nref = " << i<< std::endl;
        well.PRefineElementAbove ( sqj2_refine, 3 );
        well.DivideElementsAbove ( sqj2_refine );
        well.ExecuteSimulation ( nsubsteps,factor );
        //well.PostProcess ( 0 );
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

TPZGeoMesh * SimpleGMesh ( int ref )
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;

	gmesh->SetDimension(2);
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
TPZGeoMesh * Simple3DGMesh ( int ref )
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;

//   std::vector<std::vector<double>> co = {	{0.,0.},{1., 0. },{2., 0. },{2., 2. },{1.,2.},{0.,2.}};
//   std::vector<std::vector<int>> topopology = {{0,1,4,5},{1,2,3,4}};

	gmesh->SetDimension(3);
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
	
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (0,TopolQuad,-1,*gmesh);
	
	
	TopolQuad[0] = 6;
	TopolQuad[1] = 8;
	TopolQuad[2] = 9;
	TopolQuad[3] = 11;
	
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (1,TopolQuad,-1,*gmesh);
	
///laterais
	TopolQuad[0] = 2;
	TopolQuad[1] = 3;
	TopolQuad[2] = 9;
	TopolQuad[3] = 8;
	
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (2,TopolQuad,-1,*gmesh);
	
	
	TopolQuad[0] = 0;
	TopolQuad[1] = 5;
	TopolQuad[2] = 11;
	TopolQuad[3] = 6;
	
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (3,TopolQuad,-1,*gmesh);
	
	
	TopolQuad[0] = 0;
	TopolQuad[1] = 2;
	TopolQuad[2] = 8;
	TopolQuad[3] = 6;
	
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (4,TopolQuad,-1,*gmesh);
	
	
	TopolQuad[0] = 5;
	TopolQuad[1] = 3;
	TopolQuad[2] = 9;
	TopolQuad[3] = 11;
	
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (5,TopolQuad,-1,*gmesh);
	
	

	
    TPZVec<long> TopolArc ( 3 );
    TopolArc[0] = 0;
    TopolArc[1] = 2;
    TopolArc[2] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D>  ( 6,TopolArc,-1,*gmesh );
	
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

    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (10,TopolCube,1,*gmesh );
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
    
    std::ofstream file("cubesh222.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh, file,true );
	return gmesh;
}
// TPZGeoMesh * Simple3DGMesh ( int ref )
// {
//     TPZGeoMesh *gmesh = new TPZGeoMesh;
// 
// //   std::vector<std::vector<double>> co = {	{0.,0.},{1., 0. },{2., 0. },{2., 2. },{1.,2.},{0.,2.}};
// //   std::vector<std::vector<int>> topopology = {{0,1,4,5},{1,2,3,4}};
// 
// 	gmesh->SetDimension(3);
//     gmesh->NodeVec().Resize ( 12 );
// 
//     TPZVec<REAL> coord ( 3 );
//     coord[0] = 0.;
//     coord[1] = 8.;
//     coord[2] = 0.;
//     gmesh->NodeVec() [0] = TPZGeoNode ( 0, coord, *gmesh );
// 
//     coord[0] = sqrt ( 2. ) /2.*8.;
//     coord[1] = sqrt ( 2. ) /2.*8.;
//     coord[2] = 0.;
//     gmesh->NodeVec() [1] = TPZGeoNode ( 1, coord, *gmesh );
// 
//     coord[0] = 8.;
//     coord[1] = 0.;
//     coord[2] = 0.;
//     gmesh->NodeVec() [2] = TPZGeoNode ( 2, coord, *gmesh );
// 
//     coord[0] = 10.;
//     coord[1] = 0.;
//     coord[2] = 0.;
//     gmesh->NodeVec() [3] = TPZGeoNode ( 3, coord, *gmesh );
// 
//     coord[0] = sqrt ( 2. ) /2.*10.;
//     coord[1] = sqrt ( 2. ) /2.*10.;
//     coord[2] = 0.;
//     gmesh->NodeVec() [4] = TPZGeoNode ( 4, coord, *gmesh );
// 
//     coord[0] = 0.;
//     coord[1] = 10.;
//     coord[2] = 0.;
//     gmesh->NodeVec() [5] = TPZGeoNode ( 5, coord, *gmesh );
// 
// 
// 
//     ////
// 
//     coord[0] = 0.;
//     coord[1] = 8.;
//     coord[2] = 15.;
//     gmesh->NodeVec() [6] = TPZGeoNode ( 6, coord, *gmesh );
// 
//     coord[0] = sqrt ( 2. ) /2.*8.;
//     coord[1] = sqrt ( 2. ) /2.*8.;
//     coord[2] = 15.;
//     gmesh->NodeVec() [7] = TPZGeoNode ( 7, coord, *gmesh );
// 
//     coord[0] = 8.;
//     coord[1] = 0.;
//     coord[2] = 15.;
//     gmesh->NodeVec() [8] = TPZGeoNode ( 8, coord, *gmesh );
// 
//     coord[0] = 10.;
//     coord[1] = 0.;
//     coord[2] = 15.;
//     gmesh->NodeVec() [9] = TPZGeoNode ( 9, coord, *gmesh );
// 
//     coord[0] = sqrt ( 2. ) /2.*10.;
//     coord[1] = sqrt ( 2. ) /2.*10.;
//     coord[2] = 15.;
//     gmesh->NodeVec() [10] = TPZGeoNode ( 10, coord, *gmesh );
// 
//     coord[0] = 0.;
//     coord[1] = 10.;
//     coord[2] = 15.;
//     gmesh->NodeVec() [11] = TPZGeoNode ( 11, coord, *gmesh );
// 	
// 	TPZVec<long> TopolQuad ( 4 );
// 	TopolQuad[0] = 0;
// 	TopolQuad[1] = 2;
// 	TopolQuad[2] = 3;
// 	TopolQuad[3] = 5;
// 	
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (0,TopolQuad,-1,*gmesh);
// 	
// 	
// 	TopolQuad[0] = 6;
// 	TopolQuad[1] = 8;
// 	TopolQuad[2] = 9;
// 	TopolQuad[3] = 11;
// 	
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (1,TopolQuad,-1,*gmesh);
// 	
// ///laterais
// 	TopolQuad[0] = 2;
// 	TopolQuad[1] = 3;
// 	TopolQuad[2] = 9;
// 	TopolQuad[3] = 8;
// 	
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (2,TopolQuad,-1,*gmesh);
// 	
// 	
// 	TopolQuad[0] = 0;
// 	TopolQuad[1] = 5;
// 	TopolQuad[2] = 11;
// 	TopolQuad[3] = 6;
// 	
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (3,TopolQuad,-1,*gmesh);
// 	
// 	
// 	TopolQuad[0] = 0;
// 	TopolQuad[1] = 2;
// 	TopolQuad[2] = 8;
// 	TopolQuad[3] = 6;
// 	
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (4,TopolQuad,-1,*gmesh);
// 	
// 	
// 	TopolQuad[0] = 5;
// 	TopolQuad[1] = 3;
// 	TopolQuad[2] = 9;
// 	TopolQuad[3] = 11;
// 	
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (5,TopolQuad,-1,*gmesh);
// 	
// 	
// 
// 	
//     TPZVec<long> TopolArc ( 3 );
//     TopolArc[0] = 0;
//     TopolArc[1] = 2;
//     TopolArc[2] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZArc3D>  ( 6,TopolArc,-1,*gmesh );
// 	
//     TopolArc[0] = 5;
//     TopolArc[1] = 3;
//     TopolArc[2] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 7,TopolArc,-1,*gmesh );
// 	
//     TopolArc[0] = 6;
//     TopolArc[1] = 8;
//     TopolArc[2] = 7;
//     new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 8,TopolArc,-1,*gmesh );
// 	
//     TopolArc[0] = 11;
//     TopolArc[1] = 9;
//     TopolArc[2] = 10;
//     new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 9,TopolArc,-1,*gmesh );
// 	
//     TPZVec<long> TopolCube ( 8 );
//     TopolCube[0] = 0;
//     TopolCube[1] = 2;
//     TopolCube[2] = 3;
//     TopolCube[3] = 5;
// 
//     TopolCube[4] = 6;
//     TopolCube[5] = 8;
//     TopolCube[6] = 9;
//     TopolCube[7] = 11;
// 
//     new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (10,TopolCube,1,*gmesh );
// //     TPZVec<long> TopolCube ( 8 );
// //     TopolCube[0] = 1;
// //     TopolCube[1] = 2;
// //     TopolCube[2] = 3;
// //     TopolCube[3] = 4;
// // 
// //     TopolCube[4] = 7;
// //     TopolCube[5] = 8;
// //     TopolCube[6] = 9;
// //     TopolCube[7] = 10;
// // 
// //     new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (10,TopolCube,1,*gmesh );
// 
// 	
// //     TopolCube[0] = 0;
// //     TopolCube[1] = 1;
// //     TopolCube[2] = 4;
// //     TopolCube[3] = 5;
// // 
// //     TopolCube[4] = 6;
// //     TopolCube[5] = 7;
// //     TopolCube[6] = 10;
// //     TopolCube[7] = 11;
// // 
// //     new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > ( 11,TopolCube,1,*gmesh );
// 	
// 	
// // 	TPZVec <long> TopolLine(2);
// //  	TopolLine[0] = 0;	TopolLine[1] = 2;
// //  	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (4,TopolLine,-1,*gmesh);
// // 
// //     TPZVec<long> TopolArc ( 3 );
// //     TopolArc[0] = 0;
// //     TopolArc[1] = 2;
// //     TopolArc[2] = 1;
// //     new TPZGeoElRefPattern< pzgeom::TPZArc3D>  ( 5,TopolArc,-1,*gmesh );
// // 	
// // 
// // 	
// // 	TopolLine[0] = 5;	TopolLine[1] = 3;
// //  	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (6,TopolLine,-2,*gmesh);
// // 	
// //     TopolArc[0] = 5;
// //     TopolArc[1] = 3;
// //     TopolArc[2] = 4;
// //     new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 7,TopolArc,-2,*gmesh );
// // 	
// // 	
// // 	TopolLine[0] = 6;	TopolLine[1] = 8;
// //  	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (8,TopolLine,-3,*gmesh);
// // 	
// // 	
// //     TopolArc[0] = 6;
// //     TopolArc[1] = 8;
// //     TopolArc[2] = 7;
// //     new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 9,TopolArc,-3,*gmesh );
// // 	
// // 	
// // 	
// // 	TopolLine[0] = 11;	TopolLine[1] = 9;
// //  	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (10,TopolLine,-4,*gmesh);
// // 	
// //     TopolArc[0] = 11;
// //     TopolArc[1] = 9;
// //     TopolArc[2] = 10;
// //     new TPZGeoElRefPattern< pzgeom::TPZArc3D > ( 11,TopolArc,-4,*gmesh );
// 
// 	    gmesh->BuildConnectivity();
//     for ( int d = 0; d<ref; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//     return gmesh;
//     std::ofstream file ( "cubemesh.vtk" );
// 	
//     TPZVTKGeoMesh::PrintGMeshVTK ( gMesh, file,true );
// }




TPZVec<TPZFMatrix<REAL>> CreateFields ( TPZCompMesh * cmesh,KLAnalysis * klanal )
{
    int samples = 100;
    TPZVec<REAL> mean ( 2 );
    TPZVec<REAL> cov ( 2 );
    mean[0]=10.;
    mean[1]=30.*M_PI/180.;
    cov[0]=0.2;
    cov[1]=0.3;
    REAL corsscorrelation = -0.5;
    KLRandomField * klf = new KLRandomField ( cmesh,klanal,mean,cov,samples,corsscorrelation );
    TPZVec<TPZFMatrix<REAL>> hhat;
    hhat = klf->CreateLogNormalRandomField();

    cmesh->LoadSolution ( hhat[0] );
    // Post processing
    TPZManVector<std::string> scalarnames ( 6 ), vecnames ( 0 );
    scalarnames[0] = "vec";
    scalarnames[1] = "vec1";
    scalarnames[2] = "vec2";
    scalarnames[3] = "vec3";
    scalarnames[4] = "vec4";
    scalarnames[5] = "vec5";

    TPZAnalysis * anal = new TPZAnalysis ( cmesh );
	
	int dim = cmesh->Reference()->Dimension();

    anal->DefineGraphMesh ( dim,scalarnames,vecnames,"coes.vtk" );

    anal->PostProcess ( 0 );

    return hhat;
}
TPZCompMesh* CreateKLCMesh ( TPZGeoMesh * gmesh,int porder,REAL Lx, REAL Ly,REAL Lz, int type, int expansionorder )
{

    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    int id=1;
	int dim = gmesh->Dimension();
    KLMaterial * klmat = new KLMaterial ( id,Lx,Ly,Lz,dim,type,expansionorder );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
	
	cmesh->Print(std::cout);


    return cmesh;
}
TPZGeoMesh * MisesPressure(int h, int dir)
{
	int Qnodes = 6;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	//gRefDBase.InitializeRefPatterns();
	//	gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <long> TopolArc(3);
	TPZVec <long> TopolLine(2);
	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolPoint(1);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , 100.);//coord X
	Node[0].SetCoord(1 , 0.);//coord Y
	Node[0].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,200.);//coord X
	Node[1].SetCoord(1 ,0.);//coord Y
	Node[1].SetCoord(2, 0.);//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  100*sqrt(2.));//coord X
	Node[2].SetCoord(1 ,  100*sqrt(2.));//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  0.);//coord X
	Node[3].SetCoord(1 ,  200.);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	
	///////TRES////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  0.);//coord X
	Node[4].SetCoord(1 ,  100.);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////TRES////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  50*sqrt(2.));//coord X
	Node[5].SetCoord(1 ,  50*sqrt(2.));//coord Y
	Node[5].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	//ELEMENTO 0///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 3; 	TopolQuad[3] = 4;
	
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (0,TopolQuad,1,*gMesh);

	
// 	//ELEMENTO 2///
// 	TopolLine[0] = 0;	TopolLine[1] = 1;
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (1,TopolLine,-1,*gMesh);
// 	
// 	//ELEMENTO 2///
// 	TopolLine[0] = 3;	TopolLine[1] = 4;
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (2,TopolLine,-2,*gMesh);
// 	
// 	//ELEMENTO 2///
// 	TopolLine[0] = 4;	TopolLine[1] = 0;
// 	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (3,TopolLine,-3,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 4;	TopolArc[1] =0;		TopolArc[2] =5;
	new TPZGeoElRefPattern<  TPZArc3D  > (4,TopolArc,-4,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 1;	TopolArc[1] =3;		TopolArc[2] =2;
	new TPZGeoElRefPattern< TPZArc3D  > (5,TopolArc,-5,*gMesh);

	gMesh->BuildConnectivity();
	
	for(int ref = 0; ref < h; ref++)
	{
		TPZVec<TPZGeoEl *> tatara;
		int n = gMesh->NElements();
		for(int i = 0; i < n; i++)
		{
			TPZGeoEl * gel = gMesh->ElementVec()[i];
			//if(gel->HasSubElement()) continue;
			gel->Divide(tatara);
		}
	}

	
	return gMesh;
	
}
