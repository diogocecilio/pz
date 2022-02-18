
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

#include "pzgengrid.h"

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


RunStatsTable plast_tot ( "-tpz_plast_tot", "Raw data table statistics for the main execution." );
clarg::argInt NumberOfThreads ( "-nt", "Number of threads for WellBoreAnalysis", 8 );

TPZGeoMesh  ConfigSlope();
#include "tpzgeoelrefpattern.h"
TPZGeoMesh * SimpleGMesh2()
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;

  std::vector<std::vector<double>> co = {	{1.,0.},{3., 0. },{3., 4. },{1., 4. }};
  std::vector<std::vector<int>> topopology = {{0,1,2,3}};


  gmesh->NodeVec().Resize ( co.size() );
  int id=0;
  for ( int el=0; el<topopology.size(); el++ )
    {
      for ( int elnode =0; elnode<co.size(); elnode++ )
        {
          TPZVec<REAL> coord ( 2 );
          coord[0] = co[topopology[el][elnode]][0];
          coord[1] = co[topopology[el][elnode]][1];
          gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );
          id++;
        }
    }


  TPZVec <long> TopoQuad ( 4 );
  TopoQuad[0] = 0;
  TopoQuad[1] = 1;
  TopoQuad[2] = 2;
  TopoQuad[3] = 3;
  new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 0, TopoQuad, 1,*gmesh );

  gmesh->BuildConnectivity();
gmesh->Print();
  return gmesh;
}
TPZGeoMesh * SimpleGMesh()
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;

//   std::vector<std::vector<double>> co = {	{0.,0.},{1., 0. },{2., 0. },{2., 2. },{1.,2.},{0.,2.}};
//   std::vector<std::vector<int>> topopology = {{0,1,4,5},{1,2,3,4}};

  gmesh->NodeVec().Resize ( 6 );

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
  
//     coord[0] = 1.;
//   coord[1] =2.;
//   gmesh->NodeVec() [4] = TPZGeoNode ( 4, coord, *gmesh );
//   
//     coord[0] = 0.;
//   coord[1] = 2.;
//   gmesh->NodeVec() [5] = TPZGeoNode ( 5, coord, *gmesh );
  

  TPZVec <long> TopoQuad ( 4 );
  TopoQuad[0] = 0;
  TopoQuad[1] = 1;
  TopoQuad[2] = 2;
  TopoQuad[3] = 3;
  new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 0, TopoQuad, 1,*gmesh );
  
//   TopoQuad[0] = 1;
//   TopoQuad[1] = 2;
//   TopoQuad[2] = 3;
//   TopoQuad[3] = 4;
//   new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 1, TopoQuad, 1,*gmesh );
  

  gmesh->BuildConnectivity();
	    for ( int d = 0; d<1; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
                gel->Divide ( subels );
        }
    }
  return gmesh;
}
TPZCompMesh* CreateKLCMesh(TPZGeoMesh * gmesh);


void KLConfig(TPZCompMesh * cmesh);
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzsolve.h"
int main ( int argc, char **argv )
{

//   TPZSlopeStabilityAnalysis slopeA,slopeB;
//   int porder=2;
//   slopeA.GetCurrentConfig()->CreateGeometricMeshSlope(2);
//   slopeA.GetCurrentConfig()->CreateComputationalMeshSlope ( porder );
//   TPZGeoMesh * gmesh = &slopeA.GetCurrentConfig()->fGMesh;
  TPZGeoMesh  gmesh = ConfigSlope();
  gmesh.Print();

  TPZCompMesh * cmesh = CreateKLCMesh(&gmesh);
  KLConfig(cmesh);
  //ConfigSlope();
  return 0;
}

TPZCompMesh* CreateKLCMesh(TPZGeoMesh * gmesh)
{
  int porder=2;
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  int id=1;
  REAL Lx=20.,Ly=2.;
  int type=3;
  int expansionorder=5;
  KLMaterial * klmat = new KLMaterial(id,Lx,Ly,type,expansionorder);
  cmesh->SetDefaultOrder ( porder);
  cmesh->SetDimModel ( 2 );
  cmesh->InsertMaterialObject(klmat);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}

void KLConfig(TPZCompMesh * cmesh)
{

  
   KLAnalysis * klanal = new KLAnalysis(cmesh);
   klanal->Solve();
   
   	// Post processing
	TPZManVector<std::string> scalarnames(6), vecnames(0);
	scalarnames[0] = "vec";
    scalarnames[1] = "vec1";
    scalarnames[2] = "vec2";
    scalarnames[3] = "vec3";
    scalarnames[4] = "vec4";
    scalarnames[5] = "vec5";

    TPZAnalysis * anal = new TPZAnalysis(cmesh);
    
	anal->DefineGraphMesh(2,scalarnames,vecnames,"klsol.vtk");

	anal->PostProcess(2);
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
  for ( int i=1; i<=steps; i++ )
    {
      std::cout << "fac = " << factor<< std::endl;
      well.LoadingRamp ( factor );
      well.ExecuteInitialSimulation ( 5 );
      factor+=delta;
    }

  int nref=2;
  for ( int i=1; i<=nref; i++ )
    {
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


  
std::ofstream file("dsaadfadfd.vtk");
cout << "well.GetConfigListSize()" << well.GetConfigListSize() << std::endl;
TPZGeoMesh *gmesh= &well.GetCurrentConfig()->fGMesh;

TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file,true);

TPZCompMesh c=well.GetCurrentConfig()->fCMesh;
TPZGeoMesh *g = c.Reference();

  return *g;
}

