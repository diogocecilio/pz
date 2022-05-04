
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "KLMaterial.h"
#include "KLStrMatrix.h"
#include "KLAnalysis.h"
#include <iostream>
#include <cstdlib>
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "tpzgeoelrefpattern.h"
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
#include "TPZVTKGeoMesh.h"
using namespace pzgeom;
using namespace pztopology;

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


TPZCompMesh* CreateKLCMesh ( TPZGeoMesh * gmesh,int porder,REAL Lx, REAL Ly,REAL Lz, int type, int expansionorder );
TPZVec<TPZFMatrix<REAL>> CreateFields ( TPZCompMesh * cmesh,KLAnalysis * klanal );
TPZGeoMesh * Simple3DGMesh ( int ref );
void KLConfig();

int main ( int argc, char **argv )
{
    KLConfig();
    std::cout << " SUCESS " << std::endl;
    return 0;
}
void KLConfig()
{

    int porder=2;
    int expansionorder=100;
    int ref=2;

    TPZGeoMesh * gmesh =  Simple3DGMesh (ref);



    REAL Lx=10;
    REAL Ly=0.5;
    REAL Lz=0.;
    int type=4;


    TPZCompMesh * cmesh = CreateKLCMesh ( gmesh,porder,Lx,Ly,Lz,type, expansionorder );
    KLAnalysis * klanal = new KLAnalysis ( cmesh );
    KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );
    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );
    klanal->Solve();


    int dim = gmesh->Dimension();
    ref=3;
    string file ="outexp2h2.vtk";
    klanal->Post(file,dim,ref);

}



TPZGeoMesh * Simple3DGMesh ( int ref )
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;


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

    //cmesh->Print(std::cout);


    return cmesh;
}
