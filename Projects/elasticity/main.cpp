
#include <cmath>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include <pzgmesh.h>
#include <pzcmesh.h>
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzfunction.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzfstrmatrix.h>
#include <pzstepsolver.h> //for TPZStepSolver
#include <sstream>
#include "pzmaterial.h"
#include "pzelasmat.h"
using namespace std;

TPZGeoMesh * CreateGeoMesh();

TPZCompMesh * CreateCompMesh(TPZGeoMesh * gmesh);

TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh);


int main()
{

    cout << "Hello World"<<endl;
    TPZGeoMesh * gmesh = CreateGeoMesh();
    TPZCompMesh * cmesh = CreateCompMesh(gmesh);
    TPZAnalysis * analysis = CreateAnalysis(cmesh);
    analysis->Run();

    TPZVec<std::string> scalarVars ( 3 ),vectorVars ( 1 );
    vectorVars[0] = "Displacement";
    scalarVars[0] = "SigmaX";
    scalarVars[1] = "SigmaY";
    scalarVars[2] = "TauXY";
    string vtk = "postprocess.vtk";
    analysis->DefineGraphMesh ( 2,scalarVars,vectorVars,vtk);
    analysis->PostProcess ( 0 );

    cmesh->Print(cout);

    return 0;
}

TPZGeoMesh * CreateGeoMesh()
{
     REAL coordinates [4][2] = {	{0.,0.},{10.,0.},{10.,1.},
                                {0.,1.}};


    long elindex;
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->NodeVec().Resize(4);
    TPZVec<double> co(2);
    for(int inode=0;inode<4;inode++)
    {
        co[0]=coordinates[inode][0];
        co[1]=coordinates[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, co, *gmesh );
    }

    int matid=1;
    TPZVec<long> top(4);
    top[0]=0;
    top[1]=1;
    top[2]=2;
    top[3]=3;
    gmesh->CreateGeoElement(EQuadrilateral,top,matid,elindex);

    matid=-1;
    TPZVec<long> topline(2);
    topline[0]=1;
    topline[1]=2;
    gmesh->CreateGeoElement(EOned,topline,matid,elindex);

    matid=-2;
    TPZVec<long> toppoint(1);
    toppoint[0]=3;
    gmesh->CreateGeoElement(EPoint,toppoint,matid,elindex);

    gmesh->BuildConnectivity();
    gmesh->Print();



    for ( int d = 0; d<0; d++ ) {
    int nel = gmesh->NElements();
    TPZVec<TPZGeoEl *> subels;
    for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    std::ofstream files ( "gesimple.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;

}

TPZCompMesh * CreateCompMesh(TPZGeoMesh * gmesh)
{

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(2);
    cmesh->SetDimModel(2);

    int id=1;
    REAL E=8000.;
    REAL nu=0.;
    REAL fx=0.;
    REAL fy=0.;
    int planestress =1;
    TPZElasticityMaterial*elasticmat = new TPZElasticityMaterial(id,E,nu,fx,fy,planestress);

    cmesh->InsertMaterialObject(elasticmat);


    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
    int dirichlet =0;
    auto * clamp = elasticmat->CreateBC ( elasticmat, -1,dirichlet, val1, val2 );//bottom

    int newman =1;
    val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = -2;
    auto * load = elasticmat->CreateBC ( elasticmat, -2, newman, val1, val2 );//rigth

    cmesh->InsertMaterialObject ( clamp );
    cmesh->InsertMaterialObject ( load );

    cmesh->SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();

    return cmesh;
}

TPZAnalysis * CreateAnalysis ( TPZCompMesh *cmesh )
{
    int numthreads=0;

    TPZAnalysis * analysis =  new TPZAnalysis ( cmesh ); // Create analysis

    TPZFStructMatrix matskl ( cmesh );

    matskl.SetNumThreads ( numthreads );

    analysis->SetStructuralMatrix ( matskl );

    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );

    analysis->SetSolver ( step );

    analysis->SetSolver ( step );


    return analysis;
}
