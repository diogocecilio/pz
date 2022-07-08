
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
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzelastoplastic.h"
#include "pzelastoplasticanalysis.h"
#include <sstream>
#include "pzmaterial.h"

#include "TPZPlasticStepPV.h"
#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplastic.h"
#include "pzgeotetrahedra.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapetetra.h"
#include "tpzgeoelrefpattern.h"

#include "TPZDarcyFlowIsoPerm.h"


using namespace std;


typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;
typedef   TPZMatElastoPlastic < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat3d;


#include "readgidmesh.h"
#include <TPZParFrontStructMatrix.h>


TPZGeoMesh * CreateGMesh ( int ref );

TPZGeoMesh * CreateGMeshGid ( int ref );

TPZCompMesh * CreateCMesh( TPZGeoMesh *gmesh, int pOrder );

void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void ForcingSource(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void SolveDarcyProlem(TPZCompMesh *cmesh,string vtk);

void PostDarcy(TPZAnalysis * analysis,string vtk);

REAL waterheigth=10.;//Water heigth
REAL slopeheigth = 10.;//Slope Heigth
REAL gammawater =10.;//water specific weigth
int beta=45;

int main()
{
    int porder=4;
    
    TPZGeoMesh* gmesh = CreateGMeshGid(1);
    
    std::ofstream files ( "gmesh.vtk" );
    
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    TPZCompMesh * cmesh = CreateCMesh(gmesh,porder);
    
    std::string vtkdarcy ="darcykx2ky1.vtk";
    
    SolveDarcyProlem(cmesh,vtkdarcy);
    
    cout << "End of simulation."<<endl;
    
    return 0;
}

TPZGeoMesh *  CreateGMesh ( int ref )
{
    const std::string name ( "Slope Problem " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

    vector<vector<double>> co= {{0., 0.}, {75., 0.}, {75., 30.}, {45., 30.}, {35., 40.}, {0.,40.},
	
	{35./3., 40.},{2 * 35/3., 40.}, {30., 40.},
	
	{30., 30.}, {60.,30.},
	
	{2* 35./3.,2* 35/3.}, {45., 2* 35/3.},
	
	 {35./3., 35/3.}, {60., 35./3.}

    };
    vector<vector<int>> topol = {
		{0,  1,  14, 13},
		{1,  2,  10, 14}, 
		{14, 10, 3,  12}, 
		{13, 14, 12, 11},
		{11, 12, 3,  9}, 
		{9,  3,  4,  8},
		{11, 9,  8,  7},
		{13, 11, 7, 6},
		{0, 13,  6, 5}
    };


    gmesh->NodeVec().Resize ( co.size() );

    for ( int inode=0; inode<co.size(); inode++ ) {
        coord[0] = co[inode][0];
        coord[1] = co[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }
    TPZVec <long> TopoQuad ( 4 );
    for ( int iel=0; iel<topol.size(); iel++ ) {
        TopoQuad[0] = topol[iel][0];
        TopoQuad[1] = topol[iel][1];
        TopoQuad[2] =	topol[iel][2];
        TopoQuad[3] = topol[iel][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
    }
    int id = topol.size();
    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom

    id++;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth
	
	{
        id++;
        TopoLine[0] = 3;
        TopoLine[1] = 10;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
        
        id++;
        TopoLine[0] = 10;
        TopoLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
    }

	id++;
    TopoLine[0] = 3;
    TopoLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh );//ramp
	


    {
            id++;
            TopoLine[0] = 5;
            TopoLine[1] = 6;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
            
            id++;
            TopoLine[0] = 6;
            TopoLine[1] = 7;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
            
            id++;
            TopoLine[0] = 7;
            TopoLine[1] = 8;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
            
            id++;
            TopoLine[0] = 8;
            TopoLine[1] = 4;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
    }




    id++;
    TopoLine[0] = 0;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//left


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

TPZGeoMesh * CreateGMeshGid ( int ref )
{

    string file;
 
    if(beta==30)
    {
         file ="/home/diogo/projects/pz/data/h10-beta30.msh";
    }
    if(beta==45)
    {
         file ="/home/diogo/projects/pz/data/h10-beta45.msh";
    }
    if(beta==60)
    {
         file ="/home/diogo/projects/pz/data/h10-beta60.msh";
    }
    if(beta==90)
    {
         file ="/home/diogo/projects/pz/data/h10-beta90.msh";
    }



    readgidmesh read = readgidmesh ( file );
    read.ReadMesh();
    TPZFMatrix<int> meshtopology = read.GetTopology();
    TPZFMatrix<REAL> meshcoords = read.GetCoords();
    std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();

    int ndivs = 100000;
    TPZFMatrix<REAL> pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, idstoprigth,idsramp,idstopleft;

    std::vector<std::vector<int>> idsvec;


    TPZManVector<REAL,2> a ( 2 ), b ( 2 );


    a[0] = 0.;
    a[1] = 0.;
    b[0] = 75.;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsbottom );
    idsvec.push_back ( idsbottom );

    a = b;
    b[0] = 75.;
    b[1] = 30.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsright );
    idsvec.push_back ( idsright );

//para demais
    if(beta==30 || beta==45 || beta ==60)
    {
        a = b;
        b[0] = 45.;
        b[1] = 30.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idstoprigth );
        idsvec.push_back ( idstoprigth );
    }else{
        
        a = b;
        b[0] = 37.5;
        b[1] = 30.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idstoprigth );
        idsvec.push_back ( idstoprigth );
        
    }

    
    
    if(beta==30)
    {
        a = b;
        b[0] = 27.675;
        b[1] = 40.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idsramp );
        idsvec.push_back ( idsramp );
    }
    if(beta==45)
    {
        a = b;
        b[0] = 35.;
        b[1] = 40.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idsramp );
        idsvec.push_back ( idsramp );
    }
    
    if(beta==60)
    {
        a = b;
        b[0] = 39.2265;
        b[1] = 40.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idsramp );
        idsvec.push_back ( idsramp );
    }
    
    if(beta==90)
    {
        a = b;
        b[0] = 37.5;
        b[1] = 40.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idsramp );
        idsvec.push_back ( idsramp );
    } 

    
    a = b;
    b[0] = 0.;
    b[1] = 40.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstopleft );
    idsvec.push_back ( idstopleft );

    a = b;
    b[0] = 0.;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsleft );
    idsvec.push_back ( idsleft );


    // for(int i=0;i<idsbottom.size();i++)cout << idsvec[3][i] << endl;

    const std::string name ( "Slope Problem " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

    vector<vector<double>> co;

    int ncoords = meshcoords.Rows();
    co.resize ( ncoords );
    for ( int i=0; i<ncoords; i++ ) {
        co[i].resize ( 2 );
        co[i][0]=meshcoords ( i,0 );
        co[i][1]=meshcoords ( i,1 );
    }
    vector<vector<int>> topol;

    int ntopol = meshtopology.Rows();
    topol.resize ( ntopol );

    for ( int i=0; i<ntopol; i++ ) {
        topol[i].resize ( meshtopology.Cols() );
        for ( int j=0; j<meshtopology.Cols(); j++ ) {
            topol[i][j]=meshtopology ( i,j );
        }
    }

    gmesh->NodeVec().Resize ( co.size() );

    for ( int inode=0; inode<co.size(); inode++ ) {
        coord[0] = co[inode][0];
        coord[1] = co[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }
    if ( meshtopology.Cols() ==4 ) {
        TPZVec <long> TopoQuad ( 4 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
            TopoQuad[0] = topol[iel][0];
            TopoQuad[1] = topol[iel][1];
            TopoQuad[2] = topol[iel][2];
            TopoQuad[3] = topol[iel][3];
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
        }
    }

    if ( meshtopology.Cols() ==3 ) {
        TPZVec <long> TopoTri ( 3 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
            TopoTri[0] =topol[iel][0];
            TopoTri[1] =topol[iel][1];
            TopoTri[2] =topol[iel][2];
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
        }
    }

    if ( meshtopology.Cols() !=3 && meshtopology.Cols() !=4 ) {
        DebugStop();
    }
    TPZVec <long> TopoLine ( 2 );
    int id = topol.size();
    id++;

    for ( int ivec=0; ivec<idsvec.size(); ivec++ ) {
        int nodes = idsvec[ivec].size();
        for ( int inode=0; inode<nodes-1; inode++ ) {
            TopoLine[0] = idsvec[ivec][inode];
            TopoLine[1] = idsvec[ivec][inode+1];
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - ( ivec+1 ), *gmesh );
        }
    }


    gmesh->BuildConnectivity();
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }
    // gmesh->Print(std::cout);
    std::ofstream files ( "ge.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}

TPZCompMesh * CreateCMesh( TPZGeoMesh *gmesh, int pOrder )
{
    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    cmesh->SetDefaultOrder ( pOrder );
    cmesh->SetDimModel ( 2 );
    cmesh->SetAllCreateFunctionsContinuous();

    // Create the material:
    // TPZDarcyFlow *material = new TPZDarcyFlow(m_matID,m_dim);
	int matid=1,dim=2;
    auto *material = new TPZDarcyFlowIsoPerm ( matid,dim );
    //Bet Degan loamy sand 
	//STATE permeability = 6.3e-5;//m/s
	TPZManVector<REAL,3> permeability(3);//m/s

    permeability[0]=2.;permeability[1]=1;permeability[2]=1.;
    material->SetConstantPermeability ( permeability );
	material->SetId(1);
    
    TPZAutoPointer<TPZFunction<STATE> > sourceterm  = new TPZDummyFunction<STATE>(ForcingSource);

    material->SetForcingFunction(sourceterm);
    
    cmesh->InsertMaterialObject ( material );

    // Condition of contours
    TPZFMatrix<STATE>  val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
	
	TPZAutoPointer<TPZFunction<STATE> > pressure = new TPZDummyFunction<STATE>(ForcingBCPressao);
    
    
    //material->SetForcingFunction(rhs);
    

    TPZMaterial * BCondn1 = material->CreateBC ( material, -1, 0, val1, val2 );//tr
    
    TPZMaterial * BCondn2 = material->CreateBC ( material, -2, 0, val1, val2 );//tr
    
    TPZMaterial * BCondn6 = material->CreateBC ( material, -6, 0, val1, val2 );//tr
    
	TPZMaterial * BCond0 = material->CreateBC ( material, -3, 0, val1, val2 );//tr
	//BCond0->SetForcingFunction(pressure);

	TPZMaterial * BCond1 = material->CreateBC ( material, -4, 0, val1, val2 );//ramp
	//BCond1->SetForcingFunction(pressure);
     
	TPZMaterial * BCond2 = material->CreateBC ( material, -5, 0, val1, val2 );//tl
    //BCond2->SetForcingFunction(pressure);

    cmesh->InsertMaterialObject(BCondn1);
    cmesh->InsertMaterialObject(BCondn2);
    cmesh->InsertMaterialObject(BCondn6);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);

    //Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}

void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    const auto &x=pt[0];
    const auto &y=pt[1];
    const auto &z=pt[2];
    REAL yy=40-y;
    REAL hw = waterheigth;//Water heigth
    REAL H = slopeheigth;//Slope Heigth
    REAL gamma = gammawater;//water specific weigth

    if(yy<hw && yy<H)
    {
        disp[0]  =   gamma * yy ;
    }
    else {
        disp[0]  = hw*gamma;
    }
}

void ForcingSource(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    const auto &x=pt[0];
    const auto &y=pt[1];
    const auto &z=pt[2];
    disp[0] =-1;

}

void SolveDarcyProlem(TPZCompMesh *cmesh,string vtk)
{


    int numthreads = 0;
    TPZAnalysis *analysis =  new TPZAnalysis( cmesh ); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );
    matskl.SetNumThreads ( numthreads );
    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    analysis->SetSolver ( step );


    analysis->Run();


    PostDarcy(analysis,vtk);


}

void PostDarcy(TPZAnalysis * analysis,string vtk)
{
    TPZVec<std::string> scalarVars ( 1 ),vectorVars ( 1 );
    scalarVars[0] = "Pressure";
    vectorVars[0] = "Flux";
    analysis->DefineGraphMesh ( 2,scalarVars,vectorVars,vtk);
    constexpr int resolution{0};
    analysis->PostProcess ( resolution );
}
