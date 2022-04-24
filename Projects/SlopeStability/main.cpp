
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"


#include "pzbndcond.h"
#include "pzfunction.h"
//#include "TPZDarcyMaterial.h"
#include "TPZDarcyFlow.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <pzmanvector.h>//for TPZManVector
//#include "pznondarcyflow.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver

using namespace std;
class TPZMaterial;



int pOrder          =  3; // Polynomial order of approximation
int m_matID         =  1; // Material id
int m_dim         =  2; // Material id

// Boundary condition
int dirichlet       =  0;
int neumann         =  1;
int robin     =  2;

//bc slope
int bottom_slope = -1;
int rigth_slope = -2;
int left_slope = -3;
int topleft_slope = -4;
int toprigth_slope = -5;
int ramp_slope = -6;


// brief Function to create the geometric mesh
TPZGeoMesh *CreateGMesh ( int nelx, int nely, double hx, double hy );
void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBCZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
// brief Function to refine the geometric mesh
void UniformRefine ( TPZGeoMesh* gmesh, int nDiv );

#include "readgidmesh.h"
// brief Function to create  the multi-physical computational mesh
TPZCompMesh *CMesh_m ( TPZGeoMesh *gmesh, int pOrder );

TPZGeoMesh * CreateGeometricMeshSlope ( int ref );

TPZGeoMesh * CreateGeometricMeshSlopeGid ( int ref );

TPZCompMesh * CreateComputationalMeshSlope ( TPZGeoMesh *gmesh, int pOrder );

int main()
{
	cout << "here "<<endl;
    //TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy);
    //TPZGeoMesh *gmesh = CreateGeometricMeshSlope ( 3);
	TPZGeoMesh *gmesh = CreateGeometricMeshSlopeGid ( 0 );

    std::ofstream meshfile ( "GeoMesh.txt" );
    gmesh->Print ( meshfile );

    // Create computational mesh
    //TPZCompMesh *cmesh = CMesh_m(gmesh, pOrder);
    TPZCompMesh *cmesh = CreateComputationalMeshSlope ( gmesh, pOrder );
    std::ofstream comVfile ( "CompMesh.txt" );
    cmesh->Print ( comVfile );

    //Resolvendo o Sistema:
    int numthreads = 0;
    bool optimizeBandwidth = false; // Prevents of renumbering of the equations (As the same of Oden's result)
    TPZAnalysis analysis ( cmesh,optimizeBandwidth ); // Create analysis

    //sets number of threads to be used by the solver

    //defines storage scheme to be used for the FEM matrices
    //in this case, a symmetric skyline matrix is used
    TPZSkylineStructMatrix matskl ( cmesh );
    matskl.SetNumThreads ( numthreads );
    analysis.SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    analysis.SetSolver ( step );

    //for(int i=0;i<2;i++){

    analysis.Assemble(); // Assembla the global matrix

    std::cout << "Solving Matrix " << std::endl;
    analysis.Solve();

	analysis.Solution().Print(std::cout);
	
    ///vtk export
    TPZVec<std::string> scalarVars ( 1 ),vectorVars ( 1 );
    scalarVars[0] = "Pressure";
    vectorVars[0] = "Flux";
    analysis.DefineGraphMesh ( 2,scalarVars,vectorVars,"Darcy.vtk" );
    constexpr int resolution{0};
    analysis.PostProcess ( resolution );

    //}

    std::cout << "FINISHED!" << std::endl;


    return 0;
}


TPZGeoMesh * CreateGeometricMeshSlopeGid ( int ref )
{
	
	string file ="/home/diogo/projects/pz/data/mesh-teste-pz-fromathematica.msh";

	readgidmesh read = readgidmesh(file);
	read.ReadMesh();
	TPZFMatrix<int> meshtopology = read.GetTopology();
	TPZFMatrix<REAL> meshcoords = read.GetCoords();
	std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();
	
	TPZFMatrix<REAL> poincoords(6,3);
	poincoords(0,0)=0.; poincoords(0,1)=0.; poincoords(0,2)=0.;
	poincoords(1,0)=75.;poincoords(1,1)=0.; poincoords(1,2)=0.;
	poincoords(2,0)=75.;poincoords(2,1)=30.;poincoords(2,2)=0.;
	poincoords(3,0)=45.;poincoords(3,1)=30.;poincoords(3,2)=0.;
	poincoords(4,0)=30.;poincoords(4,1)=40.;poincoords(4,2)=0.;
	poincoords(5,0)=0.; poincoords(5,1)=40.;poincoords(5,2)=0.;
	std::vector<int> idspoints;
	for(int i =0;i<poincoords.Rows();i++)
 	{
	 
		//busca ponto
		TPZVec<double> constcoorddata ( 3,0. );
    	constcoorddata[0]=poincoords(i,0);
    	constcoorddata[1]=poincoords(i,1);
    	constcoorddata[2]=poincoords(i,2);
    	std::vector<int> idsline;
		//Direcao a buscar. 0 significa que o algoritmo e livre para buscar naquela direcao e 1 que dizer que é fixo.
    	TPZVec<int>constcoord2 ( 3 );
    	constcoord2[0]=1;//fixo
   	 	constcoord2[1]=1;//fixo
    	constcoord2[2]=1;//fixo
		read.FindIds(constcoorddata,constcoord2,idsline);

		idspoints.push_back(idsline[0]);
	
	}
	
	cout << "idspoints[i] "<< endl;
	for(int i=0;i<idspoints.size();i++)cout << idspoints[i] << endl;
	
	const std::string name ( "Slope Problem " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

    vector<vector<double>> co;
	
	int ncoords = meshcoords.Rows(); 
	co.resize(ncoords);
	for(int i=0;i<ncoords;i++)
 	{
		co[i].resize(2);
		co[i][0]=meshcoords(i,0);
		co[i][1]=meshcoords(i,1);
	}
    vector<vector<int>> topol;
	
	int ntopol = meshtopology.Rows(); 
	topol.resize(ntopol);
	
	for(int i=0;i<ntopol;i++)
 	{
		topol[i].resize(meshtopology.Cols());
		for(int j=0;j<meshtopology.Cols();j++)
  		{
	  		topol[i][j]=meshtopology(i,j);
		}
	}

    gmesh->NodeVec().Resize ( co.size() );

    for ( int inode=0; inode<co.size(); inode++ ) {
        coord[0] = co[inode][0];
        coord[1] = co[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }
    if(meshtopology.Cols()==4)
	{
		TPZVec <long> TopoQuad ( 4 );
    	for ( int iel=0; iel<topol.size(); iel++ ) {
        	TopoQuad[0] = topol[iel][0];
        	TopoQuad[1] = topol[iel][1];
        	TopoQuad[2] =topol[iel][2];
        	TopoQuad[3] = topol[iel][3];
        	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
    	}
	}

	if(meshtopology.Cols()==3)
	{
		TPZVec <long> TopoTri ( 3 );
    	for ( int iel=0; iel<topol.size(); iel++ ) {
        	TopoTri[0] =topol[iel][0];
        	TopoTri[1] =topol[iel][1];
        	TopoTri[2] =topol[iel][2];
        	new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
    	}
	}
	
	if(meshtopology.Cols()!=3 && meshtopology.Cols()!=4)
 	{
		DebugStop(); 
	}
    
	int id = topol.size();
    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = idspoints[0];
    TopoLine[1] = idspoints[1];
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, bottom_slope, *gmesh );//bottom

    id++;
    TopoLine[0] = idspoints[1];
    TopoLine[1] = idspoints[2];
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, rigth_slope, *gmesh );//rigth

    id++;
    TopoLine[0] = idspoints[0];
    TopoLine[1] = idspoints[5];
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, left_slope, *gmesh );//left

    id++;
    TopoLine[0] = idspoints[4];
    TopoLine[1] = idspoints[5];
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, topleft_slope, *gmesh ); //top left

    id++;
    TopoLine[0] = idspoints[2];
    TopoLine[1] = idspoints[3];
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, toprigth_slope, *gmesh ); // top rigth


    id++;
    TopoLine[0] = idspoints[3];
    TopoLine[1] = idspoints[4];
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, ramp_slope, *gmesh );

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

TPZGeoMesh *  CreateGeometricMeshSlope ( int ref )
{
    const std::string name ( "Slope Problem " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );


//     vector<vector<double>> co= {{0., 0.}, {75., 0.}, {75., 30.}, {45., 30.}, {35., 40.}, {0.,40.},
// 	
// 	{22., 40.},{26., 40.}, {30., 40.},
// 	
// 	{30., 30.}, {50.,30.},
// 	
// 	{26., 26.}, {45., 26.},
// 	
// 	 {22., 22.}, {50., 22.}
// 
//     };
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
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, bottom_slope, *gmesh );//bottom

    id++;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, rigth_slope, *gmesh );//rigth

    id++;
    TopoLine[0] = 0;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, left_slope, *gmesh );//left

    id++;
    TopoLine[0] = 4;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, topleft_slope, *gmesh ); //top left

    id++;
    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, toprigth_slope, *gmesh ); // top rigth


    id++;
    TopoLine[0] = 3;
    TopoLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, ramp_slope, *gmesh );

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
void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        REAL gamma  = 997.;//kg/m³
        disp[0] = -( 40-y ) *gamma;
}
void ForcingBCZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        disp[0] = 0;
}

TPZCompMesh * CreateComputationalMeshSlope ( TPZGeoMesh *gmesh, int pOrder )
{
    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    cmesh->SetDefaultOrder ( pOrder );
    cmesh->SetDimModel ( 2 );
    cmesh->SetAllCreateFunctionsContinuous();

    // Create the material:
    // TPZDarcyFlow *material = new TPZDarcyFlow(m_matID,m_dim);
    auto *material = new TPZDarcyFlow ( m_matID,m_dim );
    
	STATE permeability = 0.00001;

    material->SetConstantPermeability ( permeability );
	material->SetId(1);

    cmesh->InsertMaterialObject ( material );

    // Condition of contours
    TPZFMatrix<STATE>  val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 1,1,0. );


    REAL gamma  = 997.;//kg/m³
    REAL slopeheigth1 = 40.;//metros

//     auto pressure = [] ( const TPZVec<REAL> &loc,
//                          TPZVec<STATE>&u,
//     TPZFMatrix<STATE>&gradU ) {
//         const auto &x=loc[0];
//         const auto &y=loc[1];
//         const auto &z=loc[2];
//         REAL gamma  = 997.;//kg/m³
//         u[0] = -( 40-y ) *gamma;
// 
//     };
// 
//     auto zero = [] ( const TPZVec<REAL> &loc,
//                      TPZVec<STATE>&u,
//     TPZFMatrix<STATE>&gradU ) {
//         u[0] = 0.;
//     };
	
	//TPZDummyFunction<STATE>(ForcingBCPressao);
	
	TPZAutoPointer<TPZFunction<STATE> > pressure = new TPZDummyFunction<STATE>(ForcingBCPressao);
	TPZAutoPointer<TPZFunction<STATE> > zero = new TPZDummyFunction<STATE>(ForcingBCZero);

    //val2[0]=0.;
    //val2[1]=0.;
    TPZMaterial * BCond0 = material->CreateBC ( material, bottom_slope, dirichlet, val1, val2 );
	BCond0->SetForcingFunction(pressure);
	//BCond0->SetForcingFunction(0,pressure);
    //BCond0->SetForcingFunctionBC ( pressure );


    TPZMaterial * BCond1 = material->CreateBC ( material, left_slope, dirichlet, val1, val2 );
	BCond1->SetForcingFunction(zero);
	//BCond1->SetForcingFunction(0,zero);
    //BCond1->SetForcingFunctionBC ( zero );

    TPZMaterial * BCond2 = material->CreateBC ( material, topleft_slope, dirichlet, val1, val2 );
	BCond2->SetForcingFunction(pressure);
	//BCond2->SetForcingFunction(0,pressure);
    //BCond2->SetForcingFunctionBC ( pressure );

    TPZMaterial * BCond3 = material->CreateBC ( material, toprigth_slope, dirichlet, val1, val2 );
	BCond3->SetForcingFunction(pressure);
	//BCond3->SetForcingFunction(0,pressure);
    //BCond3->SetForcingFunctionBC ( pressure );

    TPZMaterial * BCond4 = material->CreateBC ( material, ramp_slope, dirichlet, val1, val2 );
	BCond4->SetForcingFunction(pressure);
	//BCond4->SetForcingFunction(0,pressure);
    //BCond4->SetForcingFunctionBC ( pressure );

    TPZMaterial * BCond5 = material->CreateBC ( material, rigth_slope,dirichlet, val1, val2 );
	BCond5->SetForcingFunction(zero);
	//BCond5->SetForcingFunction(0,zero);
    //BCond5->SetForcingFunctionBC ( zero );

    //cmesh->InsertMaterialObject(BCond0);
    //cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject ( BCond2 );
    cmesh->InsertMaterialObject ( BCond3 );
    cmesh->InsertMaterialObject ( BCond4 );
    //cmesh->InsertMaterialObject(BCond5);

    //Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}
