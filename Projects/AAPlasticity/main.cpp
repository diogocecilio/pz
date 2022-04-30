
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

using namespace std;


typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;




#include "readgidmesh.h"

TPZGeoMesh * CreateGMesh ( int ref );

TPZGeoMesh * CreateGMeshGid ( int ref );

TPZCompMesh * CreateCMesh(TPZGeoMesh * gmesh,int porder);

void LoadingRamp (TPZCompMesh * cmesh,  REAL factor);

void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);

void  CreatePostProcessingMesh(TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh);

void Post(TPZCompMesh * cmesh);

int main()
{
	
	int porder =1;
	int samples =10;

 	TPZGeoMesh *gmesh = CreateGMeshGid ( 0 );

    TPZCompMesh *cmesh = CreateCMesh ( gmesh, porder );

    //Resolvendo o Sistema:
    int numthreads = 0;
	
    bool optimizeBandwidth = false; // Prevents of renumbering of the equations (As the same of Oden's result)
    
    TPZAnalysis analysis ( cmesh,optimizeBandwidth ); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );
	
    matskl.SetNumThreads ( numthreads );
	
    analysis.SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    analysis.SetSolver ( step );

	
	

    std::cout << "FINISHED!" << std::endl;


    return 0;
}

void Post(TPZCompMesh * cmesh)
{
	TPZPostProcAnalysis  * postproc = new TPZPostProcAnalysis();
	
	CreatePostProcessingMesh( postproc, cmesh);
	
	//postproc->SetCompMesh(cmesh);
	
	const std::string name = "out";
	
	//postproc->Print(name, std::cout);
	
	
	TPZVec<int> PostProcMatIds(1,1);
	
	TPZStack<std::string> PostProcVars, scalNames, vecNames;
	
	PostProcessVariables(scalNames, vecNames);

    std::string vtkFile = "outFile.vtk";

    postproc->DefineGraphMesh(2,scalNames,vecNames,vtkFile);

    postproc->PostProcess(0);
}


TPZGeoMesh * CreateGMeshGid ( int ref )
{
	

	string file ="/home/diogo/projects/pz/data/mesh-teste-pz-fromathematica.msh";
    //string file ="/home/diogo/projects/pz/data/quad-gid2.msh";
	//string file ="/home/diogo/projects/pz/data/quad-gid.msh";

    

	readgidmesh read = readgidmesh(file);
	read.ReadMesh();
	TPZFMatrix<int> meshtopology = read.GetTopology();
	TPZFMatrix<REAL> meshcoords = read.GetCoords();
	std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();

    int ndivs = 10000;
    TPZFMatrix<REAL> pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, idstoprigth,idsramp,idstopleft;
    
    std::vector<std::vector<int>> idsvec;
    
    
    TPZManVector<REAL,2> a ( 2 ), b ( 2 );
    
    
    a[0] = 0.; a[1] = 0.;
    b[0] = 75.; b[1] = 0.;
    read.Line( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsbottom );
    idsvec.push_back(idsbottom);
    
    a = b;
    b[0] = 75.; b[1] = 30.;
    read.Line( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsright );
    idsvec.push_back(idsright);
    
    
    a = b;
    b[0] = 45.; b[1] = 30.;
    read.Line( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstoprigth );
    idsvec.push_back(idstoprigth);
    
    a = b;
    b[0] = 35.; b[1] = 40.;
    read.Line( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsramp );
    idsvec.push_back(idsramp);
    
    a = b;
    b[0] = 0.; b[1] = 40.;
    read.Line( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstopleft );
    idsvec.push_back(idstopleft);
    
    a = b;
    b[0] = 0.; b[1] = 0.;
    read.Line( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsleft );
    idsvec.push_back(idsleft);
    
    
   // for(int i=0;i<idsbottom.size();i++)cout << idsvec[3][i] << endl;

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
        	TopoQuad[2] = topol[iel][2];
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
    TPZVec <long> TopoLine ( 2 );
    int id = topol.size();
    id++;
    
    for(int ivec=0;ivec<idsvec.size();ivec++)
    {
        int nodes = idsvec[ivec].size();
        for(int inode=0;inode<nodes-1;inode++)
        {
            TopoLine[0] = idsvec[ivec][inode];
            TopoLine[1] = idsvec[ivec][inode+1];
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -(ivec+1), *gmesh );
        }
    }
    

    gmesh->BuildConnectivity();

    gmesh->Print(std::cout);
 	std::ofstream files ( "ge.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,files,false);
    
    return gmesh;
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

    id++;
    TopoLine[0] = 0;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh );//left

    {
            id++;
            TopoLine[0] = 5;
            TopoLine[1] = 6;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
            
            id++;
            TopoLine[0] = 6;
            TopoLine[1] = 7;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
            
            id++;
            TopoLine[0] = 7;
            TopoLine[1] = 8;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
            
            id++;
            TopoLine[0] = 8;
            TopoLine[1] = 4;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
    }


    {
        id++;
        TopoLine[0] = 3;
        TopoLine[1] = 10;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); // top rigth
        
        id++;
        TopoLine[0] = 10;
        TopoLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); // top rigth
    }


    id++;
    TopoLine[0] = 3;
    TopoLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//ramp

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

/// create the computational mesh
TPZCompMesh * CreateCMesh(TPZGeoMesh * gmesh,int porder)
{
	unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

	TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    // Mohr Coulomb data
    REAL mc_cohesion    = 10.0;
    REAL mc_phi         = ( 30.0*M_PI/180 );
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000.;

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetUp( E, nu );
	LEMC.fER =ER;
   // LEMC.SetElasticResponse( ER );
    LEMC.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );
	
    int PlaneStrain = 1;
	
	plasticmat * material = new plasticmat ( 1,PlaneStrain );
    material->SetPlasticity ( LEMC );

	material->SetId(1);
    cmesh->InsertMaterialObject ( material );

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
	int directionadirichlet =3;
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 1;
    auto * bc_bottom = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_rigth = material->CreateBC ( material, -2, directionadirichlet, val1, val2 );
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_left = material->CreateBC ( material, -3, directionadirichlet, val1, val2 );

    cmesh->InsertMaterialObject ( bc_bottom );
   	cmesh->InsertMaterialObject ( bc_rigth );
    cmesh->InsertMaterialObject ( bc_left );
	
    //cmesh->InsertMaterialObject ( top );
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();
	
	return cmesh;

}

void LoadingRamp (TPZCompMesh * cmesh,  REAL factor)
{
	plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=-factor*20.;
    body->SetBodyForce ( force );

}

void  CreatePostProcessingMesh(TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh)
{

    if (PostProcess->ReferenceCompMesh() != cmesh)
    {
        
        PostProcess->SetCompMesh(cmesh);
        
        TPZVec<int> PostProcMatIds(1,1);
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        PostProcessVariables(scalNames, vecNames);
        
        for (int i=0; i<scalNames.size(); i++) {
            PostProcVars.Push(scalNames[i]);
        }
        for (int i=0; i<vecNames.size(); i++) {
            PostProcVars.Push(vecNames[i]);
        }
        //
        PostProcess->SetPostProcessVariables(PostProcMatIds, PostProcVars);
        TPZFStructMatrix structmatrix(PostProcess->Mesh());
        structmatrix.SetNumThreads(0);
        PostProcess->SetStructuralMatrix(structmatrix);
    }
    //
    //Chamar com o analysis e nao com o postanalysis pois tem o acumulo de sols
    PostProcess->TransferSolution();
	
    //fElastoplasticAnalysis.TransferSolution(fPostprocess);
}


/// Get the post processing variables
void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames)
{
	  scalNames.Push("StrainVol");
  scalNames.Push("StrainXX");
  scalNames.Push("StrainYY");
  scalNames.Push("StrainZZ");
  scalNames.Push("StrainXY");
  scalNames.Push("StrainXZ");
  scalNames.Push("StrainYZ");

  scalNames.Push("ElStrainVol");
  scalNames.Push("ElStrainXX");
  scalNames.Push("ElStrainYY");
  scalNames.Push("ElStrainZZ");
  scalNames.Push("ElStrainXY");
  scalNames.Push("ElStrainXZ");
  scalNames.Push("ElStrainYZ");

  scalNames.Push("PlStrainVol");
  scalNames.Push("PlStrainXX");
  scalNames.Push("PlStrainYY");
  scalNames.Push("PlStrainZZ");
  scalNames.Push("PlStrainXY");
  scalNames.Push("PlStrainXZ");
  scalNames.Push("PlStrainYZ");
  
  scalNames.Push("PlStrainSqJ2");
  scalNames.Push("PlStrainSqJ2El");
  scalNames.Push("PlAlpha");

  scalNames.Push("DisplacementX");
  scalNames.Push("DisplacementY");
  scalNames.Push("DisplacementZ");
  vecNames.Push("DisplacementTotal");


  scalNames.Push("YieldSurface1");
  scalNames.Push("YieldSurface2");
  scalNames.Push("YieldSurface3");
  
  scalNames.Push("POrder");
  scalNames.Push("NSteps");

  
}
