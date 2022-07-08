
#define  PZPLAST
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
//#include "TPZDarcyFlow.h"

using namespace std;


typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;
typedef   TPZMatElastoPlastic < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat3d;

//TPZVec<REAL> fPlasticDeformSqJ2;

#include "readgidmesh.h"
#include <TPZParFrontStructMatrix.h>

TPZGeoMesh * CreateGMesh ( int ref );

TPZGeoMesh * CreateGMeshGid ( int ref );

TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder );

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor );
void SetSuportPressure ( TPZCompMesh * cmesh,  REAL factor );

void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void Post ( TPZPostProcAnalysis * postproc,std::string vtkFile );

TPZElastoPlasticAnalysis * CreateAnal ( TPZCompMesh *cmesh,bool optimize );

void ShearRed ( TPZCompMesh * cmesh );

void GravityIncrease ( TPZCompMesh * cmesh );

REAL ShearRed ( TPZCompMesh * cmesh,int isample,TPZCompMesh * incmesh );

void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphix,TPZSolVec &sol, TPZGradSolVec &dsol );


void SetFlux ( TPZCompMesh * plasticCmesh,TPZCompMesh* incmesh);

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );

void PrintMat ( std::string out,TPZFMatrix<REAL> mat );

TPZManVector<TPZCompMesh *,2> CreateFieldsDummy ( TPZGeoMesh * gmesh,int porder );


void LoadingRampRainFall ( TPZCompMesh * cmesh,  REAL factor );

template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs );

void SolveDarcyProlem(TPZCompMesh *cmesh,string vtk);

void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

TPZCompMesh * CreateCMeshDarcyDummy( TPZGeoMesh *gmesh, int pOrder );

TPZCompMesh * CreateCMeshDarcy( TPZGeoMesh *gmesh, int pOrder );

void PostDarcy(TPZAnalysis * analysis,string vtk);

REAL findnodalsol(TPZCompMesh *cmesh);

void IncrementalSolution();

int main()
{
	//IncrementalSolution();
	
	

    //return 0;
    
	chrono::steady_clock sc;
	auto start = sc.now();


    int order=2;
    TPZGeoMesh * gmesh =  CreateGMeshGid ( 0 );

    TPZCompMesh*cmesh = CreateCMeshDarcy(gmesh,order);
    string vtk  = "darcy2.vtk";
    std::string vtkplasticity ="plasticity.vtk";
    //SolveDarcyProlem(cmesh,vtk);
    
    TPZCompMesh *cmeshgi = CreateCMesh ( gmesh,order );
	
	
	cout << "\n Setting flux in mechanic comp mesh... " << endl;
	//SetFlux(cmeshsrm,darcycompmesh);
	//SetFlux(cmeshgi,cmesh);
	
    
    cout << "\n Gravity Increase routine.. " << endl;
    GravityIncrease ( cmeshgi );
    
    TPZPostProcAnalysis * postprocdetergim = new TPZPostProcAnalysis();
    
    cout << "\n Post Processing gravity increase... " << endl;
    CreatePostProcessingMesh ( postprocdetergim, cmeshgi );
    Post ( postprocdetergim,vtkplasticity );
    std::cout << "FINISHED!" << std::endl;

	auto end = sc.now();
	auto time_span = static_cast<chrono::duration<double>> ( end - start );
	cout << "| total time taken =  " << time_span.count()<< std::endl;
    return 0;
}

void IncrementalSolution()
{

	int porder =2;

	cout << "aqui "<<endl;
 	TPZGeoMesh *gmesh = CreateGMeshGid (0);
	
    TPZCompMesh *cmesh = CreateCMesh ( gmesh, porder );

	
	int maxiter=100;
	REAL tol=1;
	bool linesearch=true;
	bool checkconv=false;
	int steps=4;
	REAL finalload = 0.25;
	std::string vtkFile ="tunnelbodyload.vtk";
	std::ofstream outloadu("loadvu.nb");
	outloadu << "plot = {";
	TPZPostProcAnalysis  * postproc = new TPZPostProcAnalysis();
	
	REAL uy=0;
	for(int iloadstep=0;iloadstep<=0;iloadstep++)
 	{
		
		TPZElastoPlasticAnalysis  * analysis =  CreateAnal(cmesh,porder);

		
		REAL load = iloadstep*finalload/steps;
		cout << " \n --------- iloadstep  = "<< iloadstep << endl;
		
		cout << " \n --------- load  = "<< load << endl;
		
		//LoadingRamp(cmesh,load);
		int iters;
		bool conv = analysis->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv,iters);
		uy+=findnodalsol(cmesh);
		TPZFMatrix<REAL> sol = analysis->Solution();
		
		//cmesh->Solution().Print(cout);
		cout << " sol "<<-uy << " | load = " << load <<    endl;
		analysis->AcceptSolution();
		
		outloadu << "{ "<<-uy << ", " << load << " } ," << endl;
		
					postproc->SetCompMesh(0);
		CreatePostProcessingMesh( postproc, cmesh);
		
		Post(postproc,vtkFile);


	}
	outloadu <<  " }; ListLinePlot[plot,PlotRange->All]" << endl;

}


TPZGeoMesh * CreateGMeshGid ( int ref )
{

    string file;

    //file ="/home/diogo/projects/pz/data/tunnel.msh";
	file ="/home/diogo/projects/pz/data/tunnel-medium.msh";
    //file ="/home/diogo/projects/pz/data/tunnel-fine.msh";

    readgidmesh read = readgidmesh ( file );
	
    read.ReadMesh();
	
    TPZFMatrix<int> meshtopology = read.GetTopology();
    TPZFMatrix<REAL> meshcoords = read.GetCoords();
    std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();
    TPZVec<double> constcoorddata(3,0.);
    TPZVec<int> constcoord(3);
    std::vector<int> ids1,ids2,ids3,ids4,ids5,ids6,ids7;
    std::vector<std::vector<int>> idsvec;


    //face -1 (plano x=0 zy)
    constcoorddata[0]=0;
    constcoorddata[1]=20;
    constcoorddata[2]=20;
    constcoord[0]=1;
    constcoord[1]=0;
    constcoord[2]=0;//direceo x fixa
    read.FindElements ( constcoorddata,constcoord, ids1 );
    idsvec.push_back(ids1);
    //face -2 (plano x=50 zy)
    constcoorddata[0]=50;
    constcoorddata[1]=20;
    constcoorddata[2]=20;
    constcoord[0]=1;
    constcoord[1]=0;
    constcoord[2]=0;//direceo x fixa
    read.FindElements ( constcoorddata,constcoord, ids2 );
    idsvec.push_back(ids2);
    //face -3 (plano y=0 zx)
    constcoorddata[0]=25;
    constcoorddata[1]=0;
    constcoorddata[2]=20;
    constcoord[0]=0;
    constcoord[1]=1;
    constcoord[2]=0;//direceo y fixa
    read.FindElements ( constcoorddata,constcoord, ids3 );
    idsvec.push_back(ids3);
    //face -4 (plano y=40 zx)
    constcoorddata[0]=25;
    constcoorddata[1]=40;
    constcoorddata[2]=20;
    constcoord[0]=0;
    constcoord[1]=1;
    constcoord[2]=0;//direceo y fixa
    read.FindElements ( constcoorddata,constcoord, ids4 );
    idsvec.push_back(ids4);
    //face -5 (plano z=0 xy)
    constcoorddata[0]=25;
    constcoorddata[1]=20;
    constcoorddata[2]=0;
    constcoord[0]=0;
    constcoord[1]=0;
    constcoord[2]=1;//direceo z fixa
    read.FindElements ( constcoorddata,constcoord, ids5 );
    idsvec.push_back(ids5);
    //face -6 (plano z=40 xy)
    constcoorddata[0]=25;
    constcoorddata[1]=20;
    constcoorddata[2]=40;
    constcoord[0]=0;
    constcoord[1]=0;
    constcoord[2]=1;//direceo z fixa
    read.FindElements ( constcoorddata,constcoord, ids6 );
    idsvec.push_back(ids6);
    //face -7 (plano x=10 yz)
    constcoorddata[0]=20;
    constcoorddata[1]=20;
    constcoorddata[2]=40;
    constcoord[0]=1;
    constcoord[1]=0;
    constcoord[2]=0;//direceo x fixa
    read.FindElements ( constcoorddata,constcoord, ids7 );
    idsvec.push_back(ids7);
    //for(int i=0;i<ids7.size();i++)cout<< ids7[i]  <<endl;
    //for(int i=0; i<ids7.size(); i++)cout<< ids7[i] << " " <<  meshtopology(ids7[i],0) << " "<< meshtopology(ids7[i],1) << " "<< //meshtopology(ids7[i],2)  <<endl;
    /*
    for(int i=0;i<ids7.size();i++)
    {

    	cout <<" Triangle id = " << ids7[i] << endl;

    	cout <<" Node 1  " << endl;
    	cout<<   meshcoords (meshtopology (ids7[i],0) ,0) <<endl;
    	cout<<   meshcoords (meshtopology (ids7[i],0) ,1) <<endl;
    	cout<<   meshcoords (meshtopology (ids7[i],0) ,2) <<endl;

    	cout <<" Node 2  " << endl;
    	cout<<   meshcoords (meshtopology (ids7[i],1) ,0) <<endl;
    	cout<<   meshcoords (meshtopology (ids7[i],1) ,1) <<endl;
    	cout<<   meshcoords (meshtopology (ids7[i],1) ,2) <<endl;

    	cout <<" Node 3  " << endl;
    	cout<<   meshcoords (meshtopology (ids7[i],2) ,0) <<endl;
    	cout<<   meshcoords (meshtopology (ids7[i],2) ,1) <<endl;
    	cout<<   meshcoords (meshtopology (ids7[i],2) ,2) <<endl;
    }
    */


    const std::string name ( "Tunnel Problem " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
    gmesh->SetDimension ( 3 );

    TPZVec<REAL> coord ( 3 );

    vector<vector<double>> co;

    int ncoords = meshcoords.Rows();
    co.resize ( ncoords );
    for ( int i=0; i<ncoords; i++ ) {
        co[i].resize ( 3 );
        co[i][0]=meshcoords ( i,0 );
        co[i][1]=meshcoords ( i,1 );
        co[i][2]=meshcoords ( i,2 );
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
        coord[2] = co[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }


    int el=0;
    for ( int iel=0; iel<topol.size(); iel++ ) {
        if ( meshtopology(iel,3) ==-1 ) {

        } else {
            TPZVec <long> TopoTetra ( 4 );
            TopoTetra[0] =topol[iel][0];
            TopoTetra[1] =topol[iel][1];
            TopoTetra[2] =topol[iel][2];
            TopoTetra[3] =topol[iel][3];
            new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> ( el, TopoTetra, 1,*gmesh );
            el++;
        }
    }



    for ( int ivec=0; ivec<idsvec.size(); ivec++ )
    {
        int bcid;
        for(int jvec=0; jvec<idsvec[ivec].size(); jvec++)
        {
            TPZVec <long> TopoTri ( 3 );
            bcid=-ivec;
            int elindex = idsvec[ivec][jvec];
            TopoTri[0] = meshtopology (elindex,0);
            TopoTri[1] = meshtopology (elindex,1);
            TopoTri[2] = meshtopology (elindex,2);
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( el, TopoTri, - ( ivec+1 ), *gmesh );
            //cout << idsvec[ivec][jvec]  << " ";
            el++;
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
    //gmesh->Print(std::cout);
    std::ofstream files ( "ge.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}

void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    const auto &x=pt[0];
    const auto &y=pt[1];
    const auto &z=pt[2];
    REAL yy=40-y;
    REAL hw=2.;
    REAL H = 10.;
    REAL gamma =10.;
    REAL val=0.;
    if(yy<hw && yy<H)
    {
        //cout << " yy = "<< yy << endl;
        disp[0]  =   gamma * yy ;
    }
    else {
        disp[0]  = hw*gamma;
    }
    //cout << " disp[0] = "<< disp[0] << endl;



}
// void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
// 		const auto &x=pt[0];
//         const auto &y=pt[1];
//         const auto &z=pt[2];
//         REAL yy=40-y;
//         REAL hw=10.;
//         REAL H = 10.;
//         REAL gamma =10.;
//         REAL val=0.;
//         if(yy<(H-hw))
//         {
//             //cout << " yy = "<< yy << endl;
//             disp[0]  =   gamma * hw/(H-hw)* yy ;
//             val=disp[0];
//
//         }
//         else{
//             disp[0]  = hw*gamma;
//         }
//         //cout << " disp[0] = "<< disp[0] << endl;
//
//
//
// }

REAL findnodalsol(TPZCompMesh *cmesh) {

    //b[0] = 39.2265;
    //b[1] = 40.;
    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension();
    TPZVec<REAL> xd(dim,0.);
    TPZVec<REAL> mpt(dim,0.);
// xd[0] =39.2265; xd[1] = 40.;
    xd[0] =20.;
    xd[1] = 25.;
	xd[2]=40.;
    int id;
    int nels = gmesh->NElements();
    for(int iel=0; iel<nels; iel++) {


        TPZVec<REAL> qsi(dim,0.);
        TPZGeoEl * gel = gmesh->ElementVec()[iel];

        TPZCompEl * cel = gel->Reference();



        if(gel->MaterialId()<0)continue;
        bool check = gel->ComputeXInverse(xd,qsi,1.e-5);

        if(check==true)
        {


// 	   		cout << "elemento encontrado"<<endl;
// 			cout << "index do elemento = " <<gel->Index() <<endl;
            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(3);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                if(fabs(co[0]-xd[0])<1.e-3 && fabs(co[1]-xd[1])<1.e-3)
                {
// 					cout << "node id = "<<id <<endl;
// 					cout << " Coordinates "<< endl;
// 					cout << " x = "<< co[0] << ",  y = " << co[1] << endl;
// 					cout << " qsi = "<< qsi << endl;
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

                    TPZMaterialData data;

                    data.fNeedsSol = true;

                    intel->InitMaterialData ( data );

                    intel->ComputeRequiredData ( data, qsi );

// 					cout << " data.sol  = "<<data.sol     << endl;
                    return data.sol[0][1];
                }
            }


        }
    }

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
    analysis->DefineGraphMesh ( 3,scalarVars,vectorVars,vtk);
    constexpr int resolution{0};
    analysis->PostProcess ( resolution );
}


// void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
// 		const auto &x=pt[0];
//         const auto &y=pt[1];
//         const auto &z=pt[2];
//         REAL atm  = 10.33;//10.33 mca = 1 atm
//         disp[0] = /*kn/m^3*/  10*( 40-y )/* m */ ;/* = kn/m^2 = kPa*/
// }



void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    const auto &x=pt[0];
    const auto &y=pt[1];
    const auto &z=pt[2];
    REAL gamma=10;
    disp[0] = -(40-y)*gamma ;
}


TPZCompMesh * CreateCMeshDarcy( TPZGeoMesh *gmesh, int pOrder )
{
    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    cmesh->SetDefaultOrder ( pOrder );
    cmesh->SetDimModel ( 3 );
    cmesh->SetAllCreateFunctionsContinuous();

    // Create the material:
    // TPZDarcyFlow *material = new TPZDarcyFlow(m_matID,m_dim);
    int matid=1,dim=3;
    auto *material = new TPZDarcyFlowIsoPerm ( matid,dim );
    //Bet Degan loamy sand
    //STATE permeability = 0.0063 ;//cm/s
    //STATE permeability = 0.000063;//m/s
    //STATE permeability = 0.1;//m/s
    //STATE permeability = 1.;//m/s
    TPZManVector<REAL,3> permeability(3);
    permeability[0]=1.;
    permeability[1]=1.;
    permeability[2]=1.;
    material->SetConstantPermeability ( permeability );
    material->SetId(1);
    TPZAutoPointer<TPZFunction<STATE> > rhs = new TPZDummyFunction<STATE>(Forcing);
    //material->SetForcingFunction(rhs);

    cmesh->InsertMaterialObject ( material );

    // Condition of contours
    TPZFMatrix<STATE>  val1 ( 3,3,0. );
    TPZFMatrix<STATE>  val2 ( 3,1,0. );

    TPZAutoPointer<TPZFunction<STATE> > pressure = new TPZDummyFunction<STATE>(ForcingBCPressao);

    

    //material->SetForcingFunction(rhs);

//     REAL big = TPZMaterial::gBigNumber;
    
	
	//TPZMaterial * BCond3 = material->CreateBC ( material, -3, 1, val1, val2 );//tr
    //TPZMaterial * BCond1 = material->CreateBC ( material, 1, 0, val1, val2 );//tr
    //BCond0->SetForcingFunction(pressure);


	TPZMaterial * BCond1 = material->CreateBC ( material, -1, 0, val1, val2 );//tr
	TPZMaterial * BCond2 = material->CreateBC ( material, -2, 0, val1, val2 );//tr
	TPZMaterial * BCond3 = material->CreateBC ( material, -4, 0, val1, val2 );//tr
	TPZMaterial * BCond4 = material->CreateBC ( material, -5, 0, val1, val2 );//tr
	TPZMaterial * BCond7 = material->CreateBC ( material, -7, 0, val1, val2 );//tr
	
	BCond1->SetForcingFunction(rhs);
	BCond2->SetForcingFunction(rhs);
	BCond3->SetForcingFunction(rhs);
	BCond4->SetForcingFunction(rhs);
	
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);
	cmesh->InsertMaterialObject(BCond4);
	cmesh->InsertMaterialObject(BCond7);
    //cmesh->InsertMaterialObject(BCond1);
    //cmesh->InsertMaterialObject(BCond2);

    //Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();


    //LoadingRampRainFall ( cmesh,  1. );

    return cmesh;
}

void LoadingRampRainFall ( TPZCompMesh * cmesh,  REAL factor )
{
//      REAL big = TPZMaterial::gBigNumber;
// //
//     auto * bc1 = dynamic_cast<TPZBndCond *> (cmesh->FindMaterial(-3));
//     bc1->Val1()(0,0)=factor;

    auto * bc2 = dynamic_cast<TPZBndCond *> (cmesh->FindMaterial(-4));
    bc2->Val2()(0,0)=factor;

    auto * bc3 = dynamic_cast<TPZBndCond *> (cmesh->FindMaterial(-5));
    bc3->Val2()(0,0)=factor;



}


void ComputeSolution(TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphix,TPZSolVec &sol, TPZGradSolVec &dsol) {
    const int dim = cel->Reference()->Dimension();
    const int numdof = cel->Material()->NStateVariables();
    const int ncon = cel->NConnects();
    TPZFMatrix<STATE> &MeshSol = cel->Mesh()->Solution();


    long numbersol = MeshSol.Cols();
    sol.Resize(numbersol);
    dsol.Resize(numbersol);

    for (long is=0 ; is<numbersol; is++) {
        sol[is].Resize(numdof);
        sol[is].Fill(0.);
        dsol[is].Redim(dim, numdof);
        dsol[is].Zero();

    }

    TPZBlock<STATE> &block = cel->Mesh()->Block();
    long iv = 0, d;
    for(int in=0; in<ncon; in++) {
        TPZConnect *df = &cel->Connect(in);
        long dfseq = df->SequenceNumber();
        int dfvar = block.Size(dfseq);
        long pos = block.Position(dfseq);
        for(int jn=0; jn<dfvar; jn++) {
            for (int is=0; is<numbersol; is++) {

                sol[is][iv%numdof] += (STATE)phi(iv/numdof,0)*MeshSol(pos+jn,is);

                for(d=0; d<dim; d++) {
                    dsol[is](d,iv%numdof) += (STATE)dphix(d,iv/numdof)*MeshSol(pos+jn,is);
                }
            }
            iv++;
        }
    }
    //  cout << "------ end -----  " << endl;
}//metho


void Post ( TPZPostProcAnalysis * postproc,std::string vtkFile )
{

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    postproc->DefineGraphMesh ( 3,scalNames,vecNames,vtkFile );

    postproc->PostProcess ( 0 );
}




/// create the computational mesh
TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder )
{
    unsigned int dim  = 3;
    const std::string name ( "ElastoPlastic Tunnel " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    REAL mc_cohesion    = 20;//kpa
    REAL mc_phi         = ( 30.*M_PI/180 );
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000.;

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetUp ( E, nu );
    LEMC.fER =ER;
 
    LEMC.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int id =1;

    plasticmat3d * material = new plasticmat3d ( id );
    material->SetPlasticity ( LEMC );

    material->SetId ( id );
    cmesh->InsertMaterialObject ( material );

    TPZFMatrix<STATE> val1 ( 3,3,0. );
    TPZFMatrix<STATE>  val2 ( 3,1,0. );
    int directionadirichlet =3,newmann=1;
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    val2 ( 2,0 ) = 0;
    auto * bc1 = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );//plano yz x=0
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    val2 ( 2,0 ) = 0;
    auto * bc2 = material->CreateBC ( material, -2, directionadirichlet, val1, val2 );//plano yz x=50
    val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = 1;
    val2 ( 2,0 ) = 0;
    auto * bc3 = material->CreateBC ( material, -3, directionadirichlet, val1, val2 );//plano xz y=0
    val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = 0;
    val2 ( 2,0 ) = 1;
    auto * bc5 = material->CreateBC ( material, -5, directionadirichlet, val1, val2 );//plano xy z=0
    val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = 0;
    val2 ( 2,0 ) = 1;
    auto * bc6 = material->CreateBC ( material, -6, directionadirichlet, val1, val2 );//plano xy z=40
	
	val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = -10;
    val2 ( 2,0 ) = 0;
    auto * bc4 = material->CreateBC ( material, -4, newmann, val1, val2 );//plano xy z=40
	
	val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = 0;
    val2 ( 2,0 ) = 0;
    auto * bc7 = material->CreateBC ( material, -7, newmann, val1, val2 );//plano xy z=40

    cmesh->InsertMaterialObject ( bc1 );
    cmesh->InsertMaterialObject ( bc2 );
    cmesh->InsertMaterialObject ( bc3 );
    cmesh->InsertMaterialObject ( bc4 );
    cmesh->InsertMaterialObject ( bc5 );
	cmesh->InsertMaterialObject ( bc6 );
	cmesh->InsertMaterialObject ( bc7 );

    //cmesh->InsertMaterialObject ( top );
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    return cmesh;

}

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor )
{
    plasticmat3d * body= dynamic_cast<plasticmat3d *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=-factor*20.;
    body->SetBodyForce ( force );

}




void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh )
{

    if ( PostProcess->ReferenceCompMesh() != cmesh ) {

        PostProcess->SetCompMesh ( cmesh );

        TPZVec<int> PostProcMatIds ( 1,1 );
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        PostProcessVariables ( scalNames, vecNames );

        for ( int i=0; i<scalNames.size(); i++ ) {
            PostProcVars.Push ( scalNames[i] );
        }
        for ( int i=0; i<vecNames.size(); i++ ) {
            PostProcVars.Push ( vecNames[i] );
        }
        //
        PostProcess->SetPostProcessVariables ( PostProcMatIds, PostProcVars );
        TPZFStructMatrix structmatrix ( PostProcess->Mesh() );
        structmatrix.SetNumThreads ( 0 );
        PostProcess->SetStructuralMatrix ( structmatrix );
    }
    //
    //Chamar com o analysis e nao com o postanalysis pois tem o acumulo de sols
    PostProcess->TransferSolution();

    //fElastoplasticAnalysis.TransferSolution(fPostprocess);
}


/// Get the post processing variables
void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{
    
  scalNames.Push ("StrainVol");
  scalNames.Push ("StrainXX");
  scalNames.Push ("StrainYY");
  scalNames.Push ("StrainZZ");
  scalNames.Push ("StrainXY");
  scalNames.Push ("StrainXZ");
  scalNames.Push ("StrainYZ");
  scalNames.Push ("PlStrainSqJ2");
  scalNames.Push ("Pressure");
  vecNames.Push ( "DisplacementTotal" );
  vecNames.Push ( "Flux" );
}

void SetSuportPressure ( TPZCompMesh * cmesh,  REAL factor )
{
    auto * matbc= dynamic_cast<TPZBndCond *> ( cmesh->FindMaterial ( -7 ) );
	TPZFMatrix<STATE> val1 ( 3,3,0. );
    TPZFMatrix<STATE>  val2 ( 3,1,0. );

//bc.Val2()(0,0);
	REAL sigma0=-10;//(kPa)
	matbc->Val2()(0,0) = sigma0 * factor;
	

}
void GravityIncrease ( TPZCompMesh * cmesh )
{

    REAL FS=0.1,FSmax=1000.,FSmin=0.,tol=0.01;
    int neq = cmesh->NEquations();
    int maxcount=10;
    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;
	//LoadingRamp ( cmesh,0.1 );
    REAL norm = 1000.;
    REAL tol2 = 1.e-2;
    int NumIter = 10;
    bool linesearch = true;
    bool checkconv = false;

    do {

        std::cout << "FS = " << FS  <<" | Load step = " << counterout << " | Rhs norm = " << norm  << std::endl;
        //LoadingRamp ( cmesh,FS );
		//SetSuportPressure(cmesh,FS);
        bool optimize =true;
        TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
        chrono::steady_clock sc;
        auto start = sc.now();
        int iters;
        bool conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters);
		//bool conv =anal->IterativeProcess(cout, tol2, NumIter,linesearch,checkconv);
        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << "| total time in iterative process =  " << time_span.count()<< std::endl;
        //anal->IterativeProcess ( outnewton, tol2, NumIter);

        norm = Norm ( anal->Rhs() );

        if ( conv==false) {
            cmesh->LoadSolution(displace0);
            //cmesh->Solution().Zero();
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;

        } else {
            // uy+=findnodalsol(cmesh);
            displace0 = anal->Solution();
            FSmin = FS;
            anal->AcceptSolution();
            FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }
        // cout << "|asdadadasd =  " << std::endl;
        counterout++;

    }  while ( (( FSmax - FSmin ) / FS > tol && counterout<maxcount) );
    TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,true );
    anal->AcceptSolution();
}
void ShearRed ( TPZCompMesh * cmesh )
{
    //plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    LoadingRamp ( cmesh,1. );

    REAL FS=0.1,FSmax=10000.,FSmin=0.,tol=1.e-3;
    int neq = cmesh->NEquations();

    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    REAL norm = 1000.;
    REAL tol2 = 1.e-3;
    int NumIter = 30;
    bool linesearch = true;
    bool checkconv = false;
    do {

        bool optimize =true;
        TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );

        anal->IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
        norm = Norm ( anal->Rhs() );

        //std::cout << "| Load step = " << counterout << " | Rhs norm = " << norm << "| delta displacement norm = "<< unorm << std::endl;
        if ( norm>= tol2 ) {
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
    bool optimize =false;
    TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
    anal->IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
    anal->AcceptSolution();
}
/*
void ShearRed ( TPZCompMesh * cmesh)
{
    LoadingRamp ( cmesh,1. );

    REAL FS=0.5,FSmax=10000.,FSmin=0.,tol=1.e-2;
    int neq = cmesh->NEquations();
    int maxcount=10;
    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    REAL norm = 1000.;
    REAL tol2 = 1.;
    int NumIter = 20;
    bool linesearch = true;
    bool checkconv = false;
	std::ofstream outnewton("saida-newton.txt");
		std::ofstream outloadu("srmloadvsu.nb");
	REAL uy=0.;
	outloadu << "plot = {";
    do {

        bool optimize =true;
        TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );

        int iters;
        anal->IterativeProcess ( outnewton, tol2, NumIter,linesearch,checkconv,iters );
        norm = Norm ( anal->Rhs() );

        std::cout << "FS = " << FS << " | Load step = " << counterout << " | Rhs norm = " << norm  << std::endl;
        if ( norm>= tol2 ) {
            displace = displace0;
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {
			uy+=findnodalsol(cmesh);
			outloadu << "{ "<<-uy << ", " << FS << " } ," << endl;
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
    }  while ( (( FSmax - FSmin ) / FS > tol && counterout<maxcount) || norm>tol2);
	outloadu <<  " }; ListLinePlot[plot,PlotRange->All]" << endl;
	bool optimize =false;
	TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
    anal->IterativeProcess ( outnewton, tol2, NumIter,linesearch,checkconv );
    anal->AcceptSolution();
}*/
#include "TPZFrontMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "tpzsparseblockdiagonalstructmatrix.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
//#include "TPZStepSolver.h"
#include "bicgstab.h"
//#include <bicgstab.h>
TPZElastoPlasticAnalysis * CreateAnal ( TPZCompMesh *cmesh,bool optimize )
{
    int numthreads=8;

    TPZElastoPlasticAnalysis * analysis =  new TPZElastoPlasticAnalysis ( cmesh ); // Create analysis

//     TPZSkylineStructMatrix matskl ( cmesh );
//     // TPZFStructMatrix matskl ( cmesh );
// 
//     matskl.SetNumThreads ( numthreads );
// 
//     analysis->SetStructuralMatrix ( matskl );
// 
//     ///Setting a direct solver
//     TPZStepSolver<STATE> step;


	//void SetBiCGStab(const long numiterations, const TPZMatrixSolver<TVar> &pre, const REAL tol, const long FromCurrent);
	//step.SetBiCGStab(100,matskl,1.e-3,1);
	
	//step.SetBiCGStab();
	
	//step.SetDirect ( EBICGSTAB );
	//step.SetDirect ( E );
    //step.SetDirect ( ELDLt );
    //step.SetDirect ( ENoDecompose );
    // step.SetDirect ( ECholesky );
    //analysis->SetSolver ( step );

	analysis->SetBiCGStab(200, 1.e-3);
//     long neq = cmesh->NEquations();
//     TPZVec<long> activeEquations;
//     analysis->GetActiveEquations ( activeEquations );
//     TPZEquationFilter filter ( neq );
//     filter.SetActiveEquations ( activeEquations );
//     matskl.EquationFilter() = filter;
//     analysis->SetStructuralMatrix ( matskl );

    //step.SetDirect(ECholesky);
    //analysis->SetSolver ( step );


    return analysis;
}
void SetFlux ( TPZCompMesh * plasticCmesh,TPZCompMesh* incmesh)
{

    int nels0 = incmesh->NElements();
    //int nels0 = plasticCmesh->NElements();

    //incmesh->Solution().Print(std::cout);

    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( plasticCmesh->MaterialVec() [1] );
    TPZAdmChunkVector<TPZElastoPlasticMem>  &mem = pMatWithMem2->GetMemory();

    //cout << "mem.NElements() "<< mem.NElements()  <<endl;
    //cout << " nels0 "<< nels0  <<endl;
    if ( pMatWithMem2 ) {
        pMatWithMem2->SetUpdateMem ( true );
    }

    int globpt=0;
    for ( int iel=0; iel<nels0; iel++ ) {

        //cout << "\n out1.. " << endl;
        TPZCompEl *celplastic = plasticCmesh->Element ( iel );
        TPZInterpolationSpace *intelplastic = dynamic_cast<TPZInterpolationSpace *> ( celplastic );


        TPZGeoEl *gelplastic = celplastic->Reference();

        TPZCompEl *celdarcy = incmesh->Element ( iel );
        // cout << "\n out2.. " << endl;
        //if(!celdarcy)continue;

        TPZGeoEl *geldarcy = celdarcy->Reference();
        //if(!geldarcy)continue;

        int indexdarcy = geldarcy->Index();

        int indexplastic = gelplastic->Index();


        TPZInterpolationSpace *intelrandom1 = dynamic_cast<TPZInterpolationSpace *> ( celdarcy );


        const TPZIntPoints &intpoints = intelplastic->GetIntegrationRule();




        TPZManVector<REAL,3> point ( 3,0. );

        TPZMaterialData data,data1;

        data1.fNeedsSol = true;
        intelplastic->InitMaterialData ( data );
        intelrandom1->InitMaterialData ( data1 );

        REAL weight=0;
        int nint = intpoints.NPoints();


        if(celplastic->Material()->Id()!=1 || celdarcy->Material()->Id()!=1)
        {
            // cout << "\n Boundary El "<<endl;
            // cout << "\n celplastic->Material()->Id() = "<< celplastic->Material()->Id()<<endl;
            // cout << "\n celdarcy->Material()->Id() = "<< celdarcy->Material()->Id()<<endl;
            continue;
        }

        if(indexdarcy!=indexplastic)
        {
            cout << "\n indexdarcy = "<< indexdarcy<<endl;
            cout << "\n indexplastic = "<< indexplastic<<endl;
            DebugStop();
        }

        if(celplastic->Material()->Id() != celdarcy->Material()->Id())
        {
            cout << "\n Different material IDs "<<endl;
            cout << "\n celplastic->Material()->Id() = "<< celplastic->Material()->Id()<<endl;
            cout << "\n celdarcy->Material()->Id() = "<< celdarcy->Material()->Id()<<endl;
            DebugStop();
        }


        TPZTensor<REAL> epst,epsp;
        for ( long ip =0; ip<nint; ip++ ) {

            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;

            intelplastic->ComputeRequiredData ( data, point );
            intelrandom1->ComputeRequiredData ( data1, point );

            TPZMaterial *mat= celdarcy->Material();

            TPZVec<REAL> flux,pressure;
            int varid =7;//flux
            mat->Solution(data1,varid,flux);

            varid =1;//pressure
            mat->Solution(data1,varid,pressure);

            TPZFMatrix<REAL> flow = data1.dsol[0];

            if(fabs(data1.detjac-data.detjac)>1.e-3)
            {
                cout << "\n Different det jacs "<<endl;
                cout << "\n data.detjac = "<< data.detjac<<endl;
                cout << "\n data1.detjac = "<< data1.detjac<<endl;

                DebugStop();
            }
            if(fabs(data1.x[0]-data.x[0])>1.e-3 ||fabs(data1.x[1]-data.x[1])>1.e-3)
            {
                cout << "\n Different x "<<endl;
                cout << "\n data.x = "<< data.x<<endl;
                cout << "\n data1.x = "<< data1.x<<endl;

                DebugStop();
            }


            data.intGlobPtIndex = globpt;

            mem[globpt].fPlasticState.fflux.Resize ( 3 );
            //mem[globpt].fPlasticState.fflux[0]=data1.dsol[0][0];
            //mem[globpt].fPlasticState.fflux[1]=data1.dsol[0][1];
            mem[globpt].fPlasticState.fflux[0]=flux[0];
            mem[globpt].fPlasticState.fflux[1]=flux[1];
			mem[globpt].fPlasticState.fflux[2]=flux[2];
            mem[globpt].fPlasticState.fpressure=pressure[0];

            globpt++;

        }


        //pMatWithMem2->SetUpdateMem ( false );
    }

    pMatWithMem2->SetUpdateMem ( false );

    plasticCmesh->Solution().Zero();


}

template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs )
{
    std::vector<T> ret;
    for ( std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it ) {
        istringstream iss ( *it );
        T temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out )
{
    string line,line2, temp;

    ifstream myfile ( file );

    std::vector<std::vector<double>> coords;
    int counter=0;
    while ( getline ( myfile, line ) ) {
        std::vector<string> tokens;
        istringstream iss ( line );
        while ( iss >> temp ) tokens.push_back ( temp );
        std::vector<double> input_doub_temp= str_vec<double> ( tokens );
        coords.push_back ( input_doub_temp );
        counter++;
    }


// 	cout << "coords.size() = "<< coords.size();
// 	cout << "coords.size() = "<< coords[1].size();
// 	cout<< endl;
// 	for(int iel=0;iel<coords.size();iel++){
// 		for(int elnode=0;elnode<coords[iel].size();elnode++)
// 		{
// 			cout<< coords[iel][elnode] << " ";
// 		}
// 		cout<< endl;
// 	}
//
//
    out.Resize ( coords.size(),coords[1].size() );
    for ( int iel=0; iel<coords.size(); iel++ ) {
        for ( int elnode=0; elnode<coords[iel].size(); elnode++ ) {
            out ( iel,elnode ) = coords[iel][elnode];
        }
    }


// 	int els = coords.size();
// 	int elnodes = coords[0].size();
// 	out.Resize(els,elnodes);
// 	for(int iel=0;iel<els;iel++){
// 		for(int elnode=0;elnode<elnodes;elnode++)
// 		{
// 			out(iel,elnode)=coords[iel][elnode];
// 		}
// 	}
}

void PrintMat ( std::string out,TPZFMatrix<REAL> mat )
{
    std::ofstream print ( out );
    int row=mat.Rows();
    int cols = mat.Cols();

    for ( int i=0; i<row; i++ ) {
        for ( int j=0; j<cols; j++ ) {
            print << mat ( i,j ) << " ";
        }
        print<< std::endl;
    }


}


