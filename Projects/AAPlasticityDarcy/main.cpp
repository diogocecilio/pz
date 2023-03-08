
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
#include "TPZMohrCoulombVoigt.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplastic.h"
#include "TPZMohrCoulombVoigt.h"
#include "TPZPlasticStepVoigt.h"
//#include "TPZDarcyFlow.h"
#include "TPZDarcyFlowIsoPerm.h"

using namespace std;

//typedef TPZPlasticStepVoigt<TPZMohrCoulombVoigt,TPZElasticResponse> LEMC;
typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;


typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;



#include "readgidmesh.h"
#include <TPZParFrontStructMatrix.h>

TPZGeoMesh * CreateGMesh ( int ref );

TPZGeoMesh * CreateSimpleGMesh (  );

TPZGeoMesh * CreateGMeshGid ( int ref );

TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder );

TPZCompMesh * CreateSimpleCMesh ( TPZGeoMesh * gmesh,int porder );

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor );


void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void Post ( TPZPostProcAnalysis * postproc,std::string vtkFile );

TPZElastoPlasticAnalysis * CreateAnal ( TPZCompMesh *cmesh,bool optimize );

TPZElastoPlasticAnalysis * CreateSimpleAnal ( TPZCompMesh *cmesh,bool optimize );

void ShearRed ( TPZCompMesh * cmesh );

void GravityIncrease ( TPZCompMesh * cmesh );

REAL ShearRed ( TPZCompMesh * cmesh,int isample,TPZCompMesh * incmesh );

void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphix,TPZSolVec &sol, TPZGradSolVec &dsol );


void SetFlux ( TPZCompMesh * plasticCmesh,TPZCompMesh* incmesh);

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );

void PrintMat ( std::string out,TPZFMatrix<REAL> mat );

TPZManVector<TPZCompMesh *,2> CreateFieldsDummy ( TPZGeoMesh * gmesh,int porder );

void SolveRamp(TPZCompMesh * cmesh);
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

int betax;
REAL gammasolo;
REAL gammaaugua;
bool water;
REAL coesaofiu;
REAL phifiu;

int main()
{


//     int porder=1;
//     TPZGeoMesh *gmesh0 = CreateSimpleGMesh (  );
//     TPZCompMesh *cmesh0 = CreateSimpleCMesh ( gmesh0,porder );
//     TPZElastoPlasticAnalysis *anal0 =  CreateSimpleAnal ( cmesh0,false );
// //
//      LoadingRamp ( cmesh0,1. );
//     std::ofstream outnewton("saida-newton.txt");
// //
// //    int neq = cmesh0->NEquations();
//     REAL tol = 0.001;
//     int numiter =5;
//     REAL tol2 = 0.001;
//     int numiter2 =5;
//     REAL l =1.;
// //
//      bool linesearch = true;
//     bool checkconv=false;
//     int iters;
//     cout << "a " << endl;
//      anal0->IterativeProcessArcLength(outnewton,tol,numiter,tol2,numiter2,l,linesearch);
// //     //anal0->IterativeProcess(outnewton, tol, numiter,  linesearch,  checkconv,iters);
// //     //GravityIncrease ( cmesh0 );
// //
// //
//     return 0;


    int porder=2;
    int nref=0;

    betax=45;//30,45,60,90
    gammasolo=20.;
    gammaaugua=10.;
    water=true;
    coesaofiu=10.;
    phifiu=30.*M_PI/180;

    REAL tolfs = 1.e-2;
    int numiterfs =20;
    REAL tolres = 1.e-3;
    int numiterres =20;
    REAL l =2.;
    REAL lambda0=0.1;

    TPZGeoMesh *gmesh = CreateGMeshGid ( nref);


    string vtk = "Darcyx.vtk";


    TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    TPZPostProcAnalysis * postprocdetergim = new TPZPostProcAnalysis();
    std::string vtkFiled ="vtkfolder/deterministicflux.vtk";
    std::string vtkFiledgim ="vtkfolder/outvtk.vtk";


    TPZCompMesh *  darcycompmesh =  CreateCMeshDarcy(gmesh,porder);
    TPZCompMesh *cmeshsrm = CreateCMesh ( gmesh,porder );
    TPZCompMesh *cmeshgi = CreateCMesh ( gmesh,porder );

    if(water==true)
    {
        cout << "\n Solving Darcy... " << endl;
        SolveDarcyProlem(darcycompmesh, vtk);

        cout << "\n Setting flux in mechanic comp mesh... " << endl;
        //SetFlux(cmeshsrm,darcycompmesh);
        SetFlux(cmeshgi,darcycompmesh);
    }

    cout << "\n Gravity Increase routine.. " << endl;
    TPZElastoPlasticAnalysis  * anal0 = CreateAnal ( cmeshgi,true );
    //anal0->SetBiCGStab(30,10e-6);
    //anal0->SetBiCGStab_Jacobi(30,10e-6);
    LoadingRamp ( cmeshgi,1. );
    std::ofstream outnewton("saida-newton.txt");

    chrono::steady_clock sc;
    auto start = sc.now();
    anal0->IterativeProcessArcLength(tolfs,numiterfs,tolres,numiterres,l,lambda0);
    //GravityIncrease ( cmeshgi );
    auto end = sc.now();
    auto time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "|  IterativeProcessArcLength time =  " << time_span.count()<< std::endl;


    cout << "\n Post Processing gravity increase... " << endl;
    CreatePostProcessingMesh ( postprocdetergim, cmeshgi );
    Post ( postprocdetergim,vtkFiledgim );

    return 0;
}


void SolveRamp(TPZCompMesh * cmesh)
{
    int maxiter=100;
    REAL tol=1.;
    bool linesearch=true;
    bool checkconv=false;
    int steps=40;
    REAL finalload = 5.;
    std::string vtkFile ="tcentrifuga-p2.vtk";
    std::ofstream outloadu("loadvup2.nb");
    outloadu << "plot = {";
    TPZPostProcAnalysis  * postproc = new TPZPostProcAnalysis();
    REAL g0=1.;
    vector<double> loadvec = {0.1,0.15,0.2,0.3,0.4,0.5};

    for(int iter=0;iter<loadvec.size();iter++)loadvec[iter]*=g0;

    findnodalsol(cmesh);

    for(int iloadstep=0; iloadstep<100; iloadstep++)
    {

        TPZElastoPlasticAnalysis  * analysis =  CreateAnal(cmesh,true);
        //REAL load = iloadstep*finalload/steps;
        REAL load = loadvec[iloadstep];
        //REAL load = iloadstep;
        cout << " \n --------- iloadstep  = "<< iloadstep  << " |  load  = "<< load << " gtrial/g0 = "<< load/g0 <<endl;
        LoadingRamp(cmesh,load);
        int iters;
        bool conv = analysis->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv,iters);
        if(conv==false)
        {
            cout << "newton failed to converge!";
            break;
        }

        analysis->AcceptSolution();

    }


}

TPZElastoPlasticAnalysis * CreateSimpleAnal ( TPZCompMesh *cmesh,bool optimize )
{
    int numthreads=0;


    TPZElastoPlasticAnalysis * analysis =  new TPZElastoPlasticAnalysis ( cmesh ); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );
   // TPZFStructMatrix matskl ( cmesh );

    //matskl.SetNumThreads ( numthreads );

    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
	//step.SetDirect ( ENoDecompose );
   // step.SetDirect ( ECholesky );
     //step.SetDirect(ECholesky);

    analysis->SetSolver ( step );



    return analysis;
}

TPZCompMesh * CreateSimpleCMesh ( TPZGeoMesh * gmesh,int porder )
{
    unsigned int dim  = 2;
    const std::string name ( " simple cmesh " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    REAL mc_phi         = ( 30.*M_PI/180 );
    REAL mc_psi         = mc_phi;
    REAL mc_cohesion = 10.;

    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000.;

    //TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse> LEMC;
    LEMC  plasticstep;
   //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetUp ( E, nu );
    plasticstep.fER =ER;

    plasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int PlaneStrain = 1;

    //plasticmatcrisfield * material = new plasticmatcrisfield ( 1,PlaneStrain );
    plasticmat * material = new plasticmat ( 1,PlaneStrain );

    material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

    material->SetLoadFactor(0.);
    material->SetWhichLoadVector(0);//option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    cmesh->InsertMaterialObject ( material );

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
    int directionadirichlet =3;
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 1;
    auto * bc_bottom = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );//bottom
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_rigth = material->CreateBC ( material, -2, directionadirichlet, val1, val2 );//rigth

    cmesh->InsertMaterialObject ( bc_bottom );
    cmesh->InsertMaterialObject ( bc_rigth );

    cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    return cmesh;
}
TPZGeoMesh * CreateSimpleGMesh ( )
{
    const std::string name ( "Simple GMesh " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

    vector<vector<double>> co;

    int ncoords = 4;
    int nels =1;
    int bclines=4;

    double coords[ncoords][2]={{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}};

    long topol[nels][4]={{0,1,2,3}};

   long bclinetopol[bclines][2]= {{0,1},{1,2},{2,3},{3,0}};

   long boundaryids[bclines]={0, 1, 2, 3};

     gmesh->NodeVec().Resize(ncoords);
    for(int inode=0;inode<ncoords;inode++)
    {
        TPZVec<double> co(2);
        co[0]=coords[inode][0];co[1]=coords[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, co, *gmesh );
    }


    long id=0;
    for ( int iel=0; iel<nels; iel++ ) {
        TPZVec<long> top(4);
        top[0]=topol[iel][0];top[1]=topol[iel][1];top[2]=topol[iel][2];top[3]=topol[iel][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( id,top, 1,*gmesh );
        id++;
    }

    TPZVec <long> TopoLine ( 2 );
    id++;
    for ( int ivec=0; ivec<bclines; ivec++ )
    {
        if(boundaryids[ivec]==0)
        {
            TopoLine[0] = bclinetopol[ivec][0];
            TopoLine[1] = bclinetopol[ivec][1];
            cout<< "TopoLine[0] = " << TopoLine[0] <<endl;
            cout<< "TopoLine[1] = " << TopoLine[1] <<endl;
            cout<< " sdsad"<<endl;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 1, *gmesh );//bottom
            id++;
        }
        if(boundaryids[ivec]==1)
        {
            TopoLine[0] = bclinetopol[ivec][0];
            TopoLine[1] = bclinetopol[ivec][1];
            cout<< "TopoLine[0] = " << TopoLine[0] <<endl;
            cout<< "TopoLine[1] = " << TopoLine[1] <<endl;
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 2, *gmesh );//right
            id++;
        }
        if(boundaryids[ivec]==2)
        {
            TopoLine[0] = bclinetopol[ivec][0];
            TopoLine[1] = bclinetopol[ivec][1];
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 3, *gmesh );//top
            id++;
        }
        if(boundaryids[ivec]==3)
        {
            TopoLine[0] = bclinetopol[ivec][0];
            TopoLine[1] = bclinetopol[ivec][1];
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 4, *gmesh );//left
            id++;
        }

    }


    gmesh->BuildConnectivity();

        for ( int d = 0; d<5; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    std::ofstream outmesh("gmesh.txt");
    gmesh->Print(outmesh);
    std::ofstream files ( "gesimple.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}
// TPZGeoMesh * CreateSimpleGMesh ( )
// {
//     const std::string name ( "Simple GMesh " );
//
//     TPZGeoMesh *gmesh  =  new TPZGeoMesh();
//
//     gmesh->SetName ( name );
//     gmesh->SetDimension ( 2 );
//
//     TPZVec<REAL> coord ( 2 );
//
//     vector<vector<double>> co;
//
//     int ncoords = 13;
//     int nels =16;
//     int bclines=8;
//
//     double coords[ncoords][2]={{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}, {0.5, 0.5}, {0., 0.5},
//     {0.5, 0.}, {0.5, 1.},{0.25, 0.25},{1., 0.5}, {0.25, 0.75}, {0.75, 0.75}, {0.75, 0.25}};
//
//     long topol[nels][3]={{4, 12, 9}, {8, 4, 5}, {10, 7, 3}, {3, 5, 10}, {12, 4, 6}, {4, 8,
//   6}, {0, 6, 8}, {11, 2, 7}, {5, 0, 8}, {9, 2, 11}, {4, 7, 10}, {5, 4,
//    10}, {4, 9, 11}, {7, 4, 11}, {1, 9, 12}, {6, 1, 12}};
//
//    long bclinetopol[bclines][2]= {{1, 6}, {2, 9}, {3, 7}, {0, 5}, {5, 3}, {6, 0}, {7, 2}, {9, 1}};
//
//    long boundaryids[bclines]={0, 1, 2, 3, 3, 0, 2, 1};
//
//      gmesh->NodeVec().Resize(ncoords);
//     for(int inode=0;inode<ncoords;inode++)
//     {
//         TPZVec<double> co(2);
//         co[0]=coords[inode][0];co[1]=coords[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, co, *gmesh );
//     }
//
//
//     long id=0;
//     for ( int iel=0; iel<nels; iel++ ) {
//         TPZVec<long> top(3);
//         top[0]=topol[iel][0];top[1]=topol[iel][1];top[2]=topol[iel][2];
//         new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( id,top, 1,*gmesh );
//         id++;
//     }
//
//     TPZVec <long> TopoLine ( 2 );
//     id++;
//     for ( int ivec=0; ivec<bclines; ivec++ )
//     {
//         if(boundaryids[ivec]==0)
//         {
//             TopoLine[0] = bclinetopol[ivec][0];
//             TopoLine[1] = bclinetopol[ivec][1];
//             cout<< "TopoLine[0] = " << TopoLine[0] <<endl;
//             cout<< "TopoLine[1] = " << TopoLine[1] <<endl;
//             cout<< " sdsad"<<endl;
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 1, *gmesh );//bottom
//             id++;
//         }
//         if(boundaryids[ivec]==1)
//         {
//             TopoLine[0] = bclinetopol[ivec][0];
//             TopoLine[1] = bclinetopol[ivec][1];
//             cout<< "TopoLine[0] = " << TopoLine[0] <<endl;
//             cout<< "TopoLine[1] = " << TopoLine[1] <<endl;
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 2, *gmesh );//right
//             id++;
//         }
//         if(boundaryids[ivec]==2)
//         {
//             TopoLine[0] = bclinetopol[ivec][0];
//             TopoLine[1] = bclinetopol[ivec][1];
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 3, *gmesh );//top
//             id++;
//         }
//         if(boundaryids[ivec]==3)
//         {
//             TopoLine[0] = bclinetopol[ivec][0];
//             TopoLine[1] = bclinetopol[ivec][1];
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - 4, *gmesh );//left
//             id++;
//         }
//
//     }
//
//
//     gmesh->BuildConnectivity();
//
//         for ( int d = 0; d<0; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//
//     std::ofstream outmesh("gmesh.txt");
//     gmesh->Print(outmesh);
//     std::ofstream files ( "gesimple.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
//
//     return gmesh;
// }
void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        REAL yy=40-y;
        REAL hw=10.;
        REAL H = 10.;

        REAL val=0.;
        if(yy<hw && yy<H)
        {
            //cout << " yy = "<< yy << endl;
            disp[0]  =   gammaaugua * yy ;
        }
        else{
            disp[0]  = hw*gammaaugua;
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

REAL findnodalsol(TPZCompMesh *cmesh){

    //b[0] = 39.2265;
    //b[1] = 40.;
TPZGeoMesh * gmesh =  cmesh->Reference();
	int dim = gmesh->Dimension();
 TPZVec<REAL> xd(dim,0.);
 TPZVec<REAL> mpt(dim,0.);
// xd[0] =39.2265; xd[1] = 40.;
 xd[0] =35.; xd[1] = 40.;
int id;
 int nels = gmesh->NElements();
 for(int iel=0;iel<nels;iel++){


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
			for(int inode=0;inode<nnodes;inode++)
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
    analysis->DefineGraphMesh ( 2,scalarVars,vectorVars,vtk);
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



void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        disp[0] = -1;
}


TPZCompMesh * CreateCMeshDarcy( TPZGeoMesh *gmesh, int pOrder )
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

    permeability[0]=1.;permeability[1]=1.;permeability[2]=1.;
    material->SetConstantPermeability ( permeability );
	material->SetId(1);

    cmesh->InsertMaterialObject ( material );

    // Condition of contours
    TPZFMatrix<STATE>  val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );

	TPZAutoPointer<TPZFunction<STATE> > pressure = new TPZDummyFunction<STATE>(ForcingBCPressao);

	TPZAutoPointer<TPZFunction<STATE> > rhs = new TPZDummyFunction<STATE>(Forcing);

    //material->SetForcingFunction(rhs);

//     REAL big = TPZMaterial::gBigNumber;
     TPZMaterial * BCond0 = material->CreateBC ( material, -3, 0, val1, val2 );//tr
     BCond0->SetForcingFunction(pressure);

     TPZMaterial * BCond1 = material->CreateBC ( material, -4, 0, val1, val2 );//ramp
     BCond1->SetForcingFunction(pressure);

	//val2(0,0)=-1;
	TPZMaterial * BCond2 = material->CreateBC ( material, -5, 0, val1, val2 );//tl
	BCond2->SetForcingFunction(pressure);





    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);

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


void ComputeSolution(TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphix,TPZSolVec &sol, TPZGradSolVec &dsol){
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

                for(d=0; d<dim; d++){
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

    postproc->DefineGraphMesh ( 2,scalNames,vecNames,vtkFile );

    postproc->PostProcess ( 0 );
}


TPZGeoMesh * CreateGMeshGid ( int ref )
{

    string file;



 //   file ="/home/diogo/projects/pz/data/quad-gid2.msh";


    if(betax==30)
    {
         file ="/home/diogo/projects/pz/data/h10-beta30.msh";
    }
    if(betax==45)
    {
         file ="/home/diogo/projects/pz/data/h10-beta45.msh";
    }
    if(betax==60)
    {
         file ="/home/diogo/projects/pz/data/h10-beta60.msh";
    }
    if(betax==90)
    {
         file ="/home/diogo/projects/pz/data/h10-beta90.msh";
    }



    readgidmesh read = readgidmesh ( file );
    read.ReadMesh();
    TPZFMatrix<int> meshtopology = read.GetTopology();
    TPZFMatrix<REAL> meshcoords = read.GetCoords();
    std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();

    int ndivs = 10000;
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
    if(betax==30 || betax==45 || betax ==60)
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



    if(betax==30)
    {
        a = b;
        b[0] = 27.675;
        b[1] = 40.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idsramp );
        idsvec.push_back ( idsramp );
    }
    if(betax==45)
    {
        a = b;
        b[0] = 35.;
        b[1] = 40.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idsramp );
        idsvec.push_back ( idsramp );
    }

    if(betax==60)
    {
        a = b;
        b[0] = 39.2265;
        b[1] = 40.;
        read.Line ( a, b, ndivs, pathbottom );
        read.FindIdsInPath ( pathbottom, idsramp );
        idsvec.push_back ( idsramp );
    }

    if(betax==90)
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

/// create the computational mesh
TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder )
{
    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    // Mohr Coulomb data
   // REAL mc_cohesion    = 68;//kpa
    //REAL mc_phi         = ( 50.0*M_PI/180 );
	REAL mc_cohesion    = coesaofiu;//kpa
    REAL mc_phi         = phifiu;
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000.;

    //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    //TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse> LEMC;
    LEMC plasticstep;
    ER.SetUp ( E, nu );
    plasticstep.fER =ER;
    // LEMC.SetElasticResponse( ER );
    plasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );
    //plasticmatcrisfield * material = new plasticmatcrisfield ( 1,PlaneStrain );
    material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

    material->SetLoadFactor(1.);
    material->SetWhichLoadVector(0);//option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    cmesh->InsertMaterialObject ( material );

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
    int directionadirichlet =3;
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 1;
    auto * bc_bottom = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );//bottom
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_rigth = material->CreateBC ( material, -2, directionadirichlet, val1, val2 );//rigth
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_left = material->CreateBC ( material, -6, directionadirichlet, val1, val2 );//left

    cmesh->InsertMaterialObject ( bc_bottom );
    cmesh->InsertMaterialObject ( bc_rigth );
    cmesh->InsertMaterialObject ( bc_left );

    //cmesh->InsertMaterialObject ( top );
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    return cmesh;

}

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor )
{
    plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    //plasticmatcrisfield * body= dynamic_cast<plasticmatcrisfield *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=(gammaaugua-gammasolo);
    //force[1]=(-gammasolo);
    body->SetLoadFactor(factor);
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
    scalNames.Push ( "StrainVol" );
    scalNames.Push ( "StrainXX" );
    scalNames.Push ( "StrainYY" );
    scalNames.Push ( "StrainZZ" );
    scalNames.Push ( "StrainXY" );
    scalNames.Push ( "StrainXZ" );
    scalNames.Push ( "StrainYZ" );

    scalNames.Push ( "ElStrainVol" );
    scalNames.Push ( "ElStrainXX" );
    scalNames.Push ( "ElStrainYY" );
    scalNames.Push ( "ElStrainZZ" );
    scalNames.Push ( "ElStrainXY" );
    scalNames.Push ( "ElStrainXZ" );
    scalNames.Push ( "ElStrainYZ" );

    scalNames.Push ( "PlStrainVol" );
    scalNames.Push ( "PlStrainXX" );
    scalNames.Push ( "PlStrainYY" );
    scalNames.Push ( "PlStrainZZ" );
    scalNames.Push ( "PlStrainXY" );
    scalNames.Push ( "PlStrainXZ" );
    scalNames.Push ( "PlStrainYZ" );

    scalNames.Push ( "PlStrainSqJ2" );
    scalNames.Push ( "PlStrainSqJ2El" );
    scalNames.Push ( "PlAlpha" );

    scalNames.Push ( "DisplacementX" );
    scalNames.Push ( "DisplacementY" );
    scalNames.Push ( "DisplacementZ" );
    vecNames.Push ( "DisplacementTotal" );


    scalNames.Push ( "YieldSurface1" );
    scalNames.Push ( "YieldSurface2" );
    scalNames.Push ( "YieldSurface3" );

    scalNames.Push ( "POrder" );
    scalNames.Push ( "NSteps" );
    scalNames.Push ( "Cohesion" );
    scalNames.Push ( "FrictionAngle" );
    scalNames.Push ( "FluxX" );
    scalNames.Push ( "FluxY" );
    vecNames.Push ( "Flux" );
    vecNames.Push ( "PrincipalStress" );
    scalNames.Push ( "Pressure" );

}
// void GravityIncrease ( TPZCompMesh * cmesh )
// {
//     int neq = cmesh->NEquations();
//     REAL tol = 0.001;
//     int numiter =5;
//     REAL tol2 = 0.001;
//     int numiter2 =5;
//     REAL l =5;
//     bool linesearch =false;
//     LoadingRamp ( cmesh,1. );
//     bool optimize =false;
//     TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
//     std::ofstream outnewton("saida-newton.txt");
//     anal->IterativeProcessArcLength(outnewton,tol,numiter,tol2,numiter2,l,linesearch);
//     anal->AcceptSolution();
// }
void GravityIncrease ( TPZCompMesh * cmesh )
{

    REAL FS=0.1,FSmax=1000.,FSmin=0.,tol=0.01;
    int neq = cmesh->NEquations();
    int maxcount=100;
    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;

    REAL norm = 1000.;
    REAL tol2 = 0.1;
    int NumIter = 50;
    bool linesearch = true;
    bool checkconv = false;
	std::ofstream outnewton("saida-newton.txt");
	std::ofstream outloadu("gimloadvsu-darcy.nb");
	REAL uy=0.;
	outloadu << "plot = {";
    do {

        std::cout << "FS = " << FS  <<" | Load step = " << counterout << " | Rhs norm = " << norm  << std::endl;
        LoadingRamp ( cmesh,FS );
        bool optimize =true;
        TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
		chrono::steady_clock sc;
		auto start = sc.now();
        int iters;
       bool conv = anal->IterativeProcess ( outnewton, tol2, NumIter,linesearch,checkconv ,iters);

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
            outloadu << "{ "<<-uy << ", " << FS << " } ," << endl;
			displace0 = anal->Solution();
            FSmin = FS;
            anal->AcceptSolution();

			//FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 1.5 );
            FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }
       // cout << "|asdadadasd =  " << std::endl;
        counterout++;

    }  while ( (( FSmax - FSmin ) / FS > tol && counterout<maxcount) );
outloadu <<  " }; ListLinePlot[plot,PlotRange->All]" << endl;
	TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,true );
	anal->AcceptSolution();
}
// void ShearRed ( TPZCompMesh * cmesh )
// {
// 	//plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
//     LoadingRamp ( cmesh,1. );
//
//     REAL FS=0.1,FSmax=10000.,FSmin=0.,tol=1.e-3;
//     int neq = cmesh->NEquations();
//
//     TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );
//
//     int counterout = 0;
//
//     plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
//     //plasticmatcrisfield *material = dynamic_cast<plasticmatcrisfield *> ( cmesh->MaterialVec() [1] );
//
//     //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
//     TPZPlasticStepVoigt<TPZMohrCoulombVoigt,TPZElasticResponse> LEMC = material->GetPlasticity();
//
//     TPZElasticResponse ER = LEMC.fER;
//     REAL phi0 = LEMC.fYC.Phi();
//     REAL cohesion0 = LEMC.fYC.Cohesion();
//     REAL phi,psi,c;
//     REAL norm = 1000.;
//     REAL tol2 = 1.e-3;
//     int NumIter = 30;
//     bool linesearch = true;
//     bool checkconv = false;
//     do {
//
//         bool optimize =true;
//         TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
//
//         anal->IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
//         norm = Norm ( anal->Rhs() );
//
//         //std::cout << "| Load step = " << counterout << " | Rhs norm = " << norm << "| delta displacement norm = "<< unorm << std::endl;
//         if ( norm>= tol2 ) {
//             displace = displace0;
//             FSmax = FS;
//             FS = ( FSmin + FSmax ) / 2.;
//         } else {
//             displace0 = displace;
//             FSmin = FS;
//             FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
//         }
//         std::cout << "FS = " << FS << std::endl;
//         c=cohesion0/FS;
//         phi=atan ( tan ( phi0 ) /FS );
//         psi=phi;
//         LEMC.fYC.SetUp ( phi, psi, c, ER );
//         material->SetPlasticity ( LEMC );
//         counterout++;
//     }  while ( ( FSmax - FSmin ) / FS > tol );
//     bool optimize =false;
//     TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
//     anal->IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
//     anal->AcceptSolution();
// }
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
TPZElastoPlasticAnalysis * CreateAnal ( TPZCompMesh *cmesh,bool optimize )
{
    int numthreads=2;

    TPZElastoPlasticAnalysis * analysis =  new TPZElastoPlasticAnalysis ( cmesh ); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );
    //TPZFStructMatrix matskl ( cmesh );

    matskl.SetNumThreads ( numthreads );

    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
	//step.SetDirect ( ENoDecompose );
    //step.SetDirect ( ECholesky );
    //step.SetDirect ( ELL );
    analysis->SetSolver ( step );

//     long neq = cmesh->NEquations();
//     TPZVec<long> activeEquations;
//     analysis->GetActiveEquations ( activeEquations );
//     TPZEquationFilter filter ( neq );
//     filter.SetActiveEquations ( activeEquations );
//     matskl.EquationFilter() = filter;
//     analysis->SetStructuralMatrix ( matskl );

    //step.SetDirect(ECholesky);
    analysis->SetSolver ( step );


    return analysis;
}
void SetFlux ( TPZCompMesh * plasticCmesh,TPZCompMesh* incmesh)
{

    //int nels0 = incmesh->NElements();
    int nels0 = plasticCmesh->NElements();

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

            mem[globpt].fPlasticState.fflux.Resize ( 2 );
            //mem[globpt].fPlasticState.fflux[0]=data1.dsol[0][0];
            //mem[globpt].fPlasticState.fflux[1]=data1.dsol[0][1];
            mem[globpt].fPlasticState.fflux[0]=flux[0];
            mem[globpt].fPlasticState.fflux[1]=flux[1];
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


