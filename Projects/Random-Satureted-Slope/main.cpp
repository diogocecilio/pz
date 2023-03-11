
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

#include "TPZFrontMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "tpzsparseblockdiagonalstructmatrix.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"

#include "KLMaterial.h"
#include "KLStrMatrix.h"
#include "KLAnalysis.h"
#include "KLInterpolationSpace.h"
#include "KLRandomField.h"
#include <sys/stat.h>
#include <thread>

using namespace std;

//typedef TPZPlasticStepVoigt<TPZMohrCoulombVoigt,TPZElasticResponse> LEMC;
typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;


typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;



#include "readgidmesh.h"
#include <TPZParFrontStructMatrix.h>



TPZGeoMesh * CreateGMeshGid ( int ref );

TPZGeoMesh *  CreateGMesh ( int ref );

TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder );

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor,REAL gammasolo, REAL gammaagua );

void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void Post ( TPZPostProcAnalysis * postproc,std::string vtkFile );

TPZElastoPlasticAnalysis * CreateAnal ( TPZCompMesh *cmesh,bool optimize );

void GravityIncrease ( TPZCompMesh * cmesh );

void SetFlux ( TPZCompMesh * plasticCmesh,TPZCompMesh* incmesh);

void SolveDarcyProlem(TPZCompMesh *cmesh,string vtk);

void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

TPZCompMesh * CreateCMeshDarcy( TPZGeoMesh *gmesh, int pOrder );

void PostDarcy(TPZAnalysis * analysis,string vtk);

TPZManVector<TPZCompMesh *,2> CreateFields ( TPZGeoMesh * gmesh,int porder,int samples );

TPZFMatrix<REAL> CreateLogNormalRandomField(TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples);

TPZCompMesh * CreateCMeshRF( TPZGeoMesh* gmesh,int porder,REAL lx,REAL ly, int type, int M );

TPZManVector<  TPZCompMesh *,2>   SolveKL(TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M);

void SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZManVector<TPZCompMesh*,2> vecmesh,int isol,REAL FS );

void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol );

void PrintMat ( std::string out,TPZFMatrix<REAL> mat );

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );

void ComputeDeterministic();

void DarcySolution();

void MonteCarlo(int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder);

TPZManVector<TPZCompMesh *,2> CreateFieldsDummy ( TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M);

TPZManVector<TPZCompMesh *,2> SettingCreateFilds(TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M,bool createfield);

int betax;
bool water;
REAL coesaofiu;
REAL phifiu;
REAL gammasolo;
REAL gammaaugua;
int main()
{




    int porder=1;
    int nref=3;

    TPZGeoMesh *gmesh = CreateGMesh ( nref);
   // TPZGeoMesh *gmesh =CreateGMeshGid (  nref );

    REAL lx=20.;
    REAL ly=2.;
    int type=3;
    int M=30;
    bool createfield=false;


    TPZManVector<TPZCompMesh *,2>  vecmesh(2);
    vecmesh =  SettingCreateFilds(gmesh, porder,lx,ly,   type,  M, createfield);

     MonteCarlo(0, 3,  gmesh, vecmesh, porder);




    //phi.Print(std::cout);

    //TPZFMatrix<REAL> hhat= CreateLogNormalRandomField(phi,  mean,  cov,samples);

    //hhat.Print(std::cout);


    //ComputeDeterministic();

    return 0;
}

TPZCompMesh * CreateCMeshRF ( TPZGeoMesh* gmesh,int porder,REAL lx,REAL ly, int type, int M )
{
    int expansionorder=M;
    REAL Lz=1.;
    int id=1;
    int dim = gmesh->Dimension();

    //gmesh->ResetReference();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    KLMaterial * klmat = new KLMaterial ( id,lx,ly,Lz,dim,type,expansionorder );

    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

TPZManVector<  TPZCompMesh *,2>  SolveKL(TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M)
{
     TPZCompMesh * cmesh = CreateCMeshRF (  gmesh, porder, lx, ly,  type,  M );
     TPZCompMesh * cmesh2 = CreateCMeshRF (  gmesh, porder, lx, ly,  type,  M );

     KLAnalysis * klanal = new KLAnalysis ( cmesh );

     KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );

     klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

     klanal->Solve();

     TPZFMatrix<REAL> sol = klanal->Solution();

     TPZManVector<TPZFMatrix<REAL>,2> fields(2);
     TPZManVector< TPZCompMesh *,2> meshes(2);
    REAL meancoes =10.;
    REAL meanatrito =30.*M_PI/180.;
    REAL covcoes=0.3;
    REAL covatrito=0.2;
    int samples=1000;
    fields[0]= CreateLogNormalRandomField(sol,  meancoes,  covcoes,samples);
    fields[1]= CreateLogNormalRandomField(sol,  meanatrito,  covatrito,samples);


    cmesh->LoadSolution(fields[0]);
    cmesh2->LoadSolution(fields[1]);


    meshes[0]=cmesh;
    meshes[1]=cmesh2;


     //hhat.Print(std::cout);

//      cmesh->LoadSolution(hhat);
//
//      string file ="vtkfolder/outkl.vtk";
//      klanal->Post ( file,2,0 );

     return meshes;
}


TPZFMatrix<REAL> CreateLogNormalRandomField(TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples)
{
    TPZFMatrix<REAL>  PHIt;

     PHI.Transpose ( &PHIt );

     int M = PHI.Cols();

     cout << "M = " << M <<endl;

     std::normal_distribution<double> distribution ( 0., 1. );

    TPZFMatrix<REAL>  THETA ( M, samples, 0. );
    for ( int isample = 0; isample < samples; isample++ ) {
        for ( int irdvar = 0; irdvar < M; irdvar++ ) {
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            REAL xic = distribution ( generator );
            REAL xiphi = distribution ( generator );
            THETA(irdvar,isample) = xic;
        }
    }


      TPZFMatrix<REAL>  hhat;

      PHI.Multiply( THETA, hhat );

    REAL sdev = cov * mean;
    REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean ),2 ) ) );
    REAL lambda = log ( mean ) - xi * xi / 2.;
    for ( int i = 0; i < hhat.Rows(); i++ ) {
        for ( int j = 0; j < hhat.Cols(); j++ ) {
			REAL temp =  hhat(i,j);
            hhat(i,j) = exp ( lambda + xi * temp );
        }
    }

     return hhat;
}


void DarcySolution()
{
    int porder=2;
    int nref=0;

    gammaaugua=10.;

    TPZGeoMesh *gmesh = CreateGMeshGid ( nref);

    string vtk = "Darcyx.vtk";

    std::string vtkFiled ="vtkfolder/deterministicflux.vtk";

    TPZCompMesh *  darcycompmesh =  CreateCMeshDarcy(gmesh,porder);

    cout << "\n Solving Darcy... " << endl;
    SolveDarcyProlem(darcycompmesh, vtk);


}

void ComputeDeterministic()
{
    int porder=2;
    int nref=0;

    gammasolo=20.;
    gammaaugua=10.;
    water=false;
    if(!water)gammaaugua=0.;
    cout<<"gamma agua = "<<gammaaugua <<endl;
    coesaofiu=50.;
    phifiu=20.*M_PI/180;

    REAL tolfs = 1.e-3;
    int numiterfs =20;
    REAL tolres = 1.e-3;
    int numiterres =20;
    REAL l =0.5;
    REAL lambda0=0.1;

    TPZGeoMesh *gmesh = CreateGMeshGid ( nref);


    string vtk = "Darcyx.vtk";


    TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    TPZPostProcAnalysis * postprocdetergim = new TPZPostProcAnalysis();
    std::string vtkFiled ="vtkfolder/deterministicflux.vtk";
    std::string vtkFiledgim ="vtkfolder/outvtk.vtk";


    TPZCompMesh *  darcycompmesh =  CreateCMeshDarcy(gmesh,porder);
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

    LoadingRamp ( cmeshgi,1.,gammasolo,gammaaugua );
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
}

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

TPZCompMesh * CreateCMeshDarcy( TPZGeoMesh *gmesh, int pOrder )
{
    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    cmesh->SetDefaultOrder ( pOrder );
    cmesh->SetDimModel ( 2 );
    cmesh->SetAllCreateFunctionsContinuous();

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

     TPZMaterial * BCond0 = material->CreateBC ( material, -3, 0, val1, val2 );//tr
     BCond0->SetForcingFunction(pressure);

     TPZMaterial * BCond1 = material->CreateBC ( material, -6, 0, val1, val2 );//ramp
     BCond1->SetForcingFunction(pressure);

	//val2(0,0)=-1;
	TPZMaterial * BCond2 = material->CreateBC ( material, -4, 0, val1, val2 );//tl
	BCond2->SetForcingFunction(pressure);

    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);

    //Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        disp[0] = -1;
}

void Post ( TPZPostProcAnalysis * postproc,std::string vtkFile )
{

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    postproc->DefineGraphMesh ( 2,scalNames,vecNames,vtkFile );

    postproc->PostProcess ( 0 );
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
//     for ( int d = 0; d<ref; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             if(iel==4||iel==2||iel==6||iel==5)
//             {
//                 TPZGeoEl *gel = gmesh->ElementVec() [iel];
//                 gel->Divide ( subels );
//
//             }
//         }
//     }
        for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    //gel->HasSubElement();

    std::ofstream files ( "ge.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}


TPZGeoMesh * CreateGMeshGid ( int ref )
{
    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    string file ="/home/diogo/projects/pz/data/teste.msh";

    readgidmesh read;

    std::vector<std::vector<int>> meshtopol;
    std::vector<std::vector<double>> meshcoords;

    read.ReadMesh2(meshtopol,meshcoords,file);


    cout << "a" << endl;
    int ncoords = meshcoords.size();
    gmesh->NodeVec().Resize ( ncoords );

    TPZVec<REAL> coord(2);
    for ( int inode=0; inode<ncoords; inode++ ) {
        coord[0] = meshcoords[inode][1];
        coord[1] = meshcoords[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }

        TPZVec <long> TopoTri ( 3 );
        TPZVec <long> TopoLine ( 2 );
        for ( int iel=0; iel<meshtopol.size(); iel++ ) {
            if ( meshtopol[iel].size() ==4 ) {
                TopoTri[0] =meshtopol[iel][1]-1;
                TopoTri[1] =meshtopol[iel][2]-1;
                TopoTri[2] =meshtopol[iel][3]-1;
                new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

            }else if(meshtopol[iel].size() ==3){
                TopoLine[0] = meshtopol[iel][1]-1;
                TopoLine[1] = meshtopol[iel][2]-1;
                REAL x0 = meshcoords[TopoLine[0]][1];
                REAL y0 = meshcoords[TopoLine[0]][2];
                REAL xf = meshcoords[TopoLine[1]][1];
                REAL yf = meshcoords[TopoLine[2]][2];
                REAL tol=1.e-3;
                if((fabs((y0-0))<tol && fabs((yf-0))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
                }
                else if((fabs((x0-75))<tol && fabs((xf-75))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
                }
                else if((fabs((y0-30))<tol && fabs((yf-30))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
                }
                else if((fabs((x0-40))<tol && fabs((xf-40))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
                }
                else if((fabs((x0-0))<tol && fabs((xf-0))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
                }
                else if((fabs((xf-x0))>tol && fabs((yf-y0))>tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
                }
                else{
                    DebugStop();
                }

            }

        }

    gmesh->BuildConnectivity();
    cout << "c" << endl;
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }
    gmesh->Print(std::cout);
    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
cout << "d" << endl;
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
    auto * bc_left = material->CreateBC ( material, -5, directionadirichlet, val1, val2 );//left

    cmesh->InsertMaterialObject ( bc_bottom );
    cmesh->InsertMaterialObject ( bc_rigth );
    cmesh->InsertMaterialObject ( bc_left );

    //cmesh->InsertMaterialObject ( top );
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    return cmesh;

}

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor,REAL gammasolo, REAL gammaagua)
{
    plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    //plasticmatcrisfield * body= dynamic_cast<plasticmatcrisfield *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=(gammaagua-gammasolo);
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



TPZElastoPlasticAnalysis * CreateAnal ( TPZCompMesh *cmesh,bool optimize )
{
    int numthreads=0;

    TPZElastoPlasticAnalysis * analysis =  new TPZElastoPlasticAnalysis ( cmesh ); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );
    //TPZFStructMatrix matskl ( cmesh );

    matskl.SetNumThreads ( numthreads );

    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    //step.SetDirect ( ECholesky );
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



void SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZManVector<TPZCompMesh*,2> vecmesh,int isol,REAL FS )
{

    int nels = plasticCmesh->NElements();
    int nels0 = vecmesh[0]->NElements();
    int nels1 = vecmesh[1]->NElements();
    //vecmesh[1]->Solution().Print ( "SOLUTION PHI" );
    if ( nels!=nels0 || nels!=nels1 || nels0!=nels1 ) {
        //cout << "nels "<< nels << " | nels0 = "<< nels0 <<" | nels1 = "<< nels1 <<endl;
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

        datarandom1.fNeedsSol = true;
        datarandom2.fNeedsSol = true;
        intelplastic->InitMaterialData ( dataplastic );
        intelrandom1->InitMaterialData ( datarandom1 );
        intelrandom2->InitMaterialData ( datarandom2 );

        REAL weight=0;
        int nint = intpoints.NPoints();

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

            TPZSolVec sol1,sol2;
            ComputeSolution ( celrandom1,datarandom1.phi,sol1 );
            ComputeSolution ( celrandom2,datarandom2.phi,sol2 );

            TPZVec<REAL> cohes = sol1[isol];
            TPZVec<REAL> phi = sol2[isol];

            dataplastic.intGlobPtIndex = globpt;
            mem[globpt].fPlasticState.fmatprop[0]= cohes[0];
            mem[globpt].fPlasticState.fmatprop[1]= phi[0] ;

            TPZTensor<REAL> epsTotal,sigma;

            globpt++;

        }

        //pMatWithMem2->SetUpdateMem ( false );
    }
    pMatWithMem2->SetUpdateMem ( false );

    plasticCmesh->Solution().Zero();


}

void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol )
{

    const int numdof = cel->Material()->NStateVariables();


    int counter=0;

    std::set<long> cornerconnectlist;
    cel->BuildCornerConnectList(cornerconnectlist);
    int ncon = cornerconnectlist.size();

    TPZFMatrix<STATE> &MeshSol = cel->Mesh()->Solution();

    long numbersol = MeshSol.Cols();
    sol.Resize ( numbersol );
    for ( long is=0 ; is<numbersol; is++ ) {
        sol[is].Resize ( numdof );
        sol[is].Fill ( 0. );

    }
    TPZBlock<STATE> &block = cel->Mesh()->Block();
    long iv = 0, d;
    for ( int in=0; in<ncon; in++ ) {
        TPZConnect *df = &cel->Connect ( in );
        /** @brief Returns the Sequence number of the connect object */
        /** If the \f$ sequence number == -1 \f$ this means that the node is unused */
        long dfseq = df->SequenceNumber();

        /**
        * @brief Returns block dimension
        * @param block_diagonal Inquired block_diagonal
        */
        int dfvar = block.Size ( dfseq );
        /**
        * @brief Returns the position of first element block dependent on matrix diagonal
        * @param block_diagonal Inquired block_diagonal
        */
        long pos = block.Position ( dfseq );
        for ( int jn=0; jn<dfvar; jn++ ) {
            //cout << "\n pos+jn = " << pos+jn << endl;
            //cout << " ids[in] = " << ids[in] << endl;

            for ( int is=0; is<numbersol; is++ ) {

                sol[is][iv%numdof] += ( STATE ) phi ( iv/numdof,0 ) *MeshSol ( pos+jn,is );
            }
            iv++;
        }
    }

}//metho

void MonteCarlo(int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder)
{


	std::string vtkFile0,str;


    vtkFile0="/vtkfolder/out-gim-mc";

    REAL tolfs = 1.e-3;
    int numiterfs =20;
    REAL tolres = 1.e-3;
    int numiterres =20;
    REAL l =0.5;
    REAL lambda0=0.1;


    for ( int isample=a; isample<b; isample++ ) {

		TPZPostProcAnalysis * postproc = new TPZPostProcAnalysis();
		auto s = std::to_string ( isample );

		//create the elastoplastic comp mesh
        TPZCompMesh *cmesh = CreateCMesh ( gmesh,porder );

        chrono::steady_clock sc;
        auto start = sc.now();

        TPZManVector<REAL,10> out;


		SetMaterialParamenters ( cmesh,vecmesh,isample,1 );

        TPZElastoPlasticAnalysis  * anal0 = CreateAnal ( cmesh,true );
        //for GIM set material param with FS = 1 only once.
        //The loading will increase inside the method, but the Strength will be constant
        anal0->IterativeProcessArcLength(tolfs,numiterfs,tolres,numiterres,l,lambda0);

		auto end = sc.now();
		auto time_span = static_cast<chrono::duration<double>> ( end - start );
		cout<< "\n \n *************************** ";
		cout << "\n Gravity Increase took: " << time_span.count() << " seconds !!!";

		string  filename;
		string datafile = "/information";
		string ext = ".txt";
		filename += datafile;

		filename += s;
		filename += ext;
		std::ofstream fileinfo ( filename );
		fileinfo << "Monte Carlo Sample = " << isample << std::endl;
		fileinfo << "Safety Factor = " << out[0] << std::endl;
		fileinfo << "counterout = " << out[1] << std::endl;
		fileinfo << "rnorm = " << out[2] << std::endl;

		string vtkFile=vtkFile0;
		string ext2=".vtk";
		vtkFile+=s;
		vtkFile+=ext2;
		cout << "\n vtkFile = " << vtkFile << endl;
        postproc->SetCompMesh(0);
        CreatePostProcessingMesh ( postproc, cmesh );

        Post ( postproc,vtkFile );
		delete cmesh;
    }
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

    out.Resize ( coords.size(),coords[1].size() );
    for ( int iel=0; iel<coords.size(); iel++ ) {
        for ( int elnode=0; elnode<coords[iel].size(); elnode++ ) {
            out ( iel,elnode ) = coords[iel][elnode];
        }
    }
}
TPZManVector<TPZCompMesh *,2> CreateFieldsDummy ( TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M )
{


    TPZCompMesh * cmesh =  CreateCMeshRF ( gmesh, porder, lx, ly,  type,  M );
    TPZCompMesh * cmesh2 =  CreateCMeshRF (gmesh, porder, lx, ly,  type,  M );

    // InsertMat ( cmesh, porder );
    // InsertMat ( cmesh2, porder );

    TPZManVector<TPZCompMesh *,2> vecmesh ( 2 );

    vecmesh[0]=cmesh;
    vecmesh[1]=cmesh2;

    return vecmesh;
}

TPZManVector<TPZCompMesh *,2> SettingCreateFilds(TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M,bool createfield)
{

	string outco,outphi;
    outco="/home/diogo/projects/pzbuildrelease/Projects/Random-Satureted-Slope/coes.txt";
	outphi="/home/diogo/projects/pzbuildrelease/Projects/Random-Satureted-Slope/phi.txt";

    TPZFMatrix<REAL> readco,readphi;

	TPZManVector<TPZCompMesh *,2> vecmesh;
	//if true creates random fields
	//if false read random fields
    if (createfield ) {
        vecmesh = SolveKL( gmesh, porder, lx, ly,type, M);
        PrintMat ( outco,vecmesh[0]->Solution() );
        PrintMat ( outphi,vecmesh[1]->Solution() );

    } else {
        vecmesh = CreateFieldsDummy (gmesh, porder, lx, ly,type, M );
        ReadFile ( outco,readco );
        ReadFile ( outphi,readphi );

    }

    cout << "asd" <<endl;
    //create analysis to reorder the solution to transfer correctly to the determinitic mesh
    TPZElastoPlasticAnalysis * anal = new TPZElastoPlasticAnalysis ( vecmesh[0] );
    TPZElastoPlasticAnalysis * anal1 = new TPZElastoPlasticAnalysis ( vecmesh[1] );
	//load the random fields solution in dummy meshes
    vecmesh[0]->LoadSolution ( readco );
    vecmesh[1]->LoadSolution ( readphi );

	return vecmesh;
}
