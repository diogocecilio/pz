
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

TPZFMatrix<REAL> SolveKL(TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M);

void ComputeDeterministic();

int betax;
bool water;
REAL coesaofiu;
REAL phifiu;
REAL gammasolo;
REAL gammaaugua;
int main()
{
    int porder=2;
    int nref=3;

    TPZGeoMesh *gmesh = CreateGMesh ( nref);

    REAL lx=20.;
    REAL ly=2.;
    int type=3;
    int M=30;
    REAL mean=20.;
    REAL cov=0.2;
    int samples=10;

    TPZFMatrix<REAL> phi= SolveKL( gmesh, porder, lx, ly,   type,  M);

    //phi.Print(std::cout);

    //TPZFMatrix<REAL> hhat= CreateLogNormalRandomField(phi,  mean,  cov,samples);

    //hhat.Print(std::cout);

   // TPZCompMesh * cmesh =  CreateCMeshRF ( gmesh,porder );

//      int dim = cmesh->Reference()->Dimension();
//      KLAnalysis * klanal = new KLAnalysis ( cmesh );
//
//      KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );
//      klanal->SetExpansionOrder ( mat->GetExpansioOrder() );
//      klanal->Solve();

   // ComputeDeterministic();

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

TPZFMatrix<REAL> SolveKL(TPZGeoMesh* gmesh,REAL porder,REAL lx,REAL ly,  int type, int M)
{
     TPZCompMesh * cmesh = CreateCMeshRF (  gmesh, porder, lx, ly,  type,  M );

     KLAnalysis * klanal = new KLAnalysis ( cmesh );

     KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );

     klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

     klanal->Solve();

     TPZFMatrix<REAL> sol = klanal->Solution();

     TPZFMatrix<REAL> hhat= CreateLogNormalRandomField(sol,  20,  0.2,1000);

     //hhat.Print(std::cout);

     cmesh->LoadSolution(hhat);

     string file ="vtkfolder/outkl.vtk";
     klanal->Post ( file,2,0 );

     return sol;
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



void ComputeDeterministic()
{
    int porder=1;
    int nref=2;

    betax=45;//30,45,60,90
    gammasolo=20.;
    gammaaugua=10.;
    water=true;
    if(!water)gammaaugua=0.;
    cout<<"gamma agua = "<<gammaaugua <<endl;
    coesaofiu=10.;
    phifiu=30.*M_PI/180;

    REAL tolfs = 1.e-3;
    int numiterfs =20;
    REAL tolres = 1.e-3;
    int numiterres =20;
    REAL l =0.5;
    REAL lambda0=0.1;

    TPZGeoMesh *gmesh = CreateGMesh ( nref);


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

    string file;

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
    int numthreads=12;

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



