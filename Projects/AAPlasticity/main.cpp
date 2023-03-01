
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

using namespace std;


typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;


TPZVec<REAL> fPlasticDeformSqJ2;

#include "readgidmesh.h"

TPZGeoMesh * CreateGMesh ( int ref );
TPZGeoMesh * CreateGMeshGidTwoMats ( int ref );
TPZGeoMesh * CreateGMeshGidTwoMatsTailings ( int ref );
TPZGeoMesh * CreateGMeshGid ( int ref, TPZFMatrix<REAL> pts ,string gridname);

//TPZCompMesh * CreateCMesh(TPZGeoMesh * gmesh,int porder);
TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder, REAL c, REAL phi, REAL poisson, REAL young );

TPZCompMesh * CreateCMeshTwoMats(TPZGeoMesh * gmesh,int porder);



void LoadingRamp (TPZCompMesh * cmesh,  REAL factor);

void LoadingRamptwo (TPZCompMesh * cmesh,  REAL factor);

void IntegrateBCSol(TPZCompMesh * cmesh);

void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);

void  CreatePostProcessingMesh(TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh);

void  CreatePostProcessingMeshTwoMats(TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh);

void Post(TPZPostProcAnalysis * postproc,std::string vtkFile );

TPZElastoPlasticAnalysis * CreateAnal(TPZCompMesh *cmesh);

void ShearRed ( TPZCompMesh * cmesh);

void DivideElementsAbove(TPZCompMesh* cmesh,TPZGeoMesh* gmesh,REAL sqj2, std::set<long> &elindices);

void ComputeElementDeformation(TPZCompMesh* cmesh);

REAL findnodalsol(TPZCompMesh *cmesh);

void ShearRed ( TPZCompMesh * cmesh );

void MaterialPointMohrCoulomb();

void FadTest();

STATE ComputeBCWork(TPZCompMesh &cmesh);

TPZManVector<REAL,10> GravityIncrease ( TPZCompMesh * cmesh );

void SolveRamp(TPZCompMesh * cmesh);

TPZGeoMesh * CreateGMeshCentrifuga ( int ref );

TPZCompMesh * CreateCMeshCentrifuga(TPZGeoMesh * gmesh,int porder);
std::vector<std::vector<double>> ReadMeshPoints(std::string file);
std::vector<std::vector<double>> ReadMeshPoints(std::string file)
{
    ifstream myfile ( file );
    std::vector<std::vector<double>> coords;
    string line,temp;
    double temp2;
    
    if ( myfile.is_open() ) {
        while ( getline ( myfile, line ) ) {
            std::vector<double> tokens;
            istringstream iss ( line );
            while ( iss >> temp2 )tokens.push_back ( temp2 );
            coords.push_back ( tokens );

        }
    } else std::cout << "Unable to open file";

    return coords;
}

void LoadingRamp (TPZCompMesh * cmesh,  REAL factor,REAL gamma);

void MaterialPointMohrCoulomb2() {
//left

//     REAL mc_cohesion    =490.;
//     REAL mc_phi         = ( 20.0*M_PI/180 );
//     REAL mc_psi         = mc_phi;
//
//     TPZElasticResponse ER;
//     REAL nu = 0.48;
//     REAL E = 10000000;

  //     Mohr Coulomb data
    REAL mc_cohesion    =50.;
    REAL mc_phi         = ( 20.0*M_PI/180 );
    REAL mc_psi         = mc_phi;

    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000;

    ER.SetUp( E, nu );
    TPZMohrCoulombVoigt *mc = new  TPZMohrCoulombVoigt(mc_phi,mc_phi,mc_cohesion,ER);
//plane
    REAL epsbarnew = 0.;
     TPZFMatrix<REAL> dep;
     TPZTensor<REAL> sigma_trial,sigma_proj,a,b;


    //MAIN
    sigma_trial.XX()=3891.34;
    sigma_trial.YY()=3059.69;
    sigma_trial.ZZ()=3406.01;
    sigma_trial.XZ()=0;
    sigma_trial.YZ()=0;
    sigma_trial.XY()=168.535;
    mc->ProjectSigmaDep(sigma_trial, sigma_proj,dep,  epsbarnew);
    sigma_proj.Print(std::cout);
    dep.Print(std::cout);


    //RIGTH
    sigma_trial.XX()=5226.4;
    sigma_trial.YY()=4860.25;
    sigma_trial.ZZ()=4942.46;
    sigma_trial.XZ()=0;
    sigma_trial.YZ()=0;
    sigma_trial.XY()=-219.96;
    mc->ProjectSigmaDep(sigma_trial, sigma_proj,dep,  epsbarnew);
    sigma_proj.Print(std::cout);
    dep.Print(std::cout);


    //APEX
    sigma_trial.XX()=5024.51;
    sigma_trial.YY()=4666.28;
    sigma_trial.ZZ()=4748.49;
    sigma_trial.XZ()=0;
    sigma_trial.YZ()=0;
    sigma_trial.XY()=-194.931;
    mc->ProjectSigmaDep(sigma_trial, sigma_proj,dep,  epsbarnew);
    sigma_proj.Print(std::cout);
    dep.Print(std::cout);



    mc_cohesion    =490.;
    mc_phi         = ( 20.0*M_PI/180 );
    mc_psi         = mc_phi;
    nu = 0.48;
    E = 10000000;

    ER.SetUp( E, nu );
    TPZMohrCoulombVoigt *mc2 = new  TPZMohrCoulombVoigt(mc_phi,mc_phi,mc_cohesion,ER);

    //LEFT
    sigma_trial.XX()=-1465.15;
    sigma_trial.YY()=-1353.51;
    sigma_trial.ZZ()=-210.626;
    sigma_trial.XZ()=0;
    sigma_trial.YZ()=0;
    sigma_trial.XY()=1271.35;
    mc2->ProjectSigmaDep(sigma_trial, sigma_proj,dep,  epsbarnew);
    //mc->ReturnMapLeftEdge ( sigma_trial, sigma_proj,dep,  epsbarnew);
    sigma_proj.Print(std::cout);
    dep.Print(std::cout);

// epsttr = {-7.1436465404952029 10^-5, -5.4913608546659984 10^-5,
// 1.1423294530637051 10^-4, 0, 0, (3.7632073167565092) 10^-4};

    TPZTensor<REAL> eps;
    eps.XX()=-7.1436465404952029e-5;
    eps.YY()=-5.4913608546659984e-5;
    eps.ZZ()=1.1423294530637051e-4;
    eps.XZ()=0;
    eps.YZ()=0;
    eps.XY()=0.5*3.7632073167565092e-4;

    ER.ComputeStress(eps,sigma_trial);
    sigma_trial.Print(std::cout);

}

int main()
{
        MaterialPointMohrCoulomb2();
// int porder=2;
//     string file = "/home/diogo/Dropbox/Limit-Equilibrium-Pyton/pts-h89.9299-beta45.dat";
//    std::vector<std::vector<double>> coords = ReadMeshPoints(file);
//    int rows=coords.size();
//    int cols = coords[0].size();
//    TPZFMatrix<REAL> pts(rows,cols);
//     for(int it=0;it<coords.size();it++)
//     {
//         for(int jt=0;jt<coords[it].size();jt++)
//         {
//             pts(it,jt) =coords[it][jt];
//         }
//     }
//     pts.Print(std::cout);
//     //string gridname ="/home/diogo/Dropbox/Limit-Equilibrium-Pyton/h143.888-beta45.msh";
//     string gridname ="/home/diogo/Dropbox/Limit-Equilibrium-Pyton/h75.7305-beta45.msh";
//     TPZGeoMesh *gmesh = CreateGMeshGid ( 0 ,pts, gridname);
//     REAL c =50.;
//     REAL phi = 34.1;
//     REAL poisson=0.49;
//     REAL young = 20000.;
//     TPZCompMesh *cmesh =   CreateCMesh ( gmesh, porder,  c,  phi,  poisson,  young );
//     TPZElastoPlasticAnalysis  * analysis1 =  CreateAnal(cmesh);
//
//     ShearRed(cmesh);
//     //GravityIncrease ( cmesh );
//     //SolveRamp(cmesh);
//
//
//     int maxiter=100;
//     REAL tol=1.e-3;
//     bool linesearch=true;
//     bool checkconv=false;
//     int steps=10,iters;
//     bool conv = analysis1->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv,iters);
//     analysis1->AcceptSolution();
//
//     TPZPostProcAnalysis  * postproc1 = new TPZPostProcAnalysis();
//
//     postproc1->SetCompMesh(0);
//
//     CreatePostProcessingMesh( postproc1, cmesh);
//
//     std::string vtkFile1 ="/home/diogo/Dropbox/Limit-Equilibrium-Pyton/srm-pts-h75.7305-beta45-c50.vtk";
//     Post(postproc1,vtkFile1);
//
    return 0;
}


// int main()
// {
//     string file = "/home/diogo/Dropbox/Limit-Equilibrium-Pyton/pts.dat";
//    std::vector<std::vector<double>> coords = ReadMeshPoints(file);
//     for(int it=0;it<coords.size();it++)
//     {
//         for(int jt=0;jt<coords[it].size();jt++)
//         {
//             cout << coords[it][jt] << " ";
//         }
//         cout << endl;
//     }
//     return 0;
//     int porder =2;
// 
// 
//     int var=1;
//     enum mesh { trinta = 0, quarenta = 1, quarentacinco = 2 ,sessenta=3};
//     {
//         TPZFMatrix<REAL> pts(6,2);
//         string gridname;
//         REAL delta=0.5;
// 
//         
// //         if(mesh::trinta==var)
// //         {
// //             gridname ="/home/diogo/projects/pz/data/beta-30-h-100.msh";
// //             pts(0,0)=0.;pts(1,0)=750.;pts(2,0)=750.;pts(3,0)=461.602;pts(4,0)=288.398;pts(5,0)=0.;
// //             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=300.;pts(3,1)=300.;   pts(4,1)=400.000;pts(5,1)=400.;
// //         }
// //         if(mesh::quarenta==var)
// //         {
// //             gridname ="/home/diogo/projects/pz/data/beta-40-h-100.msh";
// //             pts(0,0)=0.;pts(1,0)=750.;pts(2,0)=750.;pts(3,0)=434.588;pts(4,0)=315.412;pts(5,0)=0.;
// //             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=300.;pts(3,1)=300.;   pts(4,1)=400.000;pts(5,1)=400.;
// //         }
// //         if(mesh::quarentacinco==var)
// //         {
// //             gridname ="/home/diogo/projects/pz/data/beta-45-h-100.msh";
// //             pts(0,0)=0.;pts(1,0)=750.;pts(2,0)=750.;pts(3,0)=450.;pts(4,0)=350.;pts(5,0)=0.;
// //             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=300.;pts(3,1)=300.;pts(4,1)=400.;pts(5,1)=400.;
// //         }
// //         if(mesh::sessenta==var)
// //         {
// //             gridname ="/home/diogo/projects/pz/data/beta-60-h-100.msh";
// //             pts(0,0)=0.;pts(1,0)=750.;pts(2,0)=750.;pts(3,0)=403.868;pts(4,0)=346.133;pts(5,0)=0.;
// //             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=300.;pts(3,1)=300.;pts(4,1)=400.;pts(5,1)=400.;
// //         }
//         if(mesh::trinta==var)
//         {
//             gridname ="/home/diogo/projects/pz/data/beta-30-h-300.msh";
//             pts(0,0)=0.;pts(1,0)=1000.;pts(2,0)=1000.;pts(3,0)=759.808;pts(4,0)=240.192;pts(5,0)=0.;
//             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=400.;pts(3,1)=400.;   pts(4,1)=700.;pts(5,1)=700.;
//         }
//         if(mesh::quarenta==var)
//         {
//             gridname ="/home/diogo/projects/pz/data/beta-40-h-300.msh";
//             pts(0,0)=0.;pts(1,0)=1000.;pts(2,0)=1000.;pts(3,0)=678.763;pts(4,0)=321.237;pts(5,0)=0.;
//             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=400.;pts(3,1)=400.;   pts(4,1)=700.;pts(5,1)=700.;
//         }
//         if(mesh::quarentacinco==var)
//         {
//             gridname ="/home/diogo/projects/pz/data/beta-45-h-300.msh";
//             pts(0,0)=0.;pts(1,0)=1000.;pts(2,0)=1000.;pts(3,0)=650;pts(4,0)=350;pts(5,0)=0.;
//             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=400.;pts(3,1)=400.;   pts(4,1)=700.;pts(5,1)=700.;
//         }
//         if(mesh::sessenta==var)
//         {
//             gridname ="/home/diogo/projects/pz/data/beta-60-h-300.msh";
//             pts(0,0)=0.;pts(1,0)=1000.;pts(2,0)=1000.;pts(3,0)=586.603;pts(4,0)=413.397;pts(5,0)=0.;
//             pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=400.;pts(3,1)=400.;   pts(4,1)=700.;pts(5,1)=700.;
//         }
// 
// 
// //         delta=0.5;
// //         gridname ="/home/diogo/projects/pz/data/meshfs1.msh";
// //         pts(0,0)=0.;pts(1,0)=240.45;pts(2,0)=240.45;pts(3,0)=140.45;pts(4,0)=100;pts(5,0)=0.;
// //         pts(0,1)=0.;pts(1,1)=0.;  pts(2,1)=59.55;pts(3,1)=59.55;   pts(4,1)=100.;pts(5,1)=100.;
// 
//         TPZGeoMesh *gmesh = CreateGMeshGid ( 0 ,pts, gridname,delta);
//        // TPZGeoMesh *gmesh = CreateGMeshCentrifuga ( 0 );
//        // TPZGeoMesh *gmesh = CreateGMesh ( 5 );
//         REAL c =50.;
//         REAL phi = 34.1;
//         REAL poisson=0.49;
//         REAL young = 20000000.;
// //         REAL c =50.;
// //         REAL phi = 20;
// //         REAL poisson=0.49;
// //         REAL young = 200000.;
//          TPZCompMesh *cmesh =   CreateCMesh ( gmesh, porder,  c,  phi,  poisson,  young );
//         //TPZCompMesh *cmesh = CreateCMesh ( gmesh, porder );
//         //TPZCompMesh *cmesh = CreateCMeshCentrifuga( gmesh, porder );
// 
//         TPZElastoPlasticAnalysis  * analysis1 =  CreateAnal(cmesh);
// 
//         //ShearRed(cmesh);
//        GravityIncrease ( cmesh );
//        //SolveRamp(cmesh);
// 
//         analysis1->AcceptSolution();
// 
//         TPZPostProcAnalysis  * postproc1 = new TPZPostProcAnalysis();
// 
//         postproc1->SetCompMesh(0);
// 
//         CreatePostProcessingMesh( postproc1, cmesh);
// 
//         std::string vtkFile1 ="shearred-nilo.vtk";
//         Post(postproc1,vtkFile1);
//     }
// 
//     return 0;
// 
//     TPZGeoMesh *gmeshtwo = CreateGMeshGidTwoMatsTailings ( 0 );
// 
//     TPZCompMesh *cmeshtwo = CreateCMeshTwoMats ( gmeshtwo, porder );
// 
//     int maxiter=100;
//     REAL tol=1.e-3;
//     bool linesearch=true;
//     bool checkconv=false;
//     int steps=10;
//     REAL finalload = 0.9;
//     std::string vtkFile ="slopx.vtk";
//     std::ofstream outloadu("loadvu.nb");
//     outloadu << "plot = {";
//     TPZPostProcAnalysis  * postproc = new TPZPostProcAnalysis();
// 
//     //ShearRed (cmesh2);
// 
// 
//     REAL uy=0;
//     for(int iloadstep=0; iloadstep<=steps; iloadstep++)
//     {
// 
// 
// 
//         TPZElastoPlasticAnalysis  * analysis =  CreateAnal(cmeshtwo);
//         // TPZElastoPlasticAnalysis  * analysis =  CreateAnal(cmesh);
// 
// 
//         REAL load = iloadstep*finalload/steps;
//         cout << " \n --------- iloadstep  = "<< iloadstep << endl;
// 
//         cout << " \n --------- load  = "<< load << endl;
// 
//         LoadingRamptwo(cmeshtwo,load);
//         int iters;
//         // LoadingRamp(cmesh,load);
//         bool conv = analysis->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv,iters);
//         
//         //uy+=findnodalsol(cmeshtwo);
//         //uy+=findnodalsol(cmesh);
//         TPZFMatrix<REAL> sol = analysis->Solution();
// 
//         //cmesh->Solution().Print(cout);
//         cout << " sol "<<-uy << " | load = " << load <<    endl;
//         analysis->AcceptSolution();
// 
//         outloadu << "{ "<<-uy << ", " << load << " } ," << endl;
// 
//         postproc->SetCompMesh(0);
//         //CreatePostProcessingMesh( postproc, cmeshtwo);
//         CreatePostProcessingMeshTwoMats( postproc, cmeshtwo);
//         //CreatePostProcessingMesh( postproc, cmesh);
// 
//         Post(postproc,vtkFile);
// 
//         // IntegrateBCSol(cmeshtwo);
//         // ComputeBCWork(*cmeshtwo);
//         //cmesh->LoadSolution();
//         //ComputeElementDeformation(cmesh);
// 
// // 		std::set<long> elindices;
// // 		REAL sqrtj2=0.005;
// 
//         //DivideElementsAbove(cmesh,gmesh,sqrtj2,elindices);
// 
// 
// 
// // 		LoadingRamp(cmesh,load);
// // 		analysis->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv);
// // 		analysis->AcceptSolution();
//         //
// 
// 
//     }
//     outloadu <<  " }; ListLinePlot[plot,PlotRange->All]" << endl;
// 
// // 		TPZElastoPlasticAnalysis *anal  =  CreateAnal(cmesh);
// //
// // 		LoadingRamp(cmesh,finalload);
// //
// //  		anal->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv);
// //
// //  		anal->AcceptSolution();
// //
// // 		CreatePostProcessingMesh( postproc, cmesh);
// //
// // 		Post(postproc,vtkFile);
// //
// //
// //     std::cout << "FINISHED!" << std::endl;
// 
// 
//     return 0;
// }

void SolveRamp(TPZCompMesh * cmesh)
{
    int maxiter=100;
    REAL tol=1.e-3;
    bool linesearch=true;
    bool checkconv=false;
    int steps=40;
    REAL finalload = 5.;
    std::string vtkFile ="tcentrifuga-p2.vtk";
    std::ofstream outloadu("loadvup2.nb");
    outloadu << "plot = {";
    TPZPostProcAnalysis  * postproc = new TPZPostProcAnalysis();
    REAL g0=1.;
    vector<double> loadvec = {1,5,10.,15.,20.,25.,30.,35.,40.,45.,50.};

    for(int iter=0;iter<loadvec.size();iter++)loadvec[iter]*=g0;
    
    findnodalsol(cmesh);
    REAL uy=0.;
    REAL gamma =17.;//kN/m^3
    //for(int iloadstep=0; iloadstep<loadvec.size(); iloadstep++)
    for(int iloadstep=0; iloadstep<100; iloadstep++)
    {

        TPZElastoPlasticAnalysis  * analysis =  CreateAnal(cmesh);
        //REAL load = iloadstep*finalload/steps;
        //REAL load = loadvec[iloadstep];
        REAL load = iloadstep;
        cout << " \n --------- iloadstep  = "<< iloadstep  << " |  load  = "<< load << " gtrial/g0 = "<< load/g0 <<endl;
        LoadingRamp(cmesh,load,gamma);
        int iters;
        bool conv = analysis->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv,iters);
        if(conv==false)
        {
            outloadu <<  " }; ListLinePlot[plot,PlotRange->All]" << endl;
            cout << "newton failed to converge!";
            break;
        }
        
        uy+=findnodalsol(cmesh);
        
        TPZFMatrix<REAL> sol = analysis->Solution();

        cout << " sol "<<-uy << " | load = " << load <<    endl;

        outloadu << "{ "<<-uy << ", " << load << " } ," << endl;
        
        //TPZFMatrix<REAL> sol = analysis->Solution();
        analysis->AcceptSolution();
        postproc->SetCompMesh(0);
        //CreatePostProcessingMeshTwoMats( postproc, cmesh);
        CreatePostProcessingMesh( postproc, cmesh);
        Post(postproc,vtkFile);
    }


}

TPZManVector<REAL,10> GravityIncrease ( TPZCompMesh * cmesh )
{

	TPZManVector<REAL,10> output(10,0.);
    REAL FS=0.5,FSmax=1000.,FSmin=0.,tol=0.001;
    int neq = cmesh->NEquations();
    int maxcount=50;
    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;

    REAL norm = 1000.;
     REAL tol2 = 1.e-2;
    int NumIter = 50;
    bool linesearch = true;
    bool checkconv = false;
	std::ofstream outloadu("outloadu.nb");
    

    bool conv=false;
    do {


        cout << "\n  FS = " << FS  <<" | Load step = " << counterout;
        LoadingRamp ( cmesh,FS );
		//LoadingRamp ( body,0.0000001 );
        bool optimize =true;
        TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh );
        chrono::steady_clock sc;
        auto start = sc.now();
		int iters;
         conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv ,iters);

        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << " | total time in iterative process =  " << time_span.count()<< std::endl;
        //anal->IterativeProcess ( outnewton, tol2, NumIter);

        norm = Norm ( anal->Rhs() );

        if ( conv==false ) {
            cmesh->LoadSolution(displace0);
			//cmesh->Solution().Zero();
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;

        } else {
            
			displace0 = anal->Solution();
            FSmin = FS;
            anal->AcceptSolution();

			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );

        }
        
        counterout++;
		delete anal;
    }  while ( (( FSmax - FSmin ) / FS > tol && counterout<maxcount) || conv==false);

    
	TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh );
	anal->AcceptSolution();
	output[0]=FS;
	output[1]=counterout;
	output[2]=norm;
    return output;
}


void MaterialPointMohrCoulomb() {

    // Mohr Coulomb data
    REAL mc_cohesion    =50.;
    REAL mc_phi         = ( 20.0*M_PI/180 );
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000;

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> mc;
    ER.SetUp( E, nu );
    mc.fER =ER;
    //mc.SetElasticResponse( ER );
    mc.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );
    TPZTensor<REAL> eps,deps;
    TPZTensor<REAL> sigma;
    TPZFMatrix<REAL> Dep;


//     (*MAIN*)
// epsttr = {3.6157562487484481 10^-2, -2.5800526093809846
// 10^-2, 0, 0, 0, (2.5111726131341262) 10^-2};
// epsptr = {0, 0, 0, 0, 0, 0};
// elastmat = {20000, 0.49, 20 Pi/180, 50};
// {sig, dep, g} = SolvePlasticity[epsttr, epsptr, elastmat];
// sig // MatrixForm
// dep // MatrixForm
    eps.XX()=0.036157562487484481;
    eps.YY()=-0.025800526093809846;
    eps.ZZ()=0;
    eps.XZ()=0;
    eps.YZ()=0;
    eps.XY()=1./2.*0.025111726131341262;
    //cout << "affter" << endl;
    mc.ApplyStrainComputeDep(eps,sigma,Dep);

    //cout << "out" << endl;
    sigma.Print(std::cout);
    cout << "Dep" << endl;
    Dep.Print(std::cout);

//     epsttr = {2.1153582845979804 10^-2, -6.1244636739267791 10^-3, 0, 0,
//    0, (-3.2774022989690031) 10^-2};
// epsptr = {0, 0, 0, 0, 0, 0};
// elastmat = {20000, 0.49, 20 Pi/180, 50};
// {sig, dep, g} = SolvePlasticity[epsttr, epsptr, elastmat];
// sig // MatrixForm
// dep // MatrixForm
    eps.XX()=2.1153582845979804e-2;
    eps.YY()=-6.1244636739267791e-3;
    eps.ZZ()=0;
    eps.XZ()=0;
    eps.YZ()=0;
    eps.XY()=1./2.*(-3.2774022989690031e-2);
    //cout << "affter" << endl;
    mc.ApplyStrainComputeDep(eps,sigma,Dep);

    //cout << "out" << endl;
    sigma.Print(std::cout);
    cout << "Dep" << endl;
    Dep.Print(std::cout);


}


// void MaterialPointMohrCoulomb() {
//
//     // Mohr Coulomb data
//     REAL mc_cohesion    =50.;
//     REAL mc_phi         = ( 20.0*M_PI/180 );
//     REAL mc_psi         = mc_phi;
//
//     /// ElastoPlastic Material using Mohr Coulomb
//     // Elastic predictor
//     TPZElasticResponse ER;
//     REAL nu = 0.49;
//     REAL E = 20000;
//
//     TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> mc;
//     ER.SetUp( E, nu );
//     mc.fER =ER;
//     //mc.SetElasticResponse( ER );
//     mc.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );
//     TPZTensor<REAL> eps,deps;
//     TPZTensor<REAL> sigma;
//     TPZFMatrix<REAL> Dep;
//
//     sigma.XX()=17658.1;
//     sigma.YY()=15494.7;
//     sigma.ZZ()=16244.8;
//     sigma.XY()=-19.7463;
//     ER.ComputeDeformation(sigma,eps);
//     REAL kprev=0.;
//     TPZVec<REAL> coef;
//     //cout << "before" << endl;
//     //eps.Print(std::cout);
//     //cout << "affter" << endl;
//     mc.ApplyStrainComputeDep(eps,sigma,Dep);
//     deps.XX()=0.000001;
//     deps.YY()=0.000002;
//     deps.ZZ()=0.0000015;
//     deps.XZ()=-0.000001;
//     deps.YZ()=0.000002;
//     deps.XY()=0.000003;
//     mc.TaylorCheck(eps,deps,kprev,coef);
//     coef.Print(cout);
//     //cout << "out" << endl;
//     //sigma.Print(std::cout);
//     cout << "Dep" << endl;
//     Dep.Print(std::cout);
//
//
// }


const auto rhs3 = [](const TPZVec<REAL>&loc, TPZVec<STATE> &dx,TPZFMatrix<STATE> &d2x) {
    const REAL &x0 = loc[0];
    const REAL &y0 = loc[1];

    const int N=2;
    dx.Resize(N);
    d2x.Resize(N,N);
    TFad<N,TFad<N,float> > xsecnd2,ysecnd2,test2;
    TFad<N,float> x(x0),y(y0);

    x.fastAccessDx(0)=1;
    x.fastAccessDx(1)=0;
    y.fastAccessDx(0)=0;
    y.fastAccessDx(1)=1;
    xsecnd2=x;
    xsecnd2.fastAccessDx(0)=1;
    xsecnd2.fastAccessDx(1)=0;
    ysecnd2=y;
    ysecnd2.fastAccessDx(0)=0;
    ysecnd2.fastAccessDx(1)=1;
    test2 = (-exp(ysecnd2) + xsecnd2/ysecnd2)*(-exp(ysecnd2) + xsecnd2/ysecnd2 + sin(xsecnd2/ysecnd2));
    
    cout << " aquiiiiiiiiii f = " << test2 << endl;

    dx[0]=test2.val().fastAccessDx(0);
    dx[1]=test2.val().fastAccessDx(1);

    d2x(0,0)=test2.fastAccessDx(0).fastAccessDx(0);
    d2x(1,1)= test2.fastAccessDx(1).fastAccessDx(1);
    d2x(1,0)= test2.fastAccessDx(0).fastAccessDx(1);
    d2x(0,1)= test2.fastAccessDx(1).fastAccessDx(0);

};

const auto rhs2 = [](const TPZFMatrix<STATE>&loc, TPZFMatrix<STATE> &dx,TPZFMatrix<STATE> &d2x) {
//     const REAL &x0 = loc.Get(0,0);
//     const REAL &y0 = loc.Get(1,0);
// 
//     const int N=2;
//     //dx.Resize(N,1);
//     //d2x.Resize(N,N);
//     TFad<N,TFad<N,float> > xsecnd2,ysecnd2,xi,beta;
//     TFad<N,float> x(x0),y(y0);
// 
//     x.fastAccessDx(0)=1;
//     x.fastAccessDx(1)=0;
//     y.fastAccessDx(0)=0;
//     y.fastAccessDx(1)=1;
//     xsecnd2=x;
//     xsecnd2.fastAccessDx(0)=1;
//     xsecnd2.fastAccessDx(1)=0;
//     ysecnd2=y;
//     ysecnd2.fastAccessDx(0)=0;
//     ysecnd2.fastAccessDx(1)=1;
// 
//     xi=x;
//     xi.fastAccessDx(0)=1;
//     xi.fastAccessDx(1)=0;
//     beta=y;
//     beta.fastAccessDx(0)=0;
//     beta.fastAccessDx(1)=1;
// 
// // 	TFad<N,TFad<N,float> > sig1(79812.41),sig2(68385.69), sig3(61904.62),phi(20.*M_PI/180.f),c(490.f),nu(0.48),Pi(M_PI),um(1.f),dois(2.f),tres(3.f),seis(6.f),dezoito(18.f),sqrt3(sqrt(3.f));
//     TFad<N,TFad<N,float> > sig1(183753.31587364181),sig2(157264.29661659044), sig3(143880.63541092159),phi(20.*M_PI/180.f),c(490.f),nu(0.48),Pi(M_PI),um(1.f),dois(2.f),tres(3.f),seis(6.f),dezoito(18.f),sqrt3(sqrt(3.f));
//     TFad<N,TFad<N,float> >  test2;
// 
//     test2=((dois*(um-dois*nu))*(sqrt3*sig1 + sqrt3*sig2 +sqrt3*sig3 - tres*xi)*(sqrt3*sig1 + sqrt3*sig2 +sqrt3*sig3 - tres*xi) +(um + nu)*(tres*sig2 - tres*sig3 + (seis*sin(beta)*(-(sqrt3*c*cos(phi)) + xi*sin(phi)))/
//             (cos(beta)*(um + sin(phi)) - (-um + sin(phi))*sin(beta + Pi/seis)))*(tres*sig2 - tres*sig3 + (seis*sin(beta)*(-(sqrt3*c*cos(phi)) + xi*sin(phi)))/
//                     (cos(beta)*(um + sin(phi)) - (-um + sin(phi))*sin(beta + Pi/seis)))+tres*(um + nu)*(-dois*sig1 + sig2 + sig3 + (dois*cos(beta)*(tres*c*cos(phi) - sqrt3* xi *sin(phi)))/
//                             (cos( beta )*(um+ sin(phi)) - (-um + sin(phi))*sin( beta  + Pi/seis)))*(-dois*sig1 + sig2 + sig3 + (dois*cos(beta)*(tres*c*cos(phi) - sqrt3* xi *sin(phi)))/
//                                     (cos( beta )*(um+ sin(phi)) - (-um + sin(phi))*sin( beta  + Pi/seis))))/dezoito;
// 
//     dx(0,0)=test2.val().fastAccessDx(0);
//     dx(1,0)=test2.val().fastAccessDx(1);
// 
//     d2x(0,0)=test2.fastAccessDx(0).fastAccessDx(0);
//     d2x(1,1)= test2.fastAccessDx(1).fastAccessDx(1);
//     d2x(1,0)= test2.fastAccessDx(0).fastAccessDx(1);
//     d2x(0,1)= test2.fastAccessDx(1).fastAccessDx(0);

    
        const REAL &x0 = loc[0];
    const REAL &y0 = loc[1];

    const int N=2;
    dx.Resize(N,1);
    d2x.Resize(N,N);
    TFad<N,TFad<N,float> > xsecnd2,ysecnd2,test2;
    TFad<N,float> x(x0),y(y0);

    x.fastAccessDx(0)=1;
    x.fastAccessDx(1)=0;
    y.fastAccessDx(0)=0;
    y.fastAccessDx(1)=1;
    xsecnd2=x;
    xsecnd2.fastAccessDx(0)=1;
    xsecnd2.fastAccessDx(1)=0;
    ysecnd2=y;
    ysecnd2.fastAccessDx(0)=0;
    ysecnd2.fastAccessDx(1)=1;
    test2 = (-exp(ysecnd2) + xsecnd2/ysecnd2)*(-exp(ysecnd2) + xsecnd2/ysecnd2 + sin(xsecnd2/ysecnd2));
    
    cout << " aquiiiiiiiiii f = " << test2 << endl;

    dx(0,0)=test2.val().fastAccessDx(0);
    dx(1,0)=test2.val().fastAccessDx(1);

    d2x(0,0)=test2.fastAccessDx(0).fastAccessDx(0);
    d2x(1,1)= test2.fastAccessDx(1).fastAccessDx(1);
    d2x(1,0)= test2.fastAccessDx(0).fastAccessDx(1);
    d2x(0,1)= test2.fastAccessDx(1).fastAccessDx(0);
    
};

void FadTest()
{

    TPZFMatrix<REAL> pt(2,1);
    pt(0,0)=10000;
    pt(1,0)=0.8;
    TPZFMatrix<REAL> jac(2,2),dx(2,1),x0(2,1),xn(2,1);
    rhs2(pt,dx,jac);
    //jac.Print(cout);
    //dx.Print(cout);

   // xsecnd2 -> 1, ysecnd2 -> 0.546494
        x0(0,0)=1.;
    x0(1,0)=0.546494;
    rhs2(x0,dx,jac);
    
    
    std::cout << " dx = " << dx   << std::endl;
    
    return;
    
    for(int it=0; it<7; it++) {
        rhs2(x0,dx,jac);
// 		cout <<  "x0 = " << x0(0,0) << " "<<x0(1,0) <<  endl;
// 		cout <<  "dx = " << dx(0,0) << " "<<dx(1,0) <<  endl;
// 		cout <<  "jac = " << jac(0,0) << " "<<jac(0,1) <<  endl;
// 		cout <<  "jac = " << jac(1,0) << " "<<jac(1,1) <<  endl;
        jac.Solve_LU(&dx);
        x0-=dx;
        xn=x0;
        x0=xn;
        cout <<  "iter = " << it << "  norm = " << sqrt(dx(0,0)*dx(0,0)+dx(1,0)*dx(1,0)) << endl;
    }
    cout <<  "xn = " << xn << endl;
    xn.Print(cout);
}



void Post(TPZPostProcAnalysis * postproc,std::string vtkFile )
{

    //TPZVec<int> PostProcMatIds(2,0);
    // PostProcMatIds[0]=1;
    //PostProcMatIds[1]=2;

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables(scalNames, vecNames);

    postproc->DefineGraphMesh(2,scalNames,vecNames,vtkFile);

    postproc->PostProcess(0);
}

void IntegrateBCSol(TPZCompMesh * cmesh)
{

    //int nels0 = cmesh->NElements();
    int nels0 = cmesh->Reference()->NElements();
    REAL integralreal=0.;
    for ( int iel=0; iel<nels0; iel++ ) {

        TPZCompEl *cel = cmesh->Element ( iel );
        if(!cel)continue;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

        if(!intel)continue;
        if(cel->Material()->Id()!=1)//rigth
        {
            continue;
        }

        int var=66;//elasticenergy
        TPZManVector<REAL,3> integral;
        integral = intel->IntegrateSolution(var);

        integralreal+=integral[0];
    }
    cout << "elastic energy integral over the rigth edge boundary is = " << integralreal << endl;
}


STATE ComputeBCWork(TPZCompMesh &cmesh)
{
    STATE work=0.,nx=0.,ny=0.,radius=0.;
    TPZManVector<REAL,3> ksi(3,0.),sol(3,0.),xvec(3,0.);
    TPZFMatrix<STATE> jac,jacinv,axes,invjac;
    STATE detjac;

    int nel = cmesh.NElements();
    cmesh.Reference()->ResetReference();

    for (int el=0; el<nel; el++)
    {

        TPZCompEl * cel = cmesh.ElementVec()[el];
        if(!cel)
        {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if(!gel)
        {
            DebugStop();
        }
        int matid = gel->MaterialId();
        if(matid!=-2)
        {
            continue;
        }

//         TPZGeoElSide gelside(gel,2);
//         TPZStack<TPZCompElSide> stackcompelside;
//         gelside.EqualLevelCompElementList(stackcompelside, 0, 0);
//         if(stackcompelside.size()!=1)
//         {
//             DebugStop();
//         }
//
//         TPZTransform t1(1);
//         TPZCompElSide compneigh = stackcompelside[0];
//         TPZGeoElSide neighbour = stackcompelside[0].Reference();
//         gelside.SideTransform3(neighbour,t1);
//         int sidefrom,sideto = neighbour.Element()->NSides()-1;
//         sidefrom=neighbour.Side();
//         TPZTransform t2 = neighbour.Element()->SideToSideTransform(sidefrom, sideto);


        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        const TPZIntPoints &rule = intel->GetIntegrationRule();
        REAL weight;
        int npoints = rule.NPoints();

        TPZMaterialData data;
        data.fNeedsSol = true;
        intel->InitMaterialData ( data );

        for(int i=0; i<npoints; i++)
        {
            rule.Point(i, ksi, weight);

            TPZManVector<REAL,3> ksi1(1,0.), ksi2(2,0.);
//             t1.Apply(ksi, ksi1);
//             t2.Apply(ksi1, ksi2);
//             compneigh.Element()->Solution(ksi2, 1, sol);

            gel->X(ksi, xvec);
            intel->GetIntegrationRule().Point(0, ksi, weight);

            gel->Jacobian(xvec,jac,axes,detjac,invjac);

            int ETotStressXX = 30;
            int ETotStressYY = 31;
            int ETotStressZZ = 32;
            TPZMaterial *mat= cel->Material();
            TPZVec<REAL> sxx,syy;
            mat->Solution(data,ETotStressXX,sxx);
            mat->Solution(data,ETotStressYY,syy);
            work +=(sxx[0]*1*sol[0]+syy[0]*0*sol[1])*weight*fabs(detjac);

        }

    }
#ifdef DEBUG
    cout << "\n WORK "<< work << endl;
#endif
    return work;

}

REAL findnodalsol(TPZCompMesh *cmesh) {

    //b[0] = 39.2265;
    //b[1] = 40.;
    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension();
    TPZVec<REAL> xd(dim,0.);
    TPZVec<REAL> mpt(dim,0.);
    xd[0] =0.121;
    xd[1] = 0.3;
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


            cout << "elemento encontrado"<<endl;
            cout << "index do elemento = " <<gel->Index() <<endl;
            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(3);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                if(fabs(co[0]-xd[0])<1.e-6 && fabs(co[1]-xd[1])<1.e-6)
                {
                   // cout << "node id = "<<id <<endl;
                   // cout << " Coordinates "<< endl;
                   // cout << " x = "<< co[0] << ",  y = " << co[1] << endl;
                  //  cout << " qsi = "<< qsi << endl;
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

                    TPZMaterialData data;

                    data.fNeedsSol = true;

                    intel->InitMaterialData ( data );

                    intel->ComputeRequiredData ( data, qsi );

                //    cout << " data.sol  = "<<data.sol     << endl;
                    return data.sol[0][1];
                }
            }


        }
    }

}


// TPZGeoMesh * CreateGMeshGid ( int ref )
// {
//
//
// 	//string file ="/home/diogo/projects/pz/data/mesh-teste-pz-fromathematica.msh";
//     string file ="/home/diogo/projects/pz/data/quad-gid.msh";
// 	//string file ="/home/diogo/projects/pz/data/gid-tri.msh";
//
//
//
// 	readgidmesh read = readgidmesh(file);
// 	read.ReadMesh();
// 	TPZFMatrix<int> meshtopology = read.GetTopology();
// 	TPZFMatrix<REAL> meshcoords = read.GetCoords();
// 	std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();
//
//     int ndivs = 10000;
//     TPZFMatrix<REAL> pathbottom, pathleft, pathright, pathdisplace;
//     std::vector<int>  idsbottom, idsleft, idsright, idstoprigth,idsramp,idstopleft;
//
//     std::vector<std::vector<int>> idsvec;
//
//
//     TPZManVector<REAL,2> a ( 2 ), b ( 2 );
//
//
//     a[0] = 0.; a[1] = 0.;
//     b[0] = 75.; b[1] = 0.;
//     read.Line( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsbottom );
//     idsvec.push_back(idsbottom);
//
//     a = b;
//     b[0] = 75.; b[1] = 30.;
//     read.Line( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsright );
//     idsvec.push_back(idsright);
//
//
//     a = b;
//     b[0] = 45.; b[1] = 30.;
//     read.Line( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idstoprigth );
//     idsvec.push_back(idstoprigth);
//
//     a = b;
//     b[0] = 35.; b[1] = 40.;
//     read.Line( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsramp );
//     idsvec.push_back(idsramp);
//
//     a = b;
//     b[0] = 0.; b[1] = 40.;
//     read.Line( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idstopleft );
//     idsvec.push_back(idstopleft);
//
//     a = b;
//     b[0] = 0.; b[1] = 0.;
//     read.Line( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsleft );
//     idsvec.push_back(idsleft);
//
//
//    // for(int i=0;i<idsbottom.size();i++)cout << idsvec[3][i] << endl;
//
//     const std::string name ( "Slope Problem " );
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
// 	int ncoords = meshcoords.Rows();
// 	co.resize(ncoords);
// 	for(int i=0;i<ncoords;i++)
//  	{
// 		co[i].resize(2);
// 		co[i][0]=meshcoords(i,0);
// 		co[i][1]=meshcoords(i,1);
// 	}
//     vector<vector<int>> topol;
//
// 	int ntopol = meshtopology.Rows();
// 	topol.resize(ntopol);
//
// 	for(int i=0;i<ntopol;i++)
//  	{
// 		topol[i].resize(meshtopology.Cols());
// 		for(int j=0;j<meshtopology.Cols();j++)
//   		{
// 	  		topol[i][j]=meshtopology(i,j);
// 		}
// 	}
//
//     gmesh->NodeVec().Resize ( co.size() );
//
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     if(meshtopology.Cols()==4)
// 	{
// 		TPZVec <long> TopoQuad ( 4 );
//     	for ( int iel=0; iel<topol.size(); iel++ ) {
//         	TopoQuad[0] = topol[iel][0];
//         	TopoQuad[1] = topol[iel][1];
//         	TopoQuad[2] = topol[iel][2];
//         	TopoQuad[3] = topol[iel][3];
//         	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//     	}
// 	}
//
// 	if(meshtopology.Cols()==3)
// 	{
// 		TPZVec <long> TopoTri ( 3 );
//     	for ( int iel=0; iel<topol.size(); iel++ ) {
//         	TopoTri[0] =topol[iel][0];
//         	TopoTri[1] =topol[iel][1];
//         	TopoTri[2] =topol[iel][2];
//         	new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
//     	}
// 	}
//
// 	if(meshtopology.Cols()!=3 && meshtopology.Cols()!=4)
//  	{
// 		DebugStop();
// 	}
//     TPZVec <long> TopoLine ( 2 );
//     int id = topol.size();
//     id++;
//
//     for(int ivec=0;ivec<idsvec.size();ivec++)
//     {
//         int nodes = idsvec[ivec].size();
//         for(int inode=0;inode<nodes-1;inode++)
//         {
//             TopoLine[0] = idsvec[ivec][inode];
//             TopoLine[1] = idsvec[ivec][inode+1];
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -(ivec+1), *gmesh );
//         }
//     }
//
//
//     gmesh->BuildConnectivity();
//     for ( int d = 0; d<ref; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//     gmesh->Print(std::cout);
//  	std::ofstream files ( "ge.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK(gmesh,files,false);
//
//     return gmesh;
// }

// TPZGeoMesh *  CreateGMesh ( int ref )
// {
//     const std::string name ( "Slope Problem " );
//
//     TPZGeoMesh *gmesh  =  new TPZGeoMesh();
//
//     gmesh->SetName ( name );
//     gmesh->SetDimension ( 2 );
//
//     TPZVec<REAL> coord ( 2 );
//
//     vector<vector<double>> co= {{0., 0.}, {75., 0.}, {75., 30.}, {45., 30.}, {35., 40.}, {0.,40.},
//
// 	{35./3., 40.},{2 * 35/3., 40.}, {30., 40.},
//
// 	{30., 30.}, {60.,30.},
//
// 	{2* 35./3.,2* 35/3.}, {45., 2* 35/3.},
//
// 	 {35./3., 35/3.}, {60., 35./3.}
//
//     };
//     vector<vector<int>> topol = {
// 		{0,  1,  14, 13},
// 		{1,  2,  10, 14},
// 		{14, 10, 3,  12},
// 		{13, 14, 12, 11},
// 		{11, 12, 3,  9},
// 		{9,  3,  4,  8},
// 		{11, 9,  8,  7},
// 		{13, 11, 7, 6},
// 		{0, 13,  6, 5}
//     };
//
//
//     gmesh->NodeVec().Resize ( co.size() );
//
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     TPZVec <long> TopoQuad ( 4 );
//     for ( int iel=0; iel<topol.size(); iel++ ) {
//         TopoQuad[0] = topol[iel][0];
//         TopoQuad[1] = topol[iel][1];
//         TopoQuad[2] =	topol[iel][2];
//         TopoQuad[3] = topol[iel][3];
//         new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//     }
//     int id = topol.size();
//     TPZVec <long> TopoLine ( 2 );
//     TopoLine[0] = 0;
//     TopoLine[1] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 1;
//     TopoLine[1] = 2;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth
//
// 	{
//         id++;
//         TopoLine[0] = 3;
//         TopoLine[1] = 10;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//
//         id++;
//         TopoLine[0] = 10;
//         TopoLine[1] = 2;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//     }
//
// 	id++;
//     TopoLine[0] = 3;
//     TopoLine[1] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh );//ramp
//
//
//
//     {
//             id++;
//             TopoLine[0] = 5;
//             TopoLine[1] = 6;
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
//
//             id++;
//             TopoLine[0] = 6;
//             TopoLine[1] = 7;
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
//
//             id++;
//             TopoLine[0] = 7;
//             TopoLine[1] = 8;
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
//
//             id++;
//             TopoLine[0] = 8;
//             TopoLine[1] = 4;
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //top left
//     }
//
//
//
//
//     id++;
//     TopoLine[0] = 0;
//     TopoLine[1] = 5;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//left
//
//
//     gmesh->BuildConnectivity();
//     for ( int d = 0; d<ref; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//
//     return gmesh;
// }


TPZGeoMesh *  CreateGMesh ( int ref )
{
    const std::string name ( "Slope Problem " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

    vector<vector<double>> co= {{0., 0.}, {7500., 0.}, {7500., 3000.}, {4500., 3000.}, {3500., 4000.}, {0.,4000.},

	{3500./3., 4000.},{2 * 3500/3., 4000.}, {3000., 4000.},

	{3000., 3000.}, {6000.,3000.},

	{2* 3500./3.,2* 3500/3.}, {4500., 2* 3500/3.},

	 {3500./3., 3500/3.}, {6000., 3500./3.}

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



/// create the computational mesh
// TPZCompMesh * CreateCMesh(TPZGeoMesh * gmesh,int porder)
// {
// 	unsigned int dim  = 2;
//     const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );
//
// 	TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
//     cmesh->SetName ( name );
//     cmesh->SetDefaultOrder ( porder );
//     cmesh->SetDimModel ( dim );
//
//     // Mohr Coulomb data
//     REAL mc_cohesion    = 10.0;
//     REAL mc_phi         = ( 30.0*M_PI/180 );
//     REAL mc_psi         = mc_phi;
//
//     /// ElastoPlastic Material using Mohr Coulomb
//     // Elastic predictor
//     TPZElasticResponse ER;
//     REAL nu = 0.49;
//     REAL E = 20000.;
//
//     TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
//     ER.SetUp( E, nu );
// 	LEMC.fER =ER;
//    // LEMC.SetElasticResponse( ER );
//     LEMC.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );
//
//     int PlaneStrain = 1;
//
// 	plasticmat * material = new plasticmat ( 1,PlaneStrain );
//     material->SetPlasticity ( LEMC );
//
// 	material->SetId(1);
//     cmesh->InsertMaterialObject ( material );
//
//     TPZFMatrix<STATE> val1 ( 2,2,0. );
//     TPZFMatrix<STATE>  val2 ( 2,1,0. );
// 	int directionadirichlet =3;
//     val2 ( 0,0 ) = 1;
//     val2 ( 1,0 ) = 1;
//     auto * bc_bottom = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );//bottom
//     val2 ( 0,0 ) = 1;
//     val2 ( 1,0 ) = 0;
//     auto * bc_rigth = material->CreateBC ( material, -2, directionadirichlet, val1, val2 );//rigth
//     val2 ( 0,0 ) = 1;
//     val2 ( 1,0 ) = 0;
//     auto * bc_left = material->CreateBC ( material, -6, directionadirichlet, val1, val2 );//left
//
//     cmesh->InsertMaterialObject ( bc_bottom );
//    	cmesh->InsertMaterialObject ( bc_rigth );
//     cmesh->InsertMaterialObject ( bc_left );
//
//     //cmesh->InsertMaterialObject ( top );
//     cmesh->SetAllCreateFunctionsContinuousWithMem();
//
//     cmesh->AutoBuild();
//
// 	return cmesh;
//
// }

void LoadingRamp (TPZCompMesh * cmesh,  REAL factor,REAL gamma)
{
    plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=-factor*gamma;
    body->SetBodyForce ( force );

}

void LoadingRamp (TPZCompMesh * cmesh,  REAL factor)
{
    plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=-factor*20.;
    body->SetBodyForce ( force );

}

void LoadingRamptwo ( TPZCompMesh * cmesh,  REAL factor )
{
    plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    plasticmat * body2= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 2 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=-factor*20.;
    body->SetBodyForce ( force );
    body2->SetBodyForce ( force );

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


void  CreatePostProcessingMeshTwoMats(TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh)
{

    if (PostProcess->ReferenceCompMesh() != cmesh)
    {

        PostProcess->SetCompMesh(cmesh);

        TPZVec<int> PostProcMatIds(2);
        PostProcMatIds[0]=1;
        PostProcMatIds[1]=2;
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





// TPZManVector<REAL,10> output(10,0.);
//     REAL FS=0.1,FSmax=50.,FSmin=0.,tol=0.01;
//     int neq = cmesh->NEquations();
//     int maxcount=50;
//     TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );
// 
//     int counterout = 0;
// 
//     REAL norm = 1000.;
//      REAL tol2 = 1.e-2;
//     int NumIter = 50;
//     bool linesearch = true;
//     bool checkconv = false;
// 	std::ofstream outloadu("outloadu.nb");
//     
// 
//     do {
// 
// 
//         cout << "\n  FS = " << FS  <<" | Load step = " << counterout;
//         LoadingRamp ( cmesh,FS );
// 		//LoadingRamp ( body,0.0000001 );
//         bool optimize =true;
//         TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh );
//         chrono::steady_clock sc;
//         auto start = sc.now();
// 		int iters;
//         bool conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv ,iters);
// 



void ShearRed ( TPZCompMesh * cmesh)
{
    
    LoadingRamp(cmesh,1.);

    REAL FS=0.8,FSmax=5.,FSmin=0.,tol=0.001;
    int neq = cmesh->NEquations();

    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    bool conv=false;
    do {

        cmesh->Solution().Zero();
        std::cout << "FS "<< FS <<  "| step = " << counterout  <<std::endl;
        TPZElastoPlasticAnalysis  * anal = CreateAnal(cmesh);
                
        REAL norm = 1000.;
        REAL tol2 = 0.01;
        int NumIter = 50;
        bool linesearch = true;
        bool checkconv = false;
        int iters;
        
        
        
        chrono::steady_clock sc;
        auto start = sc.now();

         conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters );
        norm = Norm ( anal->Rhs() );


        //anal->AcceptSolution();


        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );

        std::cout <<" | Rhs norm = " << norm  << " | IterativeProcess Time: " << time_span.count() << " seconds !!! " <<std::endl;

        if ( conv==false ) {
	
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {
			
            FSmin = FS;
			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }
        
        c=cohesion0/FS;
        phi=atan ( tan ( phi0 ) /FS );
        psi=phi;
        LEMC.fYC.SetUp ( phi, psi, c, ER );
        material->SetPlasticity ( LEMC );
        counterout++;
    }  while ( ( FSmax - FSmin ) / FS > tol || conv==false);
    TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh );
	anal->AcceptSolution();
}

TPZElastoPlasticAnalysis * CreateAnal(TPZCompMesh *cmesh)
{
    int numthreads=10;

    TPZElastoPlasticAnalysis * analysis =  new TPZElastoPlasticAnalysis( cmesh); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );

    matskl.SetNumThreads ( numthreads );

    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    analysis->SetSolver ( step );

    return analysis;
}


void DivideElementsAbove(TPZCompMesh*fcmesh,TPZGeoMesh*fgmesh,REAL sqj2, std::set<long> &elindices)
{
    fgmesh->ResetReference();
    fcmesh->LoadReferences();
    TPZManVector<REAL,3> findel(3,0.),qsi(2,0.);
    findel[0] = 0.108;
    findel[1] = 0.0148;
    long elindex = 0;
    fcmesh->Reference()->FindElement(findel, qsi, elindex, 2);
    TPZGeoEl *targetel = fcmesh->Reference()->ElementVec()[elindex];
    TPZCompEl *targetcel = targetel->Reference();
    long targetindex = targetcel->Index();

    long nelem = fcmesh->NElements();
    cout << "elements antes " << nelem<<std::endl;
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fcmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        int investigate = false;
        if (el == targetindex) {
            std::cout << "I should investigate\n";
            investigate = true;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            DebugStop();
        }
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            continue;
        }
        if (fcmesh->ElementSolution()(el,0) < sqj2) {
            continue;
        }
        int porder = intel->GetPreferredOrder();
        TPZStack<long> subels;
        long index = cel->Index();

        intel->Divide(index, subels,0);
        for (int is=0; is<subels.size(); is++) {
            elindices.insert(subels[is]);
            TPZCompEl *subcel = fcmesh->ElementVec()[subels[is]];
            TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
            if (!subintel) {
                DebugStop();
            }
            subintel->SetPreferredOrder(porder);
        }
    }
    // divide elements with more than one level difference
    bool changed = true;
    while (changed) {
        changed = false;
        std::set<long> eltodivide;
        long nelem = fcmesh->NElements();
        for (long el=0; el<nelem; el++) {
            TPZCompEl *cel = fcmesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZGeoEl *gel = cel->Reference();
            if (!gel) {
                DebugStop();
            }
            int ns = gel->NSides();
            for (int is=0; is<ns; is++) {
                TPZGeoElSide gelside(gel, is);
                if (gelside.Dimension() != 1) {
                    continue;
                }
                TPZCompElSide big = gelside.LowerLevelCompElementList2(1);
                if (!big) {
                    continue;
                }
                TPZGeoElSide geobig(big.Reference());
                // boundary elements will be refined by AdjustBoundaryElements
                if (geobig.Element()->Dimension() != 2) {
                    continue;
                }
                if (gel->Level()-geobig.Element()->Level() > 1) {
                    eltodivide.insert(big.Element()->Index());
                }
            }
        }
        std::set<long>::iterator it;
        for (it = eltodivide.begin(); it != eltodivide.end(); it++) {
            changed = true;
            long el = *it;
            TPZCompEl *cel =fcmesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
            if (!pMatWithMem2) {
                continue;
            }
            int porder = intel->GetPreferredOrder();
            TPZStack<long> subels;
            long index = cel->Index();
            intel->Divide(index, subels,0);
            for (int is=0; is<subels.size(); is++) {
                elindices.insert(subels[is]);
                TPZCompEl *subcel = fcmesh->ElementVec()[subels[is]];
                TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
                if (!subintel) {
                    DebugStop();
                }
                subintel->SetPreferredOrder(porder);
            }
        }
    }
    fcmesh->AdjustBoundaryElements();
    fcmesh->InitializeBlock();
    fcmesh->Solution().Zero();
    nelem = fcmesh->NElements();
    cout << "elements depois " << nelem<<std::endl;
    std::cout << "Number of elements prefined: " << elindices.size() << std::endl;
}

// Get the vector of element plastic deformations
void ComputeElementDeformation(TPZCompMesh * fcmesh)
{
    long nelem = fcmesh->NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.);
    fcmesh->ElementSolution().Redim(nelem, 1);
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fcmesh->MaterialVec()[1]);
    if (!pMatWithMem2) {
        fPlasticDeformSqJ2.Fill(0.);
    }
    else
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fcmesh->ElementVec()[el];
            fPlasticDeformSqJ2[el] = 0.;
            if (!cel) {
                continue;
            }
            TPZManVector<long> memindices;
            cel->GetMemoryIndices(memindices);
            int numind = memindices.size();
            REAL sqj2el = 0.;
            for (int ind=0; ind<numind; ind++)
            {
                int memoryindex = memindices[ind];
                if (memoryindex < 0) {
                    continue;
                }
                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ind]);
                TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                REAL J2 = plastic.J2();
                REAL sqj2 = sqrt(J2);
                sqj2el = max(sqj2,sqj2el);
            }
            fPlasticDeformSqJ2[el] = sqj2el;
        }
    }


    fcmesh->SetElementSolution(0, fPlasticDeformSqJ2);
}



TPZGeoMesh * CreateGMeshGid ( int ref ,TPZFMatrix<REAL> pts, string gridname)
{
    
    string file = gridname;

    readgidmesh read = readgidmesh ( file );
    read.ReadMesh();
    TPZFMatrix<int> meshtopology = read.GetTopology();
    TPZFMatrix<REAL> meshcoords = read.GetCoords();
    std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();

    int ndivs = 100000;
    TPZFMatrix<REAL> pathbottom, pathleft, pathright, pathdisplace;
    std::vector<int>  idsbottom, idsleft, idsright, idstoprigth,idsramp,idstopleft;

    std::vector<std::vector<int>> idsvec;


    TPZVec<double> constcoorddata(3,0.);
    TPZVec<int> constcoord(3);



    // (linha y=0)
    constcoorddata[0]=pts(0,0);
    constcoorddata[1]=pts(0,0);
    constcoorddata[2]=0;
    constcoord[0]=0;
    constcoord[1]=1;//direceo  fixa
    constcoord[2]=1;//direceo  fixa
    read.FindIds ( constcoorddata,constcoord, idsbottom );
    
    idsvec.push_back ( idsbottom );

    cout << "pts(1,0) = " <<pts(1,0) << endl;
    cout << "pts(1,1) = " <<pts(1,1) << endl;
    // (linha x=pts(1,0) )
    constcoorddata[0]=pts(1,0);
    constcoorddata[1]=pts(1,1);
    constcoorddata[2]=0;
    constcoord[0]=1;//direceo  fixa
    constcoord[1]=0;
    constcoord[2]=1;//direceo  fixa
    read.FindIds ( constcoorddata,constcoord, idsright );
    
    idsvec.push_back ( idsright );
    
    // (linha x=0 )
    constcoorddata[0]=0;
    constcoorddata[1]=0;
    constcoorddata[2]=0;
    constcoord[0]=1;//direceo  fixa
    constcoord[1]=0;
    constcoord[2]=1;//direceo  fixa
    read.FindIds ( constcoorddata,constcoord, idsleft );
    

    idsvec.push_back ( idsleft );
    
    cout << " *** ids *** " << endl;
    for(int j=0;j<idsvec.size();j++)
    {
        cout << " *** j *** "<< j  << endl;
        for(int i=0;i<idsvec[j].size();i++)cout << idsvec[j][i] << endl;
    }
    
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
    std::ofstream files ( "ge2.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}


// TPZGeoMesh * CreateGMeshGid ( int ref ,TPZFMatrix<REAL> pts, string gridname)
// {
//     
//     string file = gridname;
// 
//     readgidmesh read = readgidmesh ( file );
//     read.ReadMesh();
//     TPZFMatrix<int> meshtopology = read.GetTopology();
//     TPZFMatrix<REAL> meshcoords = read.GetCoords();
//     std::vector<std::vector< std::vector<double > > > allcoords = read.GetAllCoords();
// 
//     int ndivs = 100000;
//     TPZFMatrix<REAL> pathbottom, pathleft, pathright, pathdisplace;
//     std::vector<int>  idsbottom, idsleft, idsright, idstoprigth,idsramp,idstopleft;
// 
//     std::vector<std::vector<int>> idsvec;
// 
// 
//     TPZManVector<REAL,2> a ( 2 ), b ( 2 );
// 
//     REAL delta=fabs((pts(0,0)-pts(1,0)))/(ndivs/10);
//    // REAL delta=fabs((pts(0,0)-pts(1,0)))/ndivs;
//     
//     cout << "delta = " << delta <<endl; 
// 
//     a[0] = pts(0,0);
//     a[1] = pts(0,1);
//     b[0] = pts(1,0);
//     b[1] = pts(1,1);
//     
//     
// 
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsbottom ,delta);
//     idsvec.push_back ( idsbottom );
//     
//     
//     TPZVec<double> constcoorddata(3,0.);
//     TPZVec<int> constcoord(3);
//     std::vector<int> ids1,ids2,ids3,ids4,ids5,ids6,ids7;
// 
// 
// 
//     // (linha y=0)
//     constcoorddata[0]=0;
//     constcoorddata[1]=0;
//     constcoorddata[2]=0;
//     constcoord[0]=0;
//     constcoord[1]=1;//direceo  fixa
//     constcoord[2]=1;//direceo  fixa
//     read.FindIds ( constcoorddata,constcoord, ids1 );
//     
//     cout << " *** ids *** " << endl;
//     for(int j=0;j<ids1.size();j++)
//     {
//         cout << ids1[j] << endl;
//     }
//     
//     a = b;
//     b[0] = pts(2,0);
//     b[1] = pts(2,1);
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsright,delta );
//     idsvec.push_back ( idsright );
// 
//     a = b;
//     b[0] = pts(3,0);
//     b[1] = pts(3,1);
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idstoprigth,delta );
//     idsvec.push_back ( idstoprigth );
// 
//     a = b;
//     b[0] = pts(4,0);
//     b[1] = pts(4,1);
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsramp,delta );
//     idsvec.push_back ( idsramp );
// 
//     a = b;
//     b[0] = pts(5,0);
//     b[1] = pts(5,1);
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idstopleft,delta );
//     idsvec.push_back ( idstopleft );
// 
//     a = b;
//     b[0] = pts(0,0);
//     b[1] = pts(0,1);
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsleft,delta );
//     idsvec.push_back ( idsleft );
// 
// 
//     cout << " *** ids *** " << endl;
//     for(int j=0;j<idsvec.size();j++)
//     {
//         cout << " *** j *** "<< j  << endl;
//         for(int i=0;i<idsvec[j].size();i++)cout << idsvec[j][i] << endl;
//     }
// 
//     const std::string name ( "Slope Problem " );
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
//     int ncoords = meshcoords.Rows();
//     co.resize ( ncoords );
//     for ( int i=0; i<ncoords; i++ ) {
//         co[i].resize ( 2 );
//         co[i][0]=meshcoords ( i,0 );
//         co[i][1]=meshcoords ( i,1 );
//     }
//     vector<vector<int>> topol;
// 
//     int ntopol = meshtopology.Rows();
//     topol.resize ( ntopol );
// 
//     for ( int i=0; i<ntopol; i++ ) {
//         topol[i].resize ( meshtopology.Cols() );
//         for ( int j=0; j<meshtopology.Cols(); j++ ) {
//             topol[i][j]=meshtopology ( i,j );
//         }
//     }
// 
//     gmesh->NodeVec().Resize ( co.size() );
// 
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     if ( meshtopology.Cols() ==4 ) {
//         TPZVec <long> TopoQuad ( 4 );
//         for ( int iel=0; iel<topol.size(); iel++ ) {
//             TopoQuad[0] = topol[iel][0];
//             TopoQuad[1] = topol[iel][1];
//             TopoQuad[2] = topol[iel][2];
//             TopoQuad[3] = topol[iel][3];
//             new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//         }
//     }
// 
//     if ( meshtopology.Cols() ==3 ) {
//         TPZVec <long> TopoTri ( 3 );
//         for ( int iel=0; iel<topol.size(); iel++ ) {
//             TopoTri[0] =topol[iel][0];
//             TopoTri[1] =topol[iel][1];
//             TopoTri[2] =topol[iel][2];
//             new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
//         }
//     }
// 
//     if ( meshtopology.Cols() !=3 && meshtopology.Cols() !=4 ) {
//         DebugStop();
//     }
//     TPZVec <long> TopoLine ( 2 );
//     int id = topol.size();
//     id++;
// 
//     for ( int ivec=0; ivec<idsvec.size(); ivec++ ) {
//         int nodes = idsvec[ivec].size();
//         for ( int inode=0; inode<nodes-1; inode++ ) {
//             TopoLine[0] = idsvec[ivec][inode];
//             TopoLine[1] = idsvec[ivec][inode+1];
//             new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, - ( ivec+1 ), *gmesh );
//         }
//     }
// 
// 
//     gmesh->BuildConnectivity();
//     for ( int d = 0; d<ref; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//     // gmesh->Print(std::cout);
//     std::ofstream files ( "ge2.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
// 
//     return gmesh;
// }

TPZGeoMesh * CreateGMeshCentrifuga ( int ref )
{

    string file ="/home/diogo/projects/pz/data/gi.msh";

    REAL delta =1.e-6;

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


    b[1] = 0.;
    a[0] = 0.;
    a[1] = 0.;
    b[0] = 0.300;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsbottom,delta );
    idsvec.push_back ( idsbottom );

    a = b;
    b[0] = 0.300;
    b[1] = 0.300;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsright,delta );
    idsvec.push_back ( idsright );

    a = b;
    b[0] = 0.121;
    b[1] = 0.300;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstoprigth,delta );
    idsvec.push_back ( idstoprigth );

    a = b;
    b[0] = 0.100;
    b[1] = 0.060;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsramp,delta );
    idsvec.push_back ( idsramp );

    a = b;
    b[0] = 0.;
    b[1] = 0.060;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstopleft,delta );
    idsvec.push_back ( idstopleft );

    a = b;
    b[0] = 0.;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsleft,delta );
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
    std::ofstream files ( "centriguga.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}



TPZGeoMesh * CreateGMeshGidTwoMats ( int ref )
{

    REAL delta =1.e-6;

    string file ="/home/diogo/projects/pz/data/gid-tri-two-mats.msh";



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
    read.FindIdsInPath ( pathbottom, idsbottom,delta );
    idsvec.push_back ( idsbottom );

    a = b;
    b[0] = 75.;
    b[1] = 30.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsright,delta );
    idsvec.push_back ( idsright );


    a = b;
    b[0] = 45.;
    b[1] = 30.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstoprigth,delta );
    idsvec.push_back ( idstoprigth );
////h10-beta 45
    a = b;
    b[0] = 35.;
    b[1] = 40.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsramp ,delta);
    idsvec.push_back ( idsramp );

//h10-beta 30
// 	a = b;
//     b[0] = 27.675;
//     b[1] = 40.;
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsramp );
//     idsvec.push_back ( idsramp );
//h10-beta 60
// 	a = b;
//     b[0] = 39.2265;
//     b[1] = 40.;
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsramp );
//     idsvec.push_back ( idsramp );


    a = b;
    b[0] = 0.;
    b[1] = 40.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstopleft,delta );
    idsvec.push_back ( idstopleft );

    a = b;
    b[0] = 0.;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsleft,delta );
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

        std::cout << "quadrilaterals are not supported/need to be implemented!" <<std::endl;
        DebugStop();
    }

    if ( meshtopology.Cols() ==3 ) {
        TPZVec <long> TopoTri ( 3 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
            bool region =false;
            for(int inode=0; inode<3; inode++)
            {
                TopoTri[inode]=topol[iel][inode];
                REAL xco = meshcoords(TopoTri[inode],0);
                REAL yco = meshcoords(TopoTri[inode],1);
                REAL zco = meshcoords(TopoTri[inode],2);
                if(yco>30.)
                {
                    region=true;
                }
            }
            if(region)
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
            } else
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 2,*gmesh );
            }

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
    std::ofstream files ( "teste-diff-mats-geomesh.vtk" );
    std::set<int> myset;
    //TPZVTKGeoMesh::PrintGMeshVTKmy_material(gmesh,files,myset,true);
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}



/// create the computational mesh
TPZCompMesh * CreateCMeshTwoMats ( TPZGeoMesh * gmesh,int porder )
{
    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    // Mohr Coulomb data
    REAL mc_cohesion    = 50.;
    REAL mc_phi         = ( 20.*M_PI/180 );
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER,ER2;
    REAL nu = 0.49;
    REAL E = 20000.;

    REAL nu2=0.49;
    REAL E2=20000;

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC,LEMC2;
    ER.SetUp ( E, nu );
    ER2.SetUp ( E2, nu2 );
    LEMC.fER =ER;
    LEMC2.fER =ER2;
    // LEMC.SetElasticResponse( ER );
    LEMC.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    mc_cohesion=50;
    mc_phi = ( 20.*M_PI/180 );
    mc_psi=mc_phi;
    LEMC2.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER2 );
    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );
    plasticmat * material2 = new plasticmat ( 2,PlaneStrain );
    material->SetPlasticity ( LEMC );
    material2->SetPlasticity(LEMC2);

    material->SetId ( 1 );
    material2->SetId( 2 );
    cmesh->InsertMaterialObject ( material );
    cmesh->InsertMaterialObject ( material2 );

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
    int directionadirichlet =3;
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 1;
    auto * bc_bottom = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );//bottom
    auto * bc_bottom2 = material2->CreateBC ( material2, -1,directionadirichlet, val1, val2 );//bottom
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_rigth = material->CreateBC ( material, -2, directionadirichlet, val1, val2 );//rigth
    auto * bc_rigth2 = material2->CreateBC ( material2, -2, directionadirichlet, val1, val2 );//rigth

    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_left = material->CreateBC ( material, -6, directionadirichlet, val1, val2 );//left
    auto * bc_left2 = material2->CreateBC ( material2, -6, directionadirichlet, val1, val2 );//left

    cmesh->InsertMaterialObject ( bc_bottom );
    cmesh->InsertMaterialObject ( bc_rigth );
    cmesh->InsertMaterialObject ( bc_left );

    cmesh->InsertMaterialObject ( bc_bottom2 );
    cmesh->InsertMaterialObject ( bc_rigth2 );
    cmesh->InsertMaterialObject ( bc_left2 );

    //cmesh->InsertMaterialObject ( top );
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();
    std::ofstream files ( "teste-diff-mats-compmesh.vtk" );
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh,  files, true);

    return cmesh;

}


/// create the computational mesh
TPZCompMesh * CreateCMeshCentrifuga ( TPZGeoMesh * gmesh,int porder )
{
    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    // Mohr Coulomb data
    REAL mc_cohesion    = 40;//MPa
    REAL mc_phi         = ( 32.*M_PI/180 );
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.32;
    REAL E = 50000.; // 50 MPa

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetUp ( E, nu );
    LEMC.fER =ER;
    // LEMC.SetElasticResponse( ER );
    LEMC.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );
    material->SetPlasticity ( LEMC );

    material->SetId ( 1 );
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



/// create the computational mesh
TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder, REAL c, REAL phi, REAL poisson, REAL young )
{
    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    // Mohr Coulomb data
    REAL mc_cohesion    = c;
    REAL mc_phi         = (phi*M_PI/180 );
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = poisson;
    REAL E = young; // 2 000 000 kPa = 2 GPa, concreto tem aprox 30 Gpa

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetUp ( E, nu );
    LEMC.fER =ER;
    // LEMC.SetElasticResponse( ER );
    LEMC.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );
    material->SetPlasticity ( LEMC );

    material->SetId ( 1 );
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


TPZGeoMesh * CreateGMeshGidTwoMatsTailings ( int ref )
{


    //string file ="/home/diogo/projects/pz/data/gid-tailings.msh";
    //string file ="/home/diogo/projects/pz/data/gid-tailings2.msh";
    string file ="/home/diogo/projects/pz/data/gid-tailings3.msh";

     REAL delta = 1.e-6;

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
    b[0] = 375.;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsbottom,delta );
    idsvec.push_back ( idsbottom );

    a = b;
    b[0] = 375.;
    b[1] = 150.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsright ,delta);
    idsvec.push_back ( idsright );


    a = b;
    b[0] = 225.;
    b[1] = 150.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstoprigth,delta );
    idsvec.push_back ( idstoprigth );
////h10-beta 45
    a = b;
    b[0] = 175.;
    b[1] = 200.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsramp ,delta);
    idsvec.push_back ( idsramp );

//h10-beta 30
// 	a = b;
//     b[0] = 27.675;
//     b[1] = 40.;
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsramp );
//     idsvec.push_back ( idsramp );
//h10-beta 60
// 	a = b;
//     b[0] = 39.2265;
//     b[1] = 40.;
//     read.Line ( a, b, ndivs, pathbottom );
//     read.FindIdsInPath ( pathbottom, idsramp );
//     idsvec.push_back ( idsramp );

   

    a = b;
    b[0] = 0.;
    b[1] = 200.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstopleft,delta );
    idsvec.push_back ( idstopleft );

    a = b;
    b[0] = 0.;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsleft,delta );
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

        std::cout << "quadrilaterals are not supported/need to be implemented!" <<std::endl;
        DebugStop();
    }

    if ( meshtopology.Cols() ==3 ) {
        TPZVec <long> TopoTri ( 3 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
            bool region =false;
            for(int inode=0; inode<3; inode++)
            {
                TopoTri[inode]=topol[iel][inode];
                REAL xco = meshcoords(TopoTri[inode],0);
                REAL yco = meshcoords(TopoTri[inode],1);
                REAL zco = meshcoords(TopoTri[inode],2);
                if(yco>150.)
                {
                    region=true;
                }
            }
            if(region)
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
            } else
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 2,*gmesh );
            }

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
    std::ofstream files ( "teste-diff-mats-geomesh.vtk" );
    std::set<int> myset;
    //TPZVTKGeoMesh::PrintGMeshVTKmy_material(gmesh,files,myset,true);
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;
}

