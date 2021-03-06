
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


TPZVec<REAL> fPlasticDeformSqJ2;

#include "readgidmesh.h"

TPZGeoMesh * CreateGMesh ( int ref );

TPZGeoMesh * CreateGMeshGid ( int ref );


TPZCompMesh * CreateCMesh(TPZGeoMesh * gmesh,int porder);

void LoadingRamp (TPZCompMesh * cmesh,  REAL factor);

void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);

void  CreatePostProcessingMesh(TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh);

void Post(TPZPostProcAnalysis * postproc,std::string vtkFile );

TPZElastoPlasticAnalysis * CreateAnal(TPZCompMesh *cmesh);

void ShearRed ( TPZCompMesh * cmesh);

void DivideElementsAbove(TPZCompMesh* cmesh,TPZGeoMesh* gmesh,REAL sqj2, std::set<long> &elindices);

void ComputeElementDeformation(TPZCompMesh* cmesh);

REAL findnodalsol(TPZCompMesh *cmesh);

void MaterialPointMohrCoulomb();

void FadTest();

int main()
{
	//FadTest();
	MaterialPointMohrCoulomb();
	
	return 0;
	
	
	int porder =2;

 	TPZGeoMesh *gmesh = CreateGMeshGid (0);
	
	
	//TPZGeoMesh *gmesh2 = CreateGMeshGid (0);

    TPZCompMesh *cmesh = CreateCMesh ( gmesh, porder );
	TPZCompMesh *cmesh2 = CreateCMesh ( gmesh, porder );
	
	int maxiter=100;
	REAL tol=1.e-3;
	bool linesearch=true;
	bool checkconv=false;
	int steps=30;
	REAL finalload = 4.;
	std::string vtkFile ="slopx.vtk";
	std::ofstream outloadu("loadvu.nb");
	outloadu << "plot = {";
	TPZPostProcAnalysis  * postproc = new TPZPostProcAnalysis();
	
	//ShearRed (cmesh2);
	

	REAL uy=0;
	for(int iloadstep=20;iloadstep<=steps;iloadstep++)
 	{
		

		
		TPZElastoPlasticAnalysis  * analysis =  CreateAnal(cmesh);
		
		
		REAL load = iloadstep*finalload/steps;
		cout << " \n --------- iloadstep  = "<< iloadstep << endl;
		
		cout << " \n --------- load  = "<< load << endl;
		
		LoadingRamp(cmesh,load);
		analysis->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv);
		uy+=findnodalsol(cmesh);
		TPZFMatrix<REAL> sol = analysis->Solution();
		
		//cmesh->Solution().Print(cout);
		cout << " sol "<<-uy << " | load = " << load <<    endl;
		analysis->AcceptSolution();
		
		outloadu << "{ "<<-uy << ", " << load << " } ," << endl;
		
					postproc->SetCompMesh(0);
		CreatePostProcessingMesh( postproc, cmesh);
		
		Post(postproc,vtkFile);
		//cmesh->LoadSolution();
		//ComputeElementDeformation(cmesh);
		
// 		std::set<long> elindices;
// 		REAL sqrtj2=0.005;

		//DivideElementsAbove(cmesh,gmesh,sqrtj2,elindices);
		

		
// 		LoadingRamp(cmesh,load);
// 		analysis->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv);
// 		analysis->AcceptSolution();
		//


	}
	outloadu <<  " }; ListLinePlot[plot,PlotRange->All]" << endl;

// 		TPZElastoPlasticAnalysis *anal  =  CreateAnal(cmesh);
// 		
// 		LoadingRamp(cmesh,finalload);
// 		
//  		anal->IterativeProcess(std::cout,tol,maxiter,linesearch,checkconv);
// 		
//  		anal->AcceptSolution();
// 		
// 		CreatePostProcessingMesh( postproc, cmesh);
// 		
// 		Post(postproc,vtkFile);
// 
// 
//     std::cout << "FINISHED!" << std::endl;


    return 0;
}

void MaterialPointMohrCoulomb(){

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

	sigma.XX()=17658.1;sigma.YY()=15494.7;sigma.ZZ()=16244.8;sigma.XY()=-19.7463;
	ER.ComputeDeformation(sigma,eps);
	REAL kprev=0.;
	TPZVec<REAL> coef;
	//cout << "before" << endl;
	//eps.Print(std::cout);
	//cout << "affter" << endl;
	mc.ApplyStrainComputeDep(eps,sigma,Dep);
	deps.XX()=0.000001;deps.YY()=0.000002;deps.ZZ()=0.0000015;deps.XZ()=-0.000001;deps.YZ()=0.000002;deps.XY()=0.000003;
	mc.TaylorCheck(eps,deps,kprev,coef);
	coef.Print(cout);
	//cout << "out" << endl;
	//sigma.Print(std::cout);
	cout << "Dep" << endl;
	Dep.Print(std::cout);
	
	
}
  const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &dx,TPZFMatrix<STATE> &d2x){
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
	
	dx[0]=test2.val().fastAccessDx(0);
	dx[1]=test2.val().fastAccessDx(1);
	
	d2x(0,0)=test2.fastAccessDx(0).fastAccessDx(0);
	d2x(1,1)= test2.fastAccessDx(1).fastAccessDx(1);
	d2x(1,0)= test2.fastAccessDx(0).fastAccessDx(1);
	d2x(0,1)= test2.fastAccessDx(1).fastAccessDx(0);
	
  };
  
  const auto rhs2 = [](const TPZFMatrix<STATE>&loc, TPZFMatrix<STATE> &dx,TPZFMatrix<STATE> &d2x){
        const REAL &x0 = loc.Get(0,0);
        const REAL &y0 = loc.Get(1,0);
	
	const int N=2;
	//dx.Resize(N,1);
	//d2x.Resize(N,N);
	TFad<N,TFad<N,float> > xsecnd2,ysecnd2,xi,beta;
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
	
	xi=x;
	xi.fastAccessDx(0)=1;
	xi.fastAccessDx(1)=0;
	beta=y;
	beta.fastAccessDx(0)=0;
	beta.fastAccessDx(1)=1;

// 	TFad<N,TFad<N,float> > sig1(79812.41),sig2(68385.69), sig3(61904.62),phi(20.*M_PI/180.f),c(490.f),nu(0.48),Pi(M_PI),um(1.f),dois(2.f),tres(3.f),seis(6.f),dezoito(18.f),sqrt3(sqrt(3.f));
		TFad<N,TFad<N,float> > sig1(183753.31587364181),sig2(157264.29661659044), sig3(143880.63541092159),phi(20.*M_PI/180.f),c(490.f),nu(0.48),Pi(M_PI),um(1.f),dois(2.f),tres(3.f),seis(6.f),dezoito(18.f),sqrt3(sqrt(3.f));
	TFad<N,TFad<N,float> >  test2;
	
	test2=((dois*(um-dois*nu))*(sqrt3*sig1 + sqrt3*sig2 +sqrt3*sig3 - tres*xi)*(sqrt3*sig1 + sqrt3*sig2 +sqrt3*sig3 - tres*xi) +(um + nu)*(tres*sig2 - tres*sig3 + (seis*sin(beta)*(-(sqrt3*c*cos(phi)) + xi*sin(phi)))/
          (cos(beta)*(um + sin(phi)) - (-um + sin(phi))*sin(beta + Pi/seis)))*(tres*sig2 - tres*sig3 + (seis*sin(beta)*(-(sqrt3*c*cos(phi)) + xi*sin(phi)))/
          (cos(beta)*(um + sin(phi)) - (-um + sin(phi))*sin(beta + Pi/seis)))+tres*(um + nu)*(-dois*sig1 + sig2 + sig3 + (dois*cos(beta)*(tres*c*cos(phi) - sqrt3* xi *sin(phi)))/
          (cos( beta )*(um+ sin(phi)) - (-um + sin(phi))*sin( beta  + Pi/seis)))*(-dois*sig1 + sig2 + sig3 + (dois*cos(beta)*(tres*c*cos(phi) - sqrt3* xi *sin(phi)))/
          (cos( beta )*(um+ sin(phi)) - (-um + sin(phi))*sin( beta  + Pi/seis))))/dezoito;

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
	x0(0,0)=0.;
	x0(1,0)=1.0472;
	for(int it=0;it<7;it++){
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
	cout <<  "xn = "  << endl;
	xn.Print(cout);
}



void Post(TPZPostProcAnalysis * postproc,std::string vtkFile )
{
	
	TPZVec<int> PostProcMatIds(1,1);
	
	TPZStack<std::string> PostProcVars, scalNames, vecNames;
	
	PostProcessVariables(scalNames, vecNames);

    postproc->DefineGraphMesh(2,scalNames,vecNames,vtkFile);

    postproc->PostProcess(0);
}
REAL findnodalsol(TPZCompMesh *cmesh){
 
    //b[0] = 39.2265;
    //b[1] = 40.;
TPZGeoMesh * gmesh =  cmesh->Reference();
	int dim = gmesh->Dimension();
 TPZVec<REAL> xd(dim,0.);
 TPZVec<REAL> mpt(dim,0.);
 xd[0] =35; xd[1] = 40.;
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
	
        	
	   		cout << "elemento encontrado"<<endl;
			cout << "index do elemento = " <<gel->Index() <<endl;
			int nnodes=gel->NNodes();
			for(int inode=0;inode<nnodes;inode++)
			{
				TPZVec<REAL> co(3);
				TPZGeoNode* node = gel->NodePtr(inode);
				node->GetCoordinates(co);
				id = node->Id();

				if(fabs(co[0]-xd[0])<1.e-3 && fabs(co[1]-xd[1])<1.e-3)
				{
					cout << "node id = "<<id <<endl;
					cout << " Coordinates "<< endl;
					cout << " x = "<< co[0] << ",  y = " << co[1] << endl;
					cout << " qsi = "<< qsi << endl;
					TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

        			TPZMaterialData data;

        			data.fNeedsSol = true;
		
        			intel->InitMaterialData ( data );
					
					intel->ComputeRequiredData ( data, qsi );
					
					cout << " data.sol  = "<<data.sol     << endl;
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

void ShearRed ( TPZCompMesh * cmesh)
{
	LoadingRamp(cmesh,1.);
	
    REAL FS=1,FSmax=5.,FSmin=0.,tol=1.e-3;
	int neq = cmesh->NEquations();

    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);

	int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    do {

		TPZElastoPlasticAnalysis  * anal = CreateAnal(cmesh);
		REAL norm = 1000.;
        REAL tol2 = 1.e-3;
        int NumIter = 30;
        bool linesearch = true;
        bool checkconv = false;
        anal->IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
		norm = Norm(anal->Rhs());

		//std::cout << "| Load step = " << counterout << " | Rhs norm = " << norm << "| delta displacement norm = "<< unorm << std::endl;
        if ( norm>= tol2) {
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
TPZGeoMesh * CreateGMeshGid ( int ref )
{



    // string file ="/home/diogo/projects/pz/data/mesh-teste-pz-fromathematica2.msh";
    //string file ="/home/diogo/projects/pz/data/quad-gid.msh";
    //string file ="/home/diogo/projects/pz/data/gid-tri-2.msh";
   // string file ="/home/diogo/projects/pz/data/gid-tri-1kels.msh";
//string file ="/home/diogo/projects/pz/data/gid-tri-880-sessenta.msh";
string file ="/home/diogo/projects/pz/data/h10-beta45.msh";


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


    a = b;
    b[0] = 45.;
    b[1] = 30.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstoprigth );
    idsvec.push_back ( idstoprigth );

    a = b;
    b[0] = 35.;
    b[1] = 40.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsramp );
    idsvec.push_back ( idsramp );
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
    REAL mc_cohesion    = 50.;
    REAL mc_phi         = ( 20.*M_PI/180 );
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000.;

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
