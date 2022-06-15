




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

#include "KLMaterial.h"
#include "KLStrMatrix.h"
#include "KLAnalysis.h"
#include "KLInterpolationSpace.h"
#include "KLRandomField.h"
#include "TPZDarcyFlow.h"
#include <sys/stat.h>
#include <thread>
using namespace std;


typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;


//TPZVec<REAL> fPlasticDeformSqJ2;

#include "readgidmesh.h"
#include <TPZParFrontStructMatrix.h>

TPZGeoMesh * CreateGMesh ( int ref );

TPZGeoMesh * CreateGMeshGid ( int ref );

TPZGeoMesh * CreateGMeshGid2 ( int ref );

TPZCompMesh * CreateCMesh ( TPZGeoMesh * gmesh,int porder );

void LoadingRamp ( plasticmat * body,  REAL factor );


void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void Post ( TPZPostProcAnalysis * postproc,std::string vtkFile );

TPZElastoPlasticAnalysis * CreateAnal ( TPZCompMesh *cmesh,bool optimize );

void ShearRed ( TPZCompMesh * cmesh );

TPZManVector<REAL,10> ShearRed ( TPZCompMesh * cmesh,int isample,TPZManVector<TPZCompMesh *,2> vecmesh );

void  InsertMat ( TPZCompMesh * cmesh,int porder );

void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol );

void ComputeSolution2 ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol );

TPZManVector<TPZCompMesh *,2> CreateFields ( TPZGeoMesh * gmesh,int porder,int samples );

TPZCompMesh * CreateCMeshRF ( TPZGeoMesh* gmesh,int porder );

void SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZManVector<TPZCompMesh*,2> vecmesh,int isol,REAL FS );

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );

void PrintMat ( std::string out,TPZFMatrix<REAL> mat );

TPZManVector<TPZCompMesh *,2> CreateFieldsDummy ( TPZGeoMesh * gmesh,int porder );

template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs );
#include "fadType.h"
void SetFlux ( TPZCompMesh * plasticCmesh,TPZCompMesh* incmesh);
void PostDarcy(TPZAnalysis * analysis,string vtk);
void SolveDarcyProlem(TPZCompMesh *cmesh,string vtk);
TPZCompMesh * CreateCMeshDarcy( TPZGeoMesh *gmesh, int pOrder );
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
TPZManVector<REAL,10> GravityIncrease ( TPZCompMesh * cmesh );
REAL findnodalsol(TPZCompMesh *cmesh);
void MonteCarloFlow(int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder, string namefolderx,bool gim);
void MonteCarlo(int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder, string namefolderx,bool gim);
//ofstream outglobal ( "outglobal2.txt" );
//std::ofstream outloadu("gimloadvsu-darcy2.nb");

void myTreads ( int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder, string namefolderx,bool gim )
{
	MonteCarlo(a,b, gmesh,vecmesh,porder,namefolderx,gim );
}

void myTreadsFlow ( int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder, string namefolderx,bool gim )
{
	MonteCarloFlow(a,b, gmesh,vecmesh,porder,namefolderx,gim );
}

TPZManVector<TPZCompMesh *,2> SettingCreateFilds(TPZGeoMesh* gmesh,int porder,int samples, string namefolderx,bool createfield=false);

void SolveMultiThread();

void SolveSerial(int a, int b, bool flow,bool gim,string file);
void SolveMultiThread(int a0,int b0, int nthreads,bool gim,string file);

void SolveDeterministic(int porder);
void SolveDeterministicSRM(int porder);

void ComputeMeanInAReagion();
REAL ComputeMeanInAReagion(TPZCompMesh *cmesh,int isol);

int main()
{
    

//   float x = 1.5,y=0.5;
// 
//   //TinyFad<2,float> x1(x),x2(y), y1;
//   //Fad<float> x1(x),x2(y), y1;
//   TFad<2,float> x1(x),x2(y), y1;
//   x1.diff(0,2);
//   x2.diff(1,2);
//   y1 = ( sin(x1/x2) + x1/x2 - exp(x2) )* (x1/x2 - exp(x2));
// 
//   cout << "Fad     : " << y1 << endl;
//  return 0;
// 	string namefolderx = "test-folder";
// 	char* cstr = new char[namefolderx.length() + 1];
//     strcpy ( cstr, namefolderx.c_str() );
//     int check = mkdir ( cstr, 777 );
	
// 	int samples=1000;
// 	
// 	int porder=2;
// 	
// 	string namefolder;
// 	
// 	TPZGeoMesh* gmesh = CreateGMeshGid ( 0 );
// 		
// 
// 	TPZManVector<TPZCompMesh*,2> vecmesh = SettingCreateFilds(gmesh, porder, samples,namefolder,true);
// 		
// 	return 0;
// 		
//

//	SolveDeterministic(2);
//	SolveDeterministicSRM(2);

    
   // TPZManVector<TPZCompMesh*,2> vecmesh = SettingCreateFilds(gmesh, porder, samples,file);
		
	//string file = "/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/saida-flow";
 	string file = "/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/saida";
 	bool flow=false;
 	bool gim = false;
 	SolveSerial(7,8,flow,gim,file);
    
 //   ComputeMeanInAReagion();
	
//0.898321 182
//	SolveMultiThread(0,1000,1,gim,file);//1:10 hrs
	//SolveMultiThread(500,1000,4,gim);


	
    return 0;
}

void ComputeMeanInAReagion()
{
    int samples=1000;
	
	int porder=2;
	
	TPZGeoMesh* gmesh = CreateGMeshGid ( 0 );
    string file = "/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/test";
    
    bool createfield = false;

    TPZManVector<TPZCompMesh*,2> vecmesh = SettingCreateFilds(gmesh, porder, samples,file,createfield);
    
    for(int isol=0;isol<samples;isol++)
    {
        cout << " isol " << isol << endl;
       REAL mean = ComputeMeanInAReagion(vecmesh[1],isol);
       if(mean<0.26)
       {
           cout << " isol " << isol << endl;
           cout << " mean " << mean << endl;
       }
    }
    
}


REAL ComputeMeanInAReagion(TPZCompMesh *cmesh,int isol) {


    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension();

    int nels = cmesh->NElements();
    REAL summ=0.;
    
     int count=0;
    for ( int iel=0; iel<nels; iel++ ) {

        
        TPZCompEl *cel = cmesh->Element ( iel );
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        const TPZIntPoints &intpoints = intel->GetIntegrationRule();

        TPZManVector<REAL,3> point ( 3,0. );
        TPZMaterialData data;

        data.fNeedsSol = true;
        intel->InitMaterialData ( data );

        REAL weight=0;
        int nint = intpoints.NPoints();
        
        TPZTensor<REAL> epst,epsp;
    
        
        for ( long ip =0; ip<nint; ip++ ) {
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );

            REAL x =  data.x[0];
            REAL y =  data.x[1];
            if( (( 30. <  x  &&  x < 46. ) && (25. < y && y < 40. )))
            {
                TPZSolVec sol1;
                ComputeSolution ( cel,data.phi,sol1 );
                TPZVec<REAL> var = sol1[isol];
                //cout << " x  = "<< x     << endl;
               //cout << " y  = "<< y   << endl;
               // cout << " cohes  = "<< cohes     << endl;
                summ+=var[0];
                count++;
            }

        }
        
        ///cout << "mean " << summ/count<< endl ; 

    }
    REAL  mean = summ/count ; 
    cout << " mean " << mean<< endl ; 
    return mean;
}

void SolveSerial(int a, int b, bool flow,bool gim,string file)
{
	int samples=1000;
	
	int porder=2;

	
	TPZGeoMesh* gmesh = CreateGMeshGid ( 0 );
		
	
	
	if(flow)
 	{
		TPZManVector<TPZCompMesh*,2> vecmesh = SettingCreateFilds(gmesh, porder, samples,file);
	 	MonteCarloFlow(a, b,  gmesh, vecmesh,porder,file,gim);
	}else{
		TPZManVector<TPZCompMesh*,2> vecmesh = SettingCreateFilds(gmesh, porder, samples,file);
		MonteCarlo(a, b,  gmesh, vecmesh,porder,file,gim);
	}
	
}

void SolveMultiThread(int a0,int b0, int nthreads,bool gim,string file)
{
	int samples=b0-a0;
	
	int porder=2;
	


	
	std::vector <std::thread> threadsmat1,threadsmat2;
	
	int delta = int ( samples/nthreads );
	int a=a0;
	int b=a+delta;
	//a=500;b=1000;
	
	for ( int i=0; i<nthreads; i++ )
	{
		std::cout << "a = "<< a <<std::endl;
		
        std::cout << "b = "<< b <<std::endl;
		
		
 		TPZGeoMesh* gmesh = CreateGMeshGid ( 0 );
		
 		TPZManVector<TPZCompMesh*,2> vecmesh = SettingCreateFilds(gmesh, porder, samples,file);

		std::thread threadx ( myTreads,a,b, gmesh,vecmesh,porder,file ,gim);


		
		
	//	TPZGeoMesh* gmeshflow = CreateGMeshGid ( 0 );
		
		//TPZManVector<TPZCompMesh*,2> vecmesh = SettingCreateFilds(gmeshflow, porder, samples,namefolder);
		
	//	TPZManVector<TPZCompMesh*,2> vecmeshflow = SettingCreateFilds(gmeshflow, porder, samples,file);
		
	//	std::thread threadflow(myTreadsFlow,a,b, gmeshflow,vecmeshflow,porder,file,gim);
		
	
		threadsmat1.push_back ( std::move ( threadx ) );
		
	//	threadsmat2.push_back ( std::move ( threadflow ) );
		
		a=b+1;
		
        b+=delta;
		
	}

    //for ( auto &threadx: threadsmat1 ) threadx.join();
	for ( auto &threadflow: threadsmat2 ) threadflow.join();
}

TPZManVector<TPZCompMesh *,2> SettingCreateFilds(TPZGeoMesh* gmesh,int porder,int samples, string namefolderx,bool createfield)
{

	string outco,outphi;
	//outco="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/teste-cohes-field.txt";
	//outphi="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/teste-phi-field.txt";
	
		outco="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/cohesion-p2-type3-alpha-H10-beta45.txt";
	outphi="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/friction-p2-type3-alpha-H10-beta45.txt";

    TPZFMatrix<REAL> readco,readphi;

	TPZManVector<TPZCompMesh *,2> vecmesh;
	//if true creates random fields
	//if false read random fields
    if (createfield ) {
        vecmesh = CreateFields ( gmesh,porder,samples );
        PrintMat ( outco,vecmesh[0]->Solution() );
        PrintMat ( outphi,vecmesh[1]->Solution() );

    } else {
        vecmesh = CreateFieldsDummy ( gmesh,porder );
        ReadFile ( outco,readco );
        ReadFile ( outphi,readphi );

    }

    //create analysis to reorder the solution to transfer correctly to the determinitic mesh
    TPZElastoPlasticAnalysis * anal = new TPZElastoPlasticAnalysis ( vecmesh[0] );
    TPZElastoPlasticAnalysis * anal1 = new TPZElastoPlasticAnalysis ( vecmesh[1] );
	//load the random fields solution in dummy meshes
    vecmesh[0]->LoadSolution ( readco );
    vecmesh[1]->LoadSolution ( readphi );
	
	return vecmesh;
}

void MonteCarlo(int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder, string namefolderx,bool gim)
{
	
    char* cstr = new char[namefolderx.length() + 1];
    strcpy ( cstr, namefolderx.c_str() );
    int check = mkdir ( cstr, 777 );
	
	std::string vtkFile0,str;
	
	vtkFile0=namefolderx;
	if(gim==true)
 	{
	 vtkFile0+="/vtkfolder/out-gim-mc";
	}else{
		vtkFile0+="/vtkfolder/out-srm-mc";
	}
	

    for ( int isample=a; isample<b; isample++ ) {

		TPZPostProcAnalysis * postproc = new TPZPostProcAnalysis();
		auto s = std::to_string ( isample );

		//create the elastoplastic comp mesh
        TPZCompMesh *cmesh = CreateCMesh ( gmesh,porder );

        chrono::steady_clock sc;
        auto start = sc.now();
		
        TPZManVector<REAL,10> out;


		SetMaterialParamenters ( cmesh,vecmesh,isample,1 );
		if(gim==true)
		{		
			//for GIM set material param with FS = 1 only once.
			//The loading will increase inside the method, but the Strength will be constant
			out = GravityIncrease(cmesh);
			
		}else{
			
			out = ShearRed ( cmesh, isample, vecmesh );
			
		}
		
		auto end = sc.now();
		auto time_span = static_cast<chrono::duration<double>> ( end - start );
		cout<< "\n \n *************************** ";
		cout << "\n Gravity Increase took: " << time_span.count() << " seconds !!!";
        
		string  filename = namefolderx;
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

void SolveDeterministic(int porder)
{

	string vtkFile0="deterministic-gim.vtk";

	TPZPostProcAnalysis * postproc = new TPZPostProcAnalysis();
	
	TPZGeoMesh* gmesh = CreateGMeshGid ( 0 );

	//create the elastoplastic comp mesh
	TPZCompMesh *cmesh = CreateCMesh ( gmesh,porder );

	chrono::steady_clock sc;
	auto start = sc.now();
		
	TPZManVector<REAL,10> out;

	//out = GravityIncrease(cmesh);
	ShearRed(cmesh);
	
	auto end = sc.now();
	
	auto time_span = static_cast<chrono::duration<double>> ( end - start );
	
	cout << "\n Gravity Increase took: " << time_span.count() << " seconds !!!";
	
	postproc->SetCompMesh(0);
	
	CreatePostProcessingMesh ( postproc, cmesh );

	Post ( postproc,vtkFile0 );
	
}

void SolveDeterministicSRM(int porder)
{

	string vtkFile0="/vtkfolder/deterministic-gim.vtk";

	TPZPostProcAnalysis * postproc = new TPZPostProcAnalysis();
	
	TPZGeoMesh* gmesh = CreateGMeshGid ( 0 );

	//create the elastoplastic comp mesh
	TPZCompMesh *cmesh = CreateCMesh ( gmesh,porder );

	chrono::steady_clock sc;
	auto start = sc.now();
		
	TPZManVector<REAL,10> out;

	ShearRed(cmesh);
	
	auto end = sc.now();
	
	auto time_span = static_cast<chrono::duration<double>> ( end - start );
	
	cout << "\n Gravity Increase took: " << time_span.count() << " seconds !!!";
	
	postproc->SetCompMesh(0);
	
	CreatePostProcessingMesh ( postproc, cmesh );

	Post ( postproc,vtkFile0 );
	
}


void MonteCarloFlow(int a, int b, TPZGeoMesh * gmesh, TPZManVector<TPZCompMesh *,2> vecmesh,int porder, string namefolderx,bool gim)
{
	
    char* cstr = new char[namefolderx.length() + 1];
    strcpy ( cstr, namefolderx.c_str() );
    int check = mkdir ( cstr, 777 );
	
	//solve the darcy flow problem
	string vtk = namefolderx;
	vtk += "/Darcy.vtk";
    TPZCompMesh *  darcycompmesh =  CreateCMeshDarcy(gmesh,porder);
	cout << "\n Solving Darcy... " << endl;
	SolveDarcyProlem(darcycompmesh, vtk);
	
    
	
	std::string vtkFile0,str;
	
	vtkFile0=namefolderx;
	
	vtkFile0+="/vtkfolder/out-gim-mc";

    for ( int isample=a; isample<b; isample++ ) {

		TPZPostProcAnalysis * postproc = new TPZPostProcAnalysis();
		auto s = std::to_string ( isample );

		//create the elastoplastic comp mesh
        TPZCompMesh *cmesh = CreateCMesh ( gmesh,porder );

		cout << "\n Setting flux in mechanic comp mesh... " << endl;
		SetFlux(cmesh,darcycompmesh);

        chrono::steady_clock sc;
        auto start = sc.now();
		
        TPZManVector<REAL,10> out;

		//for GIM set material param with FS = 1 only once.
		//The loading will increase inside the method, but the Strength will be constant
		SetMaterialParamenters ( cmesh,vecmesh,isample,1 );
		
		if(gim==true)
		{
			out = GravityIncrease(cmesh);
		}else{
			out = ShearRed ( cmesh, isample, vecmesh );
		}
		
		
		
		
		
		auto end = sc.now();
		auto time_span = static_cast<chrono::duration<double>> ( end - start );
		cout<< "\n \n *************************** ";
		cout << "\n Gravity Increase took: " << time_span.count() << " seconds !!!";
        
		string  filename = namefolderx;
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

void ComputeSolution2 ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol )
{
    const int numdof = cel->Material()->NStateVariables();


    set<long>::iterator itr;
    std::set<long> cornerconnectlist;
    cel->BuildCornerConnectList ( cornerconnectlist );
    int ncon = cornerconnectlist.size();
    TPZVec<int> ids ( ncon );

    int counter=0;
    for ( itr = cornerconnectlist.begin(); itr != cornerconnectlist.end(); itr++ ) {
        ids[counter]=*itr;
        counter++;
    }


    //cel->BuildCornerConnectList();
    TPZFMatrix<STATE> &MeshSol = cel->Mesh()->Solution();

    long numbersol = MeshSol.Cols();
    sol.Resize ( numbersol );
    for ( long is=0 ; is<numbersol; is++ ) {
        sol[is].Resize ( numdof );
        sol[is].Fill ( 0. );

    }

    TPZBlock<STATE> &block = cel->Mesh()->Block();
    long iv = 0;
    for ( int in=0; in<ncon; in++ ) {
        TPZConnect *df = &cel->Connect ( in );
        long dfseq = df->SequenceNumber();

        int dfvar = block.Size ( dfseq );
        long pos = block.Position ( dfseq );
        for ( int jn=0; jn<dfvar; jn++ ) {
            for ( int is=0; is<numbersol; is++ ) {
                sol[is][ids[iv]%numdof] += ( STATE ) phi ( ids[iv]/numdof,0 ) *MeshSol ( pos+jn,is );
            }
            iv++;
        }
    }
}
TPZManVector<TPZCompMesh *,2> CreateFields ( TPZGeoMesh * gmesh,int porder,int samples )
{

    TPZCompMesh * cmesh =  CreateCMeshRF ( gmesh,porder );
    TPZCompMesh * cmesh2 =  CreateCMeshRF ( gmesh,porder );

    //TPZCompMesh * cmesh (&slopeA.GetCurrentConfig()->fCMesh);
    //TPZCompMesh * cmesh2(&slopeA.GetCurrentConfig()->fCMesh);
    //cmesh->Print ( std::cout );
    //InsertMat ( cmesh, porder );
    //cmesh->Print ( std::cout );
    //InsertMat ( cmesh2, porder );
    int dim = cmesh->Reference()->Dimension();
    KLAnalysis * klanal = new KLAnalysis ( cmesh );
    //TPZAnalysis *analysis = new TPZAnalysis (cmesh);
    //analysis->SolveKL();
    KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );
    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );
    klanal->Solve();

    string file ="rffolder/out-eigenfunctions.vtk";
    klanal->Post ( file,dim,0 );

    TPZManVector<TPZCompMesh *,2> vecmesh ( 2 );

    TPZVec<REAL> mean ( 2 );
    TPZVec<REAL> cov ( 2 );
    mean[0]=10.;
    mean[1]=30.*M_PI/180.;
    cov[0]=0.3;
    cov[1]=0.2;
    REAL corsscorrelation = -0.5;
    KLRandomField * klf = new KLRandomField ( cmesh,klanal,mean,cov,samples,corsscorrelation );
    TPZVec<TPZFMatrix<REAL>> hhat;
    cout << "Creating Log-Normal Random Fields "<<endl;
    hhat = klf->CreateLogNormalRandomField();

    //hhat[0].Print ( "coes" );
    cmesh->LoadSolution ( hhat[0] );
    //TPZCompMesh * cmesh2(cmesh);
    cmesh2->LoadSolution ( hhat[1] );
    vecmesh[0]=cmesh;
    vecmesh[1]=cmesh2;

    file ="rffolder/out-coessssP2.vtk";
    klanal->Post ( file,dim,0 );
    KLAnalysis * klanal2 = new KLAnalysis ( cmesh2 );
    cmesh2->LoadSolution ( hhat[1] );
    file ="rffolder/out-phissss.vtk";
    klanal2->Post ( file,dim,0 );
    return vecmesh;
}

TPZManVector<TPZCompMesh *,2> CreateFieldsDummy ( TPZGeoMesh * gmesh,int porder )
{

    TPZCompMesh * cmesh =  CreateCMeshRF ( gmesh,porder );
    TPZCompMesh * cmesh2 =  CreateCMeshRF ( gmesh,porder );

    // InsertMat ( cmesh, porder );
    // InsertMat ( cmesh2, porder );

    TPZManVector<TPZCompMesh *,2> vecmesh ( 2 );

    vecmesh[0]=cmesh;
    vecmesh[1]=cmesh2;

    return vecmesh;
}

TPZCompMesh * CreateCMeshRF ( TPZGeoMesh* gmesh,int porder )
{
    int expansionorder=300;
    REAL Lx=20.;
    REAL Ly=2.;
    REAL Lz=1.;
    int type=3;
	//int type=1;
    int id=1;
    int dim = gmesh->Dimension();

    //gmesh->ResetReference();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    KLMaterial * klmat = new KLMaterial ( id,Lx,Ly,Lz,dim,type,expansionorder );

    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

void Post ( TPZPostProcAnalysis * postproc,std::string vtkFile )
{

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    postproc->DefineGraphMesh ( 2,scalNames,vecNames,vtkFile );

    postproc->PostProcess ( 0 );
}


TPZGeoMesh * CreateGMeshGid2 ( int ref )
{


	string file ="/home/diogo/projects/pz/data/h10-beta45-287.msh";



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
    b[0] = 50.;
    b[1] = 0.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsbottom );
    idsvec.push_back ( idsbottom );

    a = b;
    b[0] = 50.;
    b[1] = 10.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsright );
    idsvec.push_back ( idsright );


    a = b;
    b[0] = 30.;
    b[1] = 10.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idstoprigth );
    idsvec.push_back ( idstoprigth );
////h10-beta 45
    a = b;
    b[0] = 20.;
    b[1] = 20.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsramp );
    idsvec.push_back ( idsramp );
 
	

    a = b;
    b[0] = 0.;
    b[1] = 20.;
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

TPZGeoMesh * CreateGMeshGid ( int ref )
{



    // string file ="/home/diogo/projects/pz/data/mesh-teste-pz-fromathematica2.msh";
    //string file ="/home/diogo/projects/pz/data/quad-gid.msh";
    // string file ="/home/diogo/projects/pz/data/gid-tri-2.msh";
   // string file ="/home/diogo/projects/pz/data/gid-tri-1kels.msh";
//string file ="/home/diogo/projects/pz/data/gid-tri-2kels.msh";
//string file ="/home/diogo/projects/pz/data/gid-tri-4k.msh";
//string file ="/home/diogo/projects/pz/data/gid-tri-2k.msh";
string file ="/home/diogo/projects/pz/data/h10-beta45.msh";
//string file ="/home/diogo/projects/pz/data/h10-beta60.msh";


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
////h10-beta 45
    a = b;
    b[0] = 35.;
    b[1] = 40.;
    read.Line ( a, b, ndivs, pathbottom );
    read.FindIdsInPath ( pathbottom, idsramp );
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
//     REAL mc_cohesion    = 10.0;
//     REAL mc_phi         = ( 30.0*M_PI/180 );
//     REAL mc_psi         = mc_phi;
	
	REAL mc_cohesion    = 50.0;
    REAL mc_phi         = ( 20.0*M_PI/180 );
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

// void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor )
// {
//     plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
//     TPZManVector<REAL, 3> force ( 3,0. );
//     force[1]=-factor*20.;
//     body->SetBodyForce ( force );
// 
// }

void LoadingRamp ( plasticmat * body,  REAL factor )
{
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

TPZManVector<REAL,10> ShearRed ( TPZCompMesh * cmesh,int isample,TPZManVector<TPZCompMesh *,2> vecmesh )
{
	TPZManVector<REAL,10> out(10,0.);
	plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    LoadingRamp ( body,1. );

    REAL FS=1.,FSmax=10000.;
    REAL FSmin=0.;
    REAL tol=1.e-2;
    REAL norm = 1000.;
    REAL tol2 = 1.;
    int NumIter = 20;
    bool linesearch = true;
    bool checkconv = false;

    int counterout = 0;

    std::ofstream files ( "saidanewton.txt" );

    int neq = cmesh->NEquations();

    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );
	int iters;

    do {

        bool optimize =true;
        TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );

        chrono::steady_clock sc;
        auto start = sc.now();

        bool conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters );
        norm = Norm ( anal->Rhs() );


        //anal->AcceptSolution();


        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );

        files<< " FS "<< FS <<  "| step = " << counterout << " | Rhs norm = " << norm  << " | IterativeProcess Time: " << time_span.count() << " seconds !!! " <<std::endl;

        std::cout << "FS "<< FS <<  "| step = " << counterout << " | Rhs norm = " << norm  << " | IterativeProcess Time: " << time_span.count() << " seconds !!! " <<std::endl;

        if ( conv==false ) {
			cmesh->Solution().Zero();
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {
			cmesh->Solution().Zero();
            FSmin = FS;
			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }

        SetMaterialParamenters ( cmesh,vecmesh,isample,FS );
        //if(( FSmax - FSmin ) / FS > tol)anal->AcceptSolution();
        counterout++;
		delete anal;
    }  while ( ( FSmax - FSmin ) / FS > tol );

    bool optimize =true;
    TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
    anal->IterativeProcess ( files, tol2, NumIter,linesearch,checkconv );
    anal->AcceptSolution();
    //cmesh->LoadSolution(displace);

	delete anal;
	out[0]=FS;
	out[1]=counterout;
	out[2]=norm;
    return out;

}


void ShearRed ( TPZCompMesh * cmesh )
{
	plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    LoadingRamp ( body,1. );

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
#include "TPZFrontMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "tpzsparseblockdiagonalstructmatrix.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
//#include <boost/mpl/void.hpp>
//#include <boost/mpl/void.hpp>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/stat.h>
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

    long neq = cmesh->NEquations();
    TPZVec<long> activeEquations;
    analysis->GetActiveEquations ( activeEquations );
    TPZEquationFilter filter ( neq );
    filter.SetActiveEquations ( activeEquations );
    matskl.EquationFilter() = filter;
    analysis->SetStructuralMatrix ( matskl );

    //step.SetDirect(ECholesky);
    analysis->SetSolver ( step );


    return analysis;
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
        //TPZFMatrix<REAL> elsol = vecmesh[0]->ElementSolution()[iel];
        // cout << "elsol "<< elsol  <<endl;
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


            //datarandom1.fCornerNodeIds
            //celrandom1->Print();
            //TPZGeoEl *gel = celrandom1->Reference();
            //gel->NCornerNodes();
            //gel->Print();

            TPZSolVec sol1,sol2;
            ComputeSolution ( celrandom1,datarandom1.phi,sol1 );
            ComputeSolution ( celrandom2,datarandom2.phi,sol2 );

            TPZVec<REAL> cohes = sol1[isol];
            TPZVec<REAL> phi = sol2[isol];
           //  cout << "Sol2 = " << datarandom2.sol[isol][0] <<endl;
           // cout << "cohes = " << cohes <<endl;
           // cout << "phi = " << phi <<endl;
            //cout << "SolIntel = " << SolIntel[0] <<endl;
          //  cout << "datarandom1.sol[0]  = " << datarandom1.sol[0] <<endl;
            //                    hhatcopy[irow][0] = hhat0[irow][0] / FS;
            //          hhatcopy[irow][1] = atan ( tan ( hhat0[irow][1] ) / FS );
            dataplastic.intGlobPtIndex = globpt;
            mem[globpt].fPlasticState.fmatprop[0]=cohes[0]/FS;
            mem[globpt].fPlasticState.fmatprop[1]=atan ( tan ( phi[0] ) /FS );
// 			mem[globpt].fPlasticState.fmatprop[0]=10;
//             mem[globpt].fPlasticState.fmatprop[1]=30*M_PI/180.;
            TPZTensor<REAL> epsTotal,sigma;
            //pMatWithMem2->ApplyStrainComputeSigma(epsTotal,sigma);
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
    auto *material = new TPZDarcyFlow ( matid,dim );
    //Bet Degan loamy sand
    //STATE permeability = 0.0063 ;//cm/s
    //STATE permeability = 0.000063;//m/s
    //STATE permeability = 0.1;//m/s
    STATE permeability = 1.;//m/s
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

void SetFlux ( TPZCompMesh * plasticCmesh,TPZCompMesh* incmesh)
{

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


        TPZCompEl *celplastic = plasticCmesh->Element ( iel );
        TPZInterpolationSpace *intelplastic = dynamic_cast<TPZInterpolationSpace *> ( celplastic );

        TPZGeoEl *gelplastic = celplastic->Reference();

        TPZCompEl *celdarcy = incmesh->Element ( iel );

        TPZGeoEl *geldarcy = celdarcy->Reference();

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

void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    const auto &x=pt[0];
    const auto &y=pt[1];
    const auto &z=pt[2];
    REAL atm  = 10.33;//10.33 mca = 1 atm
    disp[0] = /*kn/m^3*/  ( 40-y )/* m */ ;/* = kn/m^2 = kPa*/
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    const auto &x=pt[0];
    const auto &y=pt[1];
    const auto &z=pt[2];
    disp[0] = -1;
}

REAL findnodalsol(TPZCompMesh *cmesh) {

    //b[0] = 39.2265;
    //b[1] = 40.;
    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension();
    TPZVec<REAL> xd(dim,0.);
    TPZVec<REAL> mpt(dim,0.);
// xd[0] =39.2265; xd[1] = 40.;
    xd[0] =35.;
    xd[1] = 40.;
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


TPZManVector<REAL,10> GravityIncrease ( TPZCompMesh * cmesh )
{

	TPZManVector<REAL,10> output(10,0.);
    REAL FS=0.5,FSmax=5.,FSmin=0.,tol=0.01;
    int neq = cmesh->NEquations();
    int maxcount=50;
    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;

    REAL norm = 1000.;
     REAL tol2 = 0.001;
    int NumIter = 30;
    bool linesearch = true;
    bool checkconv = false;
	std::ofstream outloadu("outloadu.nb");
    
    REAL uy=0.;
	plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    outloadu << "\n plot = {";
    do {


        LoadingRamp ( body,FS );
		//LoadingRamp ( body,0.0000001 );
        bool optimize =true;
        TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,optimize );
        chrono::steady_clock sc;
        auto start = sc.now();
		int iters;
        bool conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv ,iters);

        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << "\n last FSsss = " << FS  <<" | Load step = " << counterout << "| total time in iterative process =  " << time_span.count()<< std::endl;
        //anal->IterativeProcess ( outnewton, tol2, NumIter);

        norm = Norm ( anal->Rhs() );

        if ( conv==false ) {
            cmesh->LoadSolution(displace0);
			//cmesh->Solution().Zero();
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;

        } else {
            uy+=findnodalsol(cmesh);
            outloadu << "{ "<<-uy << ", " << FS << " } ," << endl;
			displace0 = anal->Solution();
            FSmin = FS;
            anal->AcceptSolution();

			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );

        }
        
        counterout++;
		delete anal;
    }  while ( (( FSmax - FSmin ) / FS > tol && counterout<maxcount) );
	//delete body;
    outloadu <<  " };";
	outloadu <<"\n ListLinePlot[plot,PlotRange->All]" << endl;
	TPZElastoPlasticAnalysis  * anal = CreateAnal ( cmesh,true );
	anal->AcceptSolution();
	output[0]=FS;
	output[1]=counterout;
	output[2]=norm;
    return output;
}

