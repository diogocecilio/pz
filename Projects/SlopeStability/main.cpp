
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "KLMaterial.h"
#include "KLStrMatrix.h"
#include "KLAnalysis.h"
#include "KLInterpolationSpace.h"
//#include "TPZPlasticityTest.h"
#include <iostream>
#include <cstdlib>
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzanalysis.h"
// #include "pzskylstrmatrix.h"darcy
#include "TPZTensor.h"
#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "pzfstrmatrix.h"
#include "pzbndmat.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "tpzgeoelrefpattern.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZTensor.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
#include "TPZDruckerPrager.h"
#include "KLRandomField.h"
#include "pzelastoplastic2D.h"
#include <pzmathyperelastic.h>
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"
#include "clock_timer.h"

#include "TPZVTKGeoMesh.h"
using namespace pzshape; // needed for TPZShapeCube and related classes


#include "pzlog.h"
#include "SlopeStabilityAnalysis.h"
#include "pzbfilestream.h"
#include "TPZProjectEllipse.h"
#include "arglib.h"
#include "run_stats_table.h"
#include "pzmganalysis.h"
#include "pzcompel.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pztransfer.h"
#include "tpzgeoelrefpattern.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzsolve.h"
#include "pzgengrid.h"
#include "tpzarc3d.h"

#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

#include "TPZVTKGeoMesh.h"


#include "TPZGeoLinear.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "tpzquadraticline.h"
#include "TPZWavyLine.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"

#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "tpzcompmeshreferred.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include "TPZGeoElement.h"
#include "pzgeoel.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzskylstrmatrix.h"

#include <time.h>
#include <stdio.h>
#include "pzl2projection.h"
#include "tpzgeoelmapped.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>

#include "pzgeotriangle.h"
#include "tpzarc3d.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoelrefpattern.h.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "TPZGeoCube.h"
#include <pzcompel.h>

#include <sstream>
using namespace std;
class TPZMaterial;

typedef   TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > plasticmat;



TPZCompMesh * CreateElastoPlasticComputationalMeshSlopeGid(int porder, TPZGeoMesh * gmesh);
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

void SolveShearRed();

void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol );

void SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZVec<TPZCompMesh*> vecmesh,int isol,REAL FS );

TPZElastoPlasticAnalysis  CreateAnal ( TPZCompMesh *fCMesh );

TPZCompMesh * CreteCompMesh ( TPZGeoMesh* gmesh,int porder );

TPZVec<TPZCompMesh *> CreateFields ( TPZGeoMesh * gmesh,int porder,int samples );

void  InsertMat ( TPZCompMesh * cmesh,int porder );

void LoadingRamp (TPZCompMesh * cmesh,  REAL factor);

void ShearRed (  );

void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);

void  CreatePostProcessingMesh(TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh);

void PostProcessF(TPZCompMesh * cmesh,int resolution);

void Post(TPZCompMesh * cmesh);

int main()
{
	
	
	//ShearRed (  );
	//SolveShearRed();
	//ShearRed (  );

	
	int porder =1;
	int samples =10;

 	TPZGeoMesh *gmesh = CreateGeometricMeshSlopeGid ( 0 );
// 	
    
// 	TPZVec<TPZCompMesh *> vecmesh = CreateFields (gmesh1,porder,samples );
// 	
// 	
	
	cout << "here "<<endl;
    //TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy);
   //TPZGeoMesh *gmesh = CreateGeometricMeshSlope ( 3);
	//TPZGeoMesh *gmesh = CreateGeometricMeshSlopeGid ( 0 );

    std::ofstream meshfile ( "GeoMesh.txt" );
    gmesh->Print ( meshfile );

    // Create computational mesh
    //TPZCompMesh *cmesh = CMesh_m(gmesh, pOrder);
    TPZCompMesh *cmesh = CreateComputationalMeshSlope ( gmesh, porder );
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
    analysis.DefineGraphMesh ( 2,scalarVars,vectorVars,"Darcyx.vtk" );
    constexpr int resolution{0};
    analysis.PostProcess ( resolution );

    //}

    std::cout << "FINISHED!" << std::endl;


    return 0;
}


void SolveShearRed()
{
	int ref = 0;
	int porder=2;
	TPZGeoMesh * gmesh = CreateGeometricMeshSlopeGid ( ref );
	//TPZGeoMesh * gmesh = CreateGeometricMeshSlope (  ref );
	TPZCompMesh * cmesh = CreateElastoPlasticComputationalMeshSlopeGid(porder,gmesh);
	
	REAL loadfactor =1.;
	
	LoadingRamp(cmesh,loadfactor);
	
	TPZElastoPlasticAnalysis  anal = CreateAnal(cmesh);

	REAL tol2 = 1.e-3;
	int NumIter = 30;
	bool linesearch = true;
	bool checkconv = false;
	
	anal.IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
	
	anal.AcceptSolution();
	
	Post(cmesh);
	
// 	TPZPostProcAnalysis  * postproc = new TPZPostProcAnalysis();
// 	
// 	CreatePostProcessingMesh( postproc, cmesh);
// 	
// 	//postproc->SetCompMesh(cmesh);
// 	
// 	const std::string name = "out";
// 	
// 	//postproc->Print(name, std::cout);
// 	
// 	
// 	TPZVec<int> PostProcMatIds(1,1);
// 	
// 	TPZStack<std::string> PostProcVars, scalNames, vecNames;
// 	
// 	PostProcessVariables(scalNames, vecNames);
// 
//     std::string vtkFile = "outFile.vtk";
// 
//     postproc->DefineGraphMesh(2,scalNames,vecNames,vtkFile);
// 
//     postproc->PostProcess(0);



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

void ShearRed (  )
{
	REAL norm = 1000.;
	REAL tol2 = 1.e-3;
	int NumIter = 30;
	bool linesearch = true;
	bool checkconv = false;
	
	
	int ref = 4;
	
	int porder=3;
	
	//TPZGeoMesh * gmesh = CreateGeometricMeshSlopeGid ( ref );
	TPZGeoMesh * gmesh = CreateGeometricMeshSlope (  ref );
	
	TPZCompMesh * cmesh = CreateElastoPlasticComputationalMeshSlopeGid(porder,gmesh);
	
	REAL loadfactor =1.;
	
	LoadingRamp(cmesh,loadfactor);
	
    REAL FS=1,FSmax=10000.,FSmin=0.,tol=1.e-3;
	
	int neq = cmesh->NEquations();
	
    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);

	int counterout = 0;
	
	TPZFMatrix<REAL>  deltau;
	
	plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();

    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi=phi0;
    REAL c=cohesion0;
	
	
	
    do {

		TPZElastoPlasticAnalysis  anal = CreateAnal(cmesh);

        anal.IterativeProcess ( std::cout, tol2, NumIter,linesearch,checkconv );
		
		norm = Norm(anal.Rhs());
		
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
		
		plasticmat *material = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
		
        TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();

        c=cohesion0/FS;
		
        phi=atan ( tan ( phi0 ) /FS );

        LEMC.fYC.SetUp ( phi, phi, c, ER );
		
        material->SetPlasticity ( LEMC );
		
        counterout++;
    }  while ( ( FSmax - FSmin ) / FS > tol );
}


TPZElastoPlasticAnalysis  CreateAnal ( TPZCompMesh *fCMesh )
{
    TPZElastoPlasticAnalysis  analysis ( fCMesh,std::cout );

    TPZSkylineStructMatrix full ( fCMesh );

    TPZStepSolver<REAL> step;

    step.SetDirect ( ELDLt );

    analysis.SetStructuralMatrix ( full );

    analysis.SetSolver ( step );

    return analysis;
}


TPZGeoMesh * CreateGeometricMeshSlopeGid ( int ref )
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

TPZGeoMesh *  CreateGeometricMeshSlope ( int ref )
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
void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        REAL atm  = 10.33;//10.33 mca = 1 atm
        disp[0] = ( 40-y );
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        disp[0] = -1;
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
	int matid=1,dim=2;
    auto *material = new TPZDarcyFlow ( matid,dim );
    //Bet Degan loamy sand 
	//STATE permeability = 6.3e-5;//m/s
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
    

     TPZMaterial * BCond0 = material->CreateBC ( material, -3, 0, val1, val2 );//tr
     BCond0->SetForcingFunction(pressure);

     TPZMaterial * BCond1 = material->CreateBC ( material, -4, 0, val1, val2 );//ramp
     BCond1->SetForcingFunction(pressure);
     
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

/// create the computational mesh
TPZCompMesh * CreateElastoPlasticComputationalMeshSlopeGid(int porder, TPZGeoMesh * gmesh)
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

void SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZVec<TPZCompMesh*> vecmesh,int isol,REAL FS )
{

    int nels = plasticCmesh->NElements();
    int nels0 = vecmesh[0]->NElements();
    int nels1 = vecmesh[1]->NElements();
    //vecmesh[1]->Solution().Print ( "SOLUTION PHI" );
    if ( nels!=nels0 || nels!=nels1 || nels0!=nels1 ) {
        cout << "nels "<< nels << " | nels0 = "<< nels0 <<" | nels1 = "<< nels1 <<endl;
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


            TPZSolVec sol1,sol2;
            ComputeSolution ( celrandom1,datarandom1.phi,sol1 );
            ComputeSolution ( celrandom2,datarandom2.phi,sol2 );

            TPZVec<REAL> cohes = sol1[isol];
            TPZVec<REAL> phi = sol2[isol];
            // cout << "Sol2 = " << datarandom2.sol[isol][0] <<endl;
            //cout << "cohes = " << cohes <<endl;
            //cout << "phi = " << phi <<endl;
            //cout << "SolIntel = " << SolIntel[0] <<endl;
            //cout << "datarandom1.sol[0]  = " << datarandom1.sol[0] <<endl;
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
void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol )
{
 
    const int numdof = cel->Material()->NStateVariables();
    //const int ncon = cel->NConnects();
    const int ncon = 3;
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
        long dfseq = df->SequenceNumber();
        int dfvar = block.Size ( dfseq );
        long pos = block.Position ( dfseq );
        for ( int jn=0; jn<dfvar; jn++ ) {
            for ( int is=0; is<numbersol; is++ ) {
                sol[is][iv%numdof] += ( STATE ) phi ( iv/numdof,0 ) *MeshSol ( pos+jn,is );
            }
            iv++;
        }
    }

}//metho
TPZVec<TPZCompMesh *> CreateFields ( TPZGeoMesh * gmesh,int porder,int samples )
{

    TPZCompMesh * cmesh =  CreteCompMesh ( gmesh,porder );
    TPZCompMesh * cmesh2 =  CreteCompMesh ( gmesh,porder );

    //TPZCompMesh * cmesh (&slopeA.GetCurrentConfig()->fCMesh);
    //TPZCompMesh * cmesh2(&slopeA.GetCurrentConfig()->fCMesh);
    cmesh->Print ( std::cout );
    InsertMat ( cmesh, porder );
    cmesh->Print ( std::cout );
    InsertMat ( cmesh2, porder );
    int dim = cmesh->Reference()->Dimension();
    KLAnalysis * klanal = new KLAnalysis ( cmesh );
    KLMaterial * mat = dynamic_cast<KLMaterial*> ( cmesh->ElementVec() [0]->Material() );
    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );
    klanal->Solve();

    string file ="out-eigenfunctions.vtk";
    klanal->Post ( file,dim,0 );

    TPZVec<TPZCompMesh *> vecmesh ( 2 );

    TPZVec<REAL> mean ( 2 );
    TPZVec<REAL> cov ( 2 );
    mean[0]=10.;
    mean[1]=30.*M_PI/180.;
    cov[0]=0.2;
    cov[1]=0.3;
    REAL corsscorrelation = -0.5;
    KLRandomField * klf = new KLRandomField ( cmesh,klanal,mean,cov,samples,corsscorrelation );
    TPZVec<TPZFMatrix<REAL>> hhat;
    cout << "Creating Log-Normal Random Fields "<<endl;
    hhat = klf->CreateLogNormalRandomField();

    hhat[0].Print ( "coes" );
    cmesh->LoadSolution ( hhat[0] );
    //TPZCompMesh * cmesh2(cmesh);
    cmesh2->LoadSolution ( hhat[1] );
    vecmesh[0]=cmesh;
    vecmesh[1]=cmesh2;

    file ="out-coessss.vtk";
    klanal->Post ( file,dim,0 );
    KLAnalysis * klanal2 = new KLAnalysis ( cmesh2 );
    cmesh2->LoadSolution ( hhat[1] );
    file ="out-phissss.vtk";
    klanal2->Post ( file,dim,0 );
    return vecmesh;
}
TPZCompMesh * CreteCompMesh ( TPZGeoMesh* gmesh,int porder )
{
    int expansionorder=20;
    REAL Lx=20;
    REAL Ly=2;
    REAL Lz=1.;
    int type=3;
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
    // cmesh->AdjustBoundaryElements();
    //  cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

void  InsertMat ( TPZCompMesh * cmesh,int porder )
{
    int expansionorder=20;
    REAL Lx=20;
    REAL Ly=2;
    REAL Lz=1.;
    int type=3;
    int id=1;
    int dim = cmesh->Reference()->Dimension();

    KLMaterial * klmat = new KLMaterial ( id,Lx,Ly,Lz,dim,type,expansionorder );

    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
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
