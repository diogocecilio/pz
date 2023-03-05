
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"
#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

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
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <chrono>
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
#include "TPZDarcyFlow.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.adaptivity"));
static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);

// bi-dimensional problem for elasticity
int main() {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMesh();

	// Creating computational mesh (approximation space and materials)
	int p = 3;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMesh(gmesh);
	// Solving linear equations
	// Initial steps
	TPZAnalysis an (cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	an.SetStructuralMatrix(strskyl);
	// Solver (is your choose) 
	TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;

	an.Run();
	
	// Post processing
	TPZManVector<std::string> scalarnames(3), vecnames(1);
	scalarnames[0] = "SigmaX";
	scalarnames[1] = "SigmaY";
	scalarnames[2] = "TauXY";
	vecnames[0] = "displacement";
	//vecnames[1] = "";
	an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions.vtk");

	an.PostProcess(0);
}




TPZGeoMesh *CreateGeoMesh() {

    REAL co[4][2] = {{0.,0.},{1,0},{1,10},{0,10}};
    long indices[1][4] = {{0,1,2,3}};
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    long nnode = 4;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }

    long el;
    long nelem = 1;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }

    TPZVec <long> TopoLine ( 2 );

    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 7, TopoLine, - 1, *gmesh );//bottom


    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 8, TopoLine, - 2, *gmesh );//top

   TPZVec <long> node(1);
   node[0]=3;
   new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 9, node, - 3, *gmesh );//top
//     long index;
//     gmesh->CreateGeoElement(EPoint,3,-3,index);
//
//     TPZGeoElBC gbc3(elvec[0],6,-3);

    gmesh->BuildConnectivity();


    for(int d=0;d<4;d++) {
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
	}

	return gmesh;
}
//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
//TPZElasticityMaterial(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1);
    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(1,10000.,0.2,0.,-1.);//selfweigth
   // TPZMaterial * mat = new TPZElasticityMaterial(1,10000.,0.2,0.,0.);
	// Creating four boundary condition
    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bctop,*bcload,*bcpoint;

    val2(1,0)=1.;
    bctop = mat->CreateBC(mat,-2,3,val1,val2);

    val2(0,0)=1.;
    bcpoint = mat->CreateBC(mat,-3,3,val1,val2);

    //val2(1,0)=-1.;
    //bcload = mat->CreateBC(mat,-1,1,val1,val2);
    cmesh->InsertMaterialObject(mat);
	// Inserting boundary conditions into computational mesh
	cmesh->InsertMaterialObject(bctop);

    cmesh->InsertMaterialObject(bcpoint);
    //cmesh->InsertMaterialObject(bcload);

	cmesh->SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}



/*
//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh() {

    REAL co[8][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.}};
    long indices[3][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7}};
    TPZGeoEl *elvec[3];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    long nnode = 8;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    long el;
    long nelem = 3;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sub;
	
    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],4,-1);
    // bc -2 -> Neumann at the bottom y==-1
    TPZGeoElBC gbc2(elvec[0],5,-2);
    // bc -3 -> Neumann at the right x==1
    TPZGeoElBC gbc3(elvec[0],6,-3);
    
    // bc -3 -> Neumann at the right x==1
    TPZGeoElBC gbc4(elvec[1],5,-3);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc5(elvec[1],6,-4);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc6(elvec[2],5,-4);
    
    // bc -5 -> Neumann at the left x==-1
    TPZGeoElBC gbc7(elvec[2],6,-5);
    
    // bc -6 -> Homogeneous Neumann
    TPZGeoElBC gbc8(elvec[2],7,-6);
    for(int d=0;d<2;d++) {
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
	}

	return gmesh;
}
//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());

    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(1,2000000000.,0.3,0.,0.);

	// Creating four boundary condition
    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bcBottom, *bcRight;
	//val1(1,1) = 1000000.;
	val2(1,0) = 10000000.;
    bcBottom = mat->CreateBC(mat,-2,1,val1,val2);
	//val1(1,1) = 0.;
	val2(1,0) = 10000000.;
    bcRight = mat->CreateBC(mat,-4,1,val1,val2);

    cmesh->InsertMaterialObject(mat);
	// Inserting boundary conditions into computational mesh
	cmesh->InsertMaterialObject(bcBottom);
	cmesh->InsertMaterialObject(bcRight);
    
	cmesh->SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    return cmesh;
}*/

