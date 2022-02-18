#include <stdio.h>
#include <iostream>
//#include "slopeanalysis.h"

#include "pzgmesh.h"
//#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "pzinterpolationspace.h"

#include "pzbndcond.h"
//#include "TPZMatElastoPlastic2D.h"
//#include "TPZMatElastoPlastic.h"
//#include "TPZElastoPlasticMem.h"

//#include "TPZElasticCriterion.h"

#include "TPZPlasticStepPV.h"
//#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzelastoplasticanalysis.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzpostprocanalysis.h"

//#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "tpzgeoelrefpattern.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif
#include "pznonlinanalysis.h"
#include "slopeanalysis.h"
#include "pzelastoplasticSest2D.h"
class TPZInterpolationSpace;
class TPZInterpolationSpace;


int main ( int argc, char *argv[] )
{

    std::cout << "asdasdasd"<< std::endl;
    //EVertical
    //ENonPenetrating
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern ( EOned );
    gRefDBase.InitializeUniformRefPattern ( ETriangle );
    gRefDBase.InitializeUniformRefPattern ( EQuadrilateral );
    std::cout << std::setprecision ( 15 );

    TPZSlopeAnalysis well;

    //std::string output = "ConfigSlope.vtk";
    //well.SetVtkOutPutName ( output );


    //well.LoadingRamp(factor);
    REAL sqj2_refine=0.01;
    int Startfrom=0;
    const int nsubsteps = 5;
    int porder=2;
    std::cout << "\n ------- 0 -------- "<<std::endl;
    well.GetCurrentConfig()->CreateGeometricMesh();

    well.GetCurrentConfig()->CreateComputationalMesh ( porder );

    well.GetCurrentConfig()->CreatePostProcessingMesh();

    well.PostProcess ( 0 );

    int steps =5;
    REAL delta=1.77/steps;
    REAL factor=delta;
    for ( int i=1; i<=steps; i++ ) {
        std::cout << "fac = " << factor<< std::endl;
        well.LoadingRamp ( factor );
        well.ExecuteInitialSimulation2 ( 5,20 );
        factor+=delta;
    }


    int nref=4;
    for ( int i=1; i<=nref; i++ ) {
        std::cout << "nref = " << i<< std::endl;
        well.PRefineElementAbove ( sqj2_refine, 3 );
        well.DivideElementsAbove ( sqj2_refine );
        well.ExecuteSimulation2 ( nsubsteps ,factor);
        well.PostProcess ( 0 );
    }





}


