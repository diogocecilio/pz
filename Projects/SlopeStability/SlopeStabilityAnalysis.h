//
//  WellBoreAnalysis.h
//  PZ
//
//  Created by phil on 1/18/13.
//  Copyright (c) 2013 LabMeC. All rights reserved.
//

#ifndef PZ_WellBoreAnalysis_h
#define PZ_WellBoreAnalysis_h

#define PlasticPQP

#include "pzcmesh.h"
#include "TPZSandlerDimaggio.h"
#include "TPZTensor.h"
#include "pzgeoel.h"
#include "pzpostprocanalysis.h"
#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"
#include "pzstack.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzelasticSest2D.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplasticSest2D.h"
#include "TPZParSkylineStructMatrix.h"
#include "tpzsparseblockdiagonalstructmatrix.h"

#include "pzelctemp.h" // TPZIntelGen<TSHAPE>
#include "pzshapecube.h" // TPZShapeCube
#include "TPZLadeKim.h"
#include "pzmat2dlin.h"

#include "pzbfilestream.h"
#include <sstream>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>
#include "pzelasmat.h"
#include "TPZVTKGeoMesh.h"

#include "TPZProjectEllipse.h"
//#include "TPZPlasticityTest.h"
#include <iostream>
#include <cstdlib>
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
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
//#include "GeoMeshClass.h"
#include <pzmathyperelastic.h>
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"
//#include "clock_timer.h"
//#include "TPBrAcidFunc.h"

#include "pzelastoplasticanalysis.h"
#include "pzelastoplasticSest2D.h"



#include "pzvisualmatrix.h"


#include "tpzchangeel.h"
//#include "poroelastoplastic.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzplasticdiagnostic.h"
#include <iostream>
#include "pzbfilestream.h"


#include "pzl2projection.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZProjectEllipse.h"
#include "pzlog.h"

#include "pzcmesh.h"
#include "TPZSandlerDimaggio.h"
#include "TPZTensor.h"
#include "pzgeoel.h"
#include "pzpostprocanalysis.h"
#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"
#include "pzstack.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzgengrid.h"
#include "pzcmesh.h"
#include "TPZSandlerDimaggio.h"
#include "TPZTensor.h"
#include "pzgeoel.h"
#include "pzpostprocanalysis.h"
#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"
#include "pzstack.h"
#include "TPZYCMohrCoulombPV.h"

#include "pzelastoplastic2D.h"
#include "pzelastoplasticSest2D.h"

/// create de standard mesh
void AddBoundaryConditions(TPZCompMesh *CMesh, TPZMaterial * mat, TPZTensor<STATE> &Confinement, STATE pressure);

/// Simulation models
enum EPlasticModel  {ESandler, EMohrCoulomb, EElastic};



/// Class which simulates the stability of a wellbore
class TPZSlopeStabilityAnalysis
{
    
public:
    
    
    /// this class represents a particular stage of the wellbore analysis
    struct TConfig
    {
        TConfig();
        
        ~TConfig();
        
        TConfig(const TConfig &copy);
        
        TConfig &operator=(const TConfig &copy);
        
        /// Write the data to the output stream
        void Write(TPZStream &out);
        
        /// Read the data from the input stream
        void Read(TPZStream &input);
        
        /// Load the solution stored in TConfig into the CompMesh and set the ZDeformation of the material
        void LoadSolution();
        
        /// Apply the deformation of the configuration to the element
        void ApplyDeformation(TPZCompEl *cel);
        
        /// Verify tangent elasto plastic relation
        void VerifyPlasticTangent(TPZCompEl *cel);
        
        /// Compute the maximum plastic element deformation associated with each element
        void ComputeElementDeformation();
        
        /// Compute the area of the domain at which sqJ2 is above a given value
        REAL ComputeAreaAboveSqJ2(REAL sqj2);
        
        /// Compute the area of the domain
        REAL ComputeTotalArea();

        /// Compute the removed area of the domain
        REAL RemovedArea();
        
        /// Delete the elements with sqj2 above the given value;
        void DeleteElementsAbove(REAL sqj2);

        /// Change the polynomial order of element using the plastic deformation as threshold
        void PRefineElementsAbove(REAL sqj2, int porder, set<long> &elindices);
        
        /// Divide the element using the plastic deformation as threshold
        void DivideElementsAbove(REAL sqj2, set<long> &elindices);
        
        /// Initialize the plastic history of the integration points of the element
        void InitializeElement(TConfig &from, TPZCompEl *cel);
        
        /// Verify the global equilibrium of the forces by boundary condition
        void VerifyGlobalEquilibrium(std::ostream &out = std::cout);
        
        /// Verify the global equilibrium of the forces by boundary condition (second version)
        void VerifyGlobalEquilibrium2(std::ostream &out = std::cout);
        
        /// Compute the Rhs for the mesh minus the elements with matid
        // this method is cumulative (sums to the rhs)
        void ComputeRhsExceptMatid(int matid, TPZFMatrix<STATE> &rhs);
        
        /// Compute the contribution of only matid
        // this method is cumulative (sums to the rhs)
        void ComputeRhsForMatid(int matid, TPZFMatrix<STATE> &rhs);
        
        /// Zera os componentes do rhs para connects diferentes do zero
        void FilterRhs(TPZFMatrix<STATE> &rhs);
        
        /// Compute the resultant x and y force
        void ComputeXYForce(TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force);
        
        /// Create the geometric and computational mesh based on the configuration parameters
        void CreateMesh();
        
        /// Return the mesh used for computations (multiphysics mesh or fCMesh)
        TPZCompMesh *CompMeshUsed();
        
        void ModifyWellElementsToQuadratic();
        
        /// Initialize the Sandler DiMaggio object and create the computational mesh
        void CreateComputationalMesh(int porder);
        
        /// Setup post processing mesh
        void CreatePostProcessingMesh();
        
        
        void CreateGeometricMeshSlope(int ref);
		void CreateGeometricMeshSlope2(int ref);
		void CreateComputationalMeshSlope(int porder);
		
		void CreateGeometricMeshSlopeGid (std::string file,int ref);
		void CreateComputationalMeshSlopeGid(int porder);

         STATE ComputeFarFieldWork();
        
        /// project the node on the boundary
        // returns true if the coordinate was changed
        bool ProjectNode(TPZVec<REAL> &co);
        

        /// Transform from physical domain to computational domain stress
        /**
         * the outcome depends on the well configuration
         */
        void FromPhysicalDomaintoComputationalDomainStress(TPZTensor<STATE> &physicalStress, TPZTensor<STATE> &computationalStress);
        
        /// print the configuration
        void Print(ostream &out);
        
        /// this method will modify the boundary condition of the computational mesh and the forcing function
        // factor is a transition parameter between the confinement tension and well pressure
        // factor = 1 corresponds to pure well pressure
        void SetWellPressure(STATE factor = 1.);
        
        /// this method will configure the forcing function and boundary condition of the computational mesh
        void ConfigureBoundaryConditions();
        
        /// Debugging method, printing the parameters of the forcing function
        void PrintForcingFunction(std::ostream &out = std::cout);
        
        /// Return gmesh
        TPZGeoMesh * GetGeoMesh();


        /// Parameters
        //TPZElasticityMaterial fEl;
        
        /// Plastic model
        EPlasticModel fModel;

        
        // Sandler Dimaggio PV
        TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> fSDPV;

        //Mohr PV
        TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> fMCPV;
        
        //object of Elastic Material
        TPZElasticityMaterialSest2D fMatEla; 
        
        /// Geometric mesh3
        TPZGeoMesh fGMesh;
        
        /// Computational mesh
        TPZCompMesh fCMesh;
        
        /// Matrix of incremental solutions
        TPZFMatrix<STATE> fSolution;
        
        /// Z Deformation associated with the solution
        STATE fZDeformation;
        
        /// Vector containing maximum element plastic deformation
        TPZVec<REAL> fPlasticDeformSqJ2;
        
        /// The post processing mesh with the transferred solution
        TPZPostProcAnalysis fPostprocess;

        std::string fHistoryLog;
        

        
        static int gNumThreads;

		void SetGMesh(TPZGeoMesh  gmesh)
		{
			fGMesh = gmesh;
		}
        
    };

    TPZSlopeStabilityAnalysis();

    TPZSlopeStabilityAnalysis(const TPZSlopeStabilityAnalysis &copy);

    TPZSlopeStabilityAnalysis &operator=(const TPZSlopeStabilityAnalysis &copy);
    
    ~TPZSlopeStabilityAnalysis();
    
    /// write the object on the stream
    void Write(TPZStream &output);
    
    /// read the object from the stream
    void Read(TPZStream &input);
    
    /// Computes the tension state transferring the geological stress state to the hidrostatic stress state
    /**
     * @param nsteps number of loading steps (the actual number of loading steps is one higher)
     * @param numnewton number of allowed newton iterations
     */
    REAL ExecuteInitialSimulation2(int nsteps);
    
    /// Computes an equilibrium state corresponding to the current boundary conditions
    REAL ExecuteSimulation2();
	
	    /// Computes the tension state transferring the geological stress state to the hidrostatic stress state
    /**
     * @param nsteps number of loading steps (the actual number of loading steps is one higher)
     * @param numnewton number of allowed newton iterations
     */
    void ExecuteInitialSimulation(int nsteps);
    
    /// Computes an equilibrium state corresponding to the current boundary conditions
    void ExecuteSimulation();
	
    
    /// Computes an equilibrium state corresponding to the current boundary conditions
    void ExecuteSimulationShearRed(TPZVec<TPZCompMesh *> vecmehs);

    /// verify the integrity of the elasto plastic material that is being used
    static void CheckDeformation(std::string filename = "deform.nb");
    
    /// verify if the stress tangent is computed correctly
    void VerifyTangentValidity();
    
    /// transfer the solution from the current configuration to the given configuration
    void TransferSolutionTo(TConfig &config);
    
    void DeleteElementsAbove(REAL sqj2)
    {
        fCurrentConfig.DeleteElementsAbove(sqj2);
    }
    
    
    /// Set the polynomial order of the elements which exceed plastic deformation
    unsigned int PRefineElementAbove(REAL sqj2, int porder)
    {
        std::set<long> elindices;
        fCurrentConfig.PRefineElementsAbove(sqj2, porder,elindices);
#ifdef DEBUG
        std::cout << "Number of elements prefined: " << elindices.size() << std::endl;
#endif
        // subject to integration points to the deformation history
        ApplyHistory(elindices);
        
        fCurrentConfig.fCMesh.Solution().Zero();
        fCurrentConfig.fSolution = fCurrentConfig.fCMesh.Solution();

        return elindices.size();
    }
    

    /// change the material id of the geometric elements of the current configuration
    void ChangeMaterialId(long idFrom, long idTo)
    {
        TPZGeoMesh *gmesh = &fCurrentConfig.fGMesh;
        long nel = gmesh->NElements();
        for (long iel=0; iel<nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel) continue;
            int matid = gel->MaterialId();
            if (matid == idFrom) {
                gel->SetMaterialId(idTo);
            }
        }
    }
    
    /// Divide the element using the plastic deformation as threshold
    unsigned int DivideElementsAbove(REAL sqj2);
    
    /// Post process the results of the current configuration
    void PostProcess(int resolution);
    
    /// Get the post processing variables
    static void PostProcessVariables(TPZStack<std::string> &scalnames, TPZStack<std::string> &vecnames);
    
    /// GetPostProcessedValues
    void PostProcessedValues(TPZVec<REAL> &x, TPZVec<std::string> &variables, TPZFMatrix<STATE> &values);
    
    /// verify the global equilibrium for the current configuration
    void VerifyGlobalEquilibrium(std::ostream &out = std::cout)
    {
        fCurrentConfig.VerifyGlobalEquilibrium(out);
    }
    
    /// Access method
    TConfig * GetCurrentConfig ()
    {
        return &fCurrentConfig;
    }
    
    void PopConfiguration()
    {
        if (this->fSequence.size() ==0) {
            return;
        }
        fSequence.pop_back();
        fCurrentConfig = *(this->fSequence.rbegin());
    }

    /// Access method
    TConfig * GetConfig (int index) {
        if (index < 0 || index >= fSequence.size())
            DebugStop();

        list<TPZSlopeStabilityAnalysis::TConfig>::iterator inte;
        int i=0;
        for (inte=fSequence.begin(); inte!=fSequence.end(); ++inte, i++)
        {
            if (i == index)
            {
                return &(*inte);
            }
        }
        DebugStop();
        return 0;
    }

    /// Return size of config list
    int GetConfigListSize () {
        return fSequence.size();
    }

    
    

   void SetSanderDiMaggioParameters(REAL poisson, REAL Elast, REAL A, REAL B, REAL C, REAL R, REAL D, REAL W)
    {

        STATE G=Elast/(2.*(1.+poisson));
        STATE K=Elast/(3.*(1.-2*poisson));
        STATE phi=0,psi=1.,N=0.;
        fCurrentConfig.fSDPV.fYC.SetUp( A,  B, C,  D, K, G, W, R, phi, N, psi);
        fCurrentConfig.fSDPV.fER.SetUp(Elast,poisson);

        fCurrentConfig.fModel = ESandler;
    }
		
	void SetMohrCoulombParameters(REAL poisson, REAL Elast, REAL c, REAL Phi, REAL Psi)
	{
		TPZElasticResponse ER;
		ER.SetUp(Elast,poisson);
		fCurrentConfig.fMCPV.fYC.SetUp(Phi, Psi, c, ER);
		fCurrentConfig.fMCPV.fER.SetUp(Elast,poisson);

        fCurrentConfig.fModel = EMohrCoulomb;
	}
    
    void SetElasticParameters(REAL poisson, REAL Elast)
    {
        fCurrentConfig.fMatEla.SetElasticParameters(Elast, poisson);
        fCurrentConfig.fModel = EElastic;
    }
	
    void Print(std::ostream &out);
    
private:
    
    /// Compute the linear elastic stiffness matrix
    void ComputeLinearMatrix(TPZVec<long> &activeEquations);
    
    /// Set the parameters of the linear material
    void ConfigureLinearMaterial(TPZElasticityMaterialSest2D &mat);
    
    /// Recompute the plastic memory of the integration points of these elements
    void ApplyHistory(std::set<long> &elindices);
    
public:
    /// Test the linear matrix with vertical compaction
    void TestLinearMaterial();
    
    int GetPostProcessNumber () {
        return fPostProcessNumber;
    }

    void SaveConfig(stringstream &strout);
    
    inline void SetVtkOutPutName(string output)
    {
        fVtkFile=output;
    }
    void LoadingRamp ( REAL factor)
{
	TConfig &LocalConfig = fCurrentConfig;
   // TPZMatElastoPlasticSest2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > *  body = NULL;
  //  body = dynamic_cast<TPZMatElastoPlasticSest2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > *> ( fCurrentConfig.fCMesh.FindMaterial ( 1 ) );
	
 	    TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > *  body = NULL;
     body = dynamic_cast<TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > *> ( fCurrentConfig.fCMesh.FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=-factor*20.;
    body->SetBodyForce ( force );

}

void  SetMaterialParamenters ( TPZCompMesh * plasticCmesh,TPZVec<TPZCompMesh*> vecmesh,int isol,REAL factor );


private:
    /// The object with the current configuration
    TConfig fCurrentConfig;
    
    /// The list of all previous configurations
    std::list<TConfig> fSequence;
    
    /// Index associated with the post processing file
    int fPostProcessNumber;
    
    /// Linear Elastic Stiffness matrix
    TPZAutoPointer<TPZMatrix<STATE> > fLinearMatrix;
    
    std::string fVtkFile;
    
    
};

// External variables
extern int startfrom;

#endif
