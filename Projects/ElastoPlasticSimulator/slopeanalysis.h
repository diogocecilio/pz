#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#include "clock_timer.h"
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
class TPZSlopeAnalysis
{

public:
    /// order in which formation stresses are stored
    enum MFormationOrder {ESh = 0, ESH = 1, ESV = 2};

    /// boundary condition numbers
    enum MWellBCs {EInner = -2, EBottom = -4, ELeft = -5, EOuter = -3};
    /// this class represents a particular stage of the wellbore analysis
    struct TConfig {
        TConfig();

        ~TConfig();

        TConfig ( const TConfig &copy );

        TConfig &operator= ( const TConfig &copy );

        /// Write the data to the output stream
        void Write ( TPZStream &out );

        /// Read the data from the input stream
        void Read ( TPZStream &input );

        /// Load the solution stored in TConfig into the CompMesh and set the ZDeformation of the material
        void LoadSolution();

        /// Apply the deformation of the configuration to the element
        void ApplyDeformation ( TPZCompEl *cel );

        /// Verify tangent elasto plastic relation
        void VerifyPlasticTangent ( TPZCompEl *cel );

        /// Compute the maximum plastic element deformation associated with each element
        void ComputeElementDeformation();

        /// Compute the area of the domain at which sqJ2 is above a given value
        REAL ComputeAreaAboveSqJ2 ( REAL sqj2 );

        /// Compute the area of the domain at which sqJ2 is above a given value
        REAL OpeningAngle ( REAL sqj2 );

        /// Compute the area of the domain
        REAL ComputeTotalArea();

        /// Compute the removed area of the domain
        REAL RemovedArea();

        /// Delete the elements with sqj2 above the given value;
        void DeleteElementsAbove ( REAL sqj2 );

        /// Diminish the spring which stabilizes the well by the given factor
        void RelaxWellSpring ( REAL factor );

        /// Change the polynomial order of element using the plastic deformation as threshold
        void PRefineElementsAbove ( REAL sqj2, int porder, std::set<long> &elindices );

        /// Divide the element using the plastic deformation as threshold
        // returns the number of elements divided
        //int DivideElementsAbove ( REAL sqj2, std::set<long> &elindices );
		
		void DivideElementsAbove(REAL sqj2, std::set<long> &elindices);

        /// Initialize the plastic history of the integration points of the element
        void InitializeElement ( TConfig &from, TPZCompEl *cel );

        /// Verify the global equilibrium of the forces by boundary condition
        void VerifyGlobalEquilibrium ( std::ostream &out = std::cout );

        /// Verify the global equilibrium of the forces by boundary condition (second version)
        void VerifyGlobalEquilibrium2 ( std::ostream &out = std::cout );

        /// Compute the Rhs for the mesh minus the elements with matid
        // this method is cumulative (sums to the rhs)
        void ComputeRhsExceptMatid ( int matid, TPZFMatrix<STATE> &rhs );

        /// Compute the contribution of only matid
        // this method is cumulative (sums to the rhs)
        void ComputeRhsForMatid ( int matid, TPZFMatrix<STATE> &rhs );

        /// Zera os componentes do rhs para connects diferentes do zero
        void FilterRhs ( TPZFMatrix<STATE> &rhs );

        /// Compute the resultant x and y force
        void ComputeXYForce ( TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force );

        /// Add elliptic breakout
        void AddEllipticBreakout ( REAL MaiorAxis, REAL MinorAxis );

        /// Create the geometric and computational mesh based on the configuration parameters
        void CreateGeometricMesh();

        /// Add the completion ring to the geometric mesh
        void AddGeometricRingElements();

        /// Return the mesh used for computations (multiphysics mesh or fCMesh)
        TPZCompMesh *CompMeshUsed();

        void ModifyWellElementsToQuadratic();

        /// Initialize the Sandler DiMaggio object and create the computational mesh
        void CreateComputationalMesh ( int porder );

        /// Setup post processing mesh
        void CreatePostProcessingMesh();


        void CreateMesh();

        STATE ComputeFarFieldWork();

        /// project the node on the boundary
        // returns true if the coordinate was changed
        bool ProjectNode ( TPZVec<REAL> &co );

        /// return the largest y-coordinate belonging to the last ellips
        // this value will be used to identify the meaningful points to match the next ellips
        REAL MaxYfromLastBreakout();

        /// Transform from physical domain to computational domain stress
        /**
         * the outcome depends on the well configuration
         */
        void FromPhysicalDomaintoComputationalDomainStress ( TPZTensor<STATE> &physicalStress, TPZTensor<STATE> &computationalStress );

        /// print the configuration
        void Print ( std::ostream &out );


        /// this method will modify the boundary condition of the computational mesh and the forcing function
        // factor is a transition parameter between the confinement tension and well pressure
        // factor = 1 corresponds to pure well pressure
        void SetWellPressure ( STATE factor = 1. );

        /// this method will configure the forcing function and boundary condition of the computational mesh
        void ConfigureBoundaryConditions();

        /// Set the Z deformation (for adapting the compaction)
        void SetZDeformation ( STATE epsZ );

        /// Set the configuration of the plastic material to use acid parameters
        void ActivateAcidification();

        /// Compute the average vertical stress of the configuration
        STATE AverageVerticalStress();

        /// Debugging method, printing the parameters of the forcing function
        void PrintForcingFunction ( std::ostream &out = std::cout );

        /// Return gmesh
        TPZGeoMesh * GetGeoMesh();

        // Geometry of mesh
        /// radius of the well
        REAL fInnerRadius;
        /// radius of the computational domain, whether it is vertical or horizontal
        REAL fOuterRadius;

        /// number of elements in the radial and circumferential direction
        TPZManVector<int,2> fNx;

        /// Size of the first element in the radial direction
        REAL fDelx;


        /// confinement stress in the physical domain
        TPZTensor<STATE> fConfinementTotal;

        TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> fMCPV;

        /// Geometric mesh3
        TPZGeoMesh fGMesh;

        /// Computational mesh
        TPZCompMesh fCMesh;

        /// Matrix of incremental solutions
        TPZFMatrix<STATE> fSolution;


        /// Vector containing maximum element plastic deformation
        TPZVec<REAL> fPlasticDeformSqJ2;

        /// The post processing mesh with the transferred solution
        TPZPostProcAnalysis fPostprocess;

        TPZElastoPlasticAnalysis fElastoplasticAnalysis;

        /// Log que apareca na tela de execucao
        std::string fHistoryLog;

        /// Log de execucao do codigo
        TPZStack<std::string > fExecutionLog;

        static int gNumThreads;
        /// Define the geometry of the simulation

    };

    TPZSlopeAnalysis();

    TPZSlopeAnalysis ( const TPZSlopeAnalysis &copy );

    TPZSlopeAnalysis &operator= ( const TPZSlopeAnalysis &copy );

    ~TPZSlopeAnalysis();

    /// write the object on the stream
    void Write ( TPZStream &output );

    /// read the object from the stream
    void Read ( TPZStream &input );

    /// Computes the tension state transferring the geological stress state to the hidrostatic stress state
    /**
     * @param nsteps number of loading steps (the actual number of loading steps is one higher)
     * @param numnewton number of allowed newton iterations
     */
    void ExecuteInitialSimulation ( REAL factor );

    /// Computes an equilibrium state corresponding to the current boundary conditions
    void ExecuteSimulation ( REAL factor, std::ostream &out );

    void LoadingRamp ( REAL factor );
    /// Set the polynomial order of the elements which exceed plastic deformation
    unsigned int PRefineElementAbove ( REAL sqj2, int porder )
    {
        std::set<long> elindices;
        fCurrentConfig.PRefineElementsAbove ( sqj2, porder,elindices );
#ifdef DEBUG
        std::cout << "Number of elements prefined: " << elindices.size() << std::endl;
#endif
        // subject to integration points to the deformation history
        ApplyHistory ( elindices );

        fCurrentConfig.fCMesh.Solution().Zero();
        fCurrentConfig.fSolution = fCurrentConfig.fCMesh.Solution();

        return elindices.size();
    }
    void ExecuteInitialSimulation2 ( int nsteps, int numnewton );
    void ExecuteSimulation2 ( int substeps, REAL factor );
    /// Divide the element using the plastic deformation as threshold
/// Divide the element using the plastic deformation as threshold
    unsigned int DivideElementsAbove ( REAL sqj2 )
    {
        std::set<long> elindices;
        fCurrentConfig.DivideElementsAbove ( sqj2,elindices );



        // subject the integration points with the deformation history
        ApplyHistory ( elindices );
        fCurrentConfig.ComputeElementDeformation();

        fCurrentConfig.fCMesh.Solution().Zero();
        fCurrentConfig.fSolution = fCurrentConfig.fCMesh.Solution();

        // invalidate the computational mesh associated with the postprocess mesh
        fCurrentConfig.fPostprocess.SetCompMesh ( 0 );

        return elindices.size();
    }


    /// change the material id of the geometric elements of the current configuration
    void ChangeMaterialId ( long idFrom, long idTo )
    {
        TPZGeoMesh *gmesh = &fCurrentConfig.fGMesh;
        long nel = gmesh->NElements();
        for ( long iel=0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            if ( !gel ) continue;
            int matid = gel->MaterialId();
            if ( matid == idFrom ) {
                gel->SetMaterialId ( idTo );
            }
        }
    }
    void ApplyHistory ( std::set<long> &elindices );
    /// Divide the element using the plastic deformation as threshold


    /// Post process the results of the current configuration
    void PostProcess ( int resolution );

    /// Get the post processing variables
    static void PostProcessVariables ( TPZStack<std::string> &scalnames, TPZStack<std::string> &vecnames );

    /// GetPostProcessedValues
    void PostProcessedValues ( TPZVec<REAL> &x, TPZVec<std::string> &variables, TPZFMatrix<STATE> &values );

    /// verify the global equilibrium for the current configuration
    void VerifyGlobalEquilibrium ( std::ostream &out = std::cout )
    {
        fCurrentConfig.VerifyGlobalEquilibrium ( out );
    }

    /// Access method
    TConfig * GetCurrentConfig ()
    {
        return &fCurrentConfig;
    }

    void PopConfiguration()
    {
        if ( this->fSequence.size() ==0 ) {
            return;
        }
        fSequence.pop_back();
        fCurrentConfig = * ( this->fSequence.rbegin() );
    }

    /// Access method
    TConfig * GetConfig ( int index )
    {
        if ( index < 0 || index >= fSequence.size() )
            DebugStop();

        std::list<TPZSlopeAnalysis::TConfig>::iterator inte;
        int i=0;
        for ( inte=fSequence.begin(); inte!=fSequence.end(); ++inte, i++ ) {
            if ( i == index ) {
                return & ( *inte );
            }
        }
        DebugStop();
        return 0;
    }

    /// Return size of config list
    int GetConfigListSize ()
    {
        return fSequence.size();
    }


    void CreatePostProcessingMesh ( TPZCompMesh * cmesh, TPZPostProcAnalysis * PostProcess );
    /// Initialize the object with standard parameters
    static void StandardConfiguration ( TPZSlopeAnalysis &obj );

    /// Configure the wellbore analysis to perform a linear analysis
    void LinearConfiguration ( int porder );


    /// Define the geometry of the simulation
    void SetInnerOuterRadius ( REAL inner, REAL outer )
    {
        fCurrentConfig.fInnerRadius = inner;
        fCurrentConfig.fOuterRadius = outer;
    }

    void SetMeshTopology ( REAL delx, TPZVec<int> &nx )
    {
        fCurrentConfig.fDelx = delx;
        fCurrentConfig.fNx = nx;
    }


    void SetMohrCoulombParameters ( REAL poisson, REAL Elast, REAL c, REAL Phi, REAL Psi )
    {
        TPZElasticResponse ER;
        ER.SetUp ( Elast,poisson );
        fCurrentConfig.fMCPV.fYC.SetUp ( Phi, Psi, c, ER );
        fCurrentConfig.fMCPV.fER.SetUp ( Elast,poisson );
    }

    void Print ( std::ostream &out );

    void PrintInitialConfiguration ( std::ostream &out );

    /// print the execution log
    void PrintExecutionLog ( std::ostream &out );

    void SaveConfig ( stringstream &strout );
private:

    /// Recompute the plastic memory of the integration points of these elements

public:

    int GetPostProcessNumber ()
    {
        return fPostProcessNumber;
    }

    /// append the string associated with the execution log
    void AppendExecutionLog ( std::stringstream &sout );

private:
    /// The object with the current configuration
    TConfig fCurrentConfig;

    /// The list of all previous configurations
    std::list<TConfig> fSequence;

    /// Index associated with the post processing file
    int fPostProcessNumber;

    std::string fVtkFile;


};

// External variables
extern int startfrom;

