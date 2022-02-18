
//
//  Created by diogo on 09/02/2022.
//
#include "slopeanalysis.h"
#include "pzelastoplasticanalysis.h"

#include "pzinterpolationspace.h"

#include "pzvisualmatrix.h"

#include "tpzquadraticquad.h"

#include "tpzchangeel.h"
//#include "poroelastoplastic.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzplasticdiagnostic.h"
#include <iostream>

#include "pzelastoplasticSest2D.h"
#include "pzl2projection.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZProjectEllipse.h"
#include "pzlog.h"

#include "arglib.h"
#include "run_stats_table.h"
#include "TPZGeoLinear.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.plasticity.wellboreanalysis"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse(Logger::getLogger("LogEllipse"));
#endif

int TPZSlopeAnalysis::TConfig::gNumThreads = 20;

TPZSlopeAnalysis::TPZSlopeAnalysis() : fCurrentConfig(), fSequence(), fPostProcessNumber(0)
{
  fVtkFile = "defaultVTKFileName.vtk";
}

TPZSlopeAnalysis::TPZSlopeAnalysis(const TPZSlopeAnalysis &copy) : fCurrentConfig(copy.fCurrentConfig), fSequence(copy.fSequence), fPostProcessNumber(copy.fPostProcessNumber)
{

}

TPZSlopeAnalysis::~TPZSlopeAnalysis()
{
    
}


TPZSlopeAnalysis::TConfig::TConfig() : fInnerRadius(0.), fOuterRadius(0.), fConfinementTotal(), fGMesh(), fCMesh(), fSolution(), fPlasticDeformSqJ2(),fMCPV()
{
    fCMesh.SetReference(&fGMesh);
}

TPZSlopeAnalysis::TConfig::TConfig(const TConfig &conf) : fInnerRadius(0.), fOuterRadius(0.), fConfinementTotal(), 
fGMesh(), fCMesh(), fSolution(), fPlasticDeformSqJ2() ,fMCPV()
{

    fGMesh.ResetReference();
    fCMesh.SetReference(&fGMesh);
    fCMesh.LoadReferences();

}

TPZSlopeAnalysis::TConfig::~TConfig()
{
//    std::cout << "Deleting config " << (void *) this << std::endl;
    fPostprocess.SetCompMesh(0);
    fCMesh.CleanUp();
    fCMesh.SetReference(0);
}

TPZSlopeAnalysis::TConfig &TPZSlopeAnalysis::TConfig::operator=(const TPZSlopeAnalysis::TConfig &copy)
{
    if (this == &copy) {
        return *this;
    }
    fInnerRadius = copy.fInnerRadius;
    fOuterRadius = copy.fOuterRadius;
    fConfinementTotal = copy.fConfinementTotal;
    fMCPV = copy.fMCPV;

    fPostprocess = copy.fPostprocess;
    fPostprocess.SetCompMesh(0);
    fCMesh = copy.fCMesh;
    fGMesh = copy.fGMesh;
    fCMesh.SetReference(&fGMesh);
    fSolution = copy.fSolution;

    fPlasticDeformSqJ2 = copy.fPlasticDeformSqJ2;
    fHistoryLog = copy.fHistoryLog;

    fExecutionLog = copy.fExecutionLog;

    
    return *this;
}

void TPZSlopeAnalysis::ExecuteInitialSimulation2(int nsteps, int numnewton)
{
    std::stringstream strout;
    strout << "Initial Simulation";
    fCurrentConfig.fSolution.Redim(fCurrentConfig.fCMesh.NEquations(), 1);
    SaveConfig(strout);
    fSequence.push_back(fCurrentConfig);

    for (int istep = 1; istep <= nsteps; istep++) {
        REAL factor = istep*1./nsteps;
        ExecuteSimulation2(istep,factor);
    }
}

void TPZSlopeAnalysis::ExecuteSimulation2(int substeps, REAL factor)
{
    
    TConfig &LocalConfig = fCurrentConfig;
    
    TPZElastoPlasticAnalysis analysis(&LocalConfig.fCMesh,std::cout);
    
	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(TPZSlopeAnalysis::TConfig::gNumThreads);


	TPZStepSolver<REAL> step;
    step.SetDirect(ELDLt);
    
    long neq = LocalConfig.fCMesh.NEquations();
    TPZVec<long> activeEquations;
    analysis.GetActiveEquations(activeEquations);
    TPZEquationFilter filter(neq);
    filter.SetActiveEquations(activeEquations);
    full.EquationFilter() = filter;
    analysis.SetStructuralMatrix(full);

    //step.SetDirect(ECholesky);
	analysis.SetSolver(step);
    
    
    LocalConfig.fSolution.Redim(neq, 1);
    int NumIter = 50;
    bool linesearch = true;
    bool checkconv = false;
    REAL tol =1.e-5;
    bool conv;
    
    int ncycles = 1;

    
    std::stringstream strout;
    


	analysis.IterativeProcess(strout, tol, NumIter,linesearch,checkconv);

        
       
    
        
    TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
    for (int ieq=0; ieq<neq; ieq++) {
        LocalConfig.fSolution(ieq,0) = sol(ieq,0);
    }
    
    analysis.AcceptSolution();

    
    // dont forget the vertical deformation
    LocalConfig.ComputeElementDeformation();
    


//    SaveConfig(strout);
//    fCurrentConfig.fHistoryLog = strout.str();

    
    fSequence.push_back(LocalConfig);
}


void TPZSlopeAnalysis::ExecuteInitialSimulation(REAL factor)
{
    std::stringstream strout;
    strout << "Initial Simulation";
	fCurrentConfig.fSolution.Redim(fCurrentConfig.fCMesh.NEquations(), 1);
    SaveConfig(strout);
    fSequence.push_back(fCurrentConfig);
	ExecuteSimulation(factor,std::cout);
    
}


void TPZSlopeAnalysis::ExecuteSimulation(REAL factor, std::ostream &out)
{

	std::stringstream strout;
    strout << "Simulation";
    TConfig &LocalConfig = fCurrentConfig;
    
    TPZElastoPlasticAnalysis analysis(&LocalConfig.fCMesh,std::cout);
	
	LoadingRamp(factor);

	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(TPZSlopeAnalysis::TConfig::gNumThreads);
    analysis.SetStructuralMatrix(full);
	
	TPZStepSolver<REAL> step;
    step.SetDirect(ECholesky);
	analysis.SetSolver(step);
    
     long neq = LocalConfig.fCMesh.NEquations();
     TPZVec<long> activeEquations;
     analysis.GetActiveEquations(activeEquations);
     TPZEquationFilter filter(neq);
     filter.SetActiveEquations(activeEquations);
     full.EquationFilter() = filter;

    
    LocalConfig.fSolution.Redim(neq, 1);
    int NumIter = 200;
    bool linesearch = true;
    bool checkconv = false;
    REAL tol =1.e-3;
    bool conv;
    
    int ncycles;
	ncycles = 1;
    
    for (int cycle = 0; cycle < ncycles; cycle++)
    {
        analysis.IterativeProcess(out, tol, NumIter,linesearch,checkconv);
    }
        
    TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
    for (int ieq=0; ieq<neq; ieq++) {
        LocalConfig.fSolution(ieq,0) = sol(ieq,0);
    }
    
    analysis.AcceptSolution();

    // dont forget the vertical deformation
    LocalConfig.ComputeElementDeformation();
	SaveConfig(strout);
	fSequence.push_back(LocalConfig);

}
void TPZSlopeAnalysis::TConfig::ApplyDeformation(TPZCompEl *cel)
{
    TPZCompEl *cel2 = cel;
    TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
    if (!intel2) {
        DebugStop();
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        TPZGeoEl *gel = cel2->Reference();
        gel->Print(sout);
        REAL co[8][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
        for (int p = 0; p<8; p++) {
            TPZManVector<REAL,3> par(2,0.),x(3,0.);
            par[0] = co[p][0];
            par[1] = co[p][1];
            gel->X(par, x);
            sout << "point " << p << "co " << x << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZCompMesh *cmesh2 = cel2->Mesh();
    TPZGeoMesh *gmesh1 = &fGMesh;
    TPZCompMesh *cmesh1 = gmesh1->Reference();
    
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh2->MaterialVec()[1]);
    if (pMatWithMem2) {
        pMatWithMem2->SetUpdateMem(true);
    }
    else {
        DebugStop();
    }
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem1 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh1->MaterialVec()[1]);
    
    if (intel2->Material() != pMatWithMem2) {
        DebugStop();
    }
    const TPZIntPoints &intpoints = intel2->GetIntegrationRule();
    int nint = intpoints.NPoints();
    TPZManVector<REAL,3> point(3,0.);
    TPZMaterialData data1,data2;
    intel2->InitMaterialData(data2);
    data2.fNeedsSol = false;
    data1.fNeedsSol = true;
    long data2phir = data2.phi.Rows();
    long data2phic = data2.phi.Cols();
    long data2dphir = data2.dphix.Rows();
    long data2dphic = data2.dphix.Cols();
    long elementid = 0;
    TPZManVector<REAL,3> qsi(2,0.);

    for (long ip =0; ip<nint; ip++) {
        REAL weight;
        intpoints.Point(ip, point, weight);
        data2.intLocPtIndex = ip;
        intel2->ComputeRequiredData(data2, point);
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            int memoryindex = data2.intGlobPtIndex;
            std::stringstream sout;
            sout << "Local integration point index " << data2.intLocPtIndex << std::endl;
            pMatWithMem2->PrintMem(sout,memoryindex);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZGeoEl *gel1 = gmesh1->FindElement(data2.x, qsi, elementid,2);
        if (!gel1) {
            DebugStop();
        }
        TPZCompEl *cel1 = gel1->Reference();
        if (!cel1) {
            DebugStop();
        }
        TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *>(cel1);
        if (!intel1) {
            DebugStop();
        }
        intel1->InitMaterialData(data1);
        data1.fNeedsSol = true;
        intel1->ComputeRequiredData(data1, qsi);
#ifdef DEBUG
        {
            REAL diff = dist(data1.x,data2.x);
            if(diff > 1.e-6)
            {
                std::cout << "Point not found " << data2.x << std::endl;
                //DebugStop();
            }
        }
#endif
        data2.sol = data1.sol;
        data2.dsol = data1.dsol;
        data2.phi.Redim(0, 1);
        data2.dphix.Redim(2, 0);
        TPZFMatrix<STATE> ek,ef;
        pMatWithMem2->Contribute(data2, weight, ek, ef);
        data2.phi.Resize(data2phir, data2phic);
        data2.dphix.Resize(data2dphir, data2dphic);
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            int memoryindex = data2.intGlobPtIndex;
            std::stringstream sout;
            sout << "Local integration point index " << data2.intLocPtIndex << std::endl;
            sout << "qsi coordinate " << qsi << std::endl;
            pMatWithMem2->PrintMem(sout,memoryindex);
            sout << "\noriginal element index " << cel2->Index() << std::endl;
            sout << "projected element index " << cel1->Index() << std::endl;
            if(cel1->Index() == cel2->Index())
            {
                sout << "Same point index in the original mesh " << std::endl;
                sout << "qsi coordinate " << qsi << std::endl;
                pMatWithMem1->PrintMem(sout,memoryindex);
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

    }
    pMatWithMem2->SetUpdateMem(false);

}

/// Load the solution stored in TConfig into the CompMesh and set the ZDeformation of the material
void TPZSlopeAnalysis::TConfig::LoadSolution()
{

    fCMesh.LoadSolution(fSolution);
}


// Get the vector of element plastic deformations
void TPZSlopeAnalysis::TConfig::ComputeElementDeformation()
{
    long nelem = fCMesh.NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.);
    fCMesh.ElementSolution().Redim(nelem, 1);
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
        fPlasticDeformSqJ2.Fill(0.);
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
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
                //TPZTensor<REAL> &plastic = mem.m_elastoplastic_state.m_eps_p;
				TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                REAL J2 = plastic.J2();
                REAL sqj2 = sqrt(J2);
                sqj2el = std::max(sqj2,sqj2el);
            }
            fPlasticDeformSqJ2[el] = sqj2el;
        }
    }
    
    fCMesh.SetElementSolution(0, fPlasticDeformSqJ2);
}
void TPZSlopeAnalysis::ApplyHistory(std::set<long> &elindices)
{
    std::list<TConfig>::iterator listit;
    for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
        listit->LoadSolution();
    }
    
#ifdef DEBUG
    std::stringstream filename;
//    filename << "applyhistory_" << startfrom << ".txt";
    std::ofstream out(filename.str().c_str());
#endif
    
    std::set<long>::iterator it;
    for (it=elindices.begin(); it != elindices.end(); it++) {
        long elindex = *it;
        TPZCompEl *cel = fCurrentConfig.fCMesh.ElementVec()[elindex];
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            DebugStop();
        }
        // Reset the memory of the integration points of the element
        TPZManVector<long> pointindices;
        cel->GetMemoryIndices(pointindices);
        long npoints = pointindices.size();
        for (long ip = 0; ip<npoints; ip++) {
            long ind = pointindices[ip];
            pMatWithMem2->ResetMemItem(ind);
        }
        std::list<TConfig>::iterator listit;
        int confindex = 0;
        for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "configure index " << confindex;
                LOGPZ_DEBUG(logger, sout.str())
                
            }
#endif
            listit->ApplyDeformation(cel);
            confindex++;
        }
#ifdef DEBUG
        for (long ip = 0; ip<npoints; ip++) {
            long ind = pointindices[ip];
            pMatWithMem2->MemItem(ind).Print(out);
        }
#endif
    }
}

/// Change the polynomial order of element using the plastic deformation as threshold
void TPZSlopeAnalysis::TConfig::PRefineElementsAbove(REAL sqj2, int porder, std::set<long> &elindices)
{
    fGMesh.ResetReference();
    fCMesh.LoadReferences();
    long nelem = fCMesh.NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
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

        if (fCMesh.ElementSolution()(el,0) < sqj2) {
            continue;
        }
        TPZStack<long> subels;
        long index = cel->Index();
        elindices.insert(index);
        intel->SetPreferredOrder(porder);
    }
    fCMesh.AdjustBoundaryElements();
    fCMesh.InitializeBlock();
    fCMesh.Solution().Zero();
    fSolution.Resize(0, 0);
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    // force the post process mesh to be regenerated
    fPostprocess.SetCompMesh(0);
}
/// Divide the element using the plastic deformation as threshold
void TPZSlopeAnalysis::TConfig::DivideElementsAbove(REAL sqj2, std::set<long> &elindices)
{
    fGMesh.ResetReference();
    fCMesh.LoadReferences();
    TPZManVector<REAL,3> findel(3,0.),qsi(2,0.);
    findel[0] = 0.108;
    findel[1] = 0.0148;
    long elindex = 0;
    fCMesh.Reference()->FindElement(findel, qsi, elindex, 2);
    TPZGeoEl *targetel = fCMesh.Reference()->ElementVec()[elindex];
    TPZCompEl *targetcel = targetel->Reference();
    long targetindex = targetcel->Index();
    
    long nelem = fCMesh.NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
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
        if (fCMesh.ElementSolution()(el,0) < sqj2) {
            continue;
        }
        int porder = intel->GetPreferredOrder();
        TPZStack<long> subels;
        long index = cel->Index();
#ifdef LOG4CXX
        if (logger->isDebugEnabled() && investigate == true) {
            std::stringstream sout;
            cel->Reference()->Print(sout);
            for (int in=0; in< cel->Reference()->NCornerNodes(); in++) {
                cel->Reference()->NodePtr(in)->Print(sout);
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        intel->Divide(index, subels,0);
        for (int is=0; is<subels.size(); is++) {
            elindices.insert(subels[is]);
            TPZCompEl *subcel = fCMesh.ElementVec()[subels[is]];
#ifdef LOG4CXX
            if (logger->isDebugEnabled() && investigate == true) {
                std::stringstream sout;
                subcel->Reference()->Print(sout);
                for (int in=0; in< subcel->Reference()->NCornerNodes(); in++) {
                    subcel->Reference()->NodePtr(in)->Print(sout);
                }

                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
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
        long nelem = fCMesh.NElements();
        for (long el=0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
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
            TPZCompEl *cel = fCMesh.ElementVec()[el];
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
                TPZCompEl *subcel = fCMesh.ElementVec()[subels[is]];
                TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
                if (!subintel) {
                    DebugStop();
                }
                subintel->SetPreferredOrder(porder);
            }
        }
    }
    fCMesh.AdjustBoundaryElements();
    fCMesh.InitializeBlock();
    fCMesh.Solution().Zero();
    fSolution.Resize(0, 0);
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}
void TPZSlopeAnalysis::SaveConfig(std::stringstream &strout) {
    //Saving TConfig
    fCurrentConfig.fHistoryLog = strout.str();
    if(fCurrentConfig.fSolution.Cols() == 0)
    {
        fCurrentConfig.fSolution = fCurrentConfig.fCMesh.Solution();
    }
//#ifdef DEBUG
//    std::cout << "Before putting the current config in the list\n";
//    fCurrentConfig.PrintForcingFunction();
//#endif
    fSequence.push_back(fCurrentConfig);
//#ifdef DEBUG
//    std::cout << "After putting the current config in the list\n";
//    fSequence.rbegin()->PrintForcingFunction();
//#endif
}
void TPZSlopeAnalysis::TConfig::CreatePostProcessingMesh()
{

    if (fPostprocess.ReferenceCompMesh() != &fCMesh)
    {
        
        fPostprocess.SetCompMesh(&fCMesh);
        
        TPZVec<int> PostProcMatIds(1,1);
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        TPZSlopeAnalysis::PostProcessVariables(scalNames, vecNames);
        
        for (int i=0; i<scalNames.size(); i++) {
            PostProcVars.Push(scalNames[i]);
        }
        for (int i=0; i<vecNames.size(); i++) {
            PostProcVars.Push(vecNames[i]);
        }
        //
        fPostprocess.SetPostProcessVariables(PostProcMatIds, PostProcVars);
        TPZFStructMatrix structmatrix(fPostprocess.Mesh());
        structmatrix.SetNumThreads(TPZSlopeAnalysis::TConfig::gNumThreads);
        fPostprocess.SetStructuralMatrix(structmatrix);
    }
    //
    //Chamar com o analysis e nao com o postanalysis pois tem o acumulo de sols
    fPostprocess.TransferSolution();
    //fElastoplasticAnalysis.TransferSolution(fPostprocess);
}

void TPZSlopeAnalysis::CreatePostProcessingMesh(TPZCompMesh * cmesh, TPZPostProcAnalysis * PostProcess)
{

    if (PostProcess->ReferenceCompMesh() != cmesh)
    {
        
        PostProcess->SetCompMesh(cmesh);
        
        TPZVec<int> PostProcMatIds(1,1);
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        TPZSlopeAnalysis::PostProcessVariables(scalNames, vecNames);
        
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

void TPZSlopeAnalysis::LoadingRamp ( REAL factor)
{
	TConfig &LocalConfig = fCurrentConfig;
    TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > *  body = NULL;
	
    body = dynamic_cast<TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > *> ( fCurrentConfig.fCMesh.FindMaterial ( 1 ) );
    TPZManVector<STATE, 3> force ( 3,0. );
    force[1]=-factor*20.;
body->SetBulkDensity(1.);
    body->SetBodyForce ( force );

}

/// Get the post processing variables
void TPZSlopeAnalysis::PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames)
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



int passCount = 0;
void TPZSlopeAnalysis::PostProcess(int resolution)
{
	fCurrentConfig.CreatePostProcessingMesh();

    std::string vtkFile = fVtkFile;

    TPZStack<std::string> scalNames,vecNames;
    PostProcessVariables(scalNames,vecNames);

    fCurrentConfig.fPostprocess.DefineGraphMesh(2,scalNames,vecNames,vtkFile);

    fCurrentConfig.fPostprocess.SetStep(fPostProcessNumber);
    fCurrentConfig.fPostprocess.PostProcess(resolution);

    fPostProcessNumber++;

}

#include "TPZVTKGeoMesh.h"

void TPZSlopeAnalysis::TConfig::CreateGeometricMesh()
{
	const std::string name ( "ElastoPlastic GEOMESH Footing Problem " );


    fGMesh.SetName ( name );
    fGMesh.SetDimension ( 2 );
   fGMesh.NodeVec().Resize ( 11 );
    TPZVec<REAL> coord ( 2 );
    coord[0] = 0.;
    coord[1] = 0.;
    fGMesh.NodeVec() [0] = TPZGeoNode ( 0, coord, fGMesh );
    coord[0] = 75.;
    coord[1] = 0.;
   fGMesh.NodeVec() [1] = TPZGeoNode ( 1, coord, fGMesh );
    coord[0] =75.;
    coord[1] = 30.;
    fGMesh.NodeVec() [2] = TPZGeoNode ( 2, coord,fGMesh );
    coord[0] = 45.;
    coord[1] = 30.;
    fGMesh.NodeVec() [3] = TPZGeoNode ( 3, coord,fGMesh );
    coord[0] = 35.;
    coord[1] = 40.;
    fGMesh.NodeVec() [4] = TPZGeoNode ( 4, coord, fGMesh );
    coord[0] = 0.;
    coord[1] = 40.;
    fGMesh.NodeVec() [5] = TPZGeoNode ( 5, coord, fGMesh );


    coord[0] = 25.;
    coord[1] = 40.;
    fGMesh.NodeVec() [6] = TPZGeoNode ( 6, coord, fGMesh );

    coord[0] = 25.;
    coord[1] = 30.;
    fGMesh.NodeVec() [7] = TPZGeoNode ( 7, coord, fGMesh );

    coord[0] = 15.;
    coord[1] = 40.;
    fGMesh.NodeVec() [8] = TPZGeoNode ( 8, coord,fGMesh );

    coord[0] = 15.;
    coord[1] = 15.;
    fGMesh.NodeVec() [9] = TPZGeoNode ( 9, coord, fGMesh);

    coord[0] = 45.;
    coord[1] = 15.;
   fGMesh.NodeVec() [10] = TPZGeoNode ( 10, coord, fGMesh );


    TPZVec <long> TopoQuad ( 4 );
    TopoQuad[0] = 0;
    TopoQuad[1] = 1;
    TopoQuad[2] = 10;
    TopoQuad[3] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 0, TopoQuad, 1,fGMesh );

    TopoQuad[0] = 1;
    TopoQuad[1] = 2;
    TopoQuad[2] = 3;
    TopoQuad[3] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 1, TopoQuad, 1, fGMesh);

    TopoQuad[0] = 10;
    TopoQuad[1] = 3;
    TopoQuad[2] = 7;
    TopoQuad[3] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 2, TopoQuad, 1, fGMesh );

    TopoQuad[0] = 3;
    TopoQuad[1] = 4;
    TopoQuad[2] = 6;
    TopoQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 3, TopoQuad, 1, fGMesh );

    TopoQuad[0] = 7;
    TopoQuad[1] = 6;
    TopoQuad[2] = 8;
    TopoQuad[3] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 4, TopoQuad, 1, fGMesh);

    TopoQuad[0] = 9;
    TopoQuad[1] = 8;
    TopoQuad[2] = 5;
    TopoQuad[3] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( 5, TopoQuad, 1,fGMesh );

    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 6, TopoLine, -1, fGMesh );

    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 7, TopoLine, -2, fGMesh );

    TopoLine[0] = 0;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 8, TopoLine, -3, fGMesh );


    TopoLine[0] = 4;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 9, TopoLine, -4, fGMesh);
	fGMesh.BuildConnectivity();
	    for ( int d = 0; d<2; d++ ) {
        int nel = fGMesh.NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = fGMesh.ElementVec() [iel];
                gel->Divide ( subels );
        }
    }
	

    


}


void TPZSlopeAnalysis::TConfig::CreateMesh()
{
    if ((fGMesh.NElements() != 0) || (fGMesh.NNodes() != 0))
    {
        DebugStop();
    }


    TPZManVector<REAL,3> x0(3,fInnerRadius),x1(3,fOuterRadius);
    x0[1] = 0;
    x1[1] = M_PI_2;
    x0[2] = 0.;
    x1[2] = 0.;
//    int numdiv = 20;
    TPZManVector<int,2> nx(fNx);
//    nx[0] = 50;
    TPZGenGrid gengrid(nx,x0,x1);
    REAL minsize = fDelx;//fCurrentConfig.fInnerRadius*M_PI_2/numdiv;
    REAL domainsize = fOuterRadius-fInnerRadius;
    REAL geoprogression = gengrid.GeometricProgression(minsize, domainsize, nx[0]);
    TPZManVector<REAL,2> geoprogressionvec(2,1.);
    geoprogressionvec[0] = geoprogression;
    gengrid.SetGeometricProgression(geoprogressionvec);
    gengrid.Read(&fGMesh);
    /// bottom side
    gengrid.SetBC(&fGMesh, 4, EBottom);
    /// outer side
    gengrid.SetBC(&fGMesh, 5, EOuter);
    /// left side
    gengrid.SetBC(&fGMesh, 6, ELeft);
    /// inner radius
    gengrid.SetBC(&fGMesh, 7, EInner);
    
    /// wrap the mesh around
    int nnodes = fGMesh.NNodes();
    for (int in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> cyl(3,0.), cart(3,0.);
        fGMesh.NodeVec()[in].GetCoordinates(cyl);
        cart[0] = cyl[0]*cos(cyl[1]);
        cart[1] = cyl[0]*sin(cyl[1]);
        fGMesh.NodeVec()[in].SetCoord(cart);
    }

    
    std::ofstream ssout("WellboreBefore.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&fGMesh, ssout,true);
    
}


static long newnode(TPZGeoMesh &mesh, TPZVec<REAL> &coord)
{
    long newnode = mesh.NodeVec().AllocateNewElement();
    mesh.NodeVec()[newnode].Initialize(coord, mesh);
    return newnode;
}

/// Initialize the Sandler DiMaggio object and create the computational mesh
void TPZSlopeAnalysis::TConfig::CreateComputationalMesh(int porder)
{
      unsigned int dim  = fGMesh.Dimension();
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    // Setting up attributes
	fCMesh =  TPZCompMesh ( &fGMesh );
    fCMesh.SetName ( name );
    fCMesh.SetDefaultOrder ( porder );
    fCMesh.SetDimModel ( dim );

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

    TPZMatElastoPlasticSest2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new TPZMatElastoPlasticSest2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > ( 1,PlaneStrain );
    material->SetPlasticity ( LEMC );

	material->SetId(1);
    fCMesh.InsertMaterialObject ( material );

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 1;
    auto * bc_bottom = material->CreateBC ( material, -1,3, val1, val2 );
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_rigth = material->CreateBC ( material, -2, 3, val1, val2 );
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_left = material->CreateBC ( material, -3, 3, val1, val2 );

    fCMesh.InsertMaterialObject ( bc_bottom );
   	fCMesh.InsertMaterialObject ( bc_rigth );
    fCMesh.InsertMaterialObject ( bc_left );
	
    //cmesh->InsertMaterialObject ( top );
    fCMesh.SetAllCreateFunctionsContinuousWithMem();

    fCMesh.AutoBuild();

}
