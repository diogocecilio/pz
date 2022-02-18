//
//  WellBoreAnalysis.cpp
//  PZ
//
//  Created by phil on 1/18/13.
//  Copyright (c) 2013 LabMeC. All rights reserved.
//
#include "SlopeStabilityAnalysis.h"
#include "pzbndcond.h"
#include "pzelastoplasticanalysis.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplasticSest2D.h"
#include "tpzgeoelrefpattern.h"
#include "pzelasticSest2D.h"
#include "pzelasmat.h"
#include "pzvisualmatrix.h"


#include "tpzchangeel.h"
//#include "poroelastoplastic.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzplasticdiagnostic.h"
#include <iostream>
#include "pzbfilestream.h"
#include "pzelasmat.h" 
#include "pzl2projection.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZProjectEllipse.h"
#include "pzlog.h"
#include "TPZGeoLinear.h"
#include "arglib.h"
#include "run_stats_table.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.plasticity.wellboreanalysis"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse(Logger::getLogger("LogEllipse"));
#endif

int TPZSlopeStabilityAnalysis::TConfig::gNumThreads = 0;



TPZSlopeStabilityAnalysis::TPZSlopeStabilityAnalysis() : fCurrentConfig(), fSequence(), fPostProcessNumber(0), fLinearMatrix()
{
}

TPZSlopeStabilityAnalysis::TPZSlopeStabilityAnalysis(const TPZSlopeStabilityAnalysis &copy) : fCurrentConfig(copy.fCurrentConfig), fSequence(copy.fSequence), fPostProcessNumber(copy.fPostProcessNumber), fLinearMatrix(copy.fLinearMatrix)
{

}

TPZSlopeStabilityAnalysis &TPZSlopeStabilityAnalysis::operator=(const TPZSlopeStabilityAnalysis &copy)
{
    if (this == &copy) {
        return *this;
    }
    fCurrentConfig = copy.fCurrentConfig;
    fSequence = copy.fSequence;
    fPostProcessNumber = copy.fPostProcessNumber;
    fLinearMatrix = copy.fLinearMatrix;
    return *this;
}

TPZSlopeStabilityAnalysis::~TPZSlopeStabilityAnalysis()
{
    
}

/// write the object on the stream
void TPZSlopeStabilityAnalysis::Write(TPZStream &output)
{
    fCurrentConfig.Write(output);
    int seqsize = fSequence.size();
    output.Write(&seqsize);
    std::list<TPZSlopeStabilityAnalysis::TConfig>::iterator it;
    for (it = fSequence.begin(); it != fSequence.end(); it++) {
        it->Write(output);
    }
    output.Write(&fPostProcessNumber);
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCurrentConfig.fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
}

/// read the object from the stream
void TPZSlopeStabilityAnalysis::Read(TPZStream &input)
{
    fCurrentConfig.Read(input);
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCurrentConfig.fGMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int seqsize;
    input.Read(&seqsize);
    for (int i=0; i<seqsize; i++) {
        TPZSlopeStabilityAnalysis::TConfig config;
        config.Read(input);
        fSequence.push_back(config);
    }
    input.Read(&fPostProcessNumber);
    
}



void TPZSlopeStabilityAnalysis::CheckDeformation(std::string filename)
{
    
    std::ofstream out(filename.c_str());
    TPZTensor<STATE> epstotal,sigma, S, pressure ;
    sigma.XX() = 1.;
    sigma.XY() = 1.;
    STATE i1 = sigma.I1();
    pressure.Identity();
    pressure *= i1/3.;
    sigma.S(S);
    STATE sqj2 = sqrt(sigma.J2());
    S *= 1./sqj2;
    sqj2 = S.J2();
    const int nangles = 100;
    out << "deform = Table[0,{" << nangles << "}];" << std::endl;
    out << "tension = Table[0,{" << nangles << "}];" << std::endl;
    
    for (int i=0; i<nangles; i++) {
        
        TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > SD;
        SD.fYC.PreSMat(SD.fYC);


        const int nincrements = 100;
        TPZManVector<std::pair<STATE, STATE>,100 > epst(nincrements),sigvec(nincrements);
        STATE angle = i*M_PI/nangles;
        STATE ca = cos(angle);
        STATE sa = sin(angle);
        TPZTensor<STATE> deformbase,p1(pressure), s1(S);
        p1 *= ca;
        s1 *= sa;
        deformbase = p1;
        deformbase += s1;
        STATE catest = deformbase.I1();
        STATE satest = sqrt(deformbase.J2());
        catest -= ca;
        satest -= fabs(sa);
        for (int j=0; j< nincrements; j++) {
            STATE scale = 0.1*j/nincrements;
            TPZTensor<STATE> epstotal(deformbase);
            epstotal *= scale;
            TPZTensor<STATE> sigma;
            STATE i1 = epstotal.I1();
            STATE sqj2deform = sqrt(epstotal.J2());
            STATE diff[2] = {i1-ca*scale,sqj2deform-fabs(sa)*scale};
            if (fabs(diff[0]) > 1.e-6 || fabs(diff[1]) > 1.e-6) {
                DebugStop();
            }
            epst[j].first = ca*scale;
            epst[j].second = sa*scale;

            SD.ApplyStrainComputeSigma(epstotal,sigma);
            // print I1 and sqrt(J2)
            REAL i1sigma = sigma.I1();
            REAL sqj2 = sqrt(sigma.J2());
            sigvec[j].first = i1sigma;
            sigvec[j].second = sqj2;
        }
        out << "deform[[" << i+1 << "]] = " <<  epst << ";" << std::endl;
        out << "tension[[" << i+1 << "]] = " << sigvec << ";" << std::endl;
    }

}



TPZSlopeStabilityAnalysis::TConfig::TConfig() :
    fGMesh(), fCMesh(), fSolution(), fZDeformation(0.), fPlasticDeformSqJ2(), fHistoryLog(), fModel(ESandler)
  , fSDPV(), fMCPV(), fMatEla()


{
    fCMesh.SetReference(&fGMesh);
}

TPZSlopeStabilityAnalysis::TConfig::TConfig(const TConfig &conf) : 
        fGMesh(conf.fGMesh), fCMesh(conf.fCMesh), fSolution(conf.fSolution), fZDeformation(conf.fZDeformation),
    fModel(conf.fModel), fPlasticDeformSqJ2(conf.fPlasticDeformSqJ2), fHistoryLog(conf.fHistoryLog), fSDPV(conf.fSDPV), fMCPV(conf.fMCPV), fMatEla(conf.fMatEla)
{
    fSDPV = conf.fSDPV;
    fGMesh.ResetReference();
    fCMesh.SetReference(&fGMesh);
    fCMesh.LoadReferences();
//    std::cout << "Copy constructor from " << (void *) &conf << " to " << (void *) this << std::endl;

}

TPZSlopeStabilityAnalysis::TConfig::~TConfig()
{
//    std::cout << "Deleting config " << (void *) this << std::endl;
    fPostprocess.SetCompMesh(0);
    fCMesh.CleanUp();
    fCMesh.SetReference(0);
}

TPZSlopeStabilityAnalysis::TConfig &TPZSlopeStabilityAnalysis::TConfig::operator=(const TPZSlopeStabilityAnalysis::TConfig &copy)
{
    if (this == &copy) {
        return *this;
    }

    fSDPV = copy.fSDPV;
    fMCPV = copy.fMCPV;
    fMatEla = copy.fMatEla;



    fPostprocess = copy.fPostprocess;
    fPostprocess.SetCompMesh(0);
    fCMesh = copy.fCMesh;
    fGMesh = copy.fGMesh;
    fCMesh.SetReference(&fGMesh);
    fSolution = copy.fSolution;
    fZDeformation = copy.fZDeformation;
    fPlasticDeformSqJ2 = copy.fPlasticDeformSqJ2;
    fHistoryLog = copy.fHistoryLog;
    fModel = copy.fModel;


    
    std::cout << "Copying " << (void *) &copy << " to " << (void *) this << std::endl;
    return *this;
}

/// Write the data to the output stream
void TPZSlopeStabilityAnalysis::TConfig::Write(TPZStream &out)
{


    fSDPV.Write(out);
    fMCPV.Write(out);
    fMatEla.Write(out,0); //AQUIPHIL

		int IntEPlasticModel = fModel;
		out.Write(&IntEPlasticModel);

    fGMesh.Write(out, 0);
    fCMesh.Write(out, 0);
    fSolution.Write(out, 0);
    out.Write(&fZDeformation);
    TPZSaveable::WriteObjects(out,fPlasticDeformSqJ2);
    out.Write(&fHistoryLog);

    int verify = 83562;
    out.Write(&verify);
}

/// Read the data from the input stream
void TPZSlopeStabilityAnalysis::TConfig::Read(TPZStream &input)
{

    int wellconf;
    input.Read(&wellconf);

    fSDPV.Read(input);
    fMCPV.Read(input);
    fMatEla.Read(input,0); //AQUIPHIL

    int IntEPlasticModel;
    input.Read(&IntEPlasticModel);
    fModel = (EPlasticModel) IntEPlasticModel;


    fGMesh.Read(input, 0);
    fCMesh.Read(input, &fGMesh);
    fSolution.Read(input, 0);
    input.Read(&fZDeformation);
    TPZSaveable::ReadObjects(input,fPlasticDeformSqJ2);
    input.Read(&fHistoryLog);
    
    int verify = 0;
    input.Read(&verify);
    if (verify != 83562)
    {
        DebugStop();
    }
}

RunStatsTable well_init("-tpz_well_init", "Raw data table statistics for TPZElastoPlasticAnalysis::IterativeProcess");




void TPZSlopeStabilityAnalysis::ExecuteInitialSimulation(int nsteps)
{
    std::stringstream strout;
    strout << "Initial Simulation";
    fCurrentConfig.fSolution.Redim(fCurrentConfig.fCMesh.NEquations(), 1);
    SaveConfig(strout);
    fSequence.push_back(fCurrentConfig);

    for (int istep = 1; istep <= nsteps; istep++) {
        REAL factor = istep*1./nsteps;
        ExecuteSimulation(istep,factor);
    }
}



void TPZSlopeStabilityAnalysis::ExecuteSimulation(int substeps, REAL factor)
{
    
    TConfig &LocalConfig = fCurrentConfig;
    
    TPZElastoPlasticAnalysis analysis(&LocalConfig.fCMesh,std::cout);
    
	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(TPZSlopeStabilityAnalysis::TConfig::gNumThreads);


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

	well_init.start();
	analysis.IterativeProcess(std::cout, tol, NumIter,linesearch,checkconv);
	well_init.stop();
        
    TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
    for (int ieq=0; ieq<neq; ieq++) {
        LocalConfig.fSolution(ieq,0) = sol(ieq,0);
    }
    
    analysis.AcceptSolution();
    
    // dont forget the vertical deformation
    LocalConfig.ComputeElementDeformation();
    
    fSequence.push_back(LocalConfig);
}




/// Load the solution stored in TConfig into the CompMesh and set the ZDeformation of the material
void TPZSlopeStabilityAnalysis::TConfig::LoadSolution()
{
//     TPZMaterial *mat = fCMesh.FindMaterial(1);
// //     typedef TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> mattype1;
// //     typedef TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem> mattype2;
// //     typedef TPZElasticityMaterialSest2D mattype3;
// 	 typedef TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> mattype1;
//      typedef TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem> mattype2;
//      typedef TPZElasticityMaterial mattype3;
//     
//     mattype1 *matposs1 = dynamic_cast<mattype1 *>(mat);
//     mattype2 *matposs2 = dynamic_cast<mattype2 *>(mat);
//     mattype3 *matposs3 = dynamic_cast<mattype3 *>(mat);
//     
//     if (matposs1) {
//         matposs1->SetZDeformation(fZDeformation);
//     }
//     if (matposs2) {
//         matposs2->SetZDeformation(fZDeformation);
//     }
//     if (matposs3) {
//         matposs3->SetZDeformation(fZDeformation);
//     }
//     if (!matposs1 && !matposs2 && !matposs3) {
//         DebugStop();
//     }
// #ifdef DEBUG
//     if (fSolution.Cols() > 1) {
//         DebugStop();
//     }
// #endif
    fCMesh.LoadSolution(fSolution);
}




/// this method will modify the boundary condition of the computational mesh and the forcing function


/// Apply the deformation of the configuration to the element
void TPZSlopeStabilityAnalysis::TConfig::ApplyDeformation(TPZCompEl *cel)
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
#ifdef LOG4CXX2
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

/// verify if the stress tangent is computed correctly
void TPZSlopeStabilityAnalysis::VerifyTangentValidity()
{
    TPZCompMesh &cmesh = fCurrentConfig.fCMesh;
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh.Solution().Print("Mesh solution ",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    int nelem = cmesh.NElements();
    for (int el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->MaterialId() != 1) {
            continue;
        }
        fCurrentConfig.VerifyPlasticTangent(cel);
    }
}


/// Verify tangent elasto plastic relation
void TPZSlopeStabilityAnalysis::TConfig::VerifyPlasticTangent(TPZCompEl *cel)
{
    TPZCompEl *cel2 = cel;
    TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
    if (!intel2) {
        DebugStop();
    }
    TPZCompMesh *cmesh2 = cel2->Mesh();
    
    
 //   TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *pMatWithMem2 =
 //   dynamic_cast<TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *>(cmesh2->MaterialVec()[1]);

    
	TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *pMatWithMem2 =
    dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *>(cmesh2->MaterialVec()[1]);

	
    if (intel2->Material() != pMatWithMem2) {
        DebugStop();
    }
    const TPZIntPoints &intpoints = intel2->GetIntegrationRule();
    int nint = intpoints.NPoints();
    TPZManVector<REAL,3> point(3,0.);
    TPZMaterialData data2;
    intel2->InitMaterialData(data2);
    data2.fNeedsSol = true;
    TPZManVector<REAL,3> qsi(2,0.);
    std::stringstream sout;
    sout << "Diagnostic for element " << cel->Index() << std::endl;
    
    
    for (int ip =0; ip<nint; ip++) {
        REAL weight;
        intpoints.Point(ip, point, weight);
        data2.intLocPtIndex = ip;
        intel2->ComputeRequiredData(data2, point);
        TPZFNMatrix<36,REAL> deltastrain(3,1),stress(3,1),dstressdstrain(3,3);
        pMatWithMem2->ComputeDeltaStrainVector(data2, deltastrain);
        pMatWithMem2->ApplyDeltaStrainComputeDep(data2, deltastrain, stress, dstressdstrain);
        TPZFMatrix<REAL> variation(deltastrain);
        variation *= -1/10.;
        TPZManVector<STATE,20> errors(9),convrate(8,0.);
        REAL minerror = 0.;
        for (int ist=1; ist<10; ist++) {
            TPZFNMatrix<6,REAL> deltastrainnext(3,1),stressnext(3,1),stressestimate(3,1),stresserror(3,1);
            for (int i=0; i<3; i++) deltastrainnext(i,0) = deltastrain(i,0)+ist*variation(i,0);
            pMatWithMem2->ApplyDeltaStrain(data2, deltastrainnext, stressnext);
            for (int i=0; i<3; i++) {
                stressestimate(i,0) = stress(i,0);
                for (int j=0; j<3; j++) {
                    stressestimate(i,0) += dstressdstrain(i,j)*ist*variation(j,0);
                }
            }
            for (int i=0; i<3; i++) stresserror(i,0) = stressnext(i,0)-stressestimate(i,0);
            errors[ist-1] = Norm(stresserror);
            if(ist == 1) minerror = errors[ist-1];
            if (minerror > errors[ist-1]) {
                minerror = errors[ist-1];
            }
        }
        if(minerror > 1.e-12)
        {
            for (int ist=2; ist<10; ist++) {
                convrate[ist-2] = log(errors[ist-1]/errors[0])/log(ist*1.);
            }
        }
        sout << "Integration point " << ip << " coordinate " << data2.x << std::endl;
        sout << "Errors " << errors << std::endl;
        sout << "Convergence rates " << convrate << std::endl;
    }
#ifdef LOG4CXX
    LOGPZ_DEBUG(logger, sout.str())
#else
    std::cout << sout.str();
#endif
}



void TPZSlopeStabilityAnalysis::TransferSolutionTo(TConfig &config)
{
    TPZCompMesh &cmesh2 = config.fCMesh;
    TPZCompMesh &cmesh1 = fCurrentConfig.fCMesh;
    
    cmesh1.Solution() = fCurrentConfig.fSolution;
    
    long nel = cmesh2.NElements();
//    TPZGeoMesh *gmesh1 = cmesh1.Reference();
    TPZMaterial *mat1 = cmesh1.FindMaterial(1);
    if (!mat1) {
        DebugStop();
    }
    
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh2.MaterialVec()[1]);
    if (pMatWithMem2) {
        pMatWithMem2->SetUpdateMem(true);
    }
    else {
        DebugStop();
    }

    
    int varindex1 = mat1->VariableIndex("Stress");
    if (varindex1 == -1) {
        DebugStop();
    }
//    int elementid = 0;
    TPZManVector<REAL,3> qsi(3,0.);
    
    for (long el =0; el<nel; el++) {
        TPZCompEl *cel2 = cmesh2.ElementVec()[el];
        if (!cel2) {
            continue;
        }
        TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
        if (!intel2) {
            continue;
        }
        if (intel2->Material() != pMatWithMem2) {
            continue;
        }
        fCurrentConfig.ApplyDeformation(intel2);
    }
    pMatWithMem2->SetUpdateMem(false);
}

/// Compute the area of the domain at which sqJ2 is above a given value
REAL TPZSlopeStabilityAnalysis::TConfig::ComputeAreaAboveSqJ2(REAL sqj2)
{
    REAL area = 0.;
    long nelem = fCMesh.NElements();
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel) {
                continue;
            }
            bool shouldanalyse = false;
            TPZManVector<long> memindices;
            cel->GetMemoryIndices(memindices);
            long numind = memindices.size();
  //          REAL sqj2el = 0.;
            for (long ind=0; ind<numind; ind++) 
            {
                int memoryindex = memindices[ind];
                if (memoryindex < 0) {
                    continue;
                }
                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ind]);
                TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                REAL J2 = plastic.J2();
                REAL sqj2pt = sqrt(J2);
                if (sqj2pt > sqj2) {
                    shouldanalyse = true;
                }
            }
            if (shouldanalyse) {
                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                TPZGeoEl *gel = cel->Reference();
                if (! intel || !gel) {
                    DebugStop();
                }
                TPZIntPoints &rule = intel->GetIntegrationRule();
                int np = rule.NPoints();
                for (int ip = 0; ip<np; ip++) {
                    TPZManVector<REAL,3> point(2,0.);
                    REAL weight;
                    rule.Point(ip, point, weight);
                    TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
                    TPZFNMatrix<9,REAL> axes(2,3);
                    REAL detjac;
                    gel->Jacobian(point, jac, axes, detjac, jacinv);
                    TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ip]);
                    TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                    REAL J2 = plastic.J2();
                    REAL sqj2pt = sqrt(J2);
                    if (sqj2pt > sqj2) {
                        area += weight*fabs(detjac);
                    }
                }
            }
        }
    }
    return area;
}

/// Compute the area of the domain
REAL TPZSlopeStabilityAnalysis::TConfig::ComputeTotalArea()
{
    REAL area = 0.;
    long nelem = fCMesh.NElements();
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if (gel->MaterialId() != 1) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZIntPoints &rule = intel->GetIntegrationRule();
            int np = rule.NPoints();
            for (int ip = 0; ip<np; ip++) {
                TPZManVector<REAL,3> point(2,0.);
                REAL weight;
                rule.Point(ip, point, weight);
                TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
                TPZFNMatrix<9,REAL> axes(2,3);
                REAL detjac;
                gel->Jacobian(point, jac, axes, detjac, jacinv);
                area += weight*fabs(detjac);
            }
        }
    }
    return area;    
}


// Get the vector of element plastic deformations
void TPZSlopeStabilityAnalysis::TConfig::ComputeElementDeformation()
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
                TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                REAL J2 = plastic.J2();
                REAL sqj2 = sqrt(J2);
                sqj2el = max(sqj2,sqj2el);
            }
            fPlasticDeformSqJ2[el] = sqj2el;
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Element deformation " << fPlasticDeformSqJ2;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    fCMesh.SetElementSolution(0, fPlasticDeformSqJ2);
}

/// Verify the global equilibrium of the forces by boundary condition
void TPZSlopeStabilityAnalysis::TConfig::VerifyGlobalEquilibrium(std::ostream &out)
{
    long neq = fCMesh.NEquations();
    
    TPZFMatrix<STATE> rhsDirichlet(neq,1,0.);
    this->ComputeRhsExceptMatid(-4, rhsDirichlet);
//    this->ComputeRhsExceptMatid(-5, rhsDirichlet);
    
    TPZStack<int,10> boundaryconditionlist;
    std::set<int> allmaterials;
    std::map<int,TPZMaterial *>::iterator itmat;
    std::map<int,TPZMaterial *> &matmap = fCMesh.MaterialVec();
    for (itmat = matmap.begin(); itmat != matmap.end(); itmat++) {
        TPZMaterial *mat = itmat->second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        int matid = itmat->first;
        allmaterials.insert(matid);
        if (matid == -4 || matid == -5) {
            continue;
        }
        if (bnd) {
            boundaryconditionlist.Push(itmat->first);
        }
    }
    std::map<int, TPZManVector<REAL,2> > ForceResultantsExcept, ForceResultantsOnly;
    for(int ib=0; ib<boundaryconditionlist.size(); ib++)
    {
        int bc = boundaryconditionlist[ib];
        {
            TPZFMatrix<STATE> rhs(neq,1,0.);
            this->ComputeRhsForMatid(bc, rhs);
            this->ComputeXYForce(rhs , ForceResultantsOnly[ib]);
        }
        {
            TPZFMatrix<STATE> rhs(neq,1,0.);
            this->ComputeRhsExceptMatid(bc, rhs);
            rhs -= rhsDirichlet;
            this->ComputeXYForce(rhs , ForceResultantsExcept[ib]);
        }
        
    }
    
    TPZManVector<REAL,2> totalforce(2,0.), dirichletforce(2,0.);
    ComputeXYForce(rhsDirichlet, dirichletforce);
    totalforce[0] = -dirichletforce[0];
    totalforce[1] = -dirichletforce[1];
    out <<  "FORCE RESULTANTS\n";
    out << "Dirichlet force " << dirichletforce << std::endl;
    for (int ib=0; ib<boundaryconditionlist.size(); ib++) {
        int bc = boundaryconditionlist[ib];
        if (bc<0) {
            totalforce[0] += ForceResultantsOnly[ib][0];
            totalforce[1] += ForceResultantsOnly[ib][1];
        }
        out << "Force computing only BC = " << bc << " x-force = " << ForceResultantsOnly[ib][0] << " y-force = " << ForceResultantsOnly[ib][1] << std::endl;
        out << "Force computing except BC = " << bc << " x-force = " << ForceResultantsExcept[ib][0] << " y-force = " << ForceResultantsExcept[ib][1] << std::endl;
    }
    out << "Force resultants " << totalforce << std::endl;
}

/// Verify the global equilibrium of the forces by boundary condition
void TPZSlopeStabilityAnalysis::TConfig::VerifyGlobalEquilibrium2(std::ostream &out)
{
    TPZStack<int,10> boundaryconditionlist;
    std::set<int> allmaterials;
    std::map<int,TPZMaterial *>::iterator itmat;
    std::map<int,TPZMaterial *> &matmap = fCMesh.MaterialVec();
    for (itmat = matmap.begin(); itmat != matmap.end(); itmat++) {
        TPZMaterial *mat = itmat->second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        allmaterials.insert(itmat->first);
        if (bnd) {
            boundaryconditionlist.Push(itmat->first);
        }
    }
    if (boundaryconditionlist[0] != -6 || boundaryconditionlist[2] != -4) {
        DebugStop();
    }
    boundaryconditionlist[0] = -4;
    boundaryconditionlist[2] = -6;
    TPZManVector<TPZManVector<REAL,2>, 10> ForceResultants(boundaryconditionlist.size()+1);
    int neq = fCMesh.NEquations();
    TPZFMatrix<STATE> rhs4(neq,1,0.),rhs5(neq,1,0.);
    for(int ib=0; ib<boundaryconditionlist.size(); ib++)
    {
        int bc = boundaryconditionlist[ib];
        TPZFStructMatrix str(&fCMesh);
        std::set<int> matset(allmaterials);
        matset.erase(bc);
        str.SetMaterialIds(matset);
        TPZFMatrix<STATE> rhs(neq,1,0.);
        //        if (bc > -4) {
        //            matset.erase(-5);
        //            matset.erase(-4);
        //        }
        str.Assemble(rhs,0);
        // correct for the Dirichlet boundary conditions
        rhs -= rhs4;
        rhs -= rhs5;
        //        if (bc == -5) {
        //            rhs5 = rhs;
        //        }
        //        if (bc == -4) {
        //            rhs4 = rhs;
        //        }
        std::set<long> xeqs,yeqs;
        for (long el=0; el<fCMesh.NElements(); el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel || !cel->Reference()) {
                continue;
            }
            int matid = cel->Reference()->MaterialId();
            if (matid != bc) {
                continue;
            }
            int nc = cel->Reference()->NCornerNodes();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                long seq = c.SequenceNumber();
                long pos = fCMesh.Block().Position(seq);
                xeqs.insert(pos);
                yeqs.insert(pos+1);
            }
        }
        ForceResultants[ib].Resize(2,0.);
        std::set<long>::iterator it;
        for(it = xeqs.begin(); it != xeqs.end(); it++)
        {
            ForceResultants[ib][0] += rhs(*it,0);
            if (bc == -5) {
                rhs5(*it,0) = rhs(*it,0);
            }
        }
        for(it = yeqs.begin(); it != yeqs.end(); it++)
        {
            ForceResultants[ib][1] += rhs(*it,0);
            if (bc == -4) {
                rhs4(*it,0) = rhs(*it,0);
            }
        }
    }
    REAL normrhs = 0.;
    {
        boundaryconditionlist.Push(1);
        int ib = boundaryconditionlist.size()-1;
        ForceResultants[ib].Resize(2,0.);
        TPZFStructMatrix str(&fCMesh);
        TPZFMatrix<STATE> rhs(neq,1,0.);
        str.Assemble(rhs, 0);
        // zerar os residuos referente a connects que nao sao de vertice
        for (long el=0; el<fCMesh.NElements(); el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel || !cel->Reference()) {
                continue;
            }
            int ncorner = cel->Reference()->NCornerNodes();
            long nc = cel->NConnects();
            for (long ic=ncorner; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                long seq = c.SequenceNumber();
                long pos = fCMesh.Block().Position(seq);
                int blsize = fCMesh.Block().Size(seq);
                for (int ibl =0; ibl<blsize; ibl++) {
                    rhs(pos+ibl,0) = 0.;
                }
            }
        }
        
        
        rhs -= rhs4;
        rhs -= rhs5;
        normrhs = Norm(rhs);
        ForceResultants[ib].Resize(2, 0.);
        for (long i=0; i<neq/2; i++) {
            ForceResultants[ib][0] += rhs(2*i,0);
            ForceResultants[ib][1] += rhs(2*i+1,0);
        }
    }
    
    TPZManVector<REAL,2> totalforce(2,0.);
    out <<  "FORCE RESULTANTS\n";
    for (int ib=0; ib<boundaryconditionlist.size(); ib++) {
        int bc = boundaryconditionlist[ib];
        if (bc<0) {
            totalforce[0] += ForceResultants[ib][0];
            totalforce[1] += ForceResultants[ib][1];
        }
        out << "Boundary condition " << bc << " x-force = " << ForceResultants[ib][0] << " y-force = " << ForceResultants[ib][1] << std::endl;
    }
    out << "Force resultants " << totalforce << std::endl;
    out << "Norm rhs " << normrhs << std::endl;
}

/// Zera os componentes do rhs para connects diferentes do zero
void TPZSlopeStabilityAnalysis::TConfig::FilterRhs(TPZFMatrix<STATE> &rhs)
{
    // zerar os residuos referente a connects que nao sao de vertice
    for (long el=0; el<fCMesh.NElements(); el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel || !cel->Reference()) {
            continue;
        }
        int ncorner = cel->Reference()->NCornerNodes();
        int nc = cel->NConnects();
        for (int ic=ncorner; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if(c.HasDependency())
            {
                continue;
            }
            long seq = c.SequenceNumber();
            long pos = fCMesh.Block().Position(seq);
            int blsize = fCMesh.Block().Size(seq);
            for (int ibl =0; ibl<blsize; ibl++) {
                rhs(pos+ibl,0) = 0.;
            }
        }
    }
    
}



/// Compute the Rhs for the mesh minus the elements with matid
// this method is cumulative (sums to the rhs)
void TPZSlopeStabilityAnalysis::TConfig::ComputeRhsExceptMatid(int matid, TPZFMatrix<STATE> &rhs)
{
    std::set<int> allmaterials;
    std::map<int,TPZMaterial *>::iterator itmat;
    std::map<int,TPZMaterial *> &matmap = fCMesh.MaterialVec();
    for (itmat = matmap.begin(); itmat != matmap.end(); itmat++) {
        allmaterials.insert(itmat->first);
    }
    allmaterials.erase(matid);
    TPZFStructMatrix str(&fCMesh);
    str.SetNumThreads(TPZSlopeStabilityAnalysis::TConfig::gNumThreads);
    str.SetMaterialIds(allmaterials);
    str.Assemble(rhs, 0);
    FilterRhs(rhs);

    
}

/// Compute the contribution of only matid
// this method is cumulative (sums to the rhs)
void TPZSlopeStabilityAnalysis::TConfig::ComputeRhsForMatid(int matid, TPZFMatrix<STATE> &rhs)
{
    std::set<int> allmaterials;
    allmaterials.insert(matid);
    TPZFStructMatrix str(&fCMesh);
    str.SetMaterialIds(allmaterials);
    str.SetNumThreads(TPZSlopeStabilityAnalysis::TConfig::gNumThreads);
    str.Assemble(rhs, 0);
    FilterRhs(rhs);
}

/// Compute the resultant x and y force
void TPZSlopeStabilityAnalysis::TConfig::ComputeXYForce(TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force)
{
    force.Resize(2, 0.);
    force.Fill(0.);
    long nel = rhs.Rows();
    for (long i=0; i<nel; i++) {
        force[i%2] += rhs(i,0);
    }
}


void TPZSlopeStabilityAnalysis::TConfig::DeleteElementsAbove(REAL sqj2)
{
    fCMesh.Reference()->ResetReference();
    fCMesh.LoadReferences();
    long ndel = 0;
    long nelem = fCMesh.NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        if (fPlasticDeformSqJ2[el] > sqj2) {
            TPZGeoEl *gel = cel->Reference();
            // delete the neighbouring boundary elements
            int ncorner = gel->NCornerNodes();
            int ns = gel->NSides();
            if(gel->Index() == 1212)
            {
                std::cout << __PRETTY_FUNCTION__ << " I should stop\n";
                gel->Print(std::cout);
            }
            for (int is=ncorner; is<ns; is++) {
                TPZGeoElSide neighbour = gel->Neighbour(is);
                TPZGeoEl *neighgel = neighbour.Element();
                TPZCompEl *neighcel = neighgel->Reference();
                long neighcelindex = 0;
                if(neighcel) neighcelindex = neighcel->Index();
                if (neighgel->Dimension() == 1 && neighcel) {
                    delete neighcel;
                    neighgel->SetMaterialId(50);
                }
                if(gel->Index() == 1212)
                {
                    std::cout << __PRETTY_FUNCTION__ << " I should stop\n";
                    neighgel->Print(std::cout);
                    std::cout << fCMesh.ElementVec()[neighcelindex] << std::endl;
                }
            }
            delete cel;
            gel->SetMaterialId(50);
            ndel++;
        }
        //gel->ResetReference();
    }
    // put boundary conditions on the sides which have no neighbours
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZStack<TPZCompElSide> cneighbours;
            TPZGeoElSide gelside(gel,is);
            gelside.ConnectedCompElementList(cneighbours, 1, 0);
            if (cneighbours.size() == 0) {
                if (gel->Dimension() == 1) {
                    delete cel;
                    gel->SetMaterialId(50);
                }
                else if (gel->Dimension() == 2) {
                    // create a boundary condition
                    TPZGeoElBC gbc(gelside, -6);
                    TPZGeoEl *bc = gbc.CreatedElement();
                    // create the corresponding computational element
                    long index;
                    fCMesh.CreateCompEl(bc, index);
                }
            }
        }
    }
    TPZMaterial *mat = fCMesh.FindMaterial(-6);
    if (!mat) {
        DebugStop();
    }
    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
    if (!bnd) {
        DebugStop();
    }
    bnd->Val1()(0,0) = 1.e8;
    bnd->Val1()(1,1) = 1.e8;
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef DEBUG
    {
        std::ofstream file("AdjustedMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fCMesh.Reference(), file, true);
    }
#endif
    
    std::cout << "Number of elements deleted " << ndel << std::endl;
}




/// Change the polynomial order of element using the plastic deformation as threshold
void TPZSlopeStabilityAnalysis::TConfig::PRefineElementsAbove(REAL sqj2, int porder, std::set<long> &elindices)
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
void TPZSlopeStabilityAnalysis::TConfig::DivideElementsAbove(REAL sqj2, std::set<long> &elindices)
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

TPZGeoMesh * TPZSlopeStabilityAnalysis::TConfig::GetGeoMesh() {
  return &this->fGMesh;
}

/// Return the mesh used for computations (multiphysics mesh or fCMesh)
TPZCompMesh *TPZSlopeStabilityAnalysis::TConfig::CompMeshUsed()
{
    TPZCompMesh *cmesh = &fCMesh;
    return cmesh;
}


/// Reset the plastic memory of the integration points of these elements
void TPZSlopeStabilityAnalysis::ApplyHistory(std::set<long> &elindices)
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




void TPZSlopeStabilityAnalysis::TConfig::CreatePostProcessingMesh()
{
    if (fPostprocess.ReferenceCompMesh() != &fCMesh)
    {

        fPostprocess.SetCompMesh(&fCMesh);
        TPZFStructMatrix structmatrix(fPostprocess.Mesh());
        structmatrix.SetNumThreads(0);
        fPostprocess.SetStructuralMatrix(structmatrix);
        
        TPZVec<int> PostProcMatIds(1,1);
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        TPZSlopeStabilityAnalysis::PostProcessVariables(scalNames, vecNames);

        for (int i=0; i<scalNames.size(); i++) {
            PostProcVars.Push(scalNames[i]);
        }
        for (int i=0; i<vecNames.size(); i++) {
            PostProcVars.Push(vecNames[i]);
        }
        //
        fPostprocess.SetPostProcessVariables(PostProcMatIds, PostProcVars);
    }
    //
    fPostprocess.TransferSolution();
    
}

/// Get the post processing variables
void TPZSlopeStabilityAnalysis::PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames)
{
  scalNames.Resize(0);
  vecNames.Resize(0);
  

vecNames.Push("DisplacementTotal");

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
scalNames.Push("YieldSurface1");
scalNames.Push("YieldSurface2");
scalNames.Push("YieldSurface3");
scalNames.Push("POrder");
scalNames.Push("NSteps");

  
}


int passCount = 0;
void TPZSlopeStabilityAnalysis::PostProcess(int resolution)
{
    fCurrentConfig.CreatePostProcessingMesh();

    std::string vtkFile = fVtkFile;

    TPZStack<std::string> scalNames,vecNames;
    PostProcessVariables(scalNames,vecNames);

    fCurrentConfig.fPostprocess.DefineGraphMesh(2,scalNames,vecNames,vtkFile);

    fCurrentConfig.fPostprocess.SetStep(fPostProcessNumber);
    fCurrentConfig.fPostprocess.PostProcess(resolution);

    fPostProcessNumber++;

      
#ifdef TESTPOLYCHAIN
    //*** DANGER (definido pelo programador) ***
    //******************************************
    REAL J2val = 0.0004;
    std::multimap<REAL,REAL> polygonalChain;

    GetJ2Isoline(J2val, polygonalChain);
    
    {
        std::stringstream nm;
        nm << "EllipDots" << passCount << ".nb";
        std::ofstream outEllips(nm.str().c_str());
        passCount++;
        std::multimap<REAL,REAL>::iterator it;
        outEllips << "ellipDots={";
        for(it = polygonalChain.begin(); it != polygonalChain.end(); it++)
        {
            outEllips << "{" << it->first << " , " << it->second << "}";
            if(it != polygonalChain.end())
            {
                outEllips << ",";
            }
        }
        outEllips << "};\n";
        outEllips << "aa=ListPlot[ellipDots,Joined->False,AspectRatio->Automatic];\n";
        outEllips << "bb=Graphics[Circle[{0,0},4.25*0.0254]];\n";
        outEllips << "Show[aa,bb]\n";
    }
#endif
}




/// Divide the element using the plastic deformation as threshold
unsigned int TPZSlopeStabilityAnalysis::DivideElementsAbove(REAL sqj2)
{
    std::set<long> elindices;
    fCurrentConfig.DivideElementsAbove(sqj2,elindices);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) 
    {
        std::stringstream sout;
        sout << "Element indices that have been created ";
        for (std::set<long>::iterator it = elindices.begin(); it!= elindices.end(); it++) {
            sout << *it << " ";
        }
        std::cout << sout.str() << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCurrentConfig.fGMesh.Print(sout);    
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCurrentConfig.fCMesh.Print(sout);    
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // subject the integration points with the deformation history
    ApplyHistory(elindices);
    fCurrentConfig.ComputeElementDeformation();
    
    fCurrentConfig.fCMesh.Solution().Zero();
    fCurrentConfig.fSolution = fCurrentConfig.fCMesh.Solution();
    
    // invalidate the computational mesh associated with the postprocess mesh
    fCurrentConfig.fPostprocess.SetCompMesh(0);

    return elindices.size();
}

/// GetPostProcessedValues at a given coordinate x
void TPZSlopeStabilityAnalysis::PostProcessedValues(TPZVec<REAL> &x, TPZVec<std::string> &variables, TPZFMatrix<STATE> &values)
{
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCurrentConfig.fCMesh.MaterialVec()[1]);
    pMatWithMem2->SetUpdateMem(true);
    
    // identify the post process indices and number of results for each index
    int varsize = variables.size();
    TPZManVector<int,10> varindices(varsize,-1),firstindex(varsize+1,-1),numvar(varsize,-1);
    int numvalid = 0;
    firstindex[0] = 0;
    for (int iv=0; iv<varsize; iv++) {
        varindices[numvalid] = pMatWithMem2->VariableIndex(variables[iv]);
        if (varindices[numvalid] != -1) {
            numvar[numvalid] = pMatWithMem2->NSolutionVariables(varindices[numvalid]);
            firstindex[numvalid+1] = firstindex[numvalid]+numvar[numvalid];
            numvalid++;
        }
    }
    varindices.Resize(numvalid);
    firstindex.Resize(numvalid+1);
    numvar.Resize(numvalid);
    values.Resize(0, firstindex[numvalid]);
    
    // compute the number of post processing variables and their first index
    std::list<TConfig>::iterator listit;
    // create a new integration point to train the deformation history
    int pointid = pMatWithMem2->PushMemItem();
    for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
        TPZGeoMesh *gmesh1 = &(listit->fGMesh);
        TPZManVector<REAL,3> qsi(2,0.);
        long elementid = 0;
        TPZMaterialData data1;
        data1.x = x;
        data1.intLocPtIndex = 0;
        // find the element which contains the coordinate
        
        TPZGeoEl *gel1 = gmesh1->FindElement(data1.x, qsi, elementid,2);
        if (!gel1) {
            DebugStop();
        }
        TPZCompEl *cel1 = gel1->Reference();
        if (!cel1) {
            DebugStop();
        }
        // compute the solutions at the point
        TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *>(cel1);
        if (!intel1) {
            DebugStop();
        }

        // expand the result matrix
        int valuesfirstindex = values.Rows();
        int nsol = listit->fSolution.Cols();
        if (nsol != 1) {
            DebugStop();
        }
        // we don't initialize the values because they will be set in the procedure
        values.Resize(valuesfirstindex+nsol, values.Cols());
        // copy the solution to a local variable
        int solsize = listit->fSolution.Rows();
        TPZFMatrix<STATE> locsol(solsize,1);
        for (int isol = 0; isol < nsol; isol++) 
        {
            // load one solution at a time
            for (int ic=0; ic<solsize; ic++) {
                locsol(ic,0) = listit->fSolution(ic,isol);
            }
            listit->LoadSolution();
            intel1->InitMaterialData(data1);
            data1.fNeedsSol = true;
            intel1->ComputeRequiredData(data1, qsi);
            data1.intGlobPtIndex = pointid;
            // initialize the material data object
            TPZFMatrix<STATE> EF(0,0);
            data1.phi.Resize(0, 1);
            data1.dphix.Resize(data1.dphix.Rows(), 0);
            pMatWithMem2->Contribute(data1, 1., EF);
            // call the post processing method
            for (int iv=0; iv<numvalid; iv++) 
            {
                TPZManVector<STATE,10> solution(numvar[iv],0.);
                pMatWithMem2->Solution(data1, varindices[iv], solution);
                // put the data in the output matrix
                for (int is=0; is<numvar[iv]; is++) {
                    values(valuesfirstindex+isol,firstindex[iv]+is) = solution[is];
                }
            }
        }
    }
    pMatWithMem2->FreeMemItem(pointid);
    pMatWithMem2->SetUpdateMem(false);
    
}



#include "pzgengrid.h"
#include "TPZVTKGeoMesh.h"







#include "pzelasmat.h"


void TPZSlopeStabilityAnalysis::Print(std::ostream &out)
{
    out << "Number of configurations " << fSequence.size() << std::endl;
    out << "fCurrentConfig \n";
    fCurrentConfig.Print(out);
    out << "Other configurations \n";
    std::list<TPZSlopeStabilityAnalysis::TConfig>::iterator it;
    for (it= fSequence.begin(); it != fSequence.end(); it++) {
        (*it).Print(out);
    }
    if (fLinearMatrix) {
        fLinearMatrix->Print("Linear Matrix",out);
    }
    else
    {
        out << "No linear matrix ";
    }
}

void TPZSlopeStabilityAnalysis::SaveConfig(std::stringstream &strout) {
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

void TPZSlopeStabilityAnalysis::TConfig::Print(ostream &out)
{


    out << "fSDPV ";
    fSDPV.Print(out);
	out << "fMCPV ";
	fMCPV.Print(out);
//    out << "fMatEla ";
//    fMatEla.Print(out); //AQUIOMAR

    out << "fGMesh ";
    fGMesh.Print(out);
    out << "fCMesh ";
    fCMesh.Print(out);
    fSolution.Print("fSolution = ",out);
    out << "fPlasticDeform " << fPlasticDeformSqJ2 << endl;
    fPostprocess.Print("fPostProcess", out);
}


#include "TPZVTKGeoMesh.h"
#include "TPZGeoLinear.h"
void TPZSlopeStabilityAnalysis::TConfig::CreateGeometricMeshSlope(int ref)
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
	    for ( int d = 0; d<ref; d++ ) {
        int nel = fGMesh.NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = fGMesh.ElementVec() [iel];
                gel->Divide ( subels );
        }
    }
	

    


}
/// Initialize the Sandler DiMaggio object and create the computational mesh
void TPZSlopeStabilityAnalysis::TConfig::CreateComputationalMeshSlope(int porder)
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

    //TPZMatElastoPlasticSest2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new //TPZMatElastoPlasticSest2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > ( 1,PlaneStrain );
	
	TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > ( 1,PlaneStrain );
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
