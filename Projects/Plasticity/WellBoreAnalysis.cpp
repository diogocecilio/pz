//
//  WellBoreAnalysis.cpp
//  PZ
//
//  Created by phil on 1/18/13.
//  Copyright (c) 2013 LabMeC. All rights reserved.
//
#include "WellBoreAnalysis.h"
#include "GeoMeshClass.h"
#include "pzelastoplasticanalysis.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplasticSest2D.h"

#include "pzelasticSest2D.h"

#include "pzvisualmatrix.h"

#include "BrazilianTestGeoMesh.h"
#include "tpzchangeel.h"
//#include "poroelastoplastic.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzplasticdiagnostic.h"
#include <iostream>
#include "pzbfilestream.h"
#include "TPBrBiotForce.h"

#include "pzl2projection.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZProjectEllipse.h"
#include "pzlog.h"

#include "arglib.h"
#include "run_stats_table.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.plasticity.wellboreanalysis"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse(Logger::getLogger("LogEllipse"));
#endif

int TPZWellBoreAnalysis::TConfig::gNumThreads = 0;

/// this method will configure the forcing function and boundary condition of the computational mesh
void TPZWellBoreAnalysis::TConfig::ConfigureBoundaryConditions()
{

    
    TPZMaterial *prev = 0;
    TPZBndCond *prevbnd = 0;
    
    TPZMaterial *mat = fCMesh.FindMaterial(1);
    
    prev = fCMesh.FindMaterial(-2);
    TPZFNMatrix<3,REAL> f2(3,1,0.);
    TPZFNMatrix<9,REAL> k2(3,3,0.);
    k2(0,0)=fConfinementTotal.XX()+fBiotCoef*fReservoirPressure;
    k2(1,1)=fConfinementTotal.YY()+fBiotCoef*fReservoirPressure;
    k2(2,2)=fConfinementTotal.ZZ()+fBiotCoef*fReservoirPressure;
    prev = fCMesh.FindMaterial(EInner);
    prevbnd = dynamic_cast<TPZBndCond *>(prev);
    if(!prevbnd || prevbnd->Material() != mat)
    {
//        DebugStop();
        prevbnd = mat->CreateBC(mat,EInner,4,k2,f2);
        fCMesh.InsertMaterialObject(prevbnd);
    }
    else
    {
        prevbnd->SetType(4);
        prevbnd->Val1() = k2;
        prevbnd->Val2() = f2;
    }

    
    // type 6 constraints in x and y
    // type 5 only normal constraint
    TPZFNMatrix<9> k6(3,3,0.),f6(3,1,fReservoirPressure);
    for (int i=0; i<3; i++) {
        k6(i,i) = 1.e5;
    }
    prev = fCMesh.FindMaterial(-6);
    prevbnd = dynamic_cast<TPZBndCond *>(prev);
    if (!prevbnd || prevbnd->Material() != mat) {
        TPZMaterial *bc6 = mat->CreateBC(mat, -6, 6, k6, f6);
        fCMesh.InsertMaterialObject(bc6);
//        DebugStop();
    }
    else
    {
        prevbnd->Val1() = k6;
        prevbnd->Val2() = f6;
        prevbnd->SetType(6);
    }

    

   
    TPZFNMatrix<9,REAL> k3(3,3,0.);
    TPZFNMatrix<3,REAL> f3(3,1,0.);
    k3(0,0)=fConfinementTotal.XX()+fBiotCoef*fReservoirPressure;
    k3(1,1)=fConfinementTotal.YY()+fBiotCoef*fReservoirPressure;
    k3(2,2)=fConfinementTotal.ZZ()+fBiotCoef*fReservoirPressure;
    prev = fCMesh.FindMaterial(EOuter);
    prevbnd = dynamic_cast<TPZBndCond *>(prev);
    if (!prevbnd || prevbnd->Material() != mat) {
//        DebugStop();
        TPZMaterial * bc3 = mat->CreateBC(mat,EOuter,4,k3,f3);
//        k3(0,0)=TPZMaterial::gBigNumber;
//        k3(1,1)=TPZMaterial::gBigNumber;
//        k3(2,2)=TPZMaterial::gBigNumber;
//        f3.Zero();
//        TPZMaterial * bc3 = mat->CreateBC(mat,EOuter,0,k3,f3);
        fCMesh.InsertMaterialObject(bc3);
    }
    else
    {
        prevbnd->SetType(4);
        prevbnd->Val1() = k3;
        prevbnd->Val2() = f3;
    }
 
    
    
    TPZFMatrix<REAL> k4(2,2,0.);
    TPZFMatrix<REAL> f4(2,1,0.);
    f4(1,0)=1.;
   // k3(1,1)=1.e12;
    prev = fCMesh.FindMaterial(EBottom);
    prevbnd = dynamic_cast<TPZBndCond *>(prev);
    if (!prevbnd || prevbnd->Material() != mat) {
//        DebugStop();
        TPZMaterial * bc4 = mat->CreateBC(mat,EBottom,3,k4,f4);
        fCMesh.InsertMaterialObject(bc4);
    }
    else
    {
        prevbnd->SetType(3);
        prevbnd->Val1() = k4;
        prevbnd->Val2() = f4;
    }
    
    
    TPZFMatrix<REAL> k5(2,2,0.);
    TPZFMatrix<REAL> f5(2,1,0.);
     f5(0,0)=1.;
    //k4(0,0)=1.e12;
    
    prev = fCMesh.FindMaterial(ELeft);
    prevbnd = dynamic_cast<TPZBndCond *>(prev);
    if (!prevbnd || prevbnd->Material() != mat)
    {
//        DebugStop();
        TPZMaterial * bc5 = mat->CreateBC(mat,ELeft,3,k5,f5);
        fCMesh.InsertMaterialObject(bc5);
    }
    else
    {
        prevbnd->SetType(3);
        prevbnd->Val1() = k5;
        prevbnd->Val2() = f5;
    }
    

#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCMesh.Print(sout);
        LOGPZ_DEBUG(logger , sout.str())
    }
#endif
}


TPZWellBoreAnalysis::TPZWellBoreAnalysis() : fCurrentConfig(), fSequence(), fPostProcessNumber(0), fLinearMatrix()
{
}

TPZWellBoreAnalysis::TPZWellBoreAnalysis(const TPZWellBoreAnalysis &copy) : fCurrentConfig(copy.fCurrentConfig), fSequence(copy.fSequence), fPostProcessNumber(copy.fPostProcessNumber), fLinearMatrix(copy.fLinearMatrix)
{

}

TPZWellBoreAnalysis &TPZWellBoreAnalysis::operator=(const TPZWellBoreAnalysis &copy)
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

TPZWellBoreAnalysis::~TPZWellBoreAnalysis()
{
    
}

/// write the object on the stream
void TPZWellBoreAnalysis::Write(TPZStream &output)
{
    fCurrentConfig.Write(output);
    int seqsize = fSequence.size();
    output.Write(&seqsize);
    std::list<TPZWellBoreAnalysis::TConfig>::iterator it;
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
void TPZWellBoreAnalysis::Read(TPZStream &input)
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
        TPZWellBoreAnalysis::TConfig config;
        config.Read(input);
        fSequence.push_back(config);
    }
    input.Read(&fPostProcessNumber);
    
}

void TPZWellBoreAnalysis::StandardConfiguration(TPZWellBoreAnalysis &obj)
{
    DebugStop(); // Am i beeing used? If so, please treat elastic material TPZElasticityMaterialSest2D
    
    TPZGeoMesh *gmesh = &obj.fCurrentConfig.fGMesh;
    //GeoMeshClass::WellBore2d(&obj.fCurrentConfig.fGMesh);
    obj.fCurrentConfig.fInnerRadius = 4.25*0.0254;//0.1;
    obj.fCurrentConfig.fOuterRadius = 3.;//1.;
    ofstream arg("wellgeomeshlog.txt");
    obj.fCurrentConfig.fGMesh.Print(arg);
    
    obj.fCurrentConfig.fCMesh.SetReference(gmesh);
    int defaultporder = 2;
    obj.fCurrentConfig.fCMesh.SetDefaultOrder(defaultporder);
    TPZCompMesh *compmesh1 = &obj.fCurrentConfig.fCMesh;
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);

    obj.fCurrentConfig.fWellConfig = EVerticalWell;
    
    TPZManVector<STATE,3> confinementEffective(3,0.), confinementTotal(3,0.);
    REAL SH,Sh,SV;
    Sh=-62.1;
    SH=-45.9;
    SV=-48.2;
    
    confinementEffective[0] = Sh;
    confinementEffective[1] = SH;
    confinementEffective[2] = SV;
    REAL effectiveWellPressure = 19.5; // 19.5 ou 23.4 ou 28.9
    STATE biotcoef = 0.659;
    obj.SetBiotCoefficient(biotcoef);
    
    STATE WellPressure = effectiveWellPressure/(1.-biotcoef);
    for (int i=0; i<3; i++) {
        confinementTotal[i] = confinementEffective[i]-biotcoef*WellPressure;
    }
    
    obj.SetConfinementTotalStresses(confinementTotal, WellPressure);



        int materialid=1;
        bool planestrain=true;
        //obj.fCurrentConfig.fSDPV.fYC.PreSMat(obj.fCurrentConfig.fSDPV.fYC);
        TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> &SD =obj.fCurrentConfig.fSDPV;
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        REAL G=elast/(2.*(1.+poisson));
        REAL K=elast/(3.*(1.-2*poisson));
        REAL phi=0,psi=1.,N=0.;
        SD.fYC.SetUp( A,  B, C,  D, K, G, W, R, phi, N, psi);
        SD.fER.SetUp(elast,poisson);
        TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> > *PlasticSD = new TPZMatElastoPlasticSest2D< TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> >(materialid,planestrain);


    PlasticSD->SetBiot(biotcoef);
        TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> &MC =obj.fCurrentConfig.fMCPV;
        REAL cohesion = A - C;
        REAL angle = B * C;

        TPZElasticResponse ER;
        ER.SetUp(elast,poisson);
        MC.fYC.SetUp(angle, angle, cohesion, ER);
        MC.fER.SetUp(elast,poisson);

        TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> > *PlasticMC = new TPZMatElastoPlasticSest2D< TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> >(materialid,planestrain);
    PlasticMC->SetBiot(biotcoef);



    
    

    TPZTensor<REAL> initstress,finalstress;
    REAL hydro = obj.fCurrentConfig.fConfinementTotal.I1()+3.*obj.fCurrentConfig.fBiotCoef*obj.fCurrentConfig.fReservoirPressure;
//    hydro -= SD.fYC.fA*SD.fYC.fR;
    hydro -= A*R;
    hydro /= 3.;
    finalstress.XX() = hydro;
    finalstress.YY() = hydro;
    finalstress.ZZ() = hydro;
    
    PrepareInitialMat(SD, initstress, finalstress, 10);
    initstress = finalstress;
    finalstress = obj.fCurrentConfig.fConfinementTotal;
    finalstress.XX() += obj.fCurrentConfig.fBiotCoef*obj.fCurrentConfig.fReservoirPressure;
    finalstress.YY() += obj.fCurrentConfig.fBiotCoef*obj.fCurrentConfig.fReservoirPressure;
    finalstress.ZZ() += obj.fCurrentConfig.fBiotCoef*obj.fCurrentConfig.fReservoirPressure;
    PrepareInitialMat(SD, initstress, finalstress, 10);
    SD.ResetPlasticMem();

    if (obj.fCurrentConfig.fModel == EMohrCoulomb) {
        TPZMaterial *plastic(PlasticMC);

        PlasticMC->SetPlasticity(MC);
        compmesh1->InsertMaterialObject(plastic);

        obj.fCurrentConfig.ConfigureBoundaryConditions();
        compmesh1->AutoBuild();
    }
    else
    {
    
        TPZMaterial *plastic(PlasticSD);

        PlasticSD->SetPlasticity(SD);
        compmesh1->InsertMaterialObject(plastic);

        obj.fCurrentConfig.ConfigureBoundaryConditions();
        compmesh1->AutoBuild();
    }
}

void TPZWellBoreAnalysis::CheckDeformation(std::string filename)
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



TPZWellBoreAnalysis::TConfig::TConfig() : fInnerRadius(0.), fOuterRadius(0.), fNx(2,0),fDelx(0.),fGreater(),fSmaller(),fConfinementTotal(),  fWellborePressure(0.),
    fGMesh(), fCMesh(), fSolution(), fZDeformation(0.), fPlasticDeformSqJ2(), fHistoryLog(), fModel(ESandler), fWellConfig(ENoConfig), fFluidModel(ENonPenetrating), fBiotCoef(-1.)
  , fSDPV(), fMCPV(), fMatEla()


{
    fCMesh.SetReference(&fGMesh);
}

TPZWellBoreAnalysis::TConfig::TConfig(const TConfig &conf) : fInnerRadius(conf.fInnerRadius), fOuterRadius(conf.fOuterRadius),fNx(conf.fNx),fDelx(conf.fDelx),
    fGreater(conf.fGreater),fSmaller(conf.fSmaller),
        fConfinementTotal(conf.fConfinementTotal), fWellborePressure(conf.fWellborePressure),
        fGMesh(conf.fGMesh), fCMesh(conf.fCMesh), fSolution(conf.fSolution), fZDeformation(conf.fZDeformation),
    fModel(conf.fModel), fPlasticDeformSqJ2(conf.fPlasticDeformSqJ2), fHistoryLog(conf.fHistoryLog), fFluidModel(conf.fFluidModel), fWellConfig(conf.fWellConfig), fBiotCoef(conf.fBiotCoef)
  , fSDPV(conf.fSDPV), fMCPV(conf.fMCPV), fMatEla(conf.fMatEla)
{
    fSDPV = conf.fSDPV;
    fGMesh.ResetReference();
    fCMesh.SetReference(&fGMesh);
    fCMesh.LoadReferences();
//    std::cout << "Copy constructor from " << (void *) &conf << " to " << (void *) this << std::endl;

}

TPZWellBoreAnalysis::TConfig::~TConfig()
{
//    std::cout << "Deleting config " << (void *) this << std::endl;
    fPostprocess.SetCompMesh(0);
    fCMesh.CleanUp();
    fCMesh.SetReference(0);
}

TPZWellBoreAnalysis::TConfig &TPZWellBoreAnalysis::TConfig::operator=(const TPZWellBoreAnalysis::TConfig &copy)
{
    if (this == &copy) {
        return *this;
    }
    fInnerRadius = copy.fInnerRadius;
    fOuterRadius = copy.fOuterRadius;
    fConfinementTotal = copy.fConfinementTotal;
    fSDPV = copy.fSDPV;
    fMCPV = copy.fMCPV;
    fMatEla = copy.fMatEla;

    fDelx = copy.fDelx;
    fNx = copy.fNx;
    fSmaller = copy.fSmaller;
    fGreater = copy.fGreater;
    fWellborePressure = copy.fWellborePressure;
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
    fFluidModel = copy.fFluidModel;
    fBiotCoef = copy.fBiotCoef;
    fWellConfig = copy.fWellConfig;
    
    std::cout << "Copying " << (void *) &copy << " to " << (void *) this << std::endl;
    return *this;
}

/// Write the data to the output stream
void TPZWellBoreAnalysis::TConfig::Write(TPZStream &out)
{
    out.Write(&fInnerRadius);
    out.Write(&fOuterRadius);
    int wellconf = fWellConfig;
    out.Write(&wellconf);
    out.Write(&fDelx);
    TPZSaveable::WriteObjects(out, fNx);
    TPZSaveable::WriteObjects(out, fGreater);
    TPZSaveable::WriteObjects(out, fSmaller);
    fConfinementTotal.Write(out);
    fSDPV.Write(out);
    fMCPV.Write(out);
    fMatEla.Write(out,0); //AQUIPHIL

		int IntEPlasticModel = fModel;
		out.Write(&IntEPlasticModel);

        int IntEFluidModel = fFluidModel;
        out.Write(&IntEFluidModel);

    out.Write(&fBiotCoef);
	
    out.Write(&fWellborePressure);
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
void TPZWellBoreAnalysis::TConfig::Read(TPZStream &input)
{
    input.Read(&fInnerRadius);
    input.Read(&fOuterRadius);
    int wellconf;
    input.Read(&wellconf);
    fWellConfig = (EWellConfiguration) wellconf;
    input.Read(&fDelx);
    TPZSaveable::ReadObjects(input, fNx);
    TPZSaveable::ReadObjects(input, fGreater);
    TPZSaveable::ReadObjects(input, fSmaller);
    
    fConfinementTotal.Read(input);
    fSDPV.Read(input);
    fMCPV.Read(input);
    fMatEla.Read(input,0); //AQUIPHIL

    int IntEPlasticModel;
    input.Read(&IntEPlasticModel);
    fModel = (EPlasticModel) IntEPlasticModel;

    int IntEFluidModel;
    input.Read(&IntEFluidModel);
    fFluidModel = (EFluidModel) IntEFluidModel;

    input.Read(&fBiotCoef);
	
    input.Read(&fWellborePressure);
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


void TPZWellBoreAnalysis::ExecuteInitialSimulation(int nsteps, int numnewton)
{
    std::stringstream strout;
    strout << "Initial Simulation";
    fCurrentConfig.fSolution.Redim(fCurrentConfig.fCMesh.NEquations(), 1);
    SaveConfig(strout);
    fSequence.push_back(fCurrentConfig);

    for (int istep = 1; istep <= nsteps; istep++) {
        REAL factor = istep*1./nsteps;
        fCurrentConfig.SetWellPressure(factor);
        ExecuteSimulation(istep,factor);
    }
}

void TPZWellBoreAnalysis::ExecuteInitialSimulation2(int nsteps, int numnewton)
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

/// Computes the given pressure in the specified steps
void TPZWellBoreAnalysis::EvolveBothPressures(int nsteps, STATE TargetWellborePressure, STATE TargetReservoirPressure)
{
    std::stringstream strout;
    strout << "Evolve";
    SaveConfig(strout);
//    fCurrentConfig.fHistoryLog = strout.str();
//    fSequence.push_back(fCurrentConfig);

    STATE InitWellPressure = fCurrentConfig.fWellborePressure;
    STATE InitReservoirPressure = fCurrentConfig.fReservoirPressure;
    for (int istep = 1; istep <= nsteps; istep++) {
        STATE WellPressure = InitWellPressure+istep*(TargetWellborePressure-InitWellPressure)/nsteps;
        fCurrentConfig.fWellborePressure = WellPressure;
        
        STATE ReservoirPressure = InitReservoirPressure+istep*(TargetReservoirPressure-InitReservoirPressure)/nsteps;
        fCurrentConfig.fReservoirPressure = ReservoirPressure;
        
        // adjust the boundary conditions and forcing function
        fCurrentConfig.SetWellPressure();
#ifdef DEBUG
        std::cout << "&&&&&&&&&&&&&&&&&&&&&EVOLVEBOTHPRESSURE&&&&&&&&&&&&&&&&&&&&&&&&&&& Running with WellborePressure = " << fCurrentConfig.fWellborePressure << "  and ReservoirPressure = " << fCurrentConfig.fReservoirPressure;
#endif
        ExecuteSimulation(istep);
    }
}

void TPZWellBoreAnalysis::ExecuteSimulation2(int substeps, REAL factor)
{
    
    TConfig &LocalConfig = fCurrentConfig;
    
    TPZElastoPlasticAnalysis analysis(&LocalConfig.fCMesh,std::cout);
    
	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(TPZWellBoreAnalysis::TConfig::gNumThreads);


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
    if (LocalConfig.fWellConfig == EVerticalWell) {
        ncycles = 3;
    }
    
    std::stringstream strout;
    

        
	well_init.start();
	analysis.IterativeProcess(strout, tol, NumIter,linesearch,checkconv);
	well_init.stop();
        
       
    
        
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


void TPZWellBoreAnalysis::ExecuteSimulation(int substeps, REAL factor)
{
    
    TConfig &LocalConfig = fCurrentConfig;
    STATE zdeformation = 0.;
    STATE targetZStress = LocalConfig.fConfinementTotal.ZZ();
    STATE InitialZStress = 0.;
    STATE InitialZDeformation = 0.;
    STATE SecondZStress = 0.;
    STATE SecondZDeformation = 0.;
    STATE elasticity = 0.;

    if (LocalConfig.fModel == ESandler) {
        elasticity = LocalConfig.fSDPV.fER.E();
    }
    else if(LocalConfig.fModel == EMohrCoulomb)
    {
        elasticity = LocalConfig.fMCPV.fER.E();
    }
    else if(LocalConfig.fModel == EElastic)
    {
#ifdef PlasticPQP
        elasticity = LocalConfig.fMCPV.fER.E();
#else
        elasticity = LocalConfig.fMatEla.GetEyoung();
#endif
    }
    else
    {
        DebugStop();
    }

    LocalConfig.SetZDeformation(zdeformation);
    
    TPZElastoPlasticAnalysis analysis(&LocalConfig.fCMesh,std::cout);
    
	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(TPZWellBoreAnalysis::TConfig::gNumThreads);
    
    analysis.AddNoPenetration(ELeft, 0);
    analysis.AddNoPenetration(EBottom, 1);
    
    if (fCurrentConfig.fWellConfig == EVerticalWell) {
        analysis.AddNoPenetration(EOuter, 0);
        analysis.AddNoPenetration(EOuter, 1);
    }
    else if (fCurrentConfig.fWellConfig == EHorizontalWellalongH || fCurrentConfig.fWellConfig == EHorizontalWellalongh)
    {
        analysis.AddNoPenetration(EOuter, 1);
    }
    else
    {
        DebugStop();
    }

    analysis.IdentifyEquationsToZero();
    
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
    if (LocalConfig.fWellConfig == EVerticalWell) {
        ncycles = 3;
    }
    
    std::stringstream strout;
    
    for (int cycle = 0; cycle < ncycles; cycle++)
    {
    
        
        well_init.start();
        analysis.IterativeProcess(strout, tol, NumIter,linesearch,checkconv);
        well_init.stop();
		
        
        if (conv==false) {
// 			std::cout << "\n here 4" <<endl;
//             // the data of the config object indicates whether there is a vertical strain to model
//             ComputeLinearMatrix(activeEquations);
//             analysis.Solver().SetMatrix(fLinearMatrix);
//             
//             well_init.start();
//             analysis.IterativeProcess(strout, fLinearMatrix, tol, NumIter, linesearch);
//             well_init.stop();
// 			std::cout << "\n here 4b" <<endl;
        }
        
        if (LocalConfig.fWellConfig == EVerticalWell && cycle == 0)
        {
			std::cout << "\n here 5" <<endl;
            InitialZStress = LocalConfig.AverageVerticalStress();
            InitialZDeformation = zdeformation;
            zdeformation = (targetZStress-InitialZStress)/elasticity;
            if (fabs(zdeformation) < 1.e-4) {
                zdeformation = (zdeformation < 0.) ? -1.e-4 : 1.e-4;
            }
            LocalConfig.SetZDeformation(zdeformation);
        }
        else if(LocalConfig.fWellConfig == EVerticalWell && cycle == 1)
        {
			std::cout << "\n here 6" <<endl;
            SecondZStress = LocalConfig.AverageVerticalStress();
            SecondZDeformation = zdeformation;
            zdeformation = InitialZDeformation + (SecondZDeformation-InitialZDeformation)*(targetZStress-InitialZStress)/(SecondZStress-InitialZStress);
            LocalConfig.SetZDeformation(zdeformation);
        }
    }
        
    TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
    for (int ieq=0; ieq<neq; ieq++) {
        LocalConfig.fSolution(ieq,0) = sol(ieq,0);
    }
    
    analysis.AcceptSolution();
    LocalConfig.SetZDeformation(0.);
    LocalConfig.fZDeformation = zdeformation;
    
    // dont forget the vertical deformation
    LocalConfig.ComputeElementDeformation();
    

    int BCId=EInner;
    TPZMaterial * mat = fCurrentConfig.fCMesh.FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    if (!pBC) {
#ifdef DEBUG
        DebugStop();
#endif
    }
    STATE SX=pBC->Val1().GetVal(0,0);//effective
    STATE SY=pBC->Val1().GetVal(1,1);//effective
    SX=SX-fCurrentConfig.fBiotCoef*fCurrentConfig.fReservoirPressure;//total
    SY=SY-fCurrentConfig.fBiotCoef*fCurrentConfig.fReservoirPressure;//total

    
    strout << " Substep " << substeps << " - (total stresses) "<< " SX  = " << SX << " SY  = " << SY;

    SaveConfig(strout);
//    fCurrentConfig.fHistoryLog = strout.str();

    
//    fSequence.push_back(LocalConfig);
}



/// Load the solution stored in TConfig into the CompMesh and set the ZDeformation of the material
void TPZWellBoreAnalysis::TConfig::LoadSolution()
{
    TPZMaterial *mat = fCMesh.FindMaterial(1);
    typedef TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> mattype1;
    typedef TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem> mattype2;
    typedef TPZElasticityMaterialSest2D mattype3;
    
    mattype1 *matposs1 = dynamic_cast<mattype1 *>(mat);
    mattype2 *matposs2 = dynamic_cast<mattype2 *>(mat);
    mattype3 *matposs3 = dynamic_cast<mattype3 *>(mat);
    
    if (matposs1) {
        matposs1->SetZDeformation(fZDeformation);
    }
    if (matposs2) {
        matposs2->SetZDeformation(fZDeformation);
    }
    if (matposs3) {
        matposs3->SetZDeformation(fZDeformation);
    }
    if (!matposs1 && !matposs2 && !matposs3) {
        DebugStop();
    }
#ifdef DEBUG
    if (fSolution.Cols() > 1) {
        DebugStop();
    }
#endif
    fCMesh.LoadSolution(fSolution);
}



/// Set the Z deformation (for adapting the compaction)
void TPZWellBoreAnalysis::TConfig::SetZDeformation(STATE epsZ)
{
    TPZMaterial *mat = fCMesh.FindMaterial(1);
    typedef TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> mattype1;
    typedef TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem> mattype2;
    typedef TPZElasticityMaterialSest2D mattype3;

    mattype1 *matposs1 = dynamic_cast<mattype1 *>(mat);
    mattype2 *matposs2 = dynamic_cast<mattype2 *>(mat);
    mattype3 *matposs3 = dynamic_cast<mattype3 *>(mat);

    if (matposs1) {
        matposs1->SetZDeformation(epsZ);
    }
    if (matposs2) {
        matposs2->SetZDeformation(epsZ);
    }
    if (matposs3) {
        matposs3->SetZDeformation(epsZ);
    }
    if (!matposs1 && !matposs2 && !matposs3) {
        DebugStop();
    }
    fZDeformation = epsZ;
}



/// this method will modify the boundary condition of the computational mesh and the forcing function
void TPZWellBoreAnalysis::TConfig::SetWellPressure(STATE factor)
{
    int BCId=EInner;
    TPZMaterial * mat = fCMesh.FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    if (!pBC) {
        DebugStop();
    }
    mat = fCMesh.FindMaterial(1);
    TPZAutoPointer<TPZFunction<STATE> > force = mat->ForcingFunction();
    if (!force) {
      force = new TPBrBiotForce();
      mat->SetForcingFunction(force);
    }
    
    TPBrBiotForce *biotforce = dynamic_cast<TPBrBiotForce *>(force.operator->());
    if (!biotforce) {
        DebugStop();
    }
    REAL inner = fInnerRadius;
    REAL outer = fOuterRadius;
    if (fBiotCoef < 0.) {
        DebugStop();
    }
    REAL biot = fBiotCoef;
    REAL wellpressure = fWellborePressure;
    REAL reservoirpressure = fReservoirPressure;
    
    TPZTensor<STATE> boundarytensor;
    FromPhysicalDomaintoComputationalDomainStress(fConfinementTotal, boundarytensor);
    boundarytensor.XX() += biot*reservoirpressure;
    boundarytensor.YY() += biot*reservoirpressure;
    
    if (fFluidModel == EPenetrating)
    {
        biotforce->SetConstants(inner, outer, wellpressure, reservoirpressure,factor*biot);
        TPZFNMatrix<9,STATE> mattemp(3,3,0.);
        mattemp(0,0) = factor*(-(1.-biot)*fWellborePressure)+(1.-factor)*boundarytensor.XX();
        mattemp(1,1) = factor*(-(1.-biot)*fWellborePressure)+(1.-factor)*boundarytensor.YY();
        pBC->Val1()=mattemp;
    }
    else{
        biotforce->SetConstants(inner, outer, reservoirpressure, reservoirpressure,biot);
        TPZFNMatrix<9,STATE> mattemp(3,3,0.);
        mattemp(0,0) = factor*(-(wellpressure-biot*reservoirpressure))+(1.-factor)*boundarytensor.XX();
        mattemp(1,1) = factor*(-(wellpressure-biot*reservoirpressure))+(1.-factor)*boundarytensor.YY();
        pBC->Val1()=mattemp;
    }
    
}


/// Apply the deformation of the configuration to the element
void TPZWellBoreAnalysis::TConfig::ApplyDeformation(TPZCompEl *cel)
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
void TPZWellBoreAnalysis::VerifyTangentValidity()
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
void TPZWellBoreAnalysis::TConfig::VerifyPlasticTangent(TPZCompEl *cel)
{
    TPZCompEl *cel2 = cel;
    TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
    if (!intel2) {
        DebugStop();
    }
    TPZCompMesh *cmesh2 = cel2->Mesh();
    
    
    TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *pMatWithMem2 =
    dynamic_cast<TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *>(cmesh2->MaterialVec()[1]);

    
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



void TPZWellBoreAnalysis::TransferSolutionTo(TConfig &config)
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
REAL TPZWellBoreAnalysis::TConfig::ComputeAreaAboveSqJ2(REAL sqj2)
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

/// Compute the area of the domain at which sqJ2 is above a given value
REAL TPZWellBoreAnalysis::TConfig::OpeningAngle(REAL sqj2)
{
    REAL angle = 0.;
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
            int numind = memindices.size();
            REAL sqj2el;
            sqj2el = 0.;
            for (int ind=0; ind<numind; ind++) 
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
                    TPZManVector<REAL,3> x(3,0.);
                    gel->X(point, x);
                    TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ip]);
                    TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                    REAL J2 = plastic.J2();
                    REAL sqj2pt = sqrt(J2);
                    if (sqj2pt > sqj2) {
                        REAL locangle = atan2(x[1],x[0]);
                        if (locangle > angle) {
                            angle = locangle;
                        }
                    }
                }
            }
        }
    }
    return angle;
}

/// Compute the average vertical stress of the configuration
STATE TPZWellBoreAnalysis::TConfig::AverageVerticalStress()
{
    std::set<int> matids;
    matids.insert(1);
    TPZVec<STATE> verticalForce = fCMesh.Integrate("TotStressZZ", matids);
    REAL MeshArea = fCMesh.Reference()->Area();
    return verticalForce[0]/MeshArea;
}



/// Compute the removed area of the domain
REAL TPZWellBoreAnalysis::TConfig::RemovedArea()
{
    return -ComputeTotalArea()+M_PI*(fOuterRadius*fOuterRadius-fInnerRadius*fInnerRadius)/4.;
}

/// Compute the area of the domain
REAL TPZWellBoreAnalysis::TConfig::ComputeTotalArea()
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
void TPZWellBoreAnalysis::TConfig::ComputeElementDeformation()
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
void TPZWellBoreAnalysis::TConfig::VerifyGlobalEquilibrium(std::ostream &out)
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
void TPZWellBoreAnalysis::TConfig::VerifyGlobalEquilibrium2(std::ostream &out)
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
void TPZWellBoreAnalysis::TConfig::FilterRhs(TPZFMatrix<STATE> &rhs)
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
void TPZWellBoreAnalysis::TConfig::ComputeRhsExceptMatid(int matid, TPZFMatrix<STATE> &rhs)
{
    std::set<int> allmaterials;
    std::map<int,TPZMaterial *>::iterator itmat;
    std::map<int,TPZMaterial *> &matmap = fCMesh.MaterialVec();
    for (itmat = matmap.begin(); itmat != matmap.end(); itmat++) {
        allmaterials.insert(itmat->first);
    }
    allmaterials.erase(matid);
    TPZFStructMatrix str(&fCMesh);
    str.SetNumThreads(TPZWellBoreAnalysis::TConfig::gNumThreads);
    str.SetMaterialIds(allmaterials);
    str.Assemble(rhs, 0);
    FilterRhs(rhs);

    
}

/// Compute the contribution of only matid
// this method is cumulative (sums to the rhs)
void TPZWellBoreAnalysis::TConfig::ComputeRhsForMatid(int matid, TPZFMatrix<STATE> &rhs)
{
    std::set<int> allmaterials;
    allmaterials.insert(matid);
    TPZFStructMatrix str(&fCMesh);
    str.SetMaterialIds(allmaterials);
    str.SetNumThreads(TPZWellBoreAnalysis::TConfig::gNumThreads);
    str.Assemble(rhs, 0);
    FilterRhs(rhs);
}

/// Compute the resultant x and y force
void TPZWellBoreAnalysis::TConfig::ComputeXYForce(TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force)
{
    force.Resize(2, 0.);
    force.Fill(0.);
    long nel = rhs.Rows();
    for (long i=0; i<nel; i++) {
        force[i%2] += rhs(i,0);
    }
}


void TPZWellBoreAnalysis::TConfig::DeleteElementsAbove(REAL sqj2)
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

// Diminish the spring which stabilizes the well by the given factor
void TPZWellBoreAnalysis::TConfig::RelaxWellSpring(REAL factor)
{
    TPZMaterial *mat = fCMesh.FindMaterial(-6);
    if (!mat) {
        DebugStop();
    }
    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
    if (!bnd) {
        DebugStop();
    }
    bnd->Val1() *= factor;
    bnd->Val1().Print("Value of the spring",std::cout);
    bnd->Val2().Print("Fluid pressure ",std::cout);
}


/// Change the polynomial order of element using the plastic deformation as threshold
void TPZWellBoreAnalysis::TConfig::PRefineElementsAbove(REAL sqj2, int porder, std::set<long> &elindices)
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
void TPZWellBoreAnalysis::TConfig::DivideElementsAbove(REAL sqj2, std::set<long> &elindices)
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

TPZGeoMesh * TPZWellBoreAnalysis::TConfig::GetGeoMesh() {
  return &this->fGMesh;
}

/// Return the mesh used for computations (multiphysics mesh or fCMesh)
TPZCompMesh *TPZWellBoreAnalysis::TConfig::CompMeshUsed()
{
    TPZCompMesh *cmesh = &fCMesh;
    return cmesh;
}


/// Reset the plastic memory of the integration points of these elements
void TPZWellBoreAnalysis::ApplyHistory(std::set<long> &elindices)
{
    std::list<TConfig>::iterator listit;
    for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
        listit->LoadSolution();
    }
    
#ifdef DEBUG
    std::stringstream filename;
    filename << "applyhistory_" << startfrom << ".txt";
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

/// Debugging method, printing the parameters of the forcing function
void TPZWellBoreAnalysis::TConfig::PrintForcingFunction(std::ostream &out)
{
    TPZMaterial *matplast = fCMesh.FindMaterial(1);
    if(!matplast)
    {
        out << __PRETTY_FUNCTION__ << " no material\n";
        return;
    }
    TPZAutoPointer<TPZFunction<STATE> > force = matplast->ForcingFunction();
    TPBrBiotForce *biot = dynamic_cast<TPBrBiotForce *>(force.operator->());
    if (!biot) {
        out << "No biot force configured\n";
    }
    else
    {
        biot->Print(out);
    }
}



void TPZWellBoreAnalysis::TConfig::CreatePostProcessingMesh()
{
    if (fPostprocess.ReferenceCompMesh() != &fCMesh)
    {

        fPostprocess.SetCompMesh(&fCMesh);
        TPZFStructMatrix structmatrix(fPostprocess.Mesh());
        structmatrix.SetNumThreads(0);
        fPostprocess.SetStructuralMatrix(structmatrix);
        
        TPZVec<int> PostProcMatIds(1,1);
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        TPZWellBoreAnalysis::PostProcessVariables(scalNames, vecNames);

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
void TPZWellBoreAnalysis::PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames)
{
  scalNames.Resize(0);
  vecNames.Resize(0);
  
scalNames.Push("PlStrainSqJ2");
scalNames.Push("StrainVol");
scalNames.Push("ElStrainVol");
scalNames.Push("StrainYY");
scalNames.Push("StrainZZ");
scalNames.Push("StrainXY");
scalNames.Push("StrainXZ");
scalNames.Push("DisplacementX");
scalNames.Push("DisplacementY");
scalNames.Push("DisplacementZ");
vecNames.Push("DisplacementTotal");



  

  
}


int passCount = 0;
void TPZWellBoreAnalysis::PostProcess(int resolution)
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

void TPZWellBoreAnalysis::GetJ2Isoline(REAL J2val, std::multimap<REAL,REAL> & polygonalChain)
{
    TPZCompMesh *cmesh = fCurrentConfig.fPostprocess.Mesh();
    REAL tol = 0.01*J2val;
    
    TPZTensor<STATE> boundarytensor;
    fCurrentConfig.FromPhysicalDomaintoComputationalDomainStress(fCurrentConfig.fConfinementTotal, boundarytensor);
    
    long nels = cmesh->NElements();
    for(long el = 0; el < nels; el++)
    {
        TPZCompEl * cel = cmesh->ElementVec()[el];
        if(!cel)
        {
            continue;
        }
        if(cel->Dimension() == 2)
        {
            int nnodes = cel->Reference()->NCornerNodes();
            TPZVec<REAL> nodeSol(nnodes,0.);
            TPZVec< TPZVec<REAL> > qsiNode(nnodes);
            
            TPZVec<STATE> sol(1);
            int var = cel->Material()->VariableIndex("PlStrainSqJ2");
            REAL Tol = 1.e-6;
			//ZeroTolerance(Tol);
            for(int n = 0; n < nnodes; n++)
            {
                qsiNode[n].Resize(2, 0.);
                TPZVec<REAL> xNode(3,0.);
                cel->Reference()->NodePtr(n)->GetCoordinates(xNode);
                cel->Reference()->ComputeXInverse(xNode, qsiNode[n],Tol);
                cel->Solution(qsiNode[n], var, sol);
                nodeSol[n] = sol[0];
            }
            for(int n = 0; n < nnodes; n++)
            {
                int pos0 = (n)%nnodes;
                int pos1 = (n+1)%nnodes;
                REAL minSol = MIN(nodeSol[pos0],nodeSol[pos1]);
                REAL maxSol = MAX(nodeSol[pos0],nodeSol[pos1]);
                
                if(J2val - minSol > tol && maxSol - J2val > tol)
                {//Bisseccao!!!
                    TPZVec<REAL> qsiMin(2,0.), qsiMax(2,0.);
                    if(fabs(minSol - nodeSol[pos0]) < tol)
                    {
                        qsiMin = qsiNode[pos0];
                        qsiMax = qsiNode[pos1];
                    }
                    else
                    {
                        qsiMin = qsiNode[pos1];
                        qsiMax = qsiNode[pos0];
                    }
                    TPZVec<REAL> qsiMiddle(2,0.);
                    qsiMiddle[0] = (qsiMin[0] + qsiMax[0])/2.;
                    qsiMiddle[1] = (qsiMin[1] + qsiMax[1])/2.;
                    
                    cel->Solution(qsiMiddle, var, sol);
                    
                    REAL difSol = fabs(sol[0] - J2val);
                    
                    while(difSol > tol)
                    {
                        if(sol[0] < J2val)
                        {
                            qsiMin = qsiMiddle;
                        }
                        else
                        {
                            qsiMax = qsiMiddle;
                        }
                        
                        qsiMiddle[0] = (qsiMin[0] + qsiMax[0])/2.;
                        qsiMiddle[1] = (qsiMin[1] + qsiMax[1])/2.;
                        
                        cel->Solution(qsiMiddle, var, sol);
                        
                        difSol = fabs(sol[0] - J2val);
                    }
                    TPZManVector<REAL,3> searchedX(3,0.);
                    cel->Reference()->X(qsiMiddle, searchedX);
                    REAL Sx,Sy;
                    Sx=boundarytensor.XX();
                    Sy=boundarytensor.YY();
                    
                    if(fabs(Sx)>fabs(Sy))
                    {
                        polygonalChain.insert( std::make_pair(searchedX[1],searchedX[0]) );
                    }
                    else
                    {
                        polygonalChain.insert( std::make_pair(searchedX[0],searchedX[1]) );
                    }
            
                    
                }
            }
        }
    }
    
    
#ifdef LOG4CXX
    if (loggerEllipse->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\n polygonalChain"<<std::endl;
        for (std::multimap<REAL, REAL>::iterator it = polygonalChain.begin(); it!= polygonalChain.end(); it++) {
            sout << it->first << " " << it->second <<std::endl;
        }
        std::cout << sout.str() << std::endl;
        
        LOGPZ_DEBUG(loggerEllipse, sout.str())
    }
#endif
    


}

/// Modify the geometry of the domain simulating an elliptic breakout
void TPZWellBoreAnalysis::AddEllipticBreakout(REAL MaiorAxis, REAL MinorAxis)
{

    fCurrentConfig.AddEllipticBreakout(MaiorAxis, MinorAxis);
    std::set<long> elindices;
    int nel = fCurrentConfig.fCMesh.NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fCurrentConfig.fCMesh.ElementVec()[el];
        TPZGeoEl *gel = 0;
        if(cel) gel = cel->Reference();
        if(cel && gel && gel->MaterialId() == 1)
        {
            elindices.insert(el);
        }
    }
    ApplyHistory(elindices);
    fCurrentConfig.fCMesh.Solution().Zero();
    fCurrentConfig.fSolution.Resize(0, 0);
    fCurrentConfig.fPostprocess.SetCompMesh(0);
}



/// Divide the element using the plastic deformation as threshold
unsigned int TPZWellBoreAnalysis::DivideElementsAbove(REAL sqj2)
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
void TPZWellBoreAnalysis::PostProcessedValues(TPZVec<REAL> &x, TPZVec<std::string> &variables, TPZFMatrix<STATE> &values)
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

void TPZWellBoreAnalysis::ComputeAandB(REAL sqj2_refine, REAL &a,REAL &b)
{
    std::multimap<REAL, REAL> polygonalChainbase, polygonalChain;
    GetJ2Isoline(sqj2_refine, polygonalChainbase);
    REAL maxy = fCurrentConfig.MaxYfromLastBreakout();
    for (std::multimap<REAL, REAL>::iterator it = polygonalChainbase.begin(); it != polygonalChainbase.end(); it++) {
        if (it->second < maxy) {
            polygonalChain.insert(std::make_pair(it->first,it->second));
        }
    }
    
#ifdef LOG4CXX
    if (loggerEllipse->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\n sqj2_refine " << sqj2_refine <<std::endl;
        sout << "\n maxy " << maxy <<std::endl;
        
        
        sout << "\n polygonalChainbase = "<<std::endl;
        for (std::multimap<REAL, REAL>::iterator it = polygonalChainbase.begin(); it!= polygonalChainbase.end(); it++) {
            
            if (it==polygonalChainbase.begin()) {
                sout << "{{"<<it->first << "," << it->second << "},"<<std::endl;   
            }
            else if(it==polygonalChainbase.end()){
                sout << "{"<<it->first << "," << it->second << "}};"<<std::endl;    
            }
            else{
                sout << "{"<<it->first << "," << it->second << "},"<<std::endl;   
            }
            
        }
        
        sout << "\n polygonalChain ="<<std::endl;
        for (std::multimap<REAL, REAL>::iterator it = polygonalChain.begin(); it!= polygonalChain.end(); it++) {
            if (it==polygonalChain.begin()) {
                sout << "{{"<<it->first << "," << it->second << "},"<<std::endl;
            }
            else if(it==polygonalChain.end()){
                sout << "{"<<it->first << "," << it->second << "}};"<<std::endl;
            }
            else{
                sout << "{"<<it->first << "," << it->second << "},"<<std::endl;
            }
        }
        
        std::cout << sout.str() << std::endl;
        
        LOGPZ_DEBUG(loggerEllipse, sout.str())
    }
#endif
    
    if(polygonalChain.size()!=0)
    {
        TPZProjectEllipse ellips(polygonalChain);
        TPZManVector<REAL,2> center(2),ratios(2),verify(2);
        ellips.StandardFormatForSimpleEllipse(center, ratios);
        verify[0] = ratios[0]/fCurrentConfig.fInnerRadius;
        verify[1] = ratios[1]/fCurrentConfig.fInnerRadius;
        a = ratios[0];
        b = ratios[1];
        
#ifdef LOG4CXX
        if (loggerEllipse->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "\n a " << a <<std::endl;
            sout << "\n b " << b <<std::endl;
            std::cout << sout.str() << std::endl;
            
            LOGPZ_DEBUG(loggerEllipse, sout.str())
        }
#endif
        
    }
    
}


/// Add elliptic breakout
void TPZWellBoreAnalysis::TConfig::AddEllipticBreakout(REAL MaiorAxis, REAL MinorAxis)
{
    if (fGreater.size() && fGreater[fGreater.size()-1] > MaiorAxis) {
        DebugStop();
    }
    if (MaiorAxis < fInnerRadius) {
        DebugStop();
    }
    fGreater.Resize(fGreater.size()+1, MaiorAxis);
    if (fSmaller.size() && fSmaller[fSmaller.size()-1] < MinorAxis) {
        DebugStop();
    }
    if (MinorAxis > fInnerRadius) {
        DebugStop();
    }
    fSmaller.Resize(fSmaller.size()+1, MinorAxis);
    
    fPostprocess.SetCompMesh(0);
    fCMesh.CleanUp();
    fGMesh.CleanUp();
    CreateMesh();

    CreateComputationalMesh(2);
    SetWellPressure();

    
}

#include "pzgengrid.h"
#include "TPZVTKGeoMesh.h"

void TPZWellBoreAnalysis::TConfig::CreateMesh()
{
    if ((fGMesh.NElements() != 0) || (fGMesh.NNodes() != 0))
    {
        DebugStop();
    }

    TPZTensor<STATE> boundarytensor;
    FromPhysicalDomaintoComputationalDomainStress(fConfinementTotal, boundarytensor);
    
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
#ifdef DEBUG
    std::ofstream out("Wellbore.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&fGMesh, out,true);
#endif
    
    if (fGreater.size() == 0) {
        return;
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fGMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    std::ofstream ssout("WellboreBefore.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&fGMesh, ssout,true);
    
    /// number of elements in the radial and circumferential direction = nx
    
    // project the nodes on the elliptic boundaries
   for (int iy=0; iy < nx[1]+1; iy++) {//loop over the circumferential
        int index = iy*(nx[0]+1);
        TPZManVector<REAL,3> co(3);
        fGMesh.NodeVec()[index].GetCoordinates(co);
        REAL xinit = co[0];
        REAL yinit=co[1];
        bool changed = ProjectNode(co);
        if (changed) {
            REAL Sh1,Sh2;
            Sh1=boundarytensor.XX();
            Sh2=boundarytensor.YY();
            if(fabs(Sh1)<fabs(Sh2))
            {
                REAL xadjust = co[0];
                int endindex = index+nx[0];
                TPZManVector<REAL,3> endco(3);
                fGMesh.NodeVec()[endindex].GetCoordinates(endco);
                REAL xend = endco[0];
                for (int ix=0; ix< nx[0]+1; ix++) {//loop over the radial
                    TPZManVector<REAL,3> nodeco(3);
                    fGMesh.NodeVec()[index+ix].GetCoordinates(nodeco);
                    REAL factor = (nodeco[0]-xinit)/(xend-xinit);
                    REAL correctx = xadjust + factor*(xend-xadjust);
                    nodeco[0] = correctx;
                    fGMesh.NodeVec()[index+ix].SetCoord(nodeco);
                }
            }else{
                REAL yadjust = co[1];
                int endindex = index+nx[0];
                TPZManVector<REAL,3> endco(3);
                fGMesh.NodeVec()[endindex].GetCoordinates(endco);
                REAL yend = endco[1];
                for (int iy=0; iy< nx[0]+1; iy++) {//loop over the radial
                    TPZManVector<REAL,3> nodeco(3);
                    fGMesh.NodeVec()[index+iy].GetCoordinates(nodeco);
                    REAL factor = (nodeco[1]-yinit)/(yend-yinit);
                    REAL correcty = yadjust + factor*(yend-yadjust);
                    nodeco[1] = correcty;
                    fGMesh.NodeVec()[index+iy].SetCoord(nodeco);
                }
            }
        }
    }
    ModifyWellElementsToQuadratic();
}

void TPZWellBoreAnalysis::TConfig::ModifyWellElementsToQuadratic()
{
//#ifdef DEBUG
//    TPZManVector<REAL,3> findel(3,0.),qsi(2,0.);
//    findel[0] = 0.108;
//    findel[1] = 0.0148;
//    long elindex = 0;
//    fCMesh.Reference()->FindElement(findel, qsi, elindex, 2);
//    TPZGeoEl *targetel = fCMesh.Reference()->ElementVec()[elindex];
//    TPZCompEl *targetcel = targetel->Reference();
//    long targetindex = targetcel->Index();
//#endif
    //fGMesh.Print();
    int nelem = fGMesh.NElements();
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * bcGel = fGMesh.ElementVec()[el];
        
#ifdef DEBUG
        if(!bcGel)
        {
            DebugStop();
        }
#endif
        
        if(bcGel->MaterialId() != -2)
        {
            continue;
        }
        
#ifdef DEBUG
        if(bcGel->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        if(bcGel->HasSubElement() || bcGel->Father())
        {

            DebugStop();//este metodo nao funciona para elementos refinados e/ou filhos
        }
        int inner1DSide = 2;
        
        bcGel = TPZChangeEl::ChangeToQuadratic(&fGMesh, el);
        
        TPZGeoElSide bcSide(bcGel,inner1DSide);
        int nodeIdlocal = 2;
        int nodeIdnodevec = -1;
        nodeIdnodevec = bcGel->NodeIndex(nodeIdlocal);
#ifdef DEBUG
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            TPZManVector<TPZManVector<REAL,3>,3> allcoord(3);
            std::stringstream sout;
            sout << "Coordinates of the quadratic element before projecting\n";
            for (int in=0; in<3; in++) {
                allcoord[in].Resize(3, 0.);
                bcGel->NodePtr(in)->GetCoordinates(allcoord[in]);
                bcGel->NodePtr(in)->Print(sout);
            }
            for (int ic=0; ic<3; ic++) {
                allcoord[2][ic] -= (allcoord[0][ic]+allcoord[1][ic])*0.5;
            }
            sout << "Relative position of the middle node " << allcoord[2];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#endif
        
        TPZManVector<REAL,3> nodeCoord(3,0.);
        fGMesh.NodeVec()[nodeIdnodevec].GetCoordinates(nodeCoord);
        bool adjusted = ProjectNode(nodeCoord);
        
#ifdef DEBUG
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            TPZManVector<TPZManVector<REAL,3>,3> allcoord(3);
            std::stringstream sout;
            sout << "Coordinates of the quadratic element after projecting\n";
            for (int in=0; in<3; in++) {
                allcoord[in].Resize(3, 0.);
                bcGel->NodePtr(in)->GetCoordinates(allcoord[in]);
                bcGel->NodePtr(in)->Print(sout);
            }
            for (int ic=0; ic<3; ic++) {
                allcoord[2][ic] = nodeCoord[ic]-(allcoord[0][ic]+allcoord[1][ic])*0.5;
            }
            sout << "Relative position of the middle node " << allcoord[2];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#endif
        
        TPZGeoElSide neighSide(bcSide.Neighbour());
        while(neighSide != bcSide)
        {
//#ifdef DEBUG
//            long neighindex = neighSide.Element()->Index();
//            if (neighindex == elindex) {
//                std::cout << "I should stop\n";
//                std::cout << "One d element index = " << el << " node index " << nodeIdnodevec << std::endl;
//                fGMesh.NodeVec()[nodeIdnodevec].GetCoordinates(nodeCoord);
//                ProjectNode(nodeCoord);
//
//            }
//#endif
            TPZGeoEl * quadraticGel = TPZChangeEl::ChangeToQuadratic(&fGMesh, neighSide.Element()->Index());
            neighSide = quadraticGel->Neighbour(neighSide.Side());
        }
        
        if(adjusted)
        {
            fGMesh.NodeVec()[nodeIdnodevec].SetCoord(nodeCoord);
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fGMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// project the node on the boundary
// returns true if the coordinate was changed
bool TPZWellBoreAnalysis::TConfig::ProjectNode(TPZVec<REAL> &co)
{
    bool wasadjusted = false;
    bool considered = false;
    int ellips = 0;
    REAL Sh1,Sh2;
    TPZTensor<STATE> boundarytensor;
    FromPhysicalDomaintoComputationalDomainStress(fConfinementTotal, boundarytensor);
    Sh1=boundarytensor.XX();
    Sh2=boundarytensor.YY();
    
    for (ellips = fGreater.size()-1; ellips >=0; ellips--) {
        
        if(fabs(Sh1)<fabs(Sh2))
        {
   
            if (co[1] > fSmaller[ellips]) {
                continue;
            }
            
            considered = true;
            REAL a = fGreater[ellips];
            REAL b = fSmaller[ellips];
            REAL xadjust = a*sqrt(1.-co[1]*co[1]/(b*b));
        
            
            if (xadjust > co[0]) {
                co[0] = xadjust;
                wasadjusted = true;
            }
            

        }
        else{
            
            if (co[0] > fSmaller[ellips]) {//sup
                continue;
            }
            considered = true;
            REAL G = fGreater[ellips];
            REAL S = fSmaller[ellips];

            REAL yadjust = G*sqrt(1.-co[0]*co[0]/(S*S));//sup
            
            if (yadjust > co[1]) {//sup
                co[1] = yadjust;
                wasadjusted = true;
            }
        }
        
        // only one adjustment can be applied
    }
    if (!wasadjusted && !considered) {
        REAL radius = std::sqrt(co[0]*co[0]+co[1]*co[1]);
        if (fabs(radius-fInnerRadius) > 1.e-7) {
            wasadjusted = true;
        }
        co[0] *= fInnerRadius/radius;
        co[1] *= fInnerRadius/radius;
    }
    
    return wasadjusted;
}

/// return the largest y-coordinate belonging to the last ellips
// this value will be used to identify the meaningful points to match the next ellips
REAL TPZWellBoreAnalysis::TConfig::MaxYfromLastBreakout()
{
   
 /*
    if(fGreater.size() == 0) return this->fInnerRadius;
    // we will use regula falsi
    int numell = fGreater.size();
    REAL lasta = fGreater[numell-1];
    REAL lastb = fSmaller[numell-1];
    REAL y = 0.;
    REAL dely = lastb/5;
    TPZManVector<REAL,2> co(2);
    while (fabs(dely) > lastb/1000.) {
        co[1] = y;
        co[0] = sqrt(fInnerRadius*fInnerRadius-y*y);
        REAL lastx = lasta*sqrt(1.-co[1]*co[1]/(lastb*lastb));
        ProjectNode(co);
        if (co[0] == lastx) {
            y += dely;
        }
        else
        {
            y -= dely;
            dely /= 10.;
        }
    }
    return y;

*/


    //if (fGreater.size()>1) {
 
    REAL r =fInnerRadius;
    REAL y,num,denom;
    if(fGreater.size() == 0)
    {
        return this->fInnerRadius;   
    }

    
    if(fGreater.size() == 1)
    {
        REAL a =fGreater[0];
        REAL b =fSmaller[0];
        num=sqrt(a*a-r*r);
        denom=sqrt(( (a*a)/(b*b) ) -1 );
        y=num/denom;
        return y;
    }
    else
    {
        int numell = fGreater.size();
        REAL a1 = fGreater[numell-1];
        REAL b1 = fSmaller[numell-1];
        REAL a2 =fGreater[numell-2];
        REAL b2 =fSmaller[numell-2];
        num = sqrt(a1*a1-a2*a2);
        denom = sqrt(((a1*a1)/(b1*b1))-((a2*a2)/(b2*b2)));
        y=num/denom;
        return y;
    }




}






/// Initialize the Sandler DiMaggio object and create the computational mesh
void TPZWellBoreAnalysis::TConfig::CreateComputationalMesh(int porder)
{
    if ((fCMesh.NElements() != 0) || (fCMesh.NConnects() != 0)) DebugStop();

    fCMesh.SetReference(&fGMesh);
    int defaultporder = porder;
    fCMesh.SetDefaultOrder(defaultporder);
    TPZCompMesh *compmesh1 = &fCMesh;
    int materialid = 1;
    
#ifdef DEBUG
    std::ofstream sout("CreateComputationalMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&fGMesh, sout,true);
#endif
    
    TPZTensor<STATE> boundarytensor;
    FromPhysicalDomaintoComputationalDomainStress(fConfinementTotal, boundarytensor);
    
    
#ifdef PlasticPQP
    if (fModel == EMohrCoulomb || fModel == EElastic) {
#else
    if (fModel == EMohrCoulomb) {
#endif

        bool planestrain = true;
        TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> &MC = fMCPV;
        TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> > *PlasticMC = new TPZMatElastoPlasticSest2D< TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> >(materialid,planestrain);
        PlasticMC->SetBiot(fBiotCoef);


        TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
        TPZTensor<REAL> initstress(0.),finalstress(0.);
        REAL hydro = boundarytensor.I1()+3.*fReservoirPressure*fBiotCoef;
//        hydro -= //SD.fYC.A()*SD.fYC.R();
        hydro /= 3.;
        finalstress.XX() = hydro;
        finalstress.YY() = hydro;
        finalstress.ZZ() = hydro;

        PrepareInitialMat(MC, initstress, finalstress, 10);
        initstress = finalstress;
        FromPhysicalDomaintoComputationalDomainStress(fConfinementTotal,finalstress);
        finalstress.XX() += fReservoirPressure*fBiotCoef;
        finalstress.YY() += fReservoirPressure*fBiotCoef;
        finalstress.ZZ() += fReservoirPressure*fBiotCoef;
        PrepareInitialMat(MC, initstress, finalstress, 10);
        //SD.ResetPlasticMem();

        TPZMaterial *plastic(PlasticMC);

        PlasticMC->SetPlasticity(MC);
        compmesh1->InsertMaterialObject(plastic);

        ConfigureBoundaryConditions();
        compmesh1->AutoBuild();
    }
#ifndef PlasticPQP
    else if (fModel == EElastic)
    {
        fMatEla.SetPlaneStrain(); // It must be planestrain because elastoplastic is only planestrain
        fMatEla.SetBiotAlpha(fBiotCoef);
        fMatEla.SetId(materialid);
        
        TPZTensor<REAL> finalstress(0.);
        FromPhysicalDomaintoComputationalDomainStress(fConfinementTotal,finalstress);
        finalstress.XX() += fReservoirPressure*fBiotCoef;
        finalstress.YY() += fReservoirPressure*fBiotCoef;
        finalstress.ZZ() += fReservoirPressure*fBiotCoef;
        fMatEla.SetPreStress(finalstress.XX(), finalstress.XY(), finalstress.YY(), finalstress.ZZ());
        TPZMaterial *mat = &fMatEla;
        compmesh1->InsertMaterialObject(mat);
        
        ConfigureBoundaryConditions();
        compmesh1->AutoBuild();
    }
#endif
    else if (fModel == ESandler)
    {
    
        bool planestrain=true;
        TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> &SD =fSDPV;
        TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> > *PlasticSD = new TPZMatElastoPlasticSest2D< TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> >(materialid,planestrain);
        PlasticSD->SetBiot(fBiotCoef);


        TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
        TPZTensor<REAL> initstress(0.),finalstress(0.);
        REAL hydro = boundarytensor.I1()+3.*fBiotCoef*fReservoirPressure;
        hydro -= SD.fYC.A()*SD.fYC.R();
        hydro /= 3.;
        finalstress.XX() = hydro;
        finalstress.YY() = hydro;
        finalstress.ZZ() = hydro;

        PrepareInitialMat(SD, initstress, finalstress, 10);
        initstress = finalstress;
        FromPhysicalDomaintoComputationalDomainStress(fConfinementTotal, finalstress);
        finalstress.XX() += fReservoirPressure*fBiotCoef;
        finalstress.YY() += fReservoirPressure*fBiotCoef;
        finalstress.ZZ() += fReservoirPressure*fBiotCoef;

        PrepareInitialMat(SD, initstress, finalstress, 10);
        SD.ResetPlasticMem();

        TPZMaterial *plastic(PlasticSD);

        PlasticSD->SetPlasticity(SD);
        compmesh1->InsertMaterialObject(plastic);
        ConfigureBoundaryConditions();
        
        {
            TPZMaterial *plastic = compmesh1->FindMaterial(materialid);
            REAL inner = fInnerRadius;
            REAL outer = fOuterRadius;
            REAL wellpress = fWellborePressure;
            REAL reservoirpress = fReservoirPressure;
            if (fBiotCoef < 0.) {
                DebugStop();
            }
            REAL biot = fBiotCoef;
            
            
            SetWellPressure(0.);
        }

        
        compmesh1->AutoBuild();


    }
    else
    {
        DebugStop();
    }

    
}

/// Transform from physical domain to computational domain stress
/**
 * the outcome depends on the well configuration
 */
void TPZWellBoreAnalysis::TConfig::FromPhysicalDomaintoComputationalDomainStress(TPZTensor<STATE> &physicalStress, TPZTensor<STATE> &computationalStress)
{
    switch (fWellConfig) {
        case EVerticalWell:
            computationalStress = physicalStress;
            break;
        case EHorizontalWellalongh:
            computationalStress[_XX_] = physicalStress[_YY_];
            computationalStress[_YY_] = physicalStress[_ZZ_];
            computationalStress[_ZZ_] = physicalStress[_XX_];
            break;
        case EHorizontalWellalongH:
            computationalStress[_XX_] = physicalStress[_XX_];
            computationalStress[_YY_] = physicalStress[_ZZ_];
            computationalStress[_ZZ_] = physicalStress[_YY_];
            break;
        default:
            DebugStop();
            break;
    }
}


#include "pzelasmat.h"

/// Configure the wellbore analysis to perform a linear analysis
void TPZWellBoreAnalysis::LinearConfiguration(int porder)
{
    fCurrentConfig.fCMesh.SetReference(&fCurrentConfig.fGMesh);
    int defaultporder = porder;
    fCurrentConfig.fCMesh.SetDefaultOrder(defaultporder);
    TPZCompMesh *compmesh1 = &fCurrentConfig.fCMesh;
    DebugStop();
    // configure a material that simulates reservoir compaction
    TPZElasticityMaterialSest2D *elasmat = new TPZElasticityMaterialSest2D(1);
    ConfigureLinearMaterial(*elasmat);
    compmesh1->InsertMaterialObject(elasmat);
    fCurrentConfig.ConfigureBoundaryConditions();
    fCurrentConfig.SetWellPressure();
    compmesh1->AutoBuild();
    
}


/// Set the parameters of the linear material
void TPZWellBoreAnalysis::ConfigureLinearMaterial(TPZElasticityMaterialSest2D &mat)
{
    REAL E,nu, lambda,G;

    if (fCurrentConfig.fModel == EMohrCoulomb) {
        G = fCurrentConfig.fMCPV.fER.fMu;
        lambda = fCurrentConfig.fMCPV.fER.fLambda;
    }
    else if (fCurrentConfig.fModel == EElastic)
    {
#ifdef PlasticPQP
        G = fCurrentConfig.fMCPV.fER.fMu;
        lambda = fCurrentConfig.fMCPV.fER.fLambda;
#else
        G = fCurrentConfig.fMatEla.GetMu();
        lambda = fCurrentConfig.fMatEla.GetLambda();
#endif
        
    }
    else if (fCurrentConfig.fModel == ESandler)
    {
        G = fCurrentConfig.fSDPV.fER.fMu;
        lambda = fCurrentConfig.fSDPV.fER.fLambda;
    }
    else
    {
        DebugStop();
    }

    
    E=G*(3.*lambda+2.*G)/(lambda+G);
    nu = lambda/(2.*(lambda+G));
    mat.SetElasticity(E, nu);
    mat.SetPlaneStrain();

    if (fCurrentConfig.fBiotCoef < 0.) {
        DebugStop();
    }
    
    
}

/// Test the linear matrix with vertical compaction
void TPZWellBoreAnalysis::TestLinearMaterial()
{
    TPZFMatrix<STATE> rhs;
    TPZCompMesh *compmesh1 = 0;
    compmesh1 = &fCurrentConfig.fCMesh;
    TPZSkylineStructMatrix skylstr(compmesh1);
//    skylstr.SetNumThreads(0);
    
    TPZMaterial *keepmat = compmesh1->MaterialVec()[1];
    // create a different material for vertical wells
    TPZElasticityMaterialSest2D *elasmat = new TPZElasticityMaterialSest2D(1);
    ConfigureLinearMaterial(*elasmat);
    elasmat->SetForcingFunction(keepmat->ForcingFunction());
    compmesh1->InsertMaterialObject(elasmat);
    
    TPZAnalysis an(compmesh1,false);
    an.SetStructuralMatrix(skylstr);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);


    
    an.Run();
    
    static int counter = 0;
    
    std::stringstream postprocfile;
    postprocfile << "../PlotFiles/Adensamento." << counter << ".vtk";
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("SigmaZ");
    an.DefineGraphMesh(2, scalnames, vecnames, postprocfile.str());
    an.PostProcess(1);
    counter++;
    
    std::set<int> matids;
    matids.insert(1);
    TPZManVector<STATE,3> ForceZ = an.Integrate("SigmaZ", matids);
    STATE StressZ = ForceZ[0]*4./(M_PI*(fCurrentConfig.fOuterRadius*fCurrentConfig.fOuterRadius-fCurrentConfig.fInnerRadius*fCurrentConfig.fInnerRadius));
    std::cout << "Total overburden computed in linear configuration " << ForceZ << std::endl;
    std::cout << "Average SigmaZ computed in linear configuration " << StressZ << std::endl;
    {
        long neq = compmesh1->NEquations();
        TPZFMatrix<STATE> rhs(neq,1,0.);
        skylstr.Assemble(rhs, NULL);
        STATE normrhs = Norm(rhs);
        std::cout << "NormRhs = " << normrhs << std::endl;
//        rhs.Print("Rhs");
    }
    an.Solution().Zero();
    an.LoadSolution();
    
    compmesh1->InsertMaterialObject(keepmat);
    delete elasmat;

}




/// Compute the linear elastic stiffness matrix
void TPZWellBoreAnalysis::ComputeLinearMatrix(TPZVec<long> &activeEquations)
{
    
    // For vertical wells, the material should be able to compute a vertical strain component as well in a multiphysics setting
    TPZCompMesh *compmesh1 = 0;
    compmesh1 = &fCurrentConfig.fCMesh;
    
    TPZSkylineStructMatrix skylstr(compmesh1);
    skylstr.SetNumThreads(TPZWellBoreAnalysis::TConfig::gNumThreads);
    skylstr.EquationFilter().SetActiveEquations(activeEquations);
    TPZFMatrix<STATE> rhs;

    TPZMaterial *keepmat = compmesh1->MaterialVec()[1];
    // create a different material for vertical wells
    TPZElasticityMaterialSest2D *elasmat = new TPZElasticityMaterialSest2D(1);
    ConfigureLinearMaterial(*elasmat);
    compmesh1->InsertMaterialObject(elasmat);
    fLinearMatrix = skylstr.CreateAssemble(rhs, NULL);
    compmesh1->InsertMaterialObject(keepmat);
    delete elasmat;
}

STATE TPZWellBoreAnalysis::TConfig::ComputeFarFieldWork()
{
    STATE work=0.,nx=0.,ny=0.,radius=0.;
    TPZManVector<REAL,3> ksi(3,0.),sol(3,0.),xvec(3,0.);
    TPZFMatrix<STATE> jac,jacinv,axes,invjac;
    STATE detjac;
    
    int nel = fCMesh.NElements();
    fCMesh.Reference()->ResetReference();
    fPostprocess.Mesh()->LoadReferences();
    
    for (int el=0; el<nel; el++)
    {
        
        TPZCompEl * cel = fCMesh.ElementVec()[el];
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
        if(matid!=-3)
        {
            continue;
        }
        
        TPZGeoElSide gelside(gel,2);
        TPZStack<TPZCompElSide> stackcompelside;
        gelside.EqualLevelCompElementList(stackcompelside, 0, 0);
        if(stackcompelside.size()!=1)
        {
            DebugStop();
        }
        
        TPZTransform t1(1);
        TPZCompElSide compneigh = stackcompelside[0];
        TPZGeoElSide neighbour = stackcompelside[0].Reference();
        gelside.SideTransform3(neighbour,t1);
        int sidefrom,sideto = neighbour.Element()->NSides()-1;
        sidefrom=neighbour.Side();
        TPZTransform t2 = neighbour.Element()->SideToSideTransform(sidefrom, sideto);
        
        
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        const TPZIntPoints &rule = intel->GetIntegrationRule();
        REAL weight;
        int npoints = rule.NPoints();
        
        for(int i=0;i<npoints;i++)
        {
            rule.Point(i, ksi, weight);
            
            TPZManVector<REAL,3> ksi1(1,0.), ksi2(2,0.);
            t1.Apply(ksi, ksi1);
            t2.Apply(ksi1, ksi2);
            compneigh.Element()->Solution(ksi2, 1, sol);
            
            gel->X(ksi, xvec);
            intel->GetIntegrationRule().Point(0, ksi, weight);
            radius=sqrt(xvec[0]*xvec[0]+xvec[1]*xvec[1]);
            nx=xvec[0]/radius;
            ny=xvec[1]/radius;
            gel->Jacobian(xvec,jac,axes,detjac,invjac);
            DebugStop();
            work +=(fConfinementTotal.XX()*nx*sol[0]+fConfinementTotal.YY()*ny*sol[1])*weight*fabs(detjac);
            
        }
        
    }
#ifdef DEBUG
    cout << "\n WORK "<< work << endl;
#endif
    return work;
    
}

void TPZWellBoreAnalysis::Print(std::ostream &out)
{
    out << "Number of configurations " << fSequence.size() << std::endl;
    out << "fCurrentConfig \n";
    fCurrentConfig.Print(out);
    out << "Other configurations \n";
    std::list<TPZWellBoreAnalysis::TConfig>::iterator it;
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

void TPZWellBoreAnalysis::SaveConfig(std::stringstream &strout) {
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

void TPZWellBoreAnalysis::TConfig::Print(ostream &out)
{
    out << "fInnerRadius " << fInnerRadius << endl;
    out << "fOuterRadius " << fOuterRadius << endl;
    out << "fNx " << fNx << endl;
    out << "fDelx " << fDelx << endl;
    out << "fGreater " << fGreater << endl;
    out << "fSmaller " << fSmaller << endl;
    out << "fConfinementTotal " << fConfinementTotal << endl;

    out << "fSDPV ";
    fSDPV.Print(out);
	out << "fMCPV ";
	fMCPV.Print(out);
//    out << "fMatEla ";
//    fMatEla.Print(out); //AQUIOMAR

    out << "fWellborePressure " << fWellborePressure << endl;
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
void TPZWellBoreAnalysis::TConfig::CreateGeometricMeshSlope()
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
/// Initialize the Sandler DiMaggio object and create the computational mesh
void TPZWellBoreAnalysis::TConfig::CreateComputationalMeshSlope(int porder)
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
