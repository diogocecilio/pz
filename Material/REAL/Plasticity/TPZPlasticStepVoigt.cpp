#include "TPZPlasticStepVoigt.h"
#include "TPZElasticResponse.h"


#include "TPZMohrCoulombVoigt.h"
#include "TPZElasticResponse.h"


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
{

}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
{

}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::TangentOperator(TPZFMatrix<REAL> & gradient,TPZTensor<REAL>::TPZDecomposed & eps_eigen_system, TPZTensor<REAL>::TPZDecomposed & sig_eigen_system, TPZFMatrix<REAL> & Tangent){
	

	
}


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv)
{

}


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{

    std::cout<< " \n this method is not implemented in PlasticStepPV. ";
    DebugStop();
    
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal)
{
 
}

template <class YC_t, class ER_t>
TPZPlasticState<STATE>  TPZPlasticStepVoigt<YC_t, ER_t>::GetState() const
{
    return fN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{

}


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::SetState(const TPZPlasticState<REAL> &state)
{
    fN=state;
}



template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::Read(TPZStream &buf)
{
	fYC.Read(buf);
	fER.Read(buf);
	buf.Read(&fResTol);
	buf.Read(&fMaxNewton);
	fN.Read(buf);
}



template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::Write(TPZStream &buf) const
{
    fYC.Write(buf);
    fER.Write(buf);
    buf.Write(&fResTol);
    buf.Write(&fMaxNewton);
    fN.Write(buf);
    
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM,TPZTensor<STATE> &copy)
{
    FNM.Resize(6,1);
    copy.XX()=FNM(0,0);
    copy.XY()=FNM(1,0);
    copy.XZ()=FNM(2,0);
    copy.YY()=FNM(3,0);
    copy.YZ()=FNM(4,0);
    copy.ZZ()=FNM(5,0);
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::CopyFromTensorToFMatrix(TPZTensor<STATE> tensor,TPZFMatrix<STATE> &copy)
{
    
    copy(0,0)=tensor.XX();
    copy(1,0)=tensor.XY();
    copy(2,0)=tensor.XZ();
    copy(3,0)=tensor.YY();
    copy(4,0)=tensor.YZ();
    copy(5,0)=tensor.ZZ();
}


template class TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse>;
