// $Id: TPZPlasticState.h,v 1.8 2009-07-19 05:47:18 erick Exp $

#ifndef TPZPLASTICSTATE_H
#define TPZPLASTICSTATE_H

#include "fadType.h" 
#include "TPZTensor.h"
#include <iostream>

/**
 * This class holds the complete set of state variables to define a valid elastoplastic strain state
 */
template <class T>
class TPZPlasticState
{

    // class member variables
	
public:	
	// Tensors to hold the total and plastic strain tensors
	TPZTensor<T> fEpsT, fEpsP;
	
	// Plastic damage variable
	T fAlpha;
	
	TPZVec<T> fmatprop;
    
    TPZVec<T> fflux;
    
    T fpressure;
    

public:
	
	/**
	 * Default constructor - all values set to zero
	 */
    TPZPlasticState(): fEpsT(), fEpsP(), fAlpha(T(0.) ),fpressure(T(0.) ),fmatprop(),fflux(){ }
		
	/**
	 * Constructor enabling predefinition of Alpha
	 */
	TPZPlasticState(const T & alpha):fEpsT(T(0.)), fEpsP(T(0.)),fpressure(T(0.)) ,fAlpha(alpha),fmatprop(T(0.)),fflux(T(0.)){ }
	
	/**
	 * Copy constructor
	 */
	TPZPlasticState(const TPZPlasticState<T> & source):
		fEpsT(source.fEpsT), fEpsP(source.fEpsP), fAlpha(source.fAlpha),fmatprop(source.fmatprop),fflux(source.fflux), fpressure(source.fpressure){ }
	
	/**
	 * Default destructor
	 */
	~TPZPlasticState(){ }
	
	/**
	 * Operator =
	 */
	const TPZPlasticState<T> & operator=(const TPZPlasticState<T> & source);
	
	/**
	 * Operator -=
	 */
	const TPZPlasticState<T> & operator-=(const TPZPlasticState<T> & source);
	
	/**
	 * Operator +=
	 */
	const TPZPlasticState<T> & operator+=(const TPZPlasticState<T> & source);
	
	/**
	 * Operator *=
	 */
	const TPZPlasticState<T> & operator*=(const TPZPlasticState<T> & source);
	
	/**
	 * Similar to Operator=, but allows copies among different template specializations.
	 * When using FAD class Types as template classes, NO DERIVATIVES are copied using this functions,
	 * enabling copies from REAL to FAD types and also vice-versa.
	 */
	template <class T1>
	void CopyTo(TPZPlasticState<T1> & target) const;

	/**
	 * Operator<<
	 */
    friend std::ostream& operator<<( std::ostream& Out, const TPZPlasticState<T> & s )
	{
		s.Print(Out);
		return Out;
	}
	
	/**
	 * More complete then Operator<< because it allows derivatives supression.
	 */
	void Print(std::ostream& Out, int fadDerivatives = 1)const;
	
	// class Access Members (needed for const PlasticState access)
	
	const TPZTensor<T> & EpsT() const 
		{ return fEpsT; }
	const TPZTensor<T> & EpsP() const 
		{ return fEpsP; }
	const T & Alpha() const 
		{ return fAlpha; }
	const TPZVec<T> & MatProp() const 
		{ return fmatprop; }
	const TPZVec<T> & Flux() const 
		{ return fflux; }
    const T & Pressure() const 
		{ return fpressure; }


    
    void Write(TPZStream &buf) const;
    
    void Read(TPZStream &buf);

//    template<>

};

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator=(const TPZPlasticState<T> & source)
{
	fEpsT = source.EpsT();
	fEpsP = source.EpsP();
	fAlpha = source.Alpha();
	fmatprop=source.MatProp();
    fflux = source.Flux();
    fpressure =source.Pressure();
		
	return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator+=(const TPZPlasticState<T> & source)
{
	fEpsT += source.EpsT();
	fEpsP += source.EpsP();
	fAlpha+= source.Alpha();
/*	fPhi += source.fPhi();*/
		
	return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator-=(const TPZPlasticState<T> & source)
{
	fEpsT -= source.EpsT();
	fEpsP -= source.EpsP();
	fAlpha-= source.Alpha();
/*	fPhi -= source.fPhi();	*/
	return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator*=(const TPZPlasticState<T> & source)
{
	fEpsT *= source.EpsT();
	fEpsP *= source.EpsP();
	fAlpha*= source.Alpha();
/*	fPhi *= source.fPhi();*/
		
	return *this;
}

template <class T>
inline void TPZPlasticState<T>::Print(std::ostream& Out, int fadDerivatives)const
{
	if(fadDerivatives)
	{
	    Out << "\tfEpsT = ";
	    for(int i = 0; i < 6; i++)Out << fEpsT[i] << " ";
	    Out << "\n\tfEpsP = ";
	    for(int i = 0; i < 6; i++)Out << fEpsP[i] << " ";
		Out << "\n\tfAlpha = " << fAlpha;
	}else{
	    Out << "\tfEpsT = ";
	    for(int i = 0; i < 6; i++)Out << shapeFAD::val(fEpsT[i]) << " ";
	    Out << "\n\tfEpsP = ";
	    for(int i = 0; i < 6; i++)Out << shapeFAD::val(fEpsP[i]) << " ";
		Out << "\n\tfAlpha = " << shapeFAD::val(fAlpha);
	}
}

template <class T>
template <class T1>
void TPZPlasticState<T>::CopyTo(TPZPlasticState<T1> & target) const
{
	EpsT().CopyTo(target.fEpsT);
	EpsP().CopyTo(target.fEpsP);
	target.fAlpha = shapeFAD::val( Alpha() );
    target.fpressure = shapeFAD::val( Pressure() );
	target.fmatprop.Resize(2);
    target.fflux.Resize(2);
	target.fmatprop[0]= shapeFAD::val(fmatprop[0]);
	target.fmatprop[1]= shapeFAD::val(fmatprop[1]);
    target.fflux[0]= shapeFAD::val(fflux[0]);
    target.fflux[1]= shapeFAD::val(fflux[1]);
    
}

template<class T>
void TPZPlasticState<T>::Write(TPZStream &buf) const
{
    DebugStop();
}

template<class T>
void TPZPlasticState<T>::Read(TPZStream &buf)
{
    DebugStop();
}

template<>
inline void TPZPlasticState<STATE>::Write(TPZStream &buf) const
{
    buf.Write(&fEpsT.fData[0],6);
    buf.Write(&fEpsP.fData[0],6);
    buf.Write(&fAlpha);
}

template<>
inline void TPZPlasticState<STATE>::Read(TPZStream &buf)
{
    buf.Read(&fEpsT.fData[0],6);
    buf.Read(&fEpsP.fData[0],6);
    buf.Read(&fAlpha);
}



#endif
