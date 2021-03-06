// Emacs will be in -*- Mode: c++ -*-
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  NumericalTraits class to illustrate TRAITS
//
//********************************************************
#ifndef _promote_h
#define _promote_h

// ANSI C++ include
#include <complex>


template <class T> class ADPromote {
  const T& x_; 
public:
  typedef typename T::value_type value_type;
  ADPromote(const T& x) : x_(x) {;}

  value_type val() { return x_.val();}
};

//Specialization
#define ADP_SPE(type)                  \
template <> class ADPromote< type > {  \
  const type& x_;                      \
public:                                \
  typedef type value_type;             \
  ADPromote(const type& x) : x_(x) {;} \
                                       \
  value_type val() { return x_;}       \
};                                 

ADP_SPE(double)
ADP_SPE(float)
ADP_SPE(long)
ADP_SPE(int)

#undef ADP_SPE


#include <TinyFadET/tfad.h>
#include <TinyFad/tinyfad.h>

template <int Num,class T> class TFad;

template <class A, class B>
class NumericalTraits
{
public:
};

//Specialization
template <class T> class NumericalTraits<T,T> {
public:
  typedef T promote;
};

template <int Num, class T> 
class NumericalTraits< TFad<Num,T>, TFad<Num,T> > {
public:
    typedef TFad< Num, T > promote;
};

template <int Num, class A, class B> class NumericalTraits< TFad<Num,A>, B> {
public:
  typedef typename TFad<Num,A>::value_type fad_type;
  typedef typename NumericalTraits<fad_type, B>::promote p_type;

  typedef TFad<Num,p_type> promote;
};

template <int Num, class A, class B> class NumericalTraits< B, TFad<Num,A> > {
public:
  typedef typename TFad<Num,A>::value_type fad_type;
  typedef typename NumericalTraits<fad_type, B>::promote p_type;

  typedef TFad<Num,p_type> promote;
};



template <int Num, class L, class R>
class NumericalTraits< TinyFad<Num,L>, TinyFad<Num,R> > {
public:
  typedef typename TinyFad<Num,L>::value_type lv;
  typedef typename TinyFad<Num,R>::value_type rv;
  typedef typename NumericalTraits<L,R>::promote value_type;

  typedef TinyFad<Num,value_type> promote;
};

template <int Num, class L, class R> class NumericalTraits< TinyFad<Num,L>, R> {
public:
  typedef typename TinyFad<Num,L>::value_type lv;
  typedef typename ADPromote<R>::value_type rv;
  typedef typename NumericalTraits<lv,rv>::promote value_type;

  typedef TinyFad<Num,value_type> promote;
};

template <int Num, class L, class R> class NumericalTraits< L, TinyFad<Num,R> > {
public:
  typedef typename ADPromote<L>::value_type lv;
  typedef typename TinyFad<Num,R>::value_type rv;
  typedef typename NumericalTraits<lv,rv>::promote value_type;

  typedef TinyFad<Num,value_type> promote;
};


#define NT_SPE(type1,type2,type3)                \
template <> class NumericalTraits<type1,type2> { \
public:                                          \
    typedef type3 promote;                       \
};                                               \
template <> class NumericalTraits<type2,type1> { \
public:                                          \
    typedef type3 promote;                       \
};                                               

NT_SPE(double,complex<float>,complex<double>)
NT_SPE(double,float,double)
NT_SPE(double,long,double)
NT_SPE(double,int,double)
NT_SPE(float,long,float)
NT_SPE(float,int,float)


#undef NT_SPE

#endif
