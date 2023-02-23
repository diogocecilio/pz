#include "TPZMohrCoulombVoigt.h"


TPZMohrCoulombVoigt::TPZMohrCoulombVoigt() : fPhi ( 0. ),fPsi ( 0. ), fc ( 0. ), fEpsPlasticBar ( 0. ), fER()
{

}

TPZMohrCoulombVoigt::TPZMohrCoulombVoigt ( REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER ) : fPhi ( Phi ), fPsi ( Psi ), fc ( c ), fEpsPlasticBar ( 0. ), fER ( ER )
{

}

TPZMohrCoulombVoigt::TPZMohrCoulombVoigt ( const TPZMohrCoulombVoigt &cp ) // : 	fPhi(cp.fPhi), fPsi(cp.fPsi), fc(cp.fc), fEpsPlasticBar(cp.fEpsPlasticBar), fER(cp.fER)
{
    TPZMohrCoulombVoigt::operator= ( cp );
}

TPZMohrCoulombVoigt & TPZMohrCoulombVoigt::operator= ( const TPZMohrCoulombVoigt &cp )
{
    fPhi = cp.fPhi;
    fPsi = cp.fPsi;
    fc = cp.fc;
    fEpsPlasticBar = cp.fEpsPlasticBar;
    fER = cp.fER;
    return *this;
}


void TPZMohrCoulombVoigt::Read ( TPZStream &buf )
{
    buf.Read ( &fPhi );
    buf.Read ( &fPsi );
    buf.Read ( &fc );
    buf.Read ( &fEpsPlasticBar );
    fER.Read ( buf );

}

void TPZMohrCoulombVoigt::Write ( TPZStream &buf ) const
{
    buf.Write ( &fPhi );
    buf.Write ( &fPsi );
    buf.Write ( &fc );
    buf.Write ( &fEpsPlasticBar );
    fER.Write ( buf );
}

template<class T>
T TPZMohrCoulombVoigt::PhiPlane ( const TPZVec<T> &sigma ) const
{

}


REAL TPZMohrCoulombVoigt::ReturnMapPlane (  TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew)
{

    TPZTensor<REAL> avec;
     REAL a =A(sigma_trial);
     REAL I1 = sigma_trial.I1();
     REAL J2 = sigma_trial.J2();
     REAL dadt = dAdt(sigma_trial);
     REAL d2adt = d2Adt(sigma_trial);

     REAL c1=C1();
     REAL c2=C2 ( sigma_trial, a, dadt);
     REAL c3 = C3 ( sigma_trial,dadt );
     REAL c4=C4 (  sigma_trial, dadt,d2adt);
     TPZTensor<REAL> di1 = sigma_trial.dI1();
     TPZTensor<REAL> dj2 = sigma_trial.dJ2();
     TPZTensor<REAL> dj3 = sigma_trial.dJ3();

     REAL f1 = -(fc*cos(fPhi)) + a*sqrt(J2) + 0.3333333333333333*I1*sin(fPhi);

     di1*=c1;
     dj2*=c2;
     dj3*=c3;

     //TPZTensor<REAL> avec;
     avec=di1;
     avec+=dj2;
     avec+=dj3;

   //  avec.Print(std::cout);

//
      REAL dlvoigt=f1;
      TPZFMatrix<REAL> cmat = GetElasticMatrix();
//
    //  cmat.Print(std::cout);

      TPZFMatrix<REAL> amat,amatT,temp,temp2,tempT;

      FromTensorToMatVoigt(avec,amat);

      amat.Transpose(&amatT);

      cmat.Multiply(amat,temp);
      temp.Transpose(&tempT);
      //amat.Print(std::cout);
      //temp.Print(std::cout);
      tempT.Multiply(amat,temp2);
     //cmat.Multiply(avec,

      dlvoigt *= 1./temp2.Get(0,0);
      //sigtr - dlvoigt cmat.a;
      FromMatToTensor(temp,sigma_proj);
      sigma_proj*=-dlvoigt;
      sigma_proj+=sigma_trial;

      //sigma_proj.Print(std::cout);

      TPZFMatrix<REAL> Q1(6,6,0.),dadsig;
      Q1(0,0)=1.;Q1(1,1)=1.;Q1(2,2)=1.;Q1(3,3)=1.;Q1(4,4)=1.;Q1(5,5)=1.;

      dadsig =  dAdsig(  sigma_trial,  a, dadt, d2adt );

      TPZFMatrix<REAL> tempx,tempy,tempz,ect;

      dadsig.Multiply(cmat,tempx);
      tempx*=dlvoigt;
      Q1-=tempx;
      tempx.Zero();

      amat.Multiply(amatT,tempx);

      cmat.Multiply(tempx,tempy);

      tempy.Multiply(cmat,tempz);

      tempz*=1./temp2.Get(0,0);

      ect=cmat;

      ect-=tempz;

      ect.Multiply(Q1,dep);

     // dep.Print(std::cout);

      return dlvoigt;
}


bool TPZMohrCoulombVoigt::ReturnMapLeftEdge ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew)
{
    TPZTensor<REAL> avec,bvec;
    TPZFMatrix<REAL> cmat = GetElasticMatrix();
     REAL I1 = sigma_trial.I1();
     REAL J2 = sigma_trial.J2();

     REAL a =A(sigma_trial);
     REAL dadt = dAdt(sigma_trial);
     REAL d2adt = d2Adt(sigma_trial);

     FlowVector(sigma_trial,a,dadt,d2adt,avec);

     REAL a3 =A3(sigma_trial);
     REAL da3dt = dA3dt(sigma_trial);
     REAL d2a3dt = d2A3dt(sigma_trial);

     FlowVector(sigma_trial,a3,da3dt,d2a3dt,bvec);


     TPZFMatrix<REAL> temp,temp2,temp3,temp4;


     FromTensorToMatVoigt(avec,temp);
     cmat.Multiply(temp,temp2);
     REAL ax = Dot(temp,temp2);


     FromTensorToMatVoigt(bvec,temp);
     cmat.Multiply(temp,temp2);
     REAL bx =Dot(temp,temp2);


     FromTensorToMatVoigt(avec,temp);
     FromTensorToMatVoigt(bvec,temp2);
     cmat.Multiply(temp2,temp3);
     REAL dx = Dot(temp,temp3);

     REAL f11,f22;

     Yield(sigma_trial,a,f11);
     Yield(sigma_trial,a3,f22);


   REAL q = ax*bx - dx*dx;
   REAL  dl11 = (bx*f11 - dx*f22)/(ax* bx - dx*dx);
   REAL dl22 = (ax*f22 - dx*f11)/(ax*bx - dx*dx);


   FromTensorToMatVoigt(avec,temp);
   FromTensorToMatVoigt(bvec,temp2);
   cmat.Multiply(temp,temp3);
   cmat.Multiply(temp2,temp4);
   temp3*=-dl11;
   temp4*=-dl22;

   temp3+=temp4;

   FromMatToTensor(temp3,sigma_proj);
   sigma_proj+=sigma_trial;

   TPZFMatrix<REAL> dadsig1 = dAdsig(sigma_trial,  a, dadt, d2adt);
   TPZFMatrix<REAL> dadsig2 = dAdsig(sigma_trial,  a3, da3dt, d2a3dt);

   //Outer[Times, avec, avec]

    TPZFMatrix<REAL> avectemp,avectempT,bvectemp,bvectempT,aaT,abT,baT,bbT;

    FromTensorToMatVoigt(avec,avectemp);
    FromTensorToMatVoigt(bvec,bvectemp);

    avectemp.Transpose(&avectempT);
    bvectemp.Transpose(&bvectempT);

    //Outer[Times, avec, avec]
    avectemp.Multiply(avectempT,aaT);
    //Outer[Times, avec, bvec]
    avectemp.Multiply(bvectempT,abT);
    //Outer[Times, bvec, avec]
    bvectemp.Multiply(avectempT,baT);
    //Outer[Times, bvec, bvec]
    bvectemp.Multiply(bvectempT,bbT);

    TPZFMatrix<REAL> cmataaT,cmataaTcmat,cmatabT,cmatabTcmat,cmatbaT,cmatbaTcmat,cmatbbT,cmatbbTcmat,tempfinal,et2;

    //cmat.Outer[Times, avec, avec]
    cmat.Multiply(aaT,cmataaT);
    //cmat.Outer[Times, avec, avec].cmat
    cmataaT.Multiply(cmat,cmataaTcmat);

    cmataaTcmat*=bx;

    cmat.Multiply(abT,cmatabT);
    cmatabT.Multiply(cmat,cmatabTcmat);

    cmatabTcmat*=dx;


    cmat.Multiply(baT,cmatbaT);
    cmatbaT.Multiply(cmat,cmatbaTcmat);

    cmatbaTcmat*=dx;

    cmat.Multiply(bbT,cmatbbT);
    cmatbbT.Multiply(cmat,cmatbbTcmat);

    cmatbbTcmat*=ax;

    tempfinal=cmataaTcmat;

    tempfinal-=cmatabTcmat;
    tempfinal-=cmatbaTcmat;
    tempfinal+=cmatbbTcmat;

    tempfinal*=1./q;

    et2=cmat;

    et2-=tempfinal;

    TPZFMatrix<REAL> T(6,6,0.);
    T(0,0)=1.;T(1,1)=1.;T(2,2)=1.;T(3,3)=1.;T(4,4)=1.;T(5,5)=1.;

    TPZFMatrix<REAL> partea,parteb;
    dadsig1.Multiply(cmat,partea);
    dadsig2.Multiply(cmat,parteb);
    partea*=-dl11;
    parteb*=-dl22;
    T+=partea;
    T+=parteb;

    et2.Multiply(T,dep);

}

bool TPZMohrCoulombVoigt::ReturnMapRightEdge ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew)
{
    TPZTensor<REAL> avec,bvec;
    TPZFMatrix<REAL> cmat = GetElasticMatrix();
     REAL I1 = sigma_trial.I1();
     REAL J2 = sigma_trial.J2();

     REAL a =A(sigma_trial);
     REAL dadt = dAdt(sigma_trial);
     REAL d2adt = d2Adt(sigma_trial);

     FlowVector(sigma_trial,a,dadt,d2adt,avec);

     REAL a2 =A2(sigma_trial);
     REAL da2dt = dA2dt(sigma_trial);
     REAL d2a2dt = d2A2dt(sigma_trial);

     FlowVector(sigma_trial,a2,da2dt,d2a2dt,bvec);


     TPZFMatrix<REAL> temp,temp2,temp3,temp4;


     FromTensorToMatVoigt(avec,temp);
     cmat.Multiply(temp,temp2);
     REAL ax = Dot(temp,temp2);


     FromTensorToMatVoigt(bvec,temp);
     cmat.Multiply(temp,temp2);
     REAL bx =Dot(temp,temp2);


     FromTensorToMatVoigt(avec,temp);
     FromTensorToMatVoigt(bvec,temp2);
     cmat.Multiply(temp2,temp3);
     REAL dx = Dot(temp,temp3);

     REAL f11,f22;

     Yield(sigma_trial,a,f11);
     Yield(sigma_trial,a2,f22);

   REAL q = ax*bx - dx*dx;
   REAL  dl11 = (bx*f11 - dx*f22)/(ax* bx - dx*dx);
   REAL dl22 = (ax*f22 - dx*f11)/(ax*bx - dx*dx);


   FromTensorToMatVoigt(avec,temp);
   FromTensorToMatVoigt(bvec,temp2);
   cmat.Multiply(temp,temp3);
   cmat.Multiply(temp2,temp4);
   temp3*=-dl11;
   temp4*=-dl22;

   temp3+=temp4;

   FromMatToTensor(temp3,sigma_proj);
   sigma_proj+=sigma_trial;

   TPZFMatrix<REAL> dadsig1 = dAdsig(sigma_trial,  a, dadt, d2adt);
   TPZFMatrix<REAL> dadsig2 = dAdsig(sigma_trial,  a2, da2dt, d2a2dt);

   //Outer[Times, avec, avec]

    TPZFMatrix<REAL> avectemp,avectempT,bvectemp,bvectempT,aaT,abT,baT,bbT;

    FromTensorToMatVoigt(avec,avectemp);
    FromTensorToMatVoigt(bvec,bvectemp);

    avectemp.Transpose(&avectempT);
    bvectemp.Transpose(&bvectempT);

    //Outer[Times, avec, avec]
    avectemp.Multiply(avectempT,aaT);
    //Outer[Times, avec, bvec]
    avectemp.Multiply(bvectempT,abT);
    //Outer[Times, bvec, avec]
    bvectemp.Multiply(avectempT,baT);
    //Outer[Times, bvec, bvec]
    bvectemp.Multiply(bvectempT,bbT);

    TPZFMatrix<REAL> cmataaT,cmataaTcmat,cmatabT,cmatabTcmat,cmatbaT,cmatbaTcmat,cmatbbT,cmatbbTcmat,tempfinal,et2;

    //cmat.Outer[Times, avec, avec]
    cmat.Multiply(aaT,cmataaT);
    //cmat.Outer[Times, avec, avec].cmat
    cmataaT.Multiply(cmat,cmataaTcmat);

    cmataaTcmat*=bx;

    cmat.Multiply(abT,cmatabT);
    cmatabT.Multiply(cmat,cmatabTcmat);

    cmatabTcmat*=dx;


    cmat.Multiply(baT,cmatbaT);
    cmatbaT.Multiply(cmat,cmatbaTcmat);

    cmatbaTcmat*=dx;

    cmat.Multiply(bbT,cmatbbT);
    cmatbbT.Multiply(cmat,cmatbbTcmat);

    cmatbbTcmat*=ax;

    tempfinal=cmataaTcmat;

    tempfinal-=cmatabTcmat;
    tempfinal-=cmatbaTcmat;
    tempfinal+=cmatbbTcmat;

    tempfinal*=1./q;

    et2=cmat;

    et2-=tempfinal;

    TPZFMatrix<REAL> T(6,6,0.);
    T(0,0)=1.;T(1,1)=1.;T(2,2)=1.;T(3,3)=1.;T(4,4)=1.;T(5,5)=1.;

    TPZFMatrix<REAL> partea,parteb;
    dadsig1.Multiply(cmat,partea);
    dadsig2.Multiply(cmat,parteb);
    partea*=-dl11;
    parteb*=-dl22;
    T+=partea;
    T+=parteb;

    et2.Multiply(T,dep);
}

bool TPZMohrCoulombVoigt::ReturnMapApex ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew )
{
    REAL ccotanphi = fc*cos(fPhi)/sin(fPhi);
    sigma_proj.XX()=ccotanphi;sigma_proj.YY()=ccotanphi;sigma_proj.ZZ()=ccotanphi;
    TPZFMatrix<REAL> dumb(6,6,10.e-6);
    dep=dumb;
}



void TPZMohrCoulombVoigt::ProjectSigmaDep ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew )
{


    REAL I1 = sigma_trial.I1();
    REAL J2 = sigma_trial.J2();
    REAL f1=0.;

    REAL a;
    REAL dadt = dAdt(sigma_trial);
    REAL d2adt = d2Adt(sigma_trial);


    a =A(sigma_trial);
    Yield(sigma_trial,a,f1);

    if(f1>0)
    {
        TPZTensor<REAL> avec,bvec;
        REAL dlamb=0.;
        REAL epsbarnew=0.;

         dlamb = ReturnMapPlane ( sigma_trial, sigma_proj,dep,  epsbarnew);
         //sigma_proj.Print(std::cout);

        FlowVector(sigma_trial,a,dadt,d2adt,avec);


        REAL aa =A(sigma_proj);
        Yield(sigma_proj,aa,f1);

        //cout <<"aa = "<< aa << std::endl;
       //sigma_proj.Print(std::cout);
       //cout <<"f1 = "<< f1 << std::endl;

        if(f1<1.)
        {
            //cout<< "projetou na f1";
            return;
        }
        else
        {
            REAL D2 = -tan(3.* theta(sigma_trial))/(2.*J2);
            REAL D3 = -sqrt(3.)/( 2.*pow(J2,1.5)* cos(3.* theta(sigma_trial ) ));
            TPZTensor<REAL> di1 = sigma_trial.dI1();
            TPZTensor<REAL> dj2 = sigma_trial.dJ2();
            TPZTensor<REAL> dj3 = sigma_trial.dJ3();
            //dj2.Print(std::cout);
            dj2*=D2;
            dj3*=D3;
            dj2+=dj3;
            TPZFMatrix<REAL> cc;
            FromTensorToMatVoigt(dj2,cc);
            TPZFMatrix<REAL>  temp,amat,temp2,cmat,avecmat,ccT;
            //ccT.Print(std::cout);
            cc.Transpose(&ccT);
            FromTensorToMatVoigt(avec,avecmat);
            cmat=GetElasticMatrix();
            cmat.Multiply(avecmat,temp);
            ccT.Multiply(temp,temp2);
            temp2*=-dlamb;
            REAL dt = temp2.Get(0,0)*180./M_PI;
            if(dt<0)
            {

                ReturnMapRightEdge(sigma_trial, sigma_proj,dep, epsbarnew);

                //cout <<"sigproj rigth = "<< std::endl;
                //sigma_proj.Print(std::cout);


                REAL a2=A2(sigma_proj);
                REAL da2dt = dA2dt(sigma_proj);
                REAL d2a2dt = d2A2dt(sigma_proj);

                FlowVector(sigma_proj,a2,da2dt,d2a2dt,bvec);

                TPZFMatrix<REAL> avecmat,bvecmat;
                FromTensorToMatVoigt(avec,avecmat);
                FromTensorToMatVoigt(bvec,bvecmat);
                //cout <<"avec "<< std::endl;
                //avecmat.Print(std::cout);
                //cout <<"bvec "<< std::endl;
                //bvecmat.Print(std::cout);

                REAL dot = Dot(avecmat,bvecmat);
                REAL pronorm=Norm(bvecmat)*Norm(bvecmat);
                REAL val=dot/pronorm;
                //if(val<-1)val=-1.;
                //if(val>1)val=1.;
                //REAL beta=acos(val)*180/M_PI;
                if(fabs(val)<0.000174533)//acos(0.000174533) = 89.99 graus
                {
                    //cout<< "projetou apex";
                    ReturnMapApex(sigma_trial, sigma_proj,dep, epsbarnew);
                }
                //cout<< "projetou na rigth";
            }
            else
            {
                ReturnMapLeftEdge(sigma_trial, sigma_proj,dep, epsbarnew);

                REAL a2=A2(theta(sigma_proj));
                REAL da2dt = dA2dt(sigma_proj);
                REAL d2a2dt = d2A2dt(sigma_proj);

                FlowVector(sigma_proj,a2,da2dt,d2a2dt,bvec);

                TPZFMatrix<REAL> avecmat,bvecmat;
                FromTensorToMatVoigt(avec,avecmat);
                FromTensorToMatVoigt(bvec,bvecmat);

                REAL dot = Dot(avecmat,bvecmat);
                REAL pronorm=Norm(bvecmat)*Norm(bvecmat);
                REAL val=dot/pronorm;
                //if(val<-1)val=-1.;
                //if(val>1)val=1.;
                //REAL beta=acos(val)*180/M_PI;
                 if(fabs(val)<0.000174533)///acos(0.000174533) = 89.99 graus
                {
                    //cout<< "projetou apex";
                    ReturnMapApex(sigma_trial, sigma_proj,dep, epsbarnew);
                }
                //cout<< "projetou na left";
            }
        }
    }
    else
    {
        sigma_proj=sigma_trial;
        dep=GetElasticMatrix();
        epsbarnew=0;

    }

}

