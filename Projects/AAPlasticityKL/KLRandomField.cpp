#include "KLRandomField.h"



KLRandomField::KLRandomField(TPZCompMesh * cmesh,KLAnalysis *klanal,TPZVec<REAL> mean, TPZVec<REAL> cov, int samples,REAL crosscorrelation)
{

    fCMesh =cmesh;
    fKlAnal =klanal;
    fMean =mean;
    fCov =cov;
    fSamples=samples;
    fCross = crosscorrelation;
}

KLRandomField::KLRandomField(KLRandomField &copy):fCMesh(copy.fCMesh),fKlAnal(copy.fKlAnal),fMean(copy.fMean),fCov(copy.fCov),fSamples(copy.fSamples),fCross(copy.fCross)
{

}

KLRandomField::~KLRandomField()
{

}

TPZVec<TPZFMatrix<REAL>> KLRandomField::CreateNormalRandomField()
{

}

TPZVec<TPZFMatrix<REAL>> KLRandomField::CreateLogNormalRandomField()
{
    TPZFMatrix<REAL>  PHIt, PHI, vect;
    TPZFMatrix<REAL> val = fKlAnal->GetEigenVal();
    TPZFMatrix<REAL> vec = fCMesh->Solution();


    vec.Transpose( &vect );

    int M= val.Rows();
    TPZFMatrix<REAL>  Identity(M,M,0.);
    for ( int i = 0; i < M; i++ ) Identity(i,i) = sqrt ( fabs ( val(i,0) ) );

    Identity.Multiply ( vect, PHI );

	//PHI = vect;
	
    PHI.Transpose ( &PHIt );

    std::normal_distribution<double> distribution ( 0., 1. );

    TPZFMatrix<REAL>  THETA ( M, fSamples, 0. ), THETA2 ( M, fSamples, 0. );
    for ( int isample = 0; isample < fSamples; isample++ ) {
        for ( int irdvar = 0; irdvar < M; irdvar++ ) {
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            REAL xic = distribution ( generator );
            REAL xiphi = distribution ( generator );
            THETA(irdvar,isample) = xic;
            THETA2(irdvar,isample) = xic * fCross + xiphi * sqrt ( 1 - fCross * fCross );
        }
    }

	string outcohesvas="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/outcohesvas.txt";
	string outphisvas="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/outphisvas.txt";
	
	
	PrintMat(outcohesvas,THETA);
    PrintMat(outphisvas,THETA2);

    TPZFMatrix<REAL>  hhatphi, hhatcoes;

    PHIt.Multiply( THETA, hhatcoes );
    PHIt.Multiply( THETA2, hhatphi );

    TPZVec<TPZFMatrix<REAL>> HHAT(2);


    REAL mean = fMean[0];
    REAL sdev = fCov[0] * mean;
    REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean ),2 ) ) );
    REAL lambda = log ( mean ) - xi * xi / 2.;
    for ( int i = 0; i < hhatcoes.Rows(); i++ ) {
        for ( int j = 0; j < hhatcoes.Cols(); j++ ) {
			REAL temp =  hhatcoes(i,j);
            hhatcoes(i,j) = exp ( lambda + xi * temp );
        }
    }

    mean = fMean[1];
    sdev = fCov[1] * mean;
    xi = sqrt ( log ( 1 + pow ( ( sdev / mean ), 2 ) ) );
    lambda = log ( mean ) - xi * xi / 2.;
    for ( int i = 0; i < hhatphi.Rows(); i++ ) {
        for ( int j = 0; j < hhatphi.Cols(); j++ ) {
			REAL temp = hhatphi(i,j);
            hhatphi(i,j) = exp ( lambda + xi *temp );
        }

    }

    HHAT[0] = hhatcoes;
    HHAT[1] = hhatphi;
	
    return HHAT;
}

TPZVec<TPZFMatrix<REAL>> KLRandomField::CreateLogNormalRandomFieldHandF()
{
    TPZFMatrix<REAL>  PHIt, PHI, vect;
    TPZFMatrix<REAL> val = fKlAnal->GetEigenVal();
    TPZFMatrix<REAL> vec = fCMesh->Solution();


    vec.Transpose( &vect );

    int M= val.Rows();
    TPZFMatrix<REAL>  Identity(M,M,0.);
    for ( int i = 0; i < M; i++ ) Identity(i,i) = sqrt ( fabs ( val(i,0) ) );

    Identity.Multiply ( vect, PHI );

	//PHI = vect;
	
    PHI.Transpose ( &PHIt );

    std::normal_distribution<double> distribution ( 0., 1. );

    TPZFMatrix<REAL>  THETA ( M, fSamples, 0. ), THETA2 ( M, fSamples, 0. );
    for ( int isample = 0; isample < fSamples; isample++ ) {
        for ( int irdvar = 0; irdvar < M; irdvar++ ) {
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            REAL xic = distribution ( generator );
            REAL xiphi = distribution ( generator );
            THETA(irdvar,isample) = xic;
            THETA2(irdvar,isample) = xic * fCross + xiphi * sqrt ( 1 - fCross * fCross );
        }
    }

	string outcohesvas="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/outcohesvas.txt";
	string outphisvas="/home/diogo/Dropbox/Projeto-Landslides/MonteCarlo/rffolder/outphisvas.txt";
	
	
	PrintMat(outcohesvas,THETA);
    PrintMat(outphisvas,THETA2);
    
    TPZFMatrix<REAL>  hhatphi, hhatcoes,fphi,fcoes;

    PHIt.Multiply( THETA, hhatcoes );
    PHIt.Multiply( THETA2, hhatphi );

    TPZVec<TPZFMatrix<REAL>> HHAT(2);


    REAL mean = fMean[0];
    REAL sdev = fCov[0] * mean;
    REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean ),2 ) ) );
    REAL lambda = log ( mean ) - xi * xi / 2.;
    for ( int i = 0; i < hhatcoes.Rows(); i++ ) {
        for ( int j = 0; j < hhatcoes.Cols(); j++ ) {
			REAL temp =  hhatcoes(i,j);
            hhatcoes(i,j) = exp ( lambda + xi * temp );
        }
    }

    mean = fMean[1];
    sdev = fCov[1] * mean;
    xi = sqrt ( log ( 1 + pow ( ( sdev / mean ), 2 ) ) );
    lambda = log ( mean ) - xi * xi / 2.;
    for ( int i = 0; i < hhatphi.Rows(); i++ ) {
        for ( int j = 0; j < hhatphi.Cols(); j++ ) {
			REAL temp = hhatphi(i,j);
            hhatphi(i,j) = exp ( lambda + xi *temp );
        }

    }

    HHAT[0] = hhatcoes;
    HHAT[1] = hhatphi;
	
    return HHAT;
}
