/**
 * @file pzmaterial.h
 * @brief Header file for abstract class KLMaterial.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef KLRandomField_H
#define KLRandomField_H

#include "pzreal.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzadmchunk.h"
#include "pzsave.h"
#include <iostream>
#include <string>
#include "pzcmesh.h"
#include "KLAnalysis.h"
#include <random>
class  KLRandomField
{
public:
    KLRandomField(TPZCompMesh * cmesh,KLAnalysis *klanal,TPZVec<REAL> mean, TPZVec<REAL> cov, int samples, REAL crosscorrelation);

    KLRandomField(KLRandomField &copy);

    ~KLRandomField();

    TPZVec<TPZFMatrix<REAL>> CreateNormalRandomField();

    TPZVec<TPZFMatrix<REAL>> CreateLogNormalRandomField();
    
    TPZVec<TPZFMatrix<REAL>>  CreateLogNormalRandomFieldHandF();
	
	void PrintMat ( std::string out,TPZFMatrix<REAL> mat )
	{
    	std::ofstream print ( out );
    	int row=mat.Rows();
    	int cols = mat.Cols();

    	for ( int i=0; i<row; i++ ) {
        	for ( int j=0; j<cols; j++ ) {
        	    print << mat ( i,j ) << " ";
        	}
        	print<< std::endl;
    	}


	}

private:


    TPZCompMesh * fCMesh;
    TPZVec<REAL> fMean;
    TPZVec<REAL> fCov;
    KLAnalysis *fKlAnal;
    int fSamples;
    REAL fCross;




};



#endif

