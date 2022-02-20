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

private:


    TPZCompMesh * fCMesh;
    TPZVec<REAL> fMean;
    TPZVec<REAL> fCov;
    KLAnalysis *fKlAnal;
    int fSamples;
    REAL fCross;




};



#endif

