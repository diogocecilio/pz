/**
 * @file
 * @brief Contains declaration of TPZElementMatrix struct which associates an element matrix with the coeficients of its contribution in the global stiffness matrix
 */

#ifndef ELMATHPP
#define ELMATHPP


#include "pzmatrix.h"
#include "pzblock.h"
#include "pzconnect.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzstack.h"


/**
 * @brief This class associates an element matrix with the coeficients of its contribution in the global stiffness matrix. \ref interpolation "Aproximation space"
 * @ingroup interpolation
 */
/**
 This class groups all information associated with an element stiffness matrix so that it can be used independent of the element object itself
 Objects of this class provide storage as well for the constrained stiffness matrix, i.e. the stiffness matrix from which the constrained connects have been eliminated
 In future versions, the computation of the contraints will be incorporated in a method of this class
 */
struct TPZElementMatrix {
	
    enum MType{Unknown = 0, EF = 1, EK = 2};
	
	MType fType;
	
	TPZCompMesh * fMesh;
	
	/** @brief Vector of pointers to TPZConnect objects*/
	TPZStack<long> fConnect;
	/** @brief Pointer to a blocked matrix object*/
	//TPZFNMatrix<1000> fMat;
	TPZFNMatrix<1000, STATE> fMat;
	/** @brief Block structure associated with fMat*/
	//TPZBlock<REAL> fBlock;
	TPZBlock<STATE> fBlock;
	/** @brief Vector of all nodes connected to the element*/
	TPZStack<long> fConstrConnect;
	/** @brief Pointer to the constrained matrix object*/
	//TPZFNMatrix<1000> fConstrMat;
	TPZFNMatrix<1000, STATE> fConstrMat;
	/** @brief Block structure associated with fConstrMat*/
	//TPZBlock<REAL> fConstrBlock;
	TPZBlock<STATE> fConstrBlock;
	
	TPZManVector<long> fDestinationIndex, fSourceIndex;
	
	int fNumStateVars;
	
	/// Reset the data structure
	void Reset(TPZCompMesh *mesh = NULL, MType type=Unknown)
	{
      fMesh = mesh;
      fType = type;
      fConnect.Resize(0);
      fMat.Resize(0,0);
      fBlock.SetNBlocks(0);
      fConstrConnect.Resize(0);
      fConstrMat.Resize(0,0);
      fConstrBlock.SetNBlocks(0);
	}
	
	TPZElementMatrix(TPZCompMesh *mesh, MType type) : fType(type), fMesh(mesh), fConnect(), fMat(0,0), fBlock(&fMat),  fConstrConnect(), fConstrMat(0,0), fConstrBlock(&fConstrMat), fNumStateVars(0)
    {
    }

    TPZElementMatrix() : fType(Unknown), fMesh(NULL), fConnect(), fMat(0,0), fBlock(&fMat), fConstrConnect(), 
      fConstrMat(0,0), fConstrBlock(&fConstrMat), fNumStateVars(0)
    {}
	
    TPZElementMatrix(const TPZElementMatrix &copy);
	
	~TPZElementMatrix(){
	}
	
	/** @brief Returns the number of nodes of TElementMatrix*/
	int NConnects(){
		return fConnect.NElements();
	}
	
	/** @brief Returns the pointer to the ith node of the element*/
	long ConnectIndex(int i)const{
		return fConnect[i];
	}
	
	void Print(std::ostream &out);
	
	void SetMatrixSize(short NumBli, short NumBlj, short BlSizei, short BlSizej);
	
	void SetMatrixMinSize(short NumBli, short NumBlj, short BlMinSizei, short BlMinSizej);
	
	void ComputeDestinationIndices();
    
    /** @brief permute the order of the connects */
    void PermuteGather(TPZVec<long> &permute);
	
	
	/** @brief Returns true if the element has at least one dependent node. Returns false otherwise */
	bool HasDependency();
	
	/** @brief Apply the constraints applied to the nodes by transforming the tangent matrix and right hand side */
	void ApplyConstraints();
	
	void ComputeDestinationIndices ( int iel, TPZManVector<long> &SourceIndex, TPZManVector<long> &DestinationIndex );
};

#endif
