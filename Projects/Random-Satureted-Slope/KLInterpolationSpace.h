// /*
// 
// #ifndef KLInterpolationSpace_H
// #define KLInterpolationSpace_H
// 
// 
// #include "pzcompel.h"
// #include "pzinterpolationspace.h"
// #include "pzmaterial.h"
// class TPZMaterialData;
// class TPZCompEl;
// class TPZMaterial;
// /**
//  * @brief Implements the interfaces for TPZCompElDisc, TPZInterfaceElement and TPZInterpolatedElement. \ref CompElement "Computational element"
//  * @since April 11, 2007
//  * @ingroup CompElement
//  */
// class KLInterpolationSpace :  TPZCompEl
// {
// public:
// 	
// 	/** @brief Default constructor */
// 	KLInterpolationSpace();
// 	
// 	/** @brief Default destructor */
// 	virtual ~KLInterpolationSpace();
// 	
// 	/** @brief Puts a copy of the element in the referred mesh */
// 	KLInterpolationSpace(TPZCompMesh &mesh, const KLInterpolationSpace &copy);
// 	
// 	/** @brief Puts a copy of the element in the patch mesh */
// 	KLInterpolationSpace(TPZCompMesh &mesh, const KLInterpolationSpace &copy, std::map<long,long> &gl2lcElMap);
// 	
// 	/** @brief Copy of the element in the new mesh whit alocated index */
// 	KLInterpolationSpace(TPZCompMesh &mesh, const KLInterpolationSpace &copy, long &index);
// 	
//    KLInterpolationSpace ( TPZCompMesh &mesh, TPZGeoEl *gel, long &index );
// 	
// 	void CalcStiffB(TPZElementMatrix &be);
// 	
// 	void CalcStiffC(TPZElementMatrix &ce);
// 	
// 	void InitializeElementMatrix(TPZElementMatrix &ef);
// 	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
// 	int MaxOrder();
// 	virtual const TPZIntPoints &GetIntegrationRule() const = 0;
// 	
// 	virtual int NShapeF() const = 0;
// 	void InitMaterialData(TPZMaterialData &data);
// 	void ComputeRequiredData(TPZMaterialData &data,TPZVec<REAL> &qsi);
// 	virtual int NConnectShapeF(int icon) const = 0;
// 	void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data);
// 	void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
//                                          TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
//                                          REAL &detjac, TPZFMatrix<REAL> &jacinv,
//                                          TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix);
// 	void ComputeNormal(TPZMaterialData & data);
// 	virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi) = 0;
// 	REAL InnerRadius();
// 	static void Convert2Axes(const TPZFMatrix<REAL> &dphidxi, const TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &dphidaxes);
// 		void VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary = false);
// 	protected:
// 	
// 	int fPreferredOrder;
// 	
// 	
// };
// 
// #endif*/
