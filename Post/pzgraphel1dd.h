/**
 * @file
 * @brief Contains the TPZGraphEl1dd class which implements the graphical one dimensional discontinuous element.
 */

#ifndef GRAFEL1DDH
#define GRAFEL1DDH


#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphNode;
template<class TVar> 
class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical one dimensional discontinuous element. \ref post "Post processing"
 */
class TPZGraphEl1dd : public TPZGraphEl
{
public:
	/** @brief Constructor for graphical element to computational one dimensional discontinuous element */
	TPZGraphEl1dd(TPZCompEl *ce, TPZGraphMesh *gg);
	
	virtual int NPoints(TPZGraphNode *n);
	
	virtual int NElements();
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle);
	
	void Print(std::ostream &out);
	
	virtual long EqNum(TPZVec<int> &co);
	
	virtual int ExportType(TPZDrawStyle st);
	
	virtual int NNodes();
	
	
protected:
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);
	
	virtual void NextIJ(int no, TPZVec<int> &co, int incr);
	
	/** @brief Graphical node (connect) to discontinuous graphical element */
	TPZGraphNode *fConnect;
	
	virtual void SetNode(long i,TPZGraphNode *gno) {
		fConnect = gno;
	}
	
	int NConnects() { return 1; }
	MElementType Type() {return EOned;}
	
	TPZGraphNode *Connect(long i) {return fConnect;}
	
};

#endif
