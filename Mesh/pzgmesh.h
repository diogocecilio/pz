/**
 * @file
 * @brief Contains declaration of TPZMesh class which defines a geometrical mesh and contains a corresponding lists of elements, nodes and conditions.
 */

#ifndef PZGEOMESHH
#define PZGEOMESHH


#include "pzsave.h"
#include "pzreal.h"
#include "pzeltype.h"
#include "pzgnode.h"
#include "pzadmchunk.h"
#include "tpzautopointer.h"

#include <iostream>
#include <string>
#include <map>
#include <list>

class TPZMaterial;
class TPZGeoNode;
struct TPZGeoNodeBC;
class TPZGeoEl;
class TPZCosys;
template<class TVar>
class TPZMatrix;
class TPZCompMesh;
class TPZRefPattern;
class TPZRefPatternDataBase;

template<class T>
class TPZVec;

template <class TGeo> class TPZGeoElRefPattern;

/**
 * @brief This class implements a geometric mesh for the pz environment. \ref geometry "Geometry"
 * @ingroup geometry
 */
/**
 * A TPZGeoMesh defines a geometrical mesh and contains a corresponding list of geometrical elements, geometrical nodes,
 * elementwise defined boundary conditions and nodal boundary conditions. \n
 * Various methods are defined to add, delete or loop over the items which are contained within the TPZGeoMesh. \n
 * Other auxiliary data structures help in the construction of the mesh
 */

class  TPZGeoMesh : public TPZSaveable {
	
protected:
	/** @brief TPZGeoMesh name for model identification */
	std::string fName;
	
	/** @brief Computational mesh associated */
	TPZCompMesh 	*fReference;
	
	/** @brief List of pointers to finite elements */
	TPZAdmChunkVector<TPZGeoEl *> fElementVec;
	
	/** @brief List of nodes */
	TPZAdmChunkVector<TPZGeoNode> fNodeVec;
	
	/** @brief Maximum id used by all nodes of this mesh */
	long fNodeMaxId;
	
	/** @brief Maximum id used by all elements of this mesh */
	long fElementMaxId;
    
    /** @brief dimension of the geometric domain */
    int fDim;
	
	typedef std::map<std::pair<int,int>, int> InterfaceMaterialsMap;
	
	/**
	 * @brief Datastructure which indicates the index of the interface material which needs to be created between two different materials 
	 * @see AddInterfaceMaterial 
	 */
	InterfaceMaterialsMap fInterfaceMaterials;
	
public:
	/** @brief Constructors and destructor */
	TPZGeoMesh();
	
	/** @brief Copy constructor */
	TPZGeoMesh(const TPZGeoMesh &cp);
	
	/** @brief Operator of copy */
	TPZGeoMesh & operator= (const TPZGeoMesh &cp );
	
	/** @brief Destructor */
	virtual ~TPZGeoMesh();
	
	/** @brief Deletes all items in the TPZGeoMesh */
	void CleanUp();
	
	/** @brief Reset all connectivities */
	void ResetConnectivities();
	
	virtual int ClassId() const;
	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Indicates that a node with id was created */
	void SetNodeIdUsed(long id) { fNodeMaxId = (id > fNodeMaxId) ? id : fNodeMaxId; }
	
	/** @brief Indicates that an element with id was created */
	void SetElementIdUsed(long id) { fElementMaxId = (id > fElementMaxId) ? id : fElementMaxId; }
	
	/** @brief Returns ++fNodeMaxId */
	long CreateUniqueNodeId() { return ++fNodeMaxId; }
	
	/** @brief Returns ++fElementMaxId */
	long CreateUniqueElementId() { return ++fElementMaxId; }
	
	/** @brief Used in patch meshes */
	void SetMaxNodeId(long id) { fNodeMaxId = (id > fNodeMaxId) ? id : fNodeMaxId; }
	
	/** @brief Used in patch meshes */
	void SetMaxElementId(long id) { fElementMaxId = (id > fElementMaxId) ? id : fElementMaxId; }
	
	/** @brief Number of nodes of the mesh */
	long NNodes() const {return fNodeVec.NElements();}
	
	/** @brief Number of elements of the mesh */
	long NElements() const {return fElementVec.NElements();}
	
	long ReallyNEl() {return (fElementVec.NElements() - fElementVec.NFreeElements()) ; }
	
	void SetName(const std::string &nm);
	
	std::string &Name() { return fName; }
	
	/** @brief Methods for handling pzlists */
	TPZAdmChunkVector<TPZGeoEl *> &ElementVec() { return fElementVec; }
    TPZGeoEl * Element(long iel) { return fElementVec[iel]; }
	TPZAdmChunkVector<TPZGeoNode> &NodeVec() { return fNodeVec; }
	const TPZAdmChunkVector<TPZGeoEl *> &ElementVec() const { return fElementVec; }
	const TPZAdmChunkVector<TPZGeoNode> &NodeVec() const { return fNodeVec; }

    /// Compute the area of the domain
    REAL Area();
    
	/** @brief Resets all load references in elements and nodes */
	void ResetReference();
	
	/** @brief Restore all reference in elements from computational mesh criated from current
     geometrical mesh previously */
	void RestoreReference(TPZCompMesh *cmesh);
	
	/** @brief Sets the reference of the geometric grid to ref */
	void SetReference(TPZCompMesh *ref)
	{
		fReference = ref;
	}
	
	/** @brief Returns the currently loaded computational grid */
	TPZCompMesh *Reference() const {return fReference;}
	
	/** @brief Print the information of the grid to an ostream */
	virtual void Print(std::ostream & out = std::cout);
	virtual void PrintTopologicalInfo(std::ostream & out = std::cout);
	
    /** @brief Returns the nearest node to the coordinate. This method is VERY INEFFICIENT */
	TPZGeoNode* FindNode(TPZVec<REAL> &co);
    
    /** by Caju 2012 */
    /** @brief Returns the element that contains the given point x and it respective point in parametric domain qsi */
    TPZGeoEl * FindElement(TPZVec<REAL> &x, TPZVec<REAL> & qsi, long & InitialElIndex, int targetDim);
    
    /** @brief find an element/parameter close to the point */
    TPZGeoEl *FindApproxElement(TPZVec<REAL> &x, TPZVec<REAL> & qsi, long & InitialElIndex, int targetDim);
    
    /** by Caju 2013 */
    /** @brief Returns the subelement that contains the given point x and it respective point in parametric domain qsi */
    TPZGeoEl * FindSubElement(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> & qsi, long & InitialElIndex);
	
    /** by Philippe 2013 */
    /** @brief Returns the element that is close to the given point x */
    TPZGeoEl * FindCloseElement(TPZVec<REAL> &x, long & InitialElIndex, int targetDim) const;
	
	/** @brief Alternative method for computing the connectivity */
	void BuildConnectivityOld();
	
	/** @brief Build the connectivity of the grid */
	void BuildConnectivity();
	
	/** @brief Fills the nodep vector with pointers to the nodes identified by their indexes */
	void GetNodePtr(TPZVec<long> &nos,TPZVec<TPZGeoNode *> &nodep);
	
	/**
	 * @brief GetBoundaryElements returns all elements beweeen NodFrom and NodTo counterclock wise this
	 * method uses the connectivity of the elements BuildConnectivity should be called to initialize
     */
	/**
	 * The connectivity information this method will only work for grid with 2-D topology the current 
	 * version will only work for a grid with only one level
	 */
	void GetBoundaryElements(long IndexNodeFrom,long IndexNodeTo,TPZStack<TPZGeoEl *> &ElementVec,TPZStack<int> &Sides);
	
    /** @brief Find the element with ids elid */
	TPZGeoEl * FindElement(long elid);
	
	/** @brief Returns the index of the given element into the fElementVec */
	/** @since 2002-05-02 (Cesar) */
	long ElementIndex(TPZGeoEl * gel);
	
	/** @brief Returns the index of the given node into the fNodeVec */
	/** @since 2002-05-02 (Cesar) */
	long NodeIndex(TPZGeoNode * nod);
	
	/**
	 * @brief Generic method for creating a geometric element. Putting this method centrally facilitates
	 * the modification of the element type all through the code
	 * @param type element topology
	 * @param matid material id
	 * @param cornerindexes indexes of the corner nodes of the element
	 * @param index index of the element in the vector of element pointers
	 * @param reftype defines the type of refinement : 0 -> uniform 1-> refinement pattern
	 */
	virtual  TPZGeoEl *CreateGeoElement(MElementType type,TPZVec<long> &cornerindexes,int matid,long &index, int reftype = 1);
	
	/** @brief Creates a geometric element in same fashion of CreateGeoElement but here the elements are blend, as Caju master thesis */
	virtual TPZGeoEl *CreateGeoBlendElement(MElementType type, TPZVec<long>& nodeindexes, int matid, long& index);
	
	/**
	 * @brief Centralized method to delete elements
	 * @param gel pointer to the element to be deleted
	 * @param index index of the element
	 */
	void DeleteElement(TPZGeoEl *gel, long index = -1);
	
	/**
	 * @brief Add an interface material associated to left and right element materials.
	 * @since Feb 05, 2004
	 */
	/**
	 * If std::pair<left, right> already exist, nothing is done and method returns 0.
	 * If material is inserted in geomesh method returns 1.
	 */
	int AddInterfaceMaterial(int leftmaterial, int rightmaterial, int interfacematerial);
	
	/**
	 * @brief Returns the interface material associated to left and right element materials.
     * If no interface material is found GMESHNOMATERIAL is returned
	 * @since Feb 05, 2004
	 */
	int InterfaceMaterial(int leftmaterial, int rightmaterial);
	
	/**
	 * @brief Delete all interface materials in map.
	 * @since Feb 05, 2004
	 */
	void ClearInterfaceMaterialsMap();
    
    /**
	 * @brief Set Dimension.
	 * @since April 17, 2004
	 */
    void SetDimension(int dim){fDim = dim;}
    
    /**
	 * @brief Get Dimension.
	 * @since April 17, 2004
	 */
    int Dimension(){return fDim;}
	
private:
	
	/** @brief Find all elements in elmap or neighbour of elements in elmap which contain a node */
	void BuildElementsAroundNode(long currentnode,std::map<long,TPZGeoEl *> &elmap);
	
	/** @brief Method which works only for two dimensional topologies! */
	/** Find, within elmap, the element which has currentnode as its first boundary side node */
 	void FindElement(std::map<long,TPZGeoEl *> &elmap,long currentnode,TPZGeoEl* &candidate,int &candidateside);
};

#endif
