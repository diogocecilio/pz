//
//  pzhdivpressurebound.cpp
//  PZ
//
//  Created by Agnaldo Farias on 25/06/13.
//
//

#include "pzhdivpressurebound.h"
#include "pzelchdivbound2.h"
#include "pzlog.h"
#include "pzmaterial.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzhdivpressure.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivPressureBound"));
#endif


template <class TSHAPE>
TPZCompElHDivPressureBound<TSHAPE>::TPZCompElHDivPressureBound(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) :
TPZCompElHDivBound2<TSHAPE>(mesh, gel, index){
    
    
    //Creating connect of the pressure's variable
	this->fConnectIndexes.Resize(NConnects());
    this->fConnectIndexes[this->NConnects()-1] = this->CreateMidSideConnect(2);
    
    int nshape = 0;
    if(gel->Dimension()!= mesh.Dimension()-1) DebugStop();
    nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, gel->Dimension(), pzshape::TPZShapeDisc::ETensorial);
    
    if(gel->Dimension()==2 && TSHAPE::Type()==EQuadrilateral){
        nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, gel->Dimension(), pzshape::TPZShapeDisc::  ETensorial);
    }
    if (gel->Dimension()==2 && TSHAPE::Type()==ETriangle){
        nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, gel->Dimension(), pzshape::TPZShapeDisc::  EOrdemTotal);
    }
    
    int nstate = 1;
    
    long newnodeindex = mesh.AllocateNewConnect(nshape,nstate,this->fPressureOrder);
    TPZConnect &newnod = mesh.ConnectVec()[newnodeindex];
    newnod.SetLagrangeMultiplier(1);
    this->fConnectIndexes[this->NConnects()-1] = this->CreateMidSideConnect(2);
    long seqnum = newnod.SequenceNumber();
    newnod.SetLagrangeMultiplier(1);
    mesh.Block().Set(seqnum,nshape);
    mesh.ConnectVec()[this->fConnectIndexes[this->NConnects()-1]].IncrementElConnected();
}

template<class TSHAPE>
TPZCompElHDivPressureBound<TSHAPE>::TPZCompElHDivPressureBound():TPZCompElHDivBound2<TSHAPE>()
{
    fPressureOrder = 0;

}

template<class TSHAPE>
TPZCompElHDivPressureBound<TSHAPE>::TPZCompElHDivPressureBound(TPZCompMesh &mesh, const TPZCompElHDivPressureBound<TSHAPE> &copy) :
TPZCompElHDivBound2<TSHAPE>(mesh,copy)
{
    fPressureOrder = copy.fPressureOrder;
}


template<class TSHAPE>
TPZCompElHDivPressureBound<TSHAPE>::TPZCompElHDivPressureBound(TPZCompMesh &mesh, const TPZCompElHDivPressureBound<TSHAPE> &copy, std::map<long,long> & gl2lcConMap, std::map<long,long> & gl2lcElMap) :
TPZCompElHDivBound2<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap)
{
    fPressureOrder = copy.fPressureOrder;
}

template<class TSHAPE>
TPZCompElHDivPressureBound<TSHAPE>::~TPZCompElHDivPressureBound(){
	
}

template<class TSHAPE>
int TPZCompElHDivPressureBound<TSHAPE>::NConnects() const
{
	
	return TPZCompElHDivBound2<TSHAPE>::NConnects()+1;
}

template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::SetPressureOrder(int order)
{
    fPressureOrder = order;
}

template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::SetConnectIndex(int i, long connectindex)
{
	TPZCompElHDivBound2<TSHAPE>::SetConnectIndex(i, connectindex);
	
}

template<class TSHAPE>
int TPZCompElHDivPressureBound<TSHAPE>::NConnectShapeF(int connect) const
{
    
    if(connect<NConnects()-2){
        
        return TPZCompElHDivBound2<TSHAPE>::NConnectShapeF(connect);
    }
    else{
        
        int nshape = 0;
        if(TPZCompElHDivBound2<TSHAPE>::Dimension() < 1) DebugStop();
        
        if(TPZCompElHDivBound2<TSHAPE>::Dimension() == 1){
            nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::ETensorial);
        }
        
        if(TPZCompElHDivBound2<TSHAPE>::Dimension()==2 && TSHAPE::Type()==EQuadrilateral){
            nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::  ETensorial);
        }
        
        if (TPZCompElHDivBound2<TSHAPE>::Dimension()==2 && TSHAPE::Type()==ETriangle){
            nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::  EOrdemTotal);
        }

        return (nshape);
    }
}


template<class TSHAPE>
int TPZCompElHDivPressureBound<TSHAPE>::SideConnectLocId(int node, int side) const
{
	if(side == TSHAPE::NSides-1 && node == 1)
	{
		return 1;
	}
	else{
        int conectId = TPZCompElHDivBound2<TSHAPE>::SideConnectLocId(node, side);
        return(conectId);
    }
}


template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord)
{
	ord.Resize(NConnects());
	int i;
	for(i=0; i<NConnects(); i++) {
		ord[i] = ConnectOrder(i);
	}
}

/** Returns the actual interpolation order of the polynomial along the side */
template<class TSHAPE>
int TPZCompElHDivPressureBound<TSHAPE>::SideOrder(int side) const
{
	DebugStop();
    std::cout<<"\nNao implementado";
    return 0;
}

template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::SetSideOrder(int side, int order)
{
	int connectaux= SideConnectLocId(0,side);
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElHDivPressureBound::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	}
	TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order);
    long seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial * mat =this-> Material();
    if(mat) nvar = mat->NStateVariables();
    int nshape = NConnectShapeF(connectaux);
    c.SetNShape(nshape);
    c.SetNState(nvar);
    this-> Mesh()->Block().Set(seqnum,nshape*nvar);
    if(connectaux == NConnects()-1)
    {
		this->SetIntegrationRule(2*order);
    }
}


/**
 return the interpolation orderof the polynomial for connect
 **/
template<class TSHAPE>
int TPZCompElHDivPressureBound<TSHAPE>::ConnectOrder(int connect) const
{
    if(connect<NConnects()-2)
    {
        int connord =  TPZCompElHDivBound2<TSHAPE>::ConnectOrder(connect);
        return (connord);
    }
    else{
        
        const TPZConnect &c = this-> Connect(connect);
        return c.Order();
    }
}

template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{

    data.fShapeType = TPZMaterialData::EVecandShape;
    TPZCompElHDivBound2<TSHAPE>::InitMaterialData(data);
    data.numberdualfunctions = NConnectShapeF(NConnects()-1);
}

template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<long> &shapeindex)
{

    TPZManVector<long> firstshapeindex;
	FirstShapeIndex(firstshapeindex);
	int nshape = TPZIntelGen<TSHAPE>::NShapeF();
	shapeindex.Resize(nshape);
	long nsides = sides.NElements();
	long is, count=0;
	for(is=0 ; is<nsides; is++)
	{
		int side = sides[is];
		int sideorder= this->SideOrder(side);
		int NShapeFace = TSHAPE::NConnectShapeF(side,sideorder);
		int ishapeface;
		for(ishapeface=0; ishapeface<NShapeFace; ishapeface++)
		{
			shapeindex[count++] = is;
		}
	}
	shapeindex.Resize(count);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "count = " << count << " nshape " << nshape;
		sout << std::endl<<"sides associated with the normals "<< sides <<
		"\nnormal associated with each shape function : shape function indexes " << shapeindex;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
}

template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::FirstShapeIndex(TPZVec<long> &Index)
{
    DebugStop();
    std::cout << "\nNao implementado";
}

template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
{
    DebugStop();
    std::cout << "\nNao implementado";
}

/** Compute the shape function at the integration point */
template<class TSHAPE>
void TPZCompElHDivPressureBound <TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
    DebugStop();
    std::cout << "\nNao implementado";
}


/** Read the element data from a stream */
template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZCompElHDivBound2<TSHAPE>::Read(buf,context);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElHDivPressureBound<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
	TPZCompElHDivBound2<TSHAPE>::Write(buf,withclassid);
}


//---------------------
#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;
/** returns the unique identifier for reading/writing objects to streams */
template<>
int TPZCompElHDivPressureBound<TPZShapePoint>::ClassId() const
{
	return TPZHDIVBOUND2POINTID;
}

template class
TPZRestoreClass< TPZCompElHDivPressureBound<TPZShapePoint>, TPZHDIVPRESSUREBOUNDPOINTID>;

template<>
int TPZCompElHDivPressureBound<TPZShapeLinear>::ClassId() const
{
	return TPZHDIVPRESSUREBOUNDLINEARID;
}

template class
TPZRestoreClass< TPZCompElHDivPressureBound<TPZShapeLinear>, TPZHDIVPRESSUREBOUNDLINEARID>;

template<>
int TPZCompElHDivPressureBound<TPZShapeTriang>::ClassId() const
{
	return TPZHDIVPRESSUREBOUNDTRIANGLEID;
}

template class
TPZRestoreClass< TPZCompElHDivPressureBound<TPZShapeTriang>, TPZHDIVPRESSUREBOUNDTRIANGLEID>;

template<>
int TPZCompElHDivPressureBound<TPZShapeQuad>::ClassId() const
{
	return TPZHDIVPRESSUREBOUNDQUADID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZCompElHDivPressureBound<TPZShapeQuad>, TPZHDIVPRESSUREBOUNDQUADID>;
#endif


template class TPZCompElHDivPressureBound<TPZShapeTriang>;
//template class TPZCompElHDivPressureBound<TPZShapePoint>;
template class TPZCompElHDivPressureBound<TPZShapeLinear>;
template class TPZCompElHDivPressureBound<TPZShapeQuad>;
