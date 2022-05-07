// #include "KLInterpolationSpace.h"
// #include "pzinterpolationspace.h"
// #include "pzmaterialdata.h"
// #include "pzbndcond.h"
// #include "pzelmat.h"
// #include "pzquad.h"
// #include "TPZCompElDisc.h"
// #include "TPZInterfaceEl.h"
// #include "pztransfer.h"
// #include "tpzchangeel.h"
// 
// #include "pzlog.h"
// 
// KLInterpolationSpace::KLInterpolationSpace()
//   : TPZCompEl()
// {
//   fPreferredOrder = -1;
// }
// 
// KLInterpolationSpace::KLInterpolationSpace ( TPZCompMesh &mesh, const KLInterpolationSpace &copy )
//   : TPZCompEl ( mesh, copy )
// {
//   fPreferredOrder = copy.fPreferredOrder;
// }
// 
// KLInterpolationSpace::KLInterpolationSpace ( TPZCompMesh &mesh, const KLInterpolationSpace &copy, std::map<long,long> &gl2lcElMap )
//   : TPZCompEl ( mesh, copy, gl2lcElMap )
// {
//   fPreferredOrder = copy.fPreferredOrder;
// }
// 
// KLInterpolationSpace::KLInterpolationSpace ( TPZCompMesh &mesh, const KLInterpolationSpace &copy, long &index )
//   : TPZCompEl ( mesh, copy, index )
// {
//   fPreferredOrder = copy.fPreferredOrder;
// }
// 
// KLInterpolationSpace::KLInterpolationSpace ( TPZCompMesh &mesh, TPZGeoEl *gel, long &index )
//   : TPZCompEl ( mesh,gel,index )
// {
//   fPreferredOrder = mesh.GetDefaultOrder();
// }
// 
// void KLInterpolationSpace::CalcStiffB ( TPZElementMatrix &be )
// {
// 
// }
// 
// void KLInterpolationSpace::CalcStiffC ( TPZElementMatrix &ce )
// {
//   TPZFMatrix<REAL> elmatt,elmat;
//   TPZMaterial * material = Material();
//   TPZElementMatrix fe,cet;
//   this->InitializeElementMatrix ( ce,fe );
//   this->InitializeElementMatrix ( cet,fe );
//   if ( this->NConnects() == 0 ) return; //boundary discontinuous elements have this characteristic
// 
//   TPZMaterialData data1,data2;
//   this->InitMaterialData ( data1 );
//   this->InitMaterialData ( data2 );
//   data1.p = this->MaxOrder();
//   data2.p = this->MaxOrder();
// 
//   int dim = Dimension();
//   TPZManVector<REAL,3> intpoint ( dim,0. );
//   REAL weight1 = 0.,weight2=0.;
// 
//   TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
//   int order = material->IntegrationRuleOrder ( data1.p );
// 
//   int intrulepoints = intrule->NPoints();
//   for ( int int_ind = 0; int_ind < intrulepoints; ++int_ind )
//     {
//       intrule->Point ( int_ind,intpoint,weight1 );
//       data1.intLocPtIndex = int_ind;
//       this->ComputeRequiredData ( data1, intpoint );
//       weight1 *= fabs ( data1.detjac );
//       for ( int int_jnd=0; int_jnd<intrulepoints; ++int_jnd )
//         {
//           intrule->Point ( int_jnd,intpoint,weight2 );
//           data2.intLocPtIndex = int_jnd;
//           this->ComputeRequiredData ( data2, intpoint );
//           weight2 *= fabs ( data2.detjac );
//           material->ContributeC ( data1,data2, weight1,weight2, cet.fMat );
//           ce.fMat+=cet.fMat;
//         }
// 
//     }//loop over integratin points
// }
// 
// void KLInterpolationSpace::InitializeElementMatrix ( TPZElementMatrix &ek, TPZElementMatrix &ef )
// {
//   TPZMaterial *mat = this->Material();
//   const int numdof = mat->NStateVariables();
//   const int ncon = this->NConnects();
//   const int numloadcases = mat->NumLoadCases();
// 
//   ek.fMesh = Mesh();
//   ek.fType = TPZElementMatrix::EK;
//   ef.fMesh = Mesh();
//   ef.fType = TPZElementMatrix::EF;
// 
//   ek.fBlock.SetNBlocks ( ncon );
//   ef.fBlock.SetNBlocks ( ncon );
//   ek.fNumStateVars = numdof;
//   ef.fNumStateVars = numdof;
//   int i;
//   int numeq=0;
//   for ( i=0; i<ncon; i++ )
//     {
//       int nshape = NConnectShapeF ( i );
//       TPZConnect &c = Connect ( i );
//       int nstate = c.NState();
// 
// #ifdef DEBUG
//       int cNShape = c.NShape();
//       if ( cNShape != nshape || nstate != numdof )
//         {
//           DebugStop();
//         }
// #endif
//       ek.fBlock.Set ( i,nshape*nstate );
//       ef.fBlock.Set ( i,nshape*nstate );
//       numeq += nshape*nstate;
//     }
//   ek.fMat.Redim ( numeq,numeq );
//   ef.fMat.Redim ( numeq,numloadcases );
//   ek.fConnect.Resize ( ncon );
//   ef.fConnect.Resize ( ncon );
//   for ( i=0; i<ncon; i++ )
//     {
//       ( ef.fConnect ) [i] = ConnectIndex ( i );
//       ( ek.fConnect ) [i] = ConnectIndex ( i );
//     }
// }//void
// 
// void KLInterpolationSpace::InitializeElementMatrix ( TPZElementMatrix &ef )
// {
//   TPZMaterial *mat = this->Material();
//   const int numdof = mat->NStateVariables();
//   const int ncon = this->NConnects();
//   const int nshape = this->NShapeF();
//   const int numeq = nshape*numdof;
//   const int numloadcases = mat->NumLoadCases();
//   ef.fMesh = Mesh();
//   ef.fType = TPZElementMatrix::EF;
//   ef.fMat.Redim ( numeq,numloadcases );
//   ef.fBlock.SetNBlocks ( ncon );
//   ef.fNumStateVars = numdof;
//   int i;
//   for ( i=0; i<ncon; i++ )
//     {
//       unsigned int nshapec = NConnectShapeF ( i );
// #ifdef DEBUG
//       TPZConnect &c = Connect ( i );
//       if ( c.NShape() != nshapec || c.NState() != numdof )
//         {
//           DebugStop();
//         }
// #endif
//       ef.fBlock.Set ( i,nshapec*numdof );
//     }
//   ef.fConnect.Resize ( ncon );
//   for ( i=0; i<ncon; i++ )
//     {
//       ( ef.fConnect ) [i] = ConnectIndex ( i );
//     }
// }//void
// 
// int KLInterpolationSpace::MaxOrder()
// {
//   const int n = this->NConnects();
//   int result = -1;
//   int side;
//   for ( int i = 0; i < n; i++ )
//     {
//       side = this->Connect ( i ).Order();
//       if ( side > result ) result = side;
//     }//i
//   return result;
// }
// 
// void KLInterpolationSpace::InitMaterialData ( TPZMaterialData &data )
// {
// 
//   this->Material()->FillDataRequirements ( data );
//   const int dim = this->Dimension();
//   const int nshape = this->NShapeF();
//   const int nstate = this->Material()->NStateVariables();
//   data.phi.Redim ( nshape,1 );
//   data.dphix.Redim ( dim,nshape );
//   data.axes.Redim ( dim,3 );
//   data.jacobian.Redim ( dim,dim );
//   data.jacinv.Redim ( dim,dim );
//   data.x.Resize ( 3 );
//   if ( data.fNeedsSol )
//     {
//       long nsol = data.sol.size();
//       for ( long is=0; is<nsol; is++ )
//         {
//           data.sol[is].Resize ( nstate );
//           data.dsol[is].Redim ( dim,nstate );
//         }
//     }
// }//void
// 
// void KLInterpolationSpace::ComputeRequiredData ( TPZMaterialData &data,
//     TPZVec<REAL> &qsi )
// {
//   data.intGlobPtIndex = -1;
//   this->ComputeShape ( qsi, data );
// 
//   if ( data.fNeedsSol )
//     {
//       if ( data.phi.Rows() ) //if shape functions are available
//         {
//           this->ComputeSolution ( qsi, data );
//         }
//       else //if shape functions are not available
//         {
//           this->ComputeSolution ( qsi, data.sol, data.dsol, data.axes );
//         }
//     }//fNeedsSol
// 
//   data.x.Resize ( 3., 0. );
//   Reference()->X ( qsi, data.x );
// 
//   if ( data.fNeedsHSize )
//     {
//       data.HSize = 2.*this->InnerRadius();
//     }//fNeedHSize
// 
//   if ( data.fNeedsNormal )
//     {
//       this->ComputeNormal ( data );
//     }//fNeedsNormal
// }//void
// void KLInterpolationSpace::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
//                                          TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
//                                          REAL &detjac, TPZFMatrix<REAL> &jacinv,
//                                          TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix){
// 	TPZGeoEl * ref = this->Reference();
// 	if (!ref){
// 		PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
// 		return;
// 	}//if
// 	TPZFNMatrix<660> dphi(dphix.Rows(), dphix.Cols(), 0.);
//     //	int dim = this->Dimension();
// 	
// 	ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
// 	this->Shape(intpoint,phi,dphi);
//     
//     Convert2Axes(dphi, jacinv, dphix);
//     
// }
// void KLInterpolationSpace::ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data){
//     
//     this->ComputeShape(intpoint,data.x,data.jacobian,data.axes,data.detjac,data.jacinv,data.phi,data.dphix);
//     
// }
// REAL KLInterpolationSpace::InnerRadius(){
// 	if (!this->Reference()){
// 		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Reference() == NULL\n";
// 		return 0.;
// 	}
// 	return this->Reference()->ElementRadius();
// }
// 
// void KLInterpolationSpace::ComputeNormal(TPZMaterialData & data)
// {
// 	data.normal.Resize(3,0.);
// 	
// 	int thisFace, neighbourFace, i, dim;
// 	TPZGeoEl * thisGeoEl, * neighbourGeoEl;
// 	TPZManVector<REAL,3> thisCenter(3,0.), neighbourCenter(3,0.), thisXVol(3,0.), neighbourXVol(3,0.), vec(3), axes1(3), axes2(3);
// 	
// 	thisGeoEl = this->Reference();
// 	thisFace = thisGeoEl->NSides() - 1;
//     TPZGeoElSide thisside(thisGeoEl,thisFace);
//     TPZCompMesh *cmesh = this->Mesh();
// 	
// 	TPZGeoElSide neighbourGeoElSide = thisGeoEl->Neighbour(thisFace);
//     int matid = neighbourGeoElSide.Element()->MaterialId();
//     while (!cmesh->FindMaterial(matid) && neighbourGeoElSide != thisside) {
//         neighbourGeoElSide = neighbourGeoElSide.Neighbour();
//         matid = neighbourGeoElSide.Element()->MaterialId();
//     }
// 	neighbourGeoEl = neighbourGeoElSide.Element();
// 	neighbourFace = neighbourGeoEl->NSides() - 1;
//     
//     
// 	if(neighbourGeoEl == thisGeoEl)
// 	{
// 		// normal evaluation makes no sense since the internal element side doesn't present a neighbour.
// 		return; // place a breakpoint here if this is an issue
// 	}
//     
// 	thisGeoEl->CenterPoint(thisFace, thisCenter);
// 	neighbourGeoEl->CenterPoint(neighbourFace, neighbourCenter);
//     
// 	
// 	thisGeoEl->X(thisCenter,thisXVol);
// 	neighbourGeoEl->X(neighbourCenter,neighbourXVol);
// 	
// 	for(i = 0; i < 3; i++)
// 		vec[i] = -neighbourXVol[i] + thisXVol[i];// vector towards the center of the neighbour element
// 	
// 	dim = thisGeoEl->Dimension();
// 	
// 	switch(dim)
// 	{
// 		case(0): // normal points towards the x-direction
// 			data.normal[0] = 1.;
// 			data.normal[1] = 0.;
// 			data.normal[2] = 0.;
// 			break;
// 		case(1):
// 			for(i = 0 ; i < 3; i ++) axes1[i] = data.axes(0,i); // rib direction
// 			this->VectorialProd(axes1, vec, axes2);
// 			this->VectorialProd(axes2, axes1, data.normal, true);
// 			break;
// 		case(2):
// 			for(i = 0; i < 3; i++)
// 			{
// 				axes1[i] = data.axes(0,i);
// 				axes2[i] = data.axes(1,i);
// 			}
// 			this->VectorialProd(axes1, axes2, data.normal, true);
// 			break;
// 		case(3):// in this case the normal becomes senseless. A null vector is passed instead
// 			break;
// 		default:
// 			PZError << "TPZInterpolationSpace::ComputeNormal - unhandled element dimension\n";
// 	}
// 	
// 	// ensuring the normal vector points towards the neighbour element
// 	
// 	REAL dot = 0.;
// 	for(i = 0; i < 3; i++) dot += data.normal[i] * vec[i];
// 	
// 	if(dot < 0.)
// 		for(i = 0; i < 3; i++) data.normal[i] *= -1.;
// 	
// }
// 
// /// convert a shapefunction derivative in xi-eta to a function derivative in x,y,z
// void KLInterpolationSpace::Convert2Axes(const TPZFMatrix<REAL> &dphidxi, const TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &dphidaxes)
// {
// 	int nshape = dphidxi.Cols();
//     int dim = dphidxi.Rows();
//     dphidaxes.Resize(dim,nshape);
// 	int ieq;
// 	switch(dim){
// 		case 0:
// 			break;
// 		case 1:
// 			dphidaxes = dphidxi;
// 			dphidaxes *= jacinv.GetVal(0,0);
// 			break;
// 		case 2:
// 			for(ieq = 0; ieq < nshape; ieq++) {
// 				dphidaxes(0,ieq) = jacinv.GetVal(0,0)*dphidxi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphidxi.GetVal(1,ieq);
// 				dphidaxes(1,ieq) = jacinv.GetVal(0,1)*dphidxi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphidxi.GetVal(1,ieq);
// 			}
// 			break;
// 		case 3:
// 			for(ieq = 0; ieq < nshape; ieq++) {
// 				dphidaxes(0,ieq) = jacinv.GetVal(0,0)*dphidxi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphidxi.GetVal(1,ieq) + jacinv.GetVal(2,0)*dphidxi.GetVal(2,ieq);
// 				dphidaxes(1,ieq) = jacinv.GetVal(0,1)*dphidxi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphidxi.GetVal(1,ieq) + jacinv.GetVal(2,1)*dphidxi.GetVal(2,ieq);
// 				dphidaxes(2,ieq) = jacinv.GetVal(0,2)*dphidxi.GetVal(0,ieq) + jacinv.GetVal(1,2)*dphidxi.GetVal(1,ieq) + jacinv.GetVal(2,2)*dphidxi.GetVal(2,ieq);
// 			}
// 			break;
// 		default:
// 			PZError << "Error at " << __PRETTY_FUNCTION__ << " please implement the " << dim << "d Jacobian and inverse\n";
// 	} //switch
//     
// }
// void KLInterpolationSpace::VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary)
// {
// 	kvec.Resize(3);
// 	kvec[0] =  ivec[1]*jvec[2] - ivec[2]*jvec[1];
// 	kvec[1] = -ivec[0]*jvec[2] + ivec[2]*jvec[0];
// 	kvec[2] =  ivec[0]*jvec[1] - ivec[1]*jvec[0];
// 	
// 	if(unitary)
// 	{
// 		REAL size = 0.;
// 		int i;
// 		for(i = 0; i < 3; i++)size += kvec[i] * kvec[i];
// 		size = sqrt(size);
// 		//if(size <= 1.e-9)PZError << "\nTPZInterpolationSpace::VectorialProd - null result\n";
// 		for(i = 0; i < 3; i++)kvec[i] /= size;
// 	}
// }
