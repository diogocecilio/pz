/**
 * @file
 * @brief Contains the implementation of the TPZPyramid methods into the pztopology scope.
 */

#include "tpzpyramid.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzpyramid"));
#endif

using namespace std;

namespace pztopology {	

	static int FaceConnectLocId[5][9] = { {0,1,2,3,5,6,7,8,13},{0,1,4,5,10,9,14,-1,-1},
		{1,2,4,6,11,10,15,-1,-1},{3,2,4,7,11,12,16,-1,-1},{0,3,4,8,12,9,17,-1,-1} };
	
	
	static int nhighdimsides[19] = {7,7,7,7,9,3,3,3,3,3,3,3,3,1,1,1,1,1,0};
	
	int TPZPyramid::SideNodes[8][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4} };
	
	int TPZPyramid::FaceNodes[5][4]  = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };
	
	int TPZPyramid::ShapeFaceId[5][4] = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };
	
	static int sidedimension[19] = {0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,3};
	
	static int highsides[19][9] = {
		{5,8,9,13,14,17,18},
		{5,6,10,13,14,15,18},
		{6,7,11,13,15,16,18},
		{7,8,12,13,16,17,18},
		{9,10,11,12,14,15,16,17,18},
		{13,14,18},
		{13,15,18},
		{13,16,18},
		{13,17,18},
		{14,17,18},
		{14,15,18},
		{15,16,18},
		{16,17,18},
		{18},
		{18},
		{18},
		{18},
		{18},
		{-999}
	};
	
	int nsidenodes[19] = {1,1,1,1,1,
		2,2,2,2,2,2,2,2,
		4,3,3,3,3,
		5};
	
	int TPZPyramid::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	static REAL sidetosidetransforms[19][9][4][3] = {
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
		},
		{
			{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{-0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,0}}
		},
		{
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{0.5,0.5,0.5},{-99,-99,-99},{-99,-99,-99},{-0.5,-0.5,0.5}}
		},
		{
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{-0.5,0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,-0.5,0.5}}
		},
		{
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{-0.5,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,0.5}}
		},
		{
			{{0,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{0.5,-0.5,0.5},{-99,-99,-99},{-99,-99,-99},{-0.5,0.5,0.5}}
		},
		{
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,0}}
		},
		{
			{{2,0,0},{1,1,1},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{0,2,0},{-1,1,1},{-99,-99,-99},{1,-1,0}}
		},
		{
			{{2,0,0},{1,-1,1},{-99,-99,-99},{-1,1,0}}
		},
		{
			{{0,2,0},{1,1,1},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	static REAL MidSideNode[19][3] = {
		/*00*/{-1.,-1.},   /*01*/{1.,-1.},   /*02*/{1.,1.},/*03*/{-1.,1.},/*04*/{0.,0.,1.},
		/*05*/{ 0.,-1.},   /*06*/{1., 0.},   /*07*/{0.,1.},/*08*/{-1.,0.},
		/*09*/{-.5,-.5,.5},/*10*/{.5,-.5,.5},/*11*/{.5,.5,.5},/*12*/{-.5,.5,.5},
		/*13*/{0.,  0. ,  0. },/*14*/{  0.  ,-2./3.,1./3.},/*15*/{2./3.,0.,1./3.},
		/*16*/{0.,2./3.,1./3.},/*17*/{-2./3.,  0.  ,1./3.},/*18*/{  0. ,0.,1./5.} };
    
    static REAL bPiram[58][3] =
    {
        {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1}, {0,-1,-1}, {1,0,-1}, {0,1,-1}, {-1,0,-1}, {0,0,-1},// face 0
        {0,-1,0}, {0,-1,0}, {-1,-1,1}, {0,-1,0}, {-1,-3,1}, {1,-3,1}, {0,-1,1},// face 1
        {1,0,0}, {1,0,0}, {1,1,1}, {1,0,0}, {3,-1,1}, {3,1,1}, {1,0,1}, // face 2
        {0,1,0}, {0,1,0}, {1,1,1}, {0,1,0}, {-1,3,1}, {1,3,1}, {0,1,1}, // face 3
        {-1,0,0}, {-1,0,0}, {-1,1,1}, {-1,0,0}, {-3,-1,1}, {-3,1,1}, {-1,0,1},// face 4
        //internos
        //faces
        {-1,0,0}, {0,1,0}, // tang da face 0
        {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)},{sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, //face 1
        {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)},{(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))}, //face 2
        {1/sqrt(3),1/sqrt(3),-1/sqrt(3)},{-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))}, //face 3
        {-1/sqrt(3),1/sqrt(3),-1/sqrt(3)},{-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))}, //face 4
        // arestas
        {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},  {1,1,1},  {-1,1,1},  {-1,-1,1},  {1,-1,1},
        //interior
        {1,0,0} ,
        {0,1,0} ,
        {0,0,1}
    };
    static REAL t1Piram[58][3] =
    {
        {-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},// face 0
        {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, // face 1
        {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, // face 2
        {1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}, // fsce 3
        {-1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),1/sqrt(3),-1/sqrt(3)},{-1/sqrt(3),1/sqrt(3),-1/sqrt(3)},{-1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, {-1/sqrt(3),1/sqrt(3),-1/sqrt(3)},{-1/sqrt(3),1/sqrt(3),-1/sqrt(3)},{-1/sqrt(3),1/sqrt(3),-1/sqrt(3)}, // face 4
        //internos
        //faces
        {0,1,0}, {1,0,0}, //face 0
        {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))},    {1/sqrt(3),1/sqrt(3),1/sqrt(3)}, //face 1
        {(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))},      {-1/sqrt(3),1/sqrt(3),1/sqrt(3)}, //face 2
        {-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))},    {-1/sqrt(3),-1/sqrt(3),1/sqrt(3)}, //face 3
        {-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))},   {1/sqrt(3),-1/sqrt(3),1/sqrt(3)}, //face 4
        // arestas
        {0,-1,0},{1,0,0},{0,1,0},{0,0,-1},  {-1,0,1},{0-1,1},{1,0,1},{0,1,1},
        //interior
        {0,1,0} ,
        {0,0,1} ,
        {1,0,0}
    };
    
    static REAL t2Piram[58][3] =
    {
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0},// face 0
        {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, {sqrt(2/3), -(1/sqrt(6)), -(1/sqrt(6))}, // face 1
        {(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))}, {(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))},  {(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))},{(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))},{(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))},{(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))},{(1/sqrt(6)), sqrt(2/3), -(1/sqrt(6))}, // face 2
        {-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))}, {-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))}, {-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))}, {-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))}, {-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))},{-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))}, {-sqrt(2/3), (1/sqrt(6)), -(1/sqrt(6))}, // face 3
        {-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))},{-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))},{-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))},{-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))},{-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))},{-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))},{-(1/sqrt(6)), -sqrt(2/3), -(1/sqrt(6))}, // face 4
        //internos
        //faces
        {0,0,-1}, {0,0,-1}, //face 0
        {0,-1,0},{0,-1,0}, //face 1
        {1,0,0},{1,0,0}, //face 2
        {0,1,0},{0,1,0}, //face 3
        {-1,0,0},{-1,0,0}, //face 4
        // arestas
        {0,0,-1},{0,0,-1},{0,0,-1},{-1,0,0},  {0,-1,1},{1,0,1},{0,1,1},{-1,0,1},
        //interior
        {0,0,1} ,
        {1,0,0} ,
        {0,1,0}
    };
    
    static int vectorsideorderPi [58] =
    {
        0,1,2,3,5,6,7,8,13, //face 0
        0,1,4,5,10,9,14,//face 1
        1,2,4,6,11,10,15,//face 2
        3,2,4,7,11,12,16,//face 3
        0,3,4,8,12,9,17,//face 4
        13,13,//tg face 0
        14,14,//tg face 1
        15,15,//tg face 2
        16,16,//tg face 3
        17,17,//tg face 3
        5,6,7,
        8,9,10,11,12,
        18,18,18
    };
    
    static int bilinearounao [58] =   {
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,1,1,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0};
    
    static int direcaoksioueta [58] = {
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,1,0,
        1,0,1,0,1,0,1,0,0,0,
        0,0,0,0,0,0,1,2};

    int TPZPyramid:: NBilinearSides()
    {
        DebugStop();
        return 2;
    }
    
	void TPZPyramid::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZPyramid::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZPyramid::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZPyramid::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	int TPZPyramid::SideNodeLocId(int side, int node)
	{
		if(side <5 && node == 0) return side;
		if(side >= 5 && side <13 && node < 2) return SideNodes[side-5][node];
		if(side == 13 && node <4) return FaceNodes[side-13][node];
		if(side >13 && side < 18) {
			if(node <3) {
                return FaceNodes[side-13][node];
            }
			else if(node==3) {
                return -1;//Previsto receber pelas faces triangulares - Cesar 2003-01-02
            }
        }
		if(side == 18 && node < 5) return node;
		PZError << "TPZPyramid::SideNodeLocId inconsistent side or node " << side << ' ' << node << endl;
		return -1;
	}
	
	void TPZPyramid::CenterPoint(int side, TPZVec<REAL> &center) {
		//center.Resize(Dimension);
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	int TPZPyramid::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZPyramid::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform TPZPyramid::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZPyramid::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
			return TPZTransform(0);
		}
		if(sidefrom == sideto) {
			return TPZTransform(sidedimension[sidefrom]);
		}
		if(sidefrom == NSides-1) {
			return TransformElementToSide(sideto);
		}
		int nhigh = nhighdimsides[sidefrom];
		int is;
		for(is=0; is<nhigh; is++) {
			if(highsides[sidefrom][is] == sideto) {
				int dfr = sidedimension[sidefrom];
				int dto = sidedimension[sideto];
				TPZTransform trans(dto,dfr);
				int i,j;
				for(i=0; i<dto; i++) {
					for(j=0; j<dfr; j++) {
						trans.Mult()(i,j) = sidetosidetransforms[sidefrom][is][j][i];
					}
					trans.Sum()(i,0) = sidetosidetransforms[sidefrom][is][3][i];
				}
				return trans;
			}
		}
		PZError << "TPZPyramid::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform(0);
	}
	
	TPZTransform TPZPyramid::TransformElementToSide(int side){
		
		if(side<0 || side>18){
			PZError << "TPZPyramid::TransformElementToSide called with side error\n";
			return TPZTransform(0,0);
		}
		
		TPZTransform t(sidedimension[side],3);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
				return t;
			case 5:
				t.Mult()(0,0) = 1.0;
				return t;
			case 6:
				t.Mult()(0,1) = 1.0;
				return t;
			case 7:
				t.Mult()(0,0) = -1.0;
				return t;
			case 8:
				t.Mult()(0,1) = -1.0;
				return t;
			case 9:
			case 12:
				t.Mult()(0,0) = 2.0;
				t.Sum()(0,0)  = 1.0;
				return t;
			case 10:
			case 11:
				t.Mult()(0,0) = -2.0;
				t.Sum()(0,0)  =  1.0;
				return t;
			case 13:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
			case 14:
				t.Mult()(0,0) =  0.5;
				t.Mult()(0,1) = -0.5;
				t.Mult()(1,2) =  1.0;
				return t;
			case 15:
			case 16:/** CONTEM ERRO AQUI */
				t.Mult()(0,0) =  0.5;
				t.Mult()(0,1) =  0.5;
				t.Mult()(1,2) =  1.0;
				return t;
			case 17:
				t.Mult()(0,0) = -0.5;
				t.Mult()(0,1) =  0.5;
				t.Mult()(1,2) =  1.0;
				return t;
			case 18:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	TPZTransform TPZPyramid::TransformSideToElement(int side){
		
		if(side<0 || side>18){
			PZError << "TPZPyramid::TransformSideToElement side out range\n";
			return TPZTransform(0,0);
		}
		TPZTransform t(3,sidedimension[side]);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 2:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 3:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  0.0;
				return t;
			case 4:
				t.Sum()(0,0) =  0.0;
				t.Sum()(1,0) =  0.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 5:
				t.Mult()(0,0) =  1.0;
				t.Sum() (1,0) = -1.0;
				return t;
			case 6:
				t.Mult()(1,0) =  1.0;
				t.Sum() (0,0) =  1.0;
				return t;
			case 7:
				t.Mult()(0,0) = -1.0;
				t.Sum() (1,0) =  1.0;
				return t;
			case 8:
				t.Mult()(1,0) = -1.0;
				t.Sum()(0,0)  = -1.0;
				return t;
			case 9:
				t.Mult()(0,0) =  0.5;
				t.Mult()(1,0) =  0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  = -0.5;
				t.Sum()(1,0)  = -0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 10:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) =  0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  =  0.5;
				t.Sum()(1,0)  = -0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 11:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) = -0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  =  0.5;
				t.Sum()(1,0)  =  0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 12:
				t.Mult()(0,0) =  0.5;
				t.Mult()(1,0) = -0.5;
				t.Mult()(2,0) =  0.5;
				t.Sum()(0,0)  = -0.5;
				t.Sum()(1,0)  =  0.5;
				t.Sum()(2,0)  =  0.5;
				return t;
			case 13:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
			case 14:
				t.Mult()(0,0) =  2.0;
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 15:
				t.Mult()(1,0) =  2.0;
				t.Mult()(0,1) = -1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 16:
				t.Mult()(0,0) =  2.0;
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,1) = -1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 17:
				t.Mult()(1,0) =  2.0;
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 18:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform(0,0);
	}
	
	
	
	TPZIntPoints * TPZPyramid::CreateSideIntegrationRule(int side, int order){
		if(side<0 || side>18) {
			PZError << "TPZPyramid::CreateSideIntegrationRule. Bad side number.\n";
			return 0;
		}
		if(side<5)   return new TPZInt1Point(order);            // sides 0 to 4 are vertices
		if(side<13)  return new TPZInt1d(order);           // sides 5 to 12 are lines
		if(side==13) return new TPZIntQuad(order,order);   // side 13 are quadrilateral (pyramid base)
		if(side<18)  {
			return new TPZIntTriang(order);                // sides 14 to 17 are triangles
		}
		if(side==18) {
			return new IntruleType(order);               // integration of the element
		}
		return 0;
	}

	MElementType TPZPyramid::Type()
	{
		return EPiramide;
	}

	MElementType TPZPyramid::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
				return EPoint;
			case 5:
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
				return EOned;
			case 13:
				return EQuadrilateral;
			case 14:
			case 15:
			case 16:
			case 17:
				return ETriangle;
			case 18:
				return EPiramide;
			default:
				return ENoType;
		}
	}
	
	
	int TPZPyramid::NumSides() {
		return 19;
	}
	
	//Tentando criar o metodo
	int TPZPyramid::NumSides(int dimension) {
		if(dimension<0 || dimension> 3) {
			PZError << "TPZPyramid::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 5;
		if(dimension==1) return 8;
		if(dimension==2) return 5;
		if(dimension==3) return 1;
		return -1;
	}
	int TPZPyramid::NContainedSides(int side) {
		if(side<0)   return -1;
		if(side<5)   return 1;//cantos : 0 a 4
		if(side<13)  return 3;//lados : 5 a 12
		if(side==13) return 9;//face : 13 , quadrilateral
		if(side<18)	 return 7;//faces : 14 a 17 , triangulares
		if(side==18) return 19;//centro : 18
		return -1;
	}
	
	int TPZPyramid::ContainedSideLocId(int side, int node) {
		if(side<0 || side>19 || node < 0) return -1;
		if(side<5) {
			if(node==0) return side;
		} 
		else if(side<9) {//5 a 8
			int s = side-5;//0,1,2
			if(!node) return s;//0,1,2,3
			if(node==1) return (s+1)%4;//1,2,0
			if(node==2) return side;//5,6,7,8
		}
		else if(side<13) {//9 a 12
			int s = side-9;//0,1,2,3
			if(!node) return s;//0,1,2,3
			if(node==1) return 4;//
			if(node==2) return side;//9,10,11,12
		} 
		else if(side<18) {//13 a 17
			int s = side-13;
			if(node<9) return FaceConnectLocId[s][node];
		} 
		else if(side==18 && node<19){
			return node;
		}
		PZError << "TPZShapePiram::ContainedSideLocId called for node = "
		<< node << " and side = " << side << "\n";
		return -1;
	}
	
	bool TPZPyramid::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		const REAL zeta = pt[2];
		
		if( (qsi < -1. - tol) || (qsi > 1.+tol) ||
		   (eta < -1. - tol) || (eta > 1.+tol) || 
		   (zeta < 0. - tol) || (zeta > 1.+tol) ||
		   (fabs(qsi) > 1.-zeta + tol) || (fabs(eta) > 1.-zeta + tol)
		   ) {
			return false;
		}
		else{
			return true;
		}  
		
		
	}//method
    
    bool TPZPyramid::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide) {
		double zero = 1.E-5;
		
		REAL qsi = InternalPar[0]; REAL eta = InternalPar[1]; REAL zeta = InternalPar[2];
		bool regularmap = true;
		switch(side)
		{
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            {
                SidePar.Resize(0); JacToSide.Resize(0,0);
                break;
            }
			case 5://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = qsi/(1.-zeta);
                    JacToSide(0,0) = 1./(1.-zeta); JacToSide(0,1) = 0.; JacToSide(0,2) = qsi/((1.-zeta)*(1.-zeta));
				}
				break;
				
			case 6://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = eta/(1.-zeta);
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 1./(1.-zeta); JacToSide(0,2) = eta/((1.-zeta)*(1.-zeta));
				}
				break;
				
			case 7://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = qsi/(zeta-1.);
                    JacToSide(0,0) = 1./(zeta-1.); JacToSide(0,1) = 0.; JacToSide(0,2) = -qsi/((zeta-1.)*(zeta-1.));
				}
				break;
				
			case 8://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = eta/(zeta-1.);
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 1./(zeta-1.); JacToSide(0,2) = -eta/((zeta-1.)*(zeta-1.));
				}
				break;
				
			case 9://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if( fabs(qsi-1.) < zero || fabs(eta-1.) < zero )
				{
                    SidePar[0] = -1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = -((eta*(-1. + qsi + zeta) + (-1. + zeta)*(-1. + qsi + 5.*zeta))/((-1. + qsi - 3.*zeta)*(-1. + zeta) + eta*(-1. + qsi + zeta)));
					REAL den = (-1. + qsi - 3.*zeta)*(-1. + zeta) + eta*(-1. + qsi + zeta);
                    JacToSide(0,0) = (8.*(-1. + zeta)*zeta*(-1. + eta + zeta))/(den*den);
                    JacToSide(0,1) = (8.*(-1. + zeta)*zeta*(-1. + qsi + zeta))/(den*den);
                    JacToSide(0,2) = (-8.*((-1. + qsi)*(-1. + zeta)*(-1. + zeta) + eta*((-1. + zeta)*(-1. + zeta) + qsi*(-1. + 2.*zeta))))/(den*den);
				}
				break;
				
			case 10://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if( fabs(qsi+1.) < zero || fabs(eta-1.) < zero )
				{
                    SidePar[0] = -1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = (-((1. + qsi - 5.*zeta)*(-1. + zeta)) + eta*(-1. - qsi + zeta))/(eta*(1. + qsi - zeta) + (-1. + zeta)*(1. + qsi + 3.*zeta));
					REAL den = eta*(1. + qsi - zeta) + (-1. + zeta)*(1. + qsi + 3.*zeta);
                    JacToSide(0,0) = (-8.*(-1. + zeta)*zeta*(-1. + eta + zeta))/(den*den);
                    JacToSide(0,1) = (-8.*(1. + qsi - zeta)*(-1. + zeta)*zeta)/(den*den);
                    JacToSide(0,2) = (8.*((1. + qsi)*(-1. + zeta)*(-1. + zeta) - eta*(qsi + (-1. + zeta)*(-1. + zeta) - 2.*qsi*zeta)))/(den*den);
				}
				break;
				
			case 11://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if( fabs(qsi+1.) < zero || fabs(eta+1.) < zero )
				{
                    SidePar[0] = -1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else if(fabs(zeta-1.) < zero )
				{
                    SidePar[0] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = -1. + (8.*(-1. + zeta)*zeta)/(eta*(-1. - qsi + zeta) + (-1. + zeta)*(1. + qsi + 3.*zeta));
					REAL den = eta*(-1. - qsi + zeta) + (-1. + zeta)*(1. + qsi + 3.*zeta);
                    JacToSide(0,0) = (-8.*(-1. + zeta)*zeta*(-1. - eta + zeta))/(den*den);
                    JacToSide(0,1) = (-8.*(-1. + zeta)*zeta*(-1. - qsi + zeta))/(den*den);
                    JacToSide(0,2) = (8.*((1. + qsi)*(-1. + zeta)*(-1. + zeta) + eta*(qsi + (-1. + zeta)*(-1. + zeta) - 2.*qsi*zeta)))/(den*den);
				}
				break;
				
			case 12://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if( fabs(qsi-1.) < zero || fabs(eta+1.) > zero)
				{
                    SidePar[0] = -1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = (-(eta*(-1. + qsi + zeta)) + (-1. + zeta)*(-1. + qsi + 5.*zeta))/(-((-1. + qsi - 3.*zeta)*(-1. + zeta)) + eta*(-1. + qsi + zeta));
					REAL den = (-1. + qsi - 3.*zeta)*(-1. + zeta) - eta*(-1. + qsi + zeta);
                    JacToSide(0,0) = (-8.*(1. + eta - zeta)*(-1. + zeta)*zeta)/(den*den);
                    JacToSide(0,1) = (-8.*(-1. + zeta)*zeta*(-1. + qsi + zeta))/(den*den);
                    JacToSide(0,2) = (8.*(-((-1. + qsi)*(-1. + zeta)*(-1. + zeta)) + eta*((-1. + zeta)*(-1. + zeta) + qsi*(-1. + 2.*zeta))))/(den*den);
				}
				break;
				
			case 13://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.; SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = qsi/(1. - zeta);
                    SidePar[1] = eta/(1. - zeta);
                    JacToSide(0,0) = 1./(1. - zeta); JacToSide(0,1) = 0.; JacToSide(0,2) = qsi/((1. - zeta)*(1. - zeta));
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 1./(1. - zeta); JacToSide(1,2) = eta/((1. - zeta)*(1. - zeta));
					regularmap = false;
				}
				break;
				
			case 14://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.; SidePar[1] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else if(fabs(eta-1.) < zero)
				{
                    SidePar[0] = qsi/2. + 0.5; SidePar[1] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = ((-1. + eta + zeta)*(-1. - qsi + zeta))/(2.*(-1. + eta - zeta)*(-1. + zeta));
                    SidePar[1] = (2.*zeta)/(1. - eta + zeta);
                    JacToSide(0,0) = -(-1. + eta + zeta)/(2.*(-1. + eta - zeta)*(-1. + zeta));
                    JacToSide(0,1) = ((1. + qsi - zeta)*zeta)/((-1. + zeta)*(1. - eta + zeta)*(1. - eta + zeta));
                    JacToSide(0,2) = (eta*eta*qsi - (2. + qsi)*(-1. + zeta)*(-1. + zeta) + 2.*eta*(1. - (2. + qsi)*zeta + zeta*zeta))/(2.*(-1. + zeta)*(-1. + zeta)*(1 - eta + zeta)*(1 - eta + zeta));
                    JacToSide(1,0) = 0.;
                    JacToSide(1,1) = (2.*zeta)/((1. - eta + zeta)*(1. - eta + zeta));
                    JacToSide(1,2) = (2. - 2.*eta)/((1. - eta + zeta)*(1. - eta + zeta));
				}
				break;
				
			case 15://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.; SidePar[1] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else if(fabs(qsi+1.) < zero)
				{
                    SidePar[0] = eta/2. + .5; SidePar[1] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.5; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = ((1. + eta - zeta)*(-1. - qsi + zeta))/(2.*(-1. + zeta)*(1. + qsi + zeta));
                    SidePar[1] = (2.*zeta)/(1. + qsi + zeta);
                    JacToSide(0,0) = (zeta*(-1. - eta + zeta))/((-1. + zeta)*(1. + qsi + zeta)*(1. + qsi + zeta));
                    JacToSide(0,1) = (-1. - qsi + zeta)/(2.*(-1. + zeta)*(1. + qsi + zeta));
                    JacToSide(0,2) = (-2.*(1. + qsi)*(-1. + zeta)*(-1. + zeta) + eta*(qsi*qsi - (-1. + zeta)*(-1. + zeta) + 2.*qsi*zeta))/(2.*(-1. + zeta)*(-1. + zeta)*(1. + qsi + zeta)*(1. + qsi + zeta));
                    JacToSide(1,0) = (-2*zeta)/((1 + qsi + zeta)*(1 + qsi + zeta));
                    JacToSide(1,1) = 0.;
                    JacToSide(1,2) = (2.*(1. + qsi))/((1. + qsi + zeta)*(1. + qsi + zeta));
					
				}
				break;
				
			case 16://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.; SidePar[1] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else if(fabs(eta+1.) < zero)
				{
                    SidePar[0] = -qsi/2. + .5; SidePar[1] = 0.;
                    JacToSide(0,0) = -0.5; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = ((1. + eta - zeta)*(-1. + qsi + zeta))/(2.*(-1. + zeta)*(1. + eta + zeta));
                    SidePar[1] = (2.*zeta)/(1. + eta + zeta);
                    JacToSide(0,0) = (1. + eta - zeta)/(2.*(-1. + zeta)*(1. + eta + zeta));
                    JacToSide(0,1) = (zeta*(-1. + qsi + zeta))/((-1. + zeta)*(1. + eta + zeta)*(1. + eta + zeta));
                    JacToSide(0,2) = (-(eta*eta*qsi) + (-2 + qsi)*(-1. + zeta)*(-1. + zeta) - 2.*eta*(1. + (-2. + qsi)*zeta + zeta*zeta))/(2.*(-1. + zeta)*(-1. + zeta)*(1. + eta + zeta)*(1. + eta + zeta));
                    JacToSide(1,0) = 0.;
                    JacToSide(1,1) = (-2.*zeta)/((1. + eta + zeta)*(1. + eta + zeta));
                    JacToSide(1,2) = (2.*(1. + eta))/((1. + eta + zeta)*(1. + eta + zeta));
				}
				break;
				
			case 17://2D
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(zeta-1.) < zero)
				{
                    SidePar[0] = 0.; SidePar[1] = 1.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else if(fabs(qsi-1.) < zero)
				{
                    SidePar[0] = 0.5 - eta/2.; SidePar[1] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = -0.5; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) =  0.; JacToSide(1,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = ((-1. + eta + zeta)*(-1. + qsi + zeta))/(2.*(-1. + qsi - zeta)*(-1. + zeta));
                    SidePar[1] = (2.*zeta)/(1. - qsi + zeta);
                    JacToSide(0,0) = -((zeta*(-1. + eta + zeta))/((-1. + zeta)*(1. - qsi + zeta)*(1. - qsi + zeta)));
                    JacToSide(0,1) = (-1. + qsi + zeta)/(2.*(-1. + qsi - zeta)*(-1. + zeta));
                    JacToSide(0,2) = (2.*(-1. + qsi)*(-1. + zeta)*(-1. + zeta) + eta*(-qsi*qsi + (-1. + zeta)*(-1. + zeta) + 2.*qsi*zeta))/(2.*(-1. + zeta)*(-1. + zeta)*(1. - qsi + zeta)*(1. - qsi + zeta));
                    JacToSide(1,0) = (2.*zeta)/((1. - qsi + zeta)*(1. - qsi + zeta));
                    JacToSide(1,1) = 0.;
                    JacToSide(1,2) = (2. - 2.*qsi)/((1. - qsi + zeta)*(1. - qsi + zeta));
				}
				break;
		}
		if(side > 17)
		{
			cout << "Cant compute MapToSide method in TPZGeoPyramid class!\nParameter (SIDE) must be between 5 and 17!\nMethod Aborted!\n";
			DebugStop();
		}
		return regularmap;
		
	}
    
    void TPZPyramid::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        if(node > NCornerNodes)
        {
            DebugStop();
        }
        nodeCoord.Resize(Dimension, 0.);
        switch (node) {
            case (0):
            {
                nodeCoord[0] = -1.;
                nodeCoord[1] = -1.;
                nodeCoord[2] =  0.;
                break;
            }
            case (1):
            {
                nodeCoord[0] =  1.;
                nodeCoord[1] = -1.;
                nodeCoord[2] =  0.;
                break;
            }
            case (2):
            {
                nodeCoord[0] = 1.;
                nodeCoord[1] = 1.;
                nodeCoord[2] = 0.;
                break;
            }
            case (3):
            {
                nodeCoord[0] = -1.;
                nodeCoord[1] =  1.;
                nodeCoord[2] =  0.;
                break;
            }
            case (4):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 0.;
                nodeCoord[2] = 1.;
                break;
            }
            default:
            {
                DebugStop();
                break;
            }
        }
    }
	
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZPyramid::GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
		std::cout << "Please implement me\n";
		DebugStop();
	}
    
    void computedirectionsPy(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);
    
    void computedirectionsPy(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions)
    {
        REAL detgrad = 0.0;
        TPZVec<REAL> u(3);
        TPZVec<REAL> v(3);
        TPZVec<REAL> uxv(3);// result
        int cont = 0;
        
        for (int ivet=inicio; ivet<=fim; ivet++)
        {
            for (int ilin=0; ilin<3; ilin++)
            {
                u[ilin] = t1vec(ilin,ivet);
                v[ilin] = t2vec(ilin,ivet);
            }
            TPZVec<REAL> e2(3);
            detgrad = 0.0;
            REAL normaX0xX1 = 0.0;
            //TPZNumeric::ProdVetorial(u,v,e2);
            e2[0] = u[1]*v[2]-u[2]*v[1];
            e2[1] = -(u[0]*v[2]-v[0]*u[2]);
            e2[2] = u[0]*v[1]-v[0]*u[1];
            
            // calc do v gradx*b
            TPZManVector<REAL,3> dxt1(3,0.),dxt2(3,0.),dxt3(3,0.),Vvec(3,0.);
            REAL be2 = 0.0, ne2 = 0.0;
            for(int i=0;i<3;i++)
            {
                ne2 += e2[i]*e2[i];
            }
            ne2 = sqrt(fabs(ne2));
            for (int il=0; il<3; il++)
            {
                for (int i = 0 ; i<3; i++)
                {
                    dxt1[il] += gradx(il,i) * t1vec(i,ivet);
                    dxt2[il] += gradx(il,i) * t2vec(i,ivet);
                    dxt3[il] += gradx(il,i) * e2[i]/ne2;
                    Vvec[il] += gradx(il,i) * bvec(i,ivet);
                }
                be2 += bvec(il,ivet)*e2[il]/ne2;
            }
            TPZManVector<REAL,3> normal(3,0.);
            //TPZNumeric::ProdVetorial(dxt1,dxt2,normal);
            normal[0] = dxt1[1]*dxt2[2]-dxt1[2]*dxt2[1];
            normal[1] = -(dxt1[0]*dxt2[2]-dxt2[0]*dxt1[2]);
            normal[2] = dxt1[0]*dxt2[1]-dxt2[0]*dxt1[1];
            
            for (int pos=0; pos<3; pos++)
            {
                detgrad += normal[pos]*dxt3[pos];//uxv[pos]*gradx.GetVal(pos, 2);
                normaX0xX1 += normal[pos]*normal[pos]; //uxv[pos]*uxv[pos];
            }
            TPZFMatrix<REAL> Wvec(3,1);
            detgrad = fabs(detgrad);
            normaX0xX1 = sqrt(normaX0xX1);
            
            for (int il=0; il<3; il++)
            {
                Wvec(il,0) = Vvec[il]*normaX0xX1/(detgrad*be2);
                directions(il,cont) = Wvec(il,0);
            }
            cont++;
        }
        
        
    }

    
    void TPZPyramid::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        if(gradx.Cols()!=3)
        { std::cout << "Gradient dimensions are not compatible with this topology" << std::endl;
            DebugStop();
        }
        TPZFMatrix<REAL> bvec(3,58);
        int numvec = bvec.Cols();
        TPZFMatrix<REAL> t1vec(3,numvec);
        TPZFMatrix<REAL> t2vec(3,numvec);
        
        directions.Redim(3, numvec);
        for (int lin = 0; lin<numvec; lin++)
        {
            for(int col = 0;col<3;col++)
            {
                bvec.PutVal(col, lin, bPiram[lin][col]);
                t1vec.PutVal(col, lin, t1Piram[lin][col]);
                t2vec.PutVal(col, lin, t2Piram[lin][col]);
            }
        }
        
        // calcula os vetores
        switch (side) {
            case 13:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 0, fim = 8;
                computedirectionsPy( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPi[ip+inicio];
                }
                
            }
                break;
            case 14:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 9, fim = 15;
                computedirectionsPy( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPi[ip+inicio];
                }
            }
                break;
            case 15:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 16, fim = 22;
                computedirectionsPy( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPi[ip+inicio];
                }
            }
                break;
            case 16:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 23, fim = 29;
                computedirectionsPy( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPi[ip+inicio];
                }
            }
                break;
            case 17:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 30, fim = 36;
                computedirectionsPy( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPi[ip+inicio];
                }
            }
                break;
            case 18:
            {
                directions.Resize(3, 21);
                sidevectors.Resize(21);
                int inicio = 37, fim = 57;
                computedirectionsPy( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPi[ip+inicio];
                }
            }
                break;
                
            default:
                break;
        }

	}
    
    void TPZPyramid::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao)
    {
        int nsides = NumSides()*3;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderPi[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
    }

    void TPZPyramid::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao, TPZVec<int> &sidevectors)
    {
        int nsides = NumSides()*3;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderPi[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
        
        for (int i=0; i<Dimension*NumSides(); i++) {
            sidevectors[i] = vectorsideorderPi[i];
        }
    }

}
