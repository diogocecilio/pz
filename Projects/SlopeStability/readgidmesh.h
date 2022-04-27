//
// Created by Diogo Cecílio on 10/12/21.
//

#pragma once

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
//#include "TPZDarcyMaterial.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <chrono>

using namespace std;

class readgidmesh
{
public:
    /**
     * @brief Default constructor
     */
    readgidmesh();

    /**
     * @brief Class constructor
     * @param [in] mesh
     */
    readgidmesh ( string file );
    ~readgidmesh();




    template <class T>
    std::vector<T> str_vec( std::vector<std::string> &vs );


    void ReadMesh ( );
	

    
    void  FindIds ( TPZVec<double> constcoorddata,TPZVec<int> constcoord, std::vector<int>& ids )
    {
		REAL tol = 1.e-12;
		//constcoorddata vector containig info about face to search id. It must contain any coodinate locate in the in the face
        TPZFMatrix<REAL> elcoords;
		std::vector<double> elcoodsvec;
        int nels = fallcoords.size();
        GetElCoords (  0, elcoords );
        int nnodes = elcoords.Rows();
        int sum=0;
        //constcoord.size() = 1 face
        //constcoord.size() = 2 linha
        //constcoord.size() = 3 pontos

        std::vector<int> dirs;
        for ( int iconst=0; iconst<constcoord.size(); iconst++ ) {
            sum+=constcoord[iconst];
            if ( constcoord[iconst]==1 ) {
                dirs.push_back ( iconst );
            }

        }
        for ( int iel = 0; iel < nels; iel++ ) {
            GetElCoords (  iel, elcoords );

            for ( int inode = 0; inode < nnodes; inode++ ) {

                if ( sum==1 ) {
                    if ( fabs ( elcoords(inode,dirs[0]) - constcoorddata[dirs[0]] ) <tol ) {
                        ids.push_back ( fmeshtopology(iel,inode) );
                    }
                } else if ( sum==2 ) {
                    if ( fabs ( elcoords(inode,dirs[0]) - constcoorddata[dirs[0]] ) <tol && abs ( elcoords(inode,dirs[1]) - constcoorddata[dirs[1]] ) <tol ) {
                        ids.push_back (  fmeshtopology(iel,inode)  );
                    }
                } else if ( sum==3 ) {
                    if ( fabs ( elcoords(inode,dirs[0]) - constcoorddata[dirs[0]] ) <tol && abs ( elcoords(inode,dirs[1]) - constcoorddata[dirs[1]] ) <tol && abs ( elcoords(inode,dirs[2]) - constcoorddata[dirs[2]] ) <tol ) {
                        ids.push_back (  fmeshtopology(iel,inode)  );
                    }
                }


            }


        }

        sort ( ids.begin(), ids.end() );
        ids.erase ( unique ( ids.begin(), ids.end() ), ids.end() );
    }
    void  GetElCoords ( int el, TPZFMatrix<double>  & elcoords )
    {
        elcoords.Resize ( fallcoords[el].size(), 3 );
		
        for ( int j = 0; j < fallcoords[el].size(); j++ ) {
            double x = fallcoords[el][j][0];
            double y = fallcoords[el][j][1];
            double z = fallcoords[el][j][2];
            elcoords(j,0) = x;
            elcoords(j,1) = y;
            elcoords(j,2) = z;
        }
    }
    
	void  GetElCoords ( int el, vector<double>  & elcoords )
    {
        elcoords.resize ( 3 );
		
        for ( int j = 0; j < fallcoords[el].size(); j++ ) {
            double x = fallcoords[el][j][0];
            double y = fallcoords[el][j][1];
            double z = fallcoords[el][j][2];
            elcoords[0] = x;
            elcoords[1] = y;
            elcoords[2] = z;
        }
    }
	
	inline TPZFMatrix<int> GetTopology()
	{
		return fmeshtopology;
	}
	inline TPZFMatrix<double> GetCoords()
	{
		return fmeshcoords;
	}
	std::vector<std::vector< std::vector<double > > > GetAllCoords()
	{
		return fallcoords;
	}

    std::vector<std::vector<int>>   LineTopology ( std::vector<int> ids )
    {
        int  k = 0;


        std::vector<std::vector<int>> vg;
        for ( int j = 0; j < ids.size() / (fOrder+1); j++ ) {
            std::vector<int> v;
            for ( int i = 0; i < fOrder + 1; i++ ) {
                v.push_back ( ids[i + k] );
            }
            vg.push_back ( v );
            k += fOrder;
        }
        return vg;
    }

void  ToMatInt ( std::vector<std::vector<int>> in, TPZFMatrix<int> & out )
{
    int  rows = in.size();
    int cols = in[0].size();
    out.Resize ( rows, cols );
    for ( int i = 0; i < rows; i++ ) for ( int j = 0; j < cols; j++ ) out(i,j)= in[i][j];
}
	

private:

string ffile;
TPZFMatrix<int> fmeshtopology;
TPZFMatrix<double> fmeshcoords;
std::vector<std::vector< std::vector<double > > > fallcoords;
int fOrder;
};


