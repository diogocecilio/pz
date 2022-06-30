#include "readgidmesh.h"


readgidmesh::readgidmesh(){
	
}

readgidmesh::~readgidmesh(){
	
}

readgidmesh::readgidmesh(string file ){
	ffile=file;
}

template <class T>
std::vector<T> readgidmesh::str_vec( std::vector<std::string> &vs )
{
    std::vector<T> ret;
    for ( std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it ) {
        istringstream iss ( *it );
        T temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}

void readgidmesh::ReadMesh (  )
{
    std::vector<std::vector<int>> topol;
    string line, temp;

    ifstream myfile ( ffile );

	
    std::vector<std::vector<double>> coords;

    if ( myfile.is_open() ) {
		getline ( myfile, line );
		std::cout << line<< std::endl;
		getline ( myfile, line );
		std::cout << line<< std::endl;
        while ( getline ( myfile, line ) ) {
			if(line=="End Coordinates")break;
            std::vector<string> tokens;
            istringstream iss ( line );
            while ( iss >> temp )
                tokens.push_back ( temp );
			std::vector<double> input_doub_temp= str_vec<double> ( tokens );
            coords.push_back ( input_doub_temp );
        }
       // myfile.close();
    } else std::cout << "Unable to open file";
	int nodes = coords.size();
	fmeshcoords.Resize(nodes,3);
	for(int inode=0;inode<nodes;inode++){
		fmeshcoords(inode,0)=coords[inode][1];
		fmeshcoords(inode,1)=coords[inode][2];
		fmeshcoords(inode,2)=coords[inode][3];
	}
    //fmeshcoords.Print();
	
    if ( myfile.is_open() ) {
		getline ( myfile, line );
		std::cout << line<< std::endl;
		getline ( myfile, line );
		std::cout << line<< std::endl;
        while ( getline ( myfile, line ) ) {
			if(line=="End Elements")break;
            std::vector<string> tokens;
            istringstream iss ( line );
            while ( iss >> temp )
                tokens.push_back ( temp );
            std::vector<int> input_int = str_vec<int> ( tokens );
            topol.push_back ( input_int );
        }
        myfile.close();
    } else std::cout << "Unable to open file";

	int els = topol.size();
	cout << "els = "<<els<<endl;
	int elnodes = topol[0].size()-1;
	cout << "elnodes = "<<elnodes<<endl;
	fmeshtopology.Resize(els,elnodes);
	for(int iel=0;iel<els;iel++){
		for(int elnode=0;elnode<elnodes;elnode++)
		{
			if((topol[iel][elnode+1]-1)>800000000){
				fmeshtopology(iel,elnode)=-1;
			}else{
			fmeshtopology(iel,elnode)=topol[iel][elnode+1]-1;
			}
		}
	}
    fmeshtopology.Print(cout);
	//cout << "sda = "<<fmeshtopology(19318,3) << endl;
	std::vector<double> temp33 ( 3 );
    for ( int i = 0; i < fmeshtopology.Rows(); i++ ) {
        std::vector< std::vector<double> > temp22;
        for ( int j = 0; j < fmeshtopology.Cols(); j++ ) {
            int top = fmeshtopology(i,j);
			//cout << "i = " << i << " j = "<<j << " top = "<<top << endl;
			if(top!=-1){
            temp33[0] = fmeshcoords(top,0);
            temp33[1] = fmeshcoords(top,1);
            temp33[2] = fmeshcoords(top,2);
            temp22.push_back ( temp33 );
			}
        }
        fallcoords.push_back ( temp22 );
    }

}
