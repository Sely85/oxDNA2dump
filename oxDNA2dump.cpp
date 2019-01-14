//Lara Querciagrossa
#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
#include "math-func.hpp"

#define THIRD 0.333333333333333
#define ROOTTHREE 1.73205080756888       
#define PI 3.141592653589793238462643383279502884197
 
using namespace std;

int TypeToInt(char base)
{
  int type = 0;

  //Assigning a numerical index to each base
  
  if (base == 'A')
    type = 2;
  else if (base == 'T')
    type = 3;
  else if (base == 'C')
    type = 4;
  else if (base == 'G')
    type = 5;
  else if (base == 'U')
    type = 6;
  else 
    {
      std::cout << "ERROR: Unknown nucleotide!" << std::endl;
    }

  return type;
}


int main(int argc, char *argv[])
{

  if (argc != 3)
    {
      std::cout << "[ERROR] Wrong number of input parameters: ./oxDNA2dump <oxDNAtop> <oxDNAconf>" << std::endl;
      return 0;
    }

  if (string(argv[1]) == "--help" || argv[1] == "-h") {
    std::cout << " " << std::endl;
    std::cout << " Welcome to the oxDNA2dump code!" << std::endl;
    std::cout << " This code will convert an oxDNA configuration into a dump file" << std::endl;
    std::cout << " that can be read by ovito" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " ./oxDNA2dump <oxDNAtop> <oxDNAconf>" << std::endl;
    std::cout << " " << std::endl;
  }
  else {

    ifstream top ( argv[1] );
    if ( !top.is_open() )
      {
	cout<<"ERROR: Could not open topology file\n";
	return 0;
      }

    ifstream conf ( argv[2] );
    if ( !conf.is_open() )
      {
	cout<<"ERROR: Could not open configuration file\n";
	return 0;
      }

    std::cout << "Your topology file is: " << argv[1] << std::endl;
    std::cout << "Your configuration file is: " << argv[2] << std::endl;
    std::cout << " " << std::endl;


    /** Read Topology **/
    int totnuc = 0; //number of nucleotides
    int totstr = 0; //number of strands
    top >> totnuc >> totstr;
    std::cout << "Your topology file states that there are " << totnuc << " nucleotides and " << totstr << " strands." << std::endl;

    // Each line of topology is made up of: number of strand, nucleotide type, previous nucleotide (3' direction), following nucleotide (5' direction)
    int strand[totnuc], prev[totnuc], suc[totnuc]; 
    char nucl[totnuc];

    for (int i=0; i<totnuc; i++)
      {
	top >> strand[i] >> nucl[i] >> prev[i] >> suc[i];
      }

    /** Read configuration(s) **/

    //Count lines in files
    string line; 
    int totlines = 0;
    int totsteps = 0;
    while ( getline(conf,line) )
      {
	char first = line[0];
	if (first == 't')
	  {
	    totsteps++;
	  }
	totlines++;
      }
    std::cout << "In your configuration file, there are " << totsteps << " steps." << std::endl;

    //Check correct number of nucleotides w.r.t. what stated in topology
    int totbead = (totlines-3*totsteps)/totsteps; 
    if ( totbead != totnuc )
      {
	std::cout << "ERROR: topology and configuration do not match!" << std::endl;
	return 0;
      }

    //Configuration file has been read till the end: this position is cleared and than we set the new position to be read next at the beginning of the file
    conf.clear();
    conf.seekg(0, std::ios::beg);

    //Write a dump file to be read in ovito
    ofstream fileout("oxDNAconv.dump", ios::out);
    //Save energies in another file
    ofstream eneout("oxDNAconv.ene", ios::out);
    eneout << "#[1]timestep [2]etot [3]upot [4]ekin" << std::endl;

    //File will be read one line at a time to save memory
    for (int i=0; i<totsteps; i++)
      {
	std::cout << "Converting step " << i+1 << "/" << totsteps << std::endl;

	int timestep;
	double sidex, sidey, sidez;
	double etot, upot, ekin;
	string null;

	conf >> null >> null >> timestep;
	conf >> null >> null >> sidex >> sidey >> sidez;
	conf >> null >> null >> etot >> upot >> ekin;

	//Saving "start" part of dump file
	fileout << "ITEM: TIMESTEP" << std::endl;
	fileout << timestep << std::endl;
	fileout << "ITEM: NUMBER OF ATOMS" << std::endl;
	fileout << 2*totbead << std::endl;
	fileout << "ITEM: BOX BOUNDS" << std::endl;
	fileout << -sidex/2.0 << " " << sidex/2.0 << std::endl;
	fileout << -sidey/2.0 << " " << sidey/2.0 << std::endl;
	fileout << -sidez/2.0 << " " << sidez/2.0 << std::endl;

	//Saving energies into energies file
	eneout << timestep << " " << etot << " " << upot << " " << ekin << std::endl;

	int cnt = 1;
	//Saving positions into dump file
	fileout << "ITEM: ATOMS id type x y z c_q[1] c_q[2] c_q[3] c_q[4] c_shape[1] c_shape[2] c_shape[3] vx vy vz angmomx angmomy angmomz mol" << std::endl;
	for (int b=0; b<totbead; b++)
	  {
	    double rx, ry, rz; //Position
	    double bx, by, bz; //Backbone-baseversor
	    double nx, ny, nz; //Normalversor 
	    double vx, vy, vz; //Velocities 
	    double wx, wy, wz; //Angular velocities

	    conf >> rx >> ry >> rz >> bx >> by >> bz >> nx >> ny >> nz >> vx >> vy >> vz >> wx >> wy >> wz;

	    //Convert character type read in topology into a number
	    int type = TypeToInt(nucl[b]);
	    //Backbone site
	    fileout << cnt << " " << " 1 " << " " << rx-0.4*bx << " " << ry-0.4*by << " " << rz-0.4*bz << " " << " 1 0 0 0   " <<  " 0.3 0.3 0.3 " << vx << " " << vy << " " << vz << " " << wx << " " << wy << " " << wz << " " << strand[b] << std::endl;
	    cnt++;

//	    //Stacking center
//	    fileout << cnt << " " << type << " " << rx+0.34*bx << " " << ry+0.34*by << " " << rz+0.34*bz << " " << " 1 0 0 0   " <<  " 0.3 0.3 0.3 " << vx << " " << vy << " " << vz << " " << wx << " " << wy << " " << wz << " " << strand[b] << std::endl;
//	    cnt++;

	    //Base site
	    double bv[3] = {bx, by, bz};
	    double nv[3] = {nx, ny, nz};
	    double* cross = crossprod(nv, bv);

	    double mat[3][3];

	    mat[0][0] = bx;
	    mat[0][1] = by;
	    mat[0][2] = bz;
	    mat[2][0] = nx;
	    mat[2][1] = ny;
	    mat[2][2] = nz;
	    mat[1][0] = cross[0];
	    mat[1][1] = cross[1];
	    mat[1][2] = cross[2];

	    bool active = false;
	    double* q = MatrixToQuaternion(mat, active);

	    fileout << cnt << " " << type << " " << rx+0.4*bx << " " << ry+0.4*by << " " << rz+0.4*bz << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " " <<  " 0.4 0.3 0.1 " << vx << " " << vy << " " << vz << " " << wx << " " << wy << " " << wz << " " << strand[b] << std::endl;
	    cnt++;
	  }

      }

    std::cout << " " << std::endl;
    std::cout << "Files oxDNAconv.dump (configurations) and oxDNAconv.ene (energies) have been written." << std::endl;
    std::cout << " " << std::endl;

  }
}
