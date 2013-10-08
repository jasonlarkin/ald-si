/*                             Lattice.h                            */
/*                            12/02/2008                            */
/*********************************************************************
*    Header file for the LATTICE class.  This class stores all       *
*  parameters needed to define the lattice.                          *
*********************************************************************/

#ifndef LATTICE_H
#define LATTICE_H


/*DECLARE HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;


/*DEFINE PREPROCESSOR VARIABLES*/
#define PERIODIC 0
#define SCATTERING 1
#define FREE 2


/*DEFINE CLASSES*/
class LATTICE {
public:
  LATTICE();
  ~LATTICE();
	string Define(ifstream &Input);
	int KeywordToInt(ifstream &Input);
	void Output();        //Writes the variables
	void Initialize();    //Initializes variables after they have been defined (non-dimensionalizes)
  double **a;           //Direct lattice vectors (a^T*v transforms v from Cartesian to reciprocal basis)
  double **b;           //Reciprocal lattice vectors/2pi (b^T*v transforms v from reciprocal to Cartesian basis)
  double V;             //Volume of unit cell defined by direct lattice vectors
  int *BC;              //Boundary conditions
  int *N;               //Number of unit cells along each vector
  int *lower;           //Lower bounds of wave vector
  int *upper;           //Upper bounds of wave vector
};

#endif
