/*						                 LDCode.h	      				            	*/
/*							              12/02/2008	              						*/
/*********************************************************************
*    Main header file to be attached to all subroutines of programs  *
*  using the LDCode style input file.  Contains the global variables *
*  and their class definitions.						                           *
*********************************************************************/

#if !defined(LDCODE_H)
#define LDCODE_H


///*COMMENT NEXT LINE TO COMPILE SERIAL VERSION
//# define PARALLEL
# if defined PARALLEL
#  include <mpi.h>
# endif
//*/

/*COMMENT NEXT LINE TO IGNORE THE FREQUENCY SHIFT*/
#define FREQ_SHIFT


/*DEFINE GLOBAL PREPROCESSOR VARIABLES*/
#define CARTESIAN   0
#define CYLINDRICAL 1
#define SPHERICAL   2
#define DIRECT      3
#define RECIPROCAL  4


/*DEFINE GLOBAL PREPROCESSOR COMMANDS*/
#define Log if(PE==0) Log1
#define DIM 3


/*DEFINE GLOBAL HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

#include "Neighbor_List.h"
#include "Potential_Derivative.h"
#include "Material.h"
#include "Lattice.h"
#include "Unit_Cell.h"
#include "Parameter.h"
#include "Symmetry.h"

#include "Potential.h"
#include "Lennard_Jones.h"
#include "Lennard_Jones_Spline.h"
#include "Inverse_Twelfth_Power.h"
#include "Buckingham.h"
#include "Stillinger_Weber.h"
#include "REBO.h"
#include "Tersoff.h"
#include "Coulomb.h"


/*DEFINE GLOBAL VARIABLES*/
extern int PE, nPE;
extern ofstream Log1;		      //Log file
//extern int DIM;               //Dimension of system
extern int UC_DOF;            //DIM*(# of atoms in unit cell)
extern int BODY;              //Largest body loop (2,3,etc...)
extern PARAMETER *Parameter;  //(Non-)dimensionalizing factors
extern MATERIAL *Material;    //MATERIAL does not depend upon other classes
extern LATTICE *Lattice;      //LATTICE does not depend upon other classes
extern UNIT_CELL *Unit_Cell;  //UNIT_CELL must be after LATTICE, MATERIAL
extern SYMMETRY *Symmetry;    //SYMMETRY must be after LATTICE


#endif
