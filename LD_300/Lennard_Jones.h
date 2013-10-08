/*					                Lennard_Jones.h		  			          		*/
/*						              	12/02/2008  					            		*/
/*********************************************************************
*    Header file that contains the LENNARD_JONES derivative class	   *
*  definition for the Lennard-Jones potential.						           *
*********************************************************************/

#if !defined(LENNARD_JONES_H)
#define LENNARD_JONES_H


/*DEFINE HEADERS*/
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class POT_LJ : public POTENTIAL {
 public:
  double *epsilon, *sigma, *cutoff, *Acutoff;//LJ potential parameters
	POT_LJ(int, int*);          //Constructor (Allocates memory)
	~POT_LJ();      		        //Destructor (Deallocates memory)
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


#endif
