/*					                   Coulomb.h		  		  	          		*/
/*						              	12/23/2009  					            		*/
/*********************************************************************
*    Header file that contains the POT_COULOMB derivative class      *
*  definition for the Coulomb electrostatic potential.		           *
*********************************************************************/

#if !defined(COULOMB_H)
#define COULOMB_H


/*DEFINE HEADERS*/
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class POT_COULOMB : public POTENTIAL {
 public:
  double *q;                  //Charge in coulombs
  double *cutoff, *Acutoff, alpha;//Cutoffs and std. dev.
  double E_const, F_const;    //Constant parameters used in force calculation
  double erfccut, a_pi;       //Constant parameters used in force calculation
	POT_COULOMB(int, int*);     //Constructor (Allocates memory)
	~POT_COULOMB();      		    //Destructor (Deallocates memory)
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


#endif
