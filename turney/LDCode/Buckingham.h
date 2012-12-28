/*					                 Buckingham.h		  		  	          		*/
/*						              	12/16/2009  					            		*/
/*********************************************************************
*    Header file that contains the POT_BUCK derivative class    	   *
*  definition for the Buckingham potential.		    				           *
*********************************************************************/

#if !defined(BUCKINGHAM_H)
#define BUCKINGHAM_H


/*DEFINE HEADERS*/
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class POT_BUCK : public POTENTIAL {
 public:
  double *A, *rho, *C, *cutoff, *Acutoff;//Buckingham potential parameters
	POT_BUCK(int, int*);          //Constructor (Allocates memory)
	~POT_BUCK();      		        //Destructor (Deallocates memory)
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


#endif
