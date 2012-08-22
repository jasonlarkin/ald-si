/*						               POT_TERSOFF.h	                 				*/
/*							              04/12/2009              							*/
/*********************************************************************
*    Header file that contains the POT_TERSOFF derivative class		 	 *
*  definition for the TERSOFF N-body potential.	                  	 *
*********************************************************************/

#if !defined(TERSOFF_H)
#define TERSOFF_H


/*DEFINE HEADERS*/
#include "Neighbor_List.h"
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class POT_TERSOFF: public POTENTIAL {
 public:
  double *cutoff;             //Cutoff radius
	double *A, *B, *lambda;	    //Tersoff constants
	double *mu, *beta, *n;      //
	double *c, *d, *h, *R;      //
	POT_TERSOFF(int, int*);     //Constructor (Allocates memory)
	~POT_TERSOFF();      		    //Destructor (Deallocates memory)
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};

#endif
