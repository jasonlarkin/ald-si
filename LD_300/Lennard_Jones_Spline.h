/*					             Lennard_Jones_Spline.h		  	          		*/
/*						              	05/04/2009  					            		*/
/*********************************************************************
*    Header file that contains the LENNARD_JONES_SPLINE derivative   *
*  class definition for the Lennard-Jones potential.						     *
*********************************************************************/

#if !defined(LENNARD_JONES_SPLINE_H)
#define LENNARD_JONES_SPLINE_H


/*DEFINE HEADERS*/
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class POT_LJ_SPLINE : public POTENTIAL {
 public:
  double *rspline;
  double *epsilon, *sigma, *cutoff, *Acutoff;//LJ potential parameters
	POT_LJ_SPLINE(int, int*);   //Constructor (Allocates memory)
	~POT_LJ_SPLINE();           //Destructor (Deallocates memory)
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


#endif
