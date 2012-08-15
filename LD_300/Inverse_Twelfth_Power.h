/*						          Inverse_Twelfth_Power.h	            				*/
/*							              12/02/2008		              					*/
/*********************************************************************
*    Header file that contains the Inverse-Twelfth-Power derivative	 *
*  class definition for the Inverse-Twelfth-Power potential.	    	 *
*********************************************************************/

#if !defined(INVERSE_TWELFTH_POWER_H)
#define INVERSE_TWELFTH_POWER_H


/*DEFINE HEADERS*/
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class POT_INV12: public POTENTIAL {
 public:
  double *epsilon, *sigma, *cutoff, *Acutoff;//Potential parameters
	POT_INV12(int, int*);       //Constructor (Assigns values)
	~POT_INV12();   		        //Destructor (Default)
	string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


#endif
