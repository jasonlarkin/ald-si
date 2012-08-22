/*						            Stillinger_Weber.h		            				*/
/*							              12/02/2008              							*/
/*********************************************************************
*    Header file that contains the POT_SW derivative class		    	 *
*  definition for the Stillinger-Weber three body potential.	    	 *
*********************************************************************/

#if !defined(STILLINGER_WEBER_H)
#define STILLINGER_WEBER_H


/*DEFINE HEADERS*/
#include "Neighbor_List.h"
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class POT_SW2: public POTENTIAL {
 public:
  double *epsilon, *sigma;    //Energy and length scales
	double *A, *B, *p, *q, *a;  //SW Two body constants
	POT_SW2(int, int*);         //Constructor
	~POT_SW2();    			        //Destructor
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


class POT_SW3: public POTENTIAL {
 public:
  double **epsilon, **lambda; //Atoms 1, 2, 3 at the center
	double *sigma, *gamma, *a;  //SW Three body constants
	POT_SW3(int, int*);         //Constructor
	~POT_SW3();   			        //Destructor
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


#endif
