/*						                  REBO.h		                   				*/
/*							              07/08/2009              							*/
/*********************************************************************
*    Header file that contains the POT_REBO derivative class		   	 *
*  definition for the REBO N-body potential.	                    	 *
*********************************************************************/

#if !defined(REBO_H)
#define REBO_H


/*DEFINE HEADERS*/
#include "Neighbor_List.h"
#include "LDCode.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
//E = 1/2 * sum[i]{sum[j]{ A(1-Q/r[ij])exp(-alpha*r[ij]) }}
class POT_REBO2: public POTENTIAL {
 public:
  double *cutoff;             //Cutoff distance
	double *Q, *A, *alpha;	    //REBO Two body constants
	POT_REBO2(int, int*);       //Constructor
	~POT_REBO2();    			      //Destructor
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


//E = 1/2 * sum[i]{sum[j]{ -b[ij](B1*exp(-beta1*r[ij])+B2*exp(-beta2*r[ij])+B3*exp(-beta3*r[ij]))}}
//b[ij] = 1/2 * (b[ijk1,ijk2]+b[ijl1,ijl2]) + b_DH
//b[ijx,ijy] = (G[theta[ijx] + G[theta[ijy]])^-0.5
//G[theta] = Q1*cos[theta]^2 + Q2*cos[theta] + Q3
//b_DH = 0.0
class POT_REBO3: public POTENTIAL {
 public:
  double *cutoff;             //Cutoff distance
	double *B1, *B2, *B3;         //REBO constants
	double *beta1, *beta2, *beta3;//REBO constants
	double *Q1, *Q2, *Q3;         //REBO angle constants
	POT_REBO3(int, int*);       //Constructor
	~POT_REBO3();   			      //Destructor
  string ReadInput(ifstream&);//Reads and stores parameters
  void GetScale(double&, double&);//Returns the energy and length scales
  void Initialize();          //Initializes parameters for calculations
  POTENTIAL *Copy();          //Returns new pointer containing a copy of the potential data
  void Print();               //Prints potential data
  double Energy(double**);    //Energy of interaction
  void ForceConstants(int, POT_DER_X*);//Force constants
};


#endif
