/*                            Potential.h                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Header file for the POTENITAL class.  This class stores all     *
*  functions and parameters needed to define the interatomic         *
*  potentials.                                                       *
*********************************************************************/

#ifndef POTENTIAL_H
#define POTENTIAL_H


/*DECLARE HEADERS*/
#include <cstring>
#include <typeinfo>
#include "Potential_Derivative.h"


/*DECLARE GLOBAL VARIABLES*/


/*DEFINE PREPROCESSOR COMMANDS*/

/*Parameters entered as AA AB AC ... BB BC ... CC ...
For parameter that do not fit the pattern you can define multiple
parameters and optionaly enter zeros*/
/*DEFINE CLASSES*/
class POTENTIAL {
 public:
  //VARAIBLES
  int n_mat;                  //Number of unique materials involed
  int *mat;                   //List of unique materials
  int *mat_pos;               //Stores material's position in mat[]
  int n_expect;               //Number of values to expect for each parameter
//  int n_consts;               //Number of parameters for potential
//  double **par;               //Values of potential parameters
//  char **consts;              //Names of potential parameters
  char *Pot_Symb;             //Potential symbol
  POTENTIAL *next;            //Pointer to next element in linked list
  //FUNCTIONS (all but Num_Der and FillPrameter need defined for each potential)
  POTENTIAL(int, int*);       //Constructor assigns materials
  virtual ~POTENTIAL();       //Destructor deallocates memory
  virtual string ReadInput(ifstream&)=0;//Reads and stores parameters
  virtual void GetScale(double&, double&)=0;//Returns the energy and length scales
  virtual void Initialize()=0;//Initializes parameters for calculations
  virtual POTENTIAL *Copy()=0;//Returns new pointer containing a copy of the potential data
  virtual void Print()=0;     //Prints potential data
  virtual double Energy(double**)=0;//Energy of interaction
  virtual void ForceConstants(int, POT_DER_X*)=0;//Force constants
  double Num_Der(int, double**, double**);//Numerical derivative of potential
  void Fill(double*, ifstream&);//Assigns values read to parameter
};

#endif
