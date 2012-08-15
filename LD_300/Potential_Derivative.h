/*                      Potential_Derivative.h                      */
/*                            06/06/2009                            */
/*********************************************************************
*    Header file for the POT_DER classes.                            *
*********************************************************************/

#ifndef POTENTIAL_DERIVATIVE_H_INCLUDED
#define POTENTIAL_DERIVATIVE_H_INCLUDED



/*DECLARE HEADERS*/
#include "LDCode.h"


/*DEFINE PREPROCESSOR VARIABLES*/


/*DEFINE CLASSES*/
class POTENTIAL;


/*LIST OF THE DERIVATIVES OF THE POTENTIALS*/
class POT_DER_LINK_LIST {
public:
  //VARIABLES
	double *fc;                 //Array of force constants
	POT_DER_LINK_LIST *next;    //Pointer to next entry in list
	//FUNCTIONS
	POT_DER_LINK_LIST(int);     //Constructor (Allocates memory)
	~POT_DER_LINK_LIST();		    //Destructor (Deallocates memory)
	void Compute(int, double**, double**, POTENTIAL*);//Computes force constants, returns true if all zero
	bool Scale(int, double);    //Scale and zero small values
	void Print(int);            //Outputs data to stdout
};


/*LIST OF IDENTIFIERS*/
class POT_DER_IDENTIFIER {
 public:
  int *b;                     //Array of atoms in derivative (ordered, primary)
  int **l;                    //Array of unit cells in derivative (ordered, secondary)
  POT_DER_IDENTIFIER *next;   //Pointer to next entry in list
  POT_DER_IDENTIFIER(int, int*, int**);//Constructor
  ~POT_DER_IDENTIFIER();      //Destructor
  int Position(int, int*, int**, bool&);//
//  void Print(int);            //Prints identifiers
};


class POT_DER_X {
 public:
  //VARIABLES
  int b0;
  int *N;
  POT_DER_IDENTIFIER **id;
  //FUNCTIONS
  POT_DER_X(int, int);        //Constructor
  ~POT_DER_X();               //Destructor
	virtual void Compute(int, int, int*, int**, double**, POTENTIAL*)=0;//
	int Position(int, int*, int**, bool&);//Determines position to insert new entry
	POT_DER_LINK_LIST *Allocate(bool, int, int, POT_DER_LINK_LIST**);//Allocates memory and returns position
	void Remove(POT_DER_LINK_LIST *pdll);//Removes entry from list
  virtual void Send()=0;      //Divides by mass and sends to other PEs
  virtual void Recv(int)=0;   //Receives from other PEs
  virtual void Print(int)=0;  //Outputs data to stdout
};


class PD_ENERGY_FORCE : public POT_DER_X {
 public:
  //VARIABLES
  POT_DER_LINK_LIST *F0;      //Energy
  POT_DER_LINK_LIST *F1_0;    //atom 0 (Forces)
	//FUNCTIONS
	PD_ENERGY_FORCE(int b0);    //Constructor (Initializes to NULL)
	~PD_ENERGY_FORCE();		      //Destructor (Deallocates memory)
	void Compute(int, int, int*, int**, double**, POTENTIAL*);//
  void Send();                //Divides by mass and sends to other PEs
  void Recv(int);             //Receives from other PEs
  void Print(int);            //Outputs data to stdout
};


class PD_HARMONIC : public POT_DER_X {
 public:
  //VARIABLES
  POT_DER_LINK_LIST *F2_00;   //atom 0, atom 0 (self term)
  POT_DER_LINK_LIST *F2_01;   //atom 0, atom 1
	//FUNCTIONS
	PD_HARMONIC(int b0);        //Constructor (Initializes to NULL)
	~PD_HARMONIC();		          //Destructor (Deallocates memory)
	void Compute(int, int, int*, int**, double**, POTENTIAL*);//
	void Send();                //Divides by mass and sends to other PEs
  void Recv(int);             //Receives from other PEs
  void Print(int);            //Outputs data to stdout
};


class PD_ANHARMONIC : public POT_DER_X {
 public:
  //VARIABLES
  /* Third Derivatives of Potential Energy */
  POT_DER_LINK_LIST *F3_000;  //atom 0, atom 0, atom 0 (self term)
  POT_DER_LINK_LIST *F3_001;  //atom 0, atom 0, atom 1
  POT_DER_LINK_LIST *F3_011;  //atom 0, atom 1, atom 1
  POT_DER_LINK_LIST *F3_012;  //atom 0, atom 1, atom 2
  /* Fourth Derivatives of Potential Energy */
  POT_DER_LINK_LIST *F4_0000; //atom 0, atom 0, atom 0, atom 0 (self term)
  POT_DER_LINK_LIST *F4_0001; //atom 0, atom 0, atom 0, atom 1
  POT_DER_LINK_LIST *F4_0011; //atom 0, atom 0, atom 1, atom 1
  POT_DER_LINK_LIST *F4_0111; //atom 0, atom 1, atom 1, atom 1
  POT_DER_LINK_LIST *F4_0012; //atom 0, atom 0, atom 1, atom 2
  POT_DER_LINK_LIST *F4_0112; //atom 0, atom 1, atom 1, atom 2
  POT_DER_LINK_LIST *F4_0122; //atom 0, atom 1, atom 2, atom 2
  POT_DER_LINK_LIST *F4_0123; //atom 0, atom 1, atom 2, atom 3
	//FUNCTIONS
	PD_ANHARMONIC(int b0);    //Constructor (Initializes to NULL)
	~PD_ANHARMONIC();		      //Destructor (Deallocates memory)
	void Compute(int, int, int*, int**, double**, POTENTIAL*);//
	void Send();                //Divides by mass and sends to other PEs
  void Recv(int);             //Receives from other PEs
  void Print(int);            //Outputs data to stdout
};


#endif // POTENTIAL_DERIVATIVE_H_INCLUDED
