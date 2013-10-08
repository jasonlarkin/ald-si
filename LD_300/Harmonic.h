/*                            Harmonic.h                            */
/*                            06/05/2009                            */
/*********************************************************************
*    Header file for the HARMONIC class.  This class is responsible  *
*  for performing harmonic LD calculations.                          *
*********************************************************************/

#ifndef HARMONIC_H
#define HARMONIC_H


/*DECLARE HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class HARMONIC {
public:
  bool master_flag;           //Turns all harmonic calculations on or off
  bool thermo;                //Turns computation of thermodynamic properties on or off (not used yet)
	HARMONIC();                 //Constructor
	~HARMONIC();                //Destructor
	string ReadInput(ifstream &Input);//
	bool Boolean(string str);   //
	void Compute(int &anh_iter, PD_HARMONIC **FC);//Performs harmonic LD calculation
	void Output();        //Writes the variables
	bool Initialize();    //Initializes variables after they have been defined (non-dimensionalizes)
};


#endif // HARMONIC
