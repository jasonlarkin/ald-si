/*                           Dispersion.h                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Header file for the DISPERSION class.  This class reads and     *
*  stores the information required to output the dispersion relation *
*  for an arbitrary number of points along any given direction.      *
*********************************************************************/

#ifndef DISPERSION_H
#define DISPERSION_H


/*DECLARE HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;


/*DEFINE PREPROCESSOR VARIABLES*/


/*DEFINE CLASSES*/
class DISPERSION {
public:
	DISPERSION();
	~DISPERSION();
	string Define(ifstream &Input);
	int KeywordToInt(ifstream &Input);
	void Compute(PD_HARMONIC **FC);
	void Output();        //Writes the variables
	bool Initialize();    //Initializes variables after they have been defined (non-dimensionalizes)
	int points;
	double *begin;
	double *end;
private:
  int coordinate;
};


#endif // DISPERSION_H
