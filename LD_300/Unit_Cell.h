/*                            Unit_Cell.h                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Header file for the UNIT_CELL class.  This class reads and      *
*  stores the atomic positions of the atoms in the unit cell.        *
*********************************************************************/


#ifndef UNIT_CELL_H
#define UNIT_CELL_H


/*DECLARE HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;


/*DEFINE PREPROCESSOR VARIABLES*/


/*DEFINE CLASSES*/
class UNIT_CELL {
public:
	UNIT_CELL();
	~UNIT_CELL();
	string Define(ifstream &Input);
	int natom;
	int *mat;
	double **X;
	int KeywordToInt(ifstream &Input);
	void Output();        //Writes the variables
	void Initialize();     //Initalizes variables after they have been defined (non-dimensionalizes)
private:
  int coordinate;
};


#endif // UNIT_CELL_H
