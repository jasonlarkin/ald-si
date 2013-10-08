/*                            Material.h                            */
/*                            12/02/2008                            */
/*********************************************************************
*    Header file for the MATERIAL class.  This class reads and       *
*  temporarily stores all material symbols defined in the input      *
*  file.                                                             *
*********************************************************************/

#ifndef MATERIAL_H
#define MATERIAL_H


/*DECLARE HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;
#include "Potential.h"


/*DEFINE PREPROCESSOR VARIABLES*/


/*DEFINE CLASSES*/
class MATERIAL {
public:
  MATERIAL();
	~MATERIAL();
	string Define(ifstream &Input, POTENTIAL ***Pot);
	int GetMaterial(string str);
	void Output();        //Writes the variables
	void Initialize();    //Initializes variables after they have been defined (non-dimensionalizes)
	int n_mat;
  string *symb;
  double *mass;
};


#endif
