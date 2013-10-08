/*                           PreProcess.h                           */
/*                            06/30/2009                            */
/*********************************************************************
*    Header file for the PREPROCESS class.  This class stores all    *
*  functions and parameters needed to define the pre-processing      *
*  routines.                                                         *
*********************************************************************/

#ifndef PREPROCESS_H_INCLUDED
#define PREPROCESS_H_INCLUDED


/*DECLARE HEADERS*/
#include <cstring>


/*DECLARE GLOBAL VARIABLES*/


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class PREPROCESS {
 public:
  bool master_flag;           //Master flag for turning PREPROCESS on/off
  int opt_itr;                //Max number of iterations for optimization
  bool struct_flag;           //Flag for outputing structure
  bool ef_flag;               //Flag for outputing energies and forces
  PREPROCESS();               //Constructor sets variables to defaults
  ~PREPROCESS();              //Destructor deallocates memory
  string ReadInput(ifstream&);//Reads and stores parameters
  void Output();              //Writes the variables
	bool Initialize();          //Initializes variables after they have been defined (non-dimensionalizes)
  void Optimize(POTENTIAL **Pot);//Optimizes the positions of the atoms in the unit cell via simple gradient method
  void Structure();           //Outputs the sturcture of the crystal
  void Energy_Force(POTENTIAL**);//Computes energies and forces on atoms
};


#endif // PREPROCESS_H_INCLUDED
