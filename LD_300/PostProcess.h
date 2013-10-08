/*                           PostProcess.h                          */
/*                            06/20/2009                            */
/*********************************************************************
*    Header file for the POSTPROCESS class.  This class stores all   *
*  functions and parameters needed to define the post processing     *
*  routines.                                                         *
*********************************************************************/

#ifndef POSTPROCESS_H_INCLUDED
#define POSTPROCESS_H_INCLUDED


/*DECLARE HEADERS*/
#include <cstring>


/*DECLARE GLOBAL VARIABLES*/


/*DEFINE PREPROCESSOR COMMANDS*/
#define DATA_NONE   0
#define DATA_NO_SYM 1
#define DATA_SYM    2


/*DEFINE CLASSES*/
class POSTPROCESS {
 public:
  bool master_flag;           //Master flag for turning POSTPROCESS on/off
  int TC_beg;                 //Beginning iteration for thermal conductivity calculation
  int TC_end;                 //Ending iteration for thermal conductivity calculation
  int data_itr;               //Iteration number for extracting data
  int *data_beg;              //Beginning wave vector for extracting data along arbitrary direction
  int *data_inc;              //Vector to increment when extracting data along arbitrary direction
  int data_num;               //Number of times to increment data along arbitrary direction
  int data_fBZ;               //Flag for outputing the data for the full BZ
  POSTPROCESS();              //Constructor sets variables to defaults
  ~POSTPROCESS();             //Destructor deallocates memory
  string ReadInput(ifstream&);//Reads and stores parameters
  void Output();              //Writes the variables
	void Initialize();          //Initializes variables after they have been defined (non-dimensionalizes)
  void Data_Dir();            //Outputs data along arbitary direction
  void Data_FBZ();            //Outputs data for the full BZ
};

#endif // POSTPROCESS_H_INCLUDED
