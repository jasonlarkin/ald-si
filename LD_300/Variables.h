/*                            Variables.h                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Header file for the VARIABLES class.  This class stores all     *
*  varaibles defined in the input file under the keyword category    *
*  'VARIABLES'.                                                      *
*********************************************************************/

#if !defined(VARIABLES_H)
#define VARIABLES_H


/*DECLARE HEADERS*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;


/*DEFINE PREPROCESSOR VARIABLES*/


/*DEFINE CLASSES*/
class VARIABLES {
public:
	string Define(ifstream &Input);
	double GetValue(string str);
	VARIABLES();
	~VARIABLES();
private:
  string symb;
  double value;
  VARIABLES *next;
};


class MATH_EXPRESSION {
public:
  MATH_EXPRESSION(const char *c_str);
  ~MATH_EXPRESSION();
  double number;
  int op;
  MATH_EXPRESSION *next;
  void Evaluate(int op1, int op2);
  void Copy();
};


#endif
