/*                           Parameter.cpp                          */
/*                            12/02/2008                            */
/*********************************************************************
*    Source code file for functions in the PARAMETER class.          *
*********************************************************************/


/*DECLARE HEADERS*/
#include "LDCode.h"
#include <cmath>


/*DEFINE PREPROCESSOR COMMANDS*/


/*DECLARE GLOBAL VARIABLES*/


/*PARAMETER CONSTRUCTOR: DEFINES (NON-)DIMENSIONALIZING
  VALUES AS AVERAGE OF ALL POTENTIAL PARAMETERS*/
PARAMETER::PARAMETER(POTENTIAL **Pot) {
  //DECLARE LOCAL VARIABLES
  int i, j;
  int n_mat=Material->n_mat;
  int n_sum = 0;
  double e_sum=0.0;
  double l_sum=0.0;
	double m_sum=0.0;
	double e, l;                //Energy and length scales for each potentail
	POTENTIAL *P_ptr;
	//COMPUTE PARAMETER SUMS
  for(i=0;i<n_mat;i++) {
    m_sum += Material->mass[i];
    P_ptr = Pot[i];
    while(P_ptr!=NULL) {
      P_ptr->GetScale(e, l);
      n_sum += 1;
      e_sum += e;
	    l_sum += l;
      P_ptr = P_ptr->next;
    }
  }
	//ASSIGN VALUES TO PARAMETERS
	mass = m_sum/double(n_mat);
	energy = e_sum/double(n_sum);
	length = l_sum/double(n_sum);
	k_B = 1.3806e-23;           // J/K
	temp = energy/k_B;
	time = sqrt(mass / energy) * length;
	pi = 4.0*atan(1.0);
	return;
}


/*PARAMETER DESTRUCTOR: DOES NOTHING*/
PARAMETER::~PARAMETER() {return;}
