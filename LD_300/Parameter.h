/*                            Parameter.h                           */
/*                            12/02/2008                            */
/*********************************************************************
*    Header file for the PARAMETER class.  This class calculates     *
*  appropriate (non-)dimensionalizing parameters and provides        *
*  functions for outputing and scaling the input data.               *
*********************************************************************/

#ifndef PARAMETER_H
#define PARAMETER_H


/*DECLARE HEADERS*/
#include "Potential.h"
#include "Dispersion.h"
#include "Anharmonic.h"


/*DEFINE PREPROCESSOR VARIABLES*/


/*DEFINE CLASSES*/
class PARAMETER {
public:
	double energy;		//Dimensional energy scale (J)
	double length;			//Dimensional length scale (m)
	double mass;			//Dimensional mass scale (kg)
	double k_B;				//Boltzmann constant (J/K)
	double temp;      //Dimensional temperature scale factor (K)
	double time;			//Dimensional time scale factor (s)
	double pi;				//Pi (3.14159...)
	PARAMETER(POTENTIAL**);//Constructor
	~PARAMETER();     //Destructor
	void Output_Scale(int scale, DISPERSION *Dispersion, ANHARMONIC *Anharmonic);//Outputs and scales data
};


#endif // PARAMETER_H
