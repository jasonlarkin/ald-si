/*                            Transform.h                           */
/*                            12/10/2008                            */
/*********************************************************************
*    Header file and code for the TRANSFORM class object.  TRANSFORM *
*  computes and stores intermediate steps and the final transform	   *
*  in the three and four phonon interaction expressions.  G4[][9]	   *
*  and G3[][9] are similar to each other.  They are to be computed	 *
*  after all wave vectors are specified and are defined as:		    	 *
*  G3[jj][DOF*a'+a"] = \sum_{i,a}\sum_{l',l"}e(kj|ia)				         *
*   *\Phi_{a,a',a"}(0i,l'i',l"i")exp[i(k'.r(l'0)+k".r(l"0))],	    	 *
*  G4[jj][DOF*a"+a~] = \sum_{i,a}\sum_{i',a'}\sum_{l',l",l~}e(kj|ia) *
*   *e(-kj|i'a')\Phi_{a,a',a",a~}(0i,l'i',l"i",l~i~)exp[-ik.r(l'0)]	 *
*   *exp[ik'.(r(l"0)-r(l~0))],									                  	 *
*  where jj = i'*neighbors+i".  H3[natom][3] is computed after the	 *
*  second eigenvector is specified (after j_in) and is defined by:	 *
*  H3[i"][a"] = \sum_{i',a'}G4[jj][DOF*a'+a"]e(k'j'|i'a').	     		 *
*  T3 and T4 are the final step and return the transformed values.	 *
*********************************************************************/

#if !defined(TRANSFORM_H)
#define TRANSFORM_H


/*DEFINE GLOBAL HEADERS*/
#include "LDCode.h"


/*DEFINE GLOBAL CLASSES*/
class TRANSFORM {
 public:
	TRANSFORM(PD_ANHARMONIC **FC);		        //Constructor (Allocates memory and defines listG and pointG)
	~TRANSFORM(void);		        //Destructor (Dealocates memory)
	void TransformG3(double **E, double **k);//Function that computes G[][3x3]
	void TransformG4(double **E, double **k);//Function that computes G4[][3x3]
	void TransformH3(double **E);//Function to compute H[natom][3]
	double TransformT3(double **E);//Function to compute the total transform 3
	double TransformT4(double **E);//Function to compute the total transform 4
	double Transform_3_0(double ***E, double **k);//Zone center optical modes
 protected:
	int *pointG;			          //Points to region of listG[]
	int *listG;				          //List of atom numbers which are pairs to atom 0
	double **G3_re, **G3_im;    //G3[G][3x3], stores intermediate value of potential transform 3
	double **G4_re, **G4_im;    //G4[G][3x3], stores intermediate value of potential transform 4
	double **H3_re, **H3_im;    //H[natom][3], stores data used directly in transform 3
	double pi2;                 //2.0*Pi
	PD_ANHARMONIC **FC;
	int Location(int, int);
};


#endif
