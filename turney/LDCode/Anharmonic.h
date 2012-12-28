/*						               Anharmonic.h		                			*/
/*							              06/06/2009	              						*/
/*********************************************************************
*    Header file for the class ANHARMONIC which stores information  *
*  pertaining to the frequency shift and linewidth calculations.     *
*********************************************************************/

#ifndef ANHARMONIC_H_INCLUDED
#define ANHARMONIC_H_INCLUDED


/*DECLARE HEADERS*/
#include "LDCode.h"
#include "Transform.h"


/*DEFINE PREPROCESSOR COMMANDS*/


/*DEFINE CLASSES*/
class ANHARMONIC {
public:
	ANHARMONIC();              //Constructor (Builds neighbor lists)
	~ANHARMONIC();		          //Destructor (Deallocates memory)
	string ReadInput(ifstream &Input);
	void Initialize(PD_ANHARMONIC**);
	void Compute(int iter);
	void Output();        //Writes the variables
	//
	int iterations;
	int continue_;
	double Temperature;
	double fs_guess;
	double lw_guess;
	string Q_Interpolate;
	string C_Interpolate;
	//
	int *UC;
	int *lower;
	int *upper;
	int Sym_F;
	int Sym_E;
	//
	double *QH_f;
	double *AH_fQ;
	double *AH_fC;
	double *AH_lwQ;
	double *AH_lwC;
	double ***E;
	TRANSFORM *Transform;
	ofstream S_W_Out;
	ofstream S_W_Out_C;
private:
  int *wv_start;
  int  nu_start;
  void LDGuess(void);
  void Diagnostics(void);
  void PV_DV(double *F, double n1, double n2, double e, double& PV, double& DV);
  void Write(int n, double F_qh, double D3Q, double D3C, double D4Q, double D4C,
    double GQ, double GC, ofstream &Quant, ofstream &Class);
  void Guess(ofstream &Out);
  void Interpolate(string In_str, ofstream &Out);
  void Progress(int &mntr, int outer_F, int begin);
};


#endif // ANHARMONIC_H_INCLUDED
