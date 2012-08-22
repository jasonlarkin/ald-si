/*						            QC_Temperature.cpp	             					*/
/*							              12/30/2008		              					*/
/*********************************************************************
*    This subroutine produces a table of classical versus quantum    *
*  temperautres due to quantum corrections (assuming harmonic        *
*  systems).  It also outputs the                                                        *
*********************************************************************/


/*DECLARE HEADERS*/
#include "LDCode.h"
#include <cmath>


/*DECLARE PREPROCESSOR DEFINITIONS*/


/*DECLARE GLOBAL VARIABLES*/


/*MAIN ROUTINE*/
int QC_Temperature(int argc, char **argv) {
	//DELCARE LOCAL VARIABLES
	int i, j;	                  //Counters
	int *wv_N=new int[DIM];			//Wave vector counters
	int *lower = Lattice->lower, *upper = Lattice->upper;	//Lower and upper bounds of 1st Brillioun zone
	double T, x, n, C_ph;       //Temp, generic variable, distribution, and heat capacity
	double factor;				      //Multiplicity factor for symmetry conditions
	double E, E_0, E_C, dT=0.0; //Total quantum, zero-point, and classical energies and dT^C/dT^Q
	double k_B = 1.3806e-23;    //Boltzmann's constant [J/K]
  double h_bar = 1.05457e-34; //Plank's constant [J-s]


	//INITIALIZE FILES AND DATA
	ofstream Output("QC_Temperature.xls");
	Output <<"T^Q\tT^C w/o ZP\tT^C w/ ZP\tdT^C/dT^Q"<<endl;
	ifstream Frequency("Frequency.txt");//Input file for frequencies
	if(!Frequency.is_open()) {Log <<"Could not open 'Freqeuncy.txt'"<<endl;exit(0);}
	Frequency.ignore(10000, '\n');		//Skip header line
	double *QH_F = new double[Symmetry->Sym_F];
	for(i=0;i<Symmetry->Sym_F;i++) {
	  Frequency >>QH_F[i];
	  QH_F[i] *= 1.0e12;
	}
	Frequency.close();


	//CALCULATE QAUNTUM CORRECTION TO TEMPERATURE
	cout <<"\n\nBeginning Quantum Correction Calculation for Temperature.\n\n";


	//CALCULATE CLASSICAL AND ZERO-POINT ENERGIES
	E_C = k_B*double(Symmetry->TotalDOF - DIM);//[J/K]
  E_0 = 0.0;
  i = 0;
  for(wv_N[0]=upper[0];wv_N[0]>=lower[0];wv_N[0]--) {
  for(wv_N[1]=upper[1];wv_N[1]>=lower[1];wv_N[1]--) {
  for(wv_N[2]=upper[2];wv_N[2]>=lower[2];wv_N[2]--) {
    //Perform symmetry operations
    if((factor=Symmetry->Sym_Negative(wv_N))<0.0) {continue;}
    Symmetry->Sym_Freq(wv_N, &i);
    for(j=0;j<UC_DOF;j++) {
      //Accumulate The Quantum Energy
      E_0 += factor * QH_F[i];
      i += 1;
    }  //for(j)
  }  //for(wv_N[2])
  }  //for(wv_N[1])
  }  //for(wv_N[0])
  E_0 *= 0.5*h_bar;           //[J]
  Output <<0.0<<'\t'<<0.0<<'\t'<<E_0/E_C<<'\t'<<0.0<<endl;



	//CALCULATE QUANTUM ENERGIES AT EACH TEMP
	i = 0;
	for(T=1.0;T<=1000.0;T+=1.0) {
	  E = 0.0;
    for(wv_N[0]=upper[0];wv_N[0]>=lower[0];wv_N[0]--) {
    for(wv_N[1]=upper[1];wv_N[1]>=lower[1];wv_N[1]--) {
    for(wv_N[2]=upper[2];wv_N[2]>=lower[2];wv_N[2]--) {

      //Perform symmetry operations
      if((factor=Symmetry->Sym_Negative(wv_N))<0.0) {continue;}
      Symmetry->Sym_Freq(wv_N, &i);

      for(j=0;j<UC_DOF;j++) {
        //Accumulate the quantum energy
        if(QH_F[i]<1.0e-6) {i+=1;continue;}
        x = h_bar*QH_F[i]/(k_B*T);
        if((C_ph=exp(x))<1.0e200) {
          n = 1.0/(C_ph-1.0);
          E += factor * QH_F[i] * n;
          C_ph = x*x*C_ph*n*n;
          dT += factor * C_ph;
          //Log <<i<<" "<<QH_F[i]<<" "<<n<<" "<<E<<endl;
        }
        i += 1;
      }  //for(j)
    }  //for(wv_N[2])
    }  //for(wv_N[1])
    }  //for(wv_N[0])


    //Output data
    E *= h_bar;
    dT *= k_B;
    Output <<T<<'\t'<<E/E_C<<'\t'<<(E+E_0)/E_C<<'\t'<<dT/E_C<<endl;

	}  //for(T)
  Output.close();

	return(0);
}
