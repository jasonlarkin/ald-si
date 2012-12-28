/*							            Anharmonic.cpp	  				          		*/
/*							              06/06/2009			              				*/
/*********************************************************************
*    Subroutine that computes the phonon frequency shift and		     *
*  linewidth due to third and fourth order anharmonic effects.	  	 *
*********************************************************************/


/*PREPROCESSOR DEFINITIONS*/


/*DEFINE HEADERS*/
#include "LDCode.h"
#include "Transform.h"
#include "Anharmonic.h"
#include <cmath>
#include <ctime>


/*SHIFT_WIDTH CONSTRUCTOR: INITIALIZES VARIABLES*/
ANHARMONIC::ANHARMONIC(void) {
  iterations = 0;
	continue_ = 0;
	wv_start = NULL;
	nu_start = 0;
	Temperature = 300;
	fs_guess = 0.005;
	lw_guess = 0.005;
	return;
}


/*SHIFT_WIDTH FUNCTION Define: DEFINES VARIABLES*/
string ANHARMONIC::ReadInput(ifstream &Input) {
  //DECLARE LOCAL VARIABLES
  int i;
  string str;
  string Read_Next(ifstream &Input);


  //PARSE THROUGH FILE
  str = Read_Next(Input);
  iterations = int(atof(str.c_str())+0.0001);
  while(!Input.eof()) {
    //identify keyword in ANHARMONIC category
    str = Read_Next(Input);
    i = 0;
    while(str[i]) {str[i]=tolower(str[i]);i++;}
    if     (str.compare(0, 2, "fs_guess",     2)==0) {fs_guess = atof(Read_Next(Input).c_str());}
    else if(str.compare(0, 2, "lw_guess",     2)==0) {lw_guess = atof(Read_Next(Input).c_str());}
    else if(str.compare(0, 1, "temperature",  1)==0) {Temperature = atof(Read_Next(Input).c_str());}
    else if(str.compare(0, 4, "q_interpolate",4)==0) {Q_Interpolate = Read_Next(Input);}
    else if(str.compare(0, 4, "c_interpolate",4)==0) {C_Interpolate = Read_Next(Input);}
    else if(str.compare(0, 4, "continue",     4)==0) {continue_ = int(atof(Read_Next(Input).c_str())+0.001);}
    else if(str.compare(0, 4, "conductivity", 4)==0) {
      str = Read_Next(Input);
      Log <<"WARNING: The 'conductivity' flag is no longer used.  The thermal conductivity is always computed."<<endl;
    }
    else{break;}
  }
  return(str);
}


/*ANHARMONIC FUNCTION Output: WRITES DATA TO LOG FILE*/
void ANHARMONIC::Output() {
  int i;
  Log <<"\nANHARMONIC = "<<iterations<<endl;
  Log <<"  continue     = "<<continue_<<endl;
  Log <<"  temperature  = "<<Temperature<<endl;
  if( (Q_Interpolate.size()==0)||(C_Interpolate.size()==0) ) {
    Log <<"  fs_guess     = "<<fs_guess<<endl;
    Log <<"  lw_guess     = "<<lw_guess<<endl;
  }
  if(Q_Interpolate.size()!=0) {
    Log <<"  Q_Interpolate = "<<Q_Interpolate<<endl;
  }
  if(C_Interpolate.size()!=0) {
    Log <<"  C_Interpolate = "<<C_Interpolate<<endl;
  }
  return;
}


/*ANHARMONIC FUNCTION Initialize: INITIALIZES VARIABLES*/
void ANHARMONIC::Initialize(PD_ANHARMONIC **FC) {
  //DECLARE LOCAL VARIABLES
  int i, j;
  double per_ps = 1.0e12*Parameter->time;	//(Non-)dimensionalizing factor (1/ps)
  int *wv = new int[3];
  int *n = new int[3];


  //INITIALIZE VARIABLES
  UC = Lattice->N;
  lower = Lattice->lower;
  upper = Lattice->upper;
  Sym_E = Symmetry->Sym_E;
  Sym_F = Symmetry->Sym_F;
	Transform = new TRANSFORM(FC);


  //ALLOCATE MEMORY FOR AND READ PHONON PROPERTIES FROM FILES
  QH_f   = new double[Sym_F];	//Quasi-harmonic phonon frequencies
  AH_fQ  = new double[Sym_F];	//Quantum anharmonic frequencies
	AH_fC  = new double[Sym_F]; //Classical anharmonic frequencies
	AH_lwQ = new double[Sym_F]; //Quantum anharmonic linewidth
	AH_lwC = new double[Sym_F]; //Classical anharmonic linewidth
	//Quasi-harmonic frequencies
	ifstream Frequency("Frequency.txt");	//Input file for frequencies (THz)
	if(!Frequency.is_open()) {Log <<"Could not open file 'Frequency.txt'"<<endl;exit(0);}
	Frequency.ignore(10000, '\n');			//Skip first line
	for(i=0;i<Sym_F;i++) {
	  Frequency >>QH_f[i];
    QH_f[i] *= per_ps;		    //Non-dimensionalize frequencies
	}
	Frequency.close();
	//Quasi-harmonic eigenvectors
	ifstream Eigenvector("Eigenvector.txt");//Input file for mode shapes (polar)
	if(!Eigenvector.is_open()) {Log <<"Could not open file 'Eigenvector.txt'"<<endl;exit(0);}
	Eigenvector.ignore(10000, '\n');			//Skip first line
	E = new double** [Sym_E];		//Eigenvectors/Mode shapes in polar
	E[0] = new double* [Sym_E*UC_DOF];		  //form: [a+b*i]
	for(i=0;i<Sym_E;i++) {
		if(i<Sym_E-1) {E[i+1] = E[i] + UC_DOF;}
		E[i][0] = new double [UC_DOF*2];
		for(j=0;j<UC_DOF;j++) {
		  if(j<UC_DOF-1) {E[i][j+1] = E[i][j] + 2;}
		  Eigenvector >>E[i][j][0]>>E[i][j][1];
    }
	}
	Eigenvector.close();


	//CREATE INITIAL GUESS FOR FREQ SHIFT AND LINEWIDTH
	if(PE==0) {
    if(continue_==0) {LDGuess();}
    else {Diagnostics();}
	}
#if defined PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);	//Ensure nodes do not get too far out of sync
#endif

  return;
}


/*SUBROUTINE Compute COMPUTES THE FREQUENCY SHIFT AND LINEWIDTH*/
void ANHARMONIC::Compute(int iter) {
	//DECLARE LOCAL VARIABLES
	int begin = time(0);	      //Begin time
	int i, mntr;			          //Generic counters
  int **wv_N;                 //Wave vector counters for outer, middle, and inner loops
	double **wv;		            //Wave vector (wv_N/N)
	int *nu = new int[3];		    //Dispersion branch counters for outer, middle, and inner loops
	int line_F, line_E;		      //Position of third wave vector in arrays
	int inner_F, inner_E, outer_F, outer_E;//Counters for inner and outer wave vector loops
	double per_ps, DV, PV, DV_C, PV_C;//Dimensionalizing factor (THz) or (1/s), Delta value
	double B_E;				          //Bose-Einstein constant factor B_E = freq*h_bar/(k_B*T)
	double n1, n2, e, T, T_C;	  //Bose-Einstein distributions, epsilon, Fourier transforms
	double Gamma, Delta3, Delta4;//Total phonon mode frequency shift and linewidth
	double n1_C, n2_C, Gamma_C, Delta3_C, Delta4_C;//Classical values
	double factor3, factor4;    //Duplication factors
	bool sym;				            //True if wave vectors has known symmetry condition
	char filename[256];


	//INITIALIZE VARIABLES
	B_E = 7.63850231e-12;					//B_E = h_bar/(k_B*T) for use with non-
	B_E /= Parameter->time*Temperature*Parameter->temp;//dimensional frequency


	//ALLOCATE MEMORY
	//Wave vectors
	wv_N = new int*[3];
	wv_N[0] = new int[3*3];
	wv = new double*[3];
	wv[0] = new double[3*3];
	for(i=1;i<3;i++) {
	  wv_N[i] = wv_N[i-1] + 3;
	  wv[i] = wv[i-1] + 3;
	}
	//Anharmonic frequencies and linewidths
	sprintf(filename, "%s%i.txt", "Shift_Width", iter-1);
	ifstream IO_Quant(filename);
	if(!IO_Quant.is_open()) {Log <<"Could not open file '"<<filename<<"' for input."<<endl;exit(0);}
	sprintf(filename, "%s%i.txt", "Classical_Shift_Width", iter-1);
	ifstream IO_Class(filename);
	if(!IO_Class.is_open()) {Log <<"Could not open file '"<<filename<<"' for input."<<endl;exit(0);}
	IO_Quant.ignore(10000, '\n');
	IO_Class.ignore(10000, '\n');
	per_ps = 1.0e12*Parameter->time;	//(Non-)dimensionalizing factor (1/ps)
	for(i=0;i<Sym_F;i++) {
	  IO_Quant >>n1>>Delta3>>Delta4>>Gamma;
    IO_Class >>n1_C>>Delta3_C>>Delta4_C>>Gamma_C;
    if((AH_fQ[i]=(n1+Delta3+Delta4)*per_ps)<1.0e-8) {AH_fQ[i]=n1*per_ps;}
    if((AH_fC[i]=(n1_C+Delta3_C+Delta4_C)*per_ps)<1.0e-8) {AH_fC[i]=n1_C*per_ps;}
    if(Gamma<1.0e-8) {AH_lwQ[i] = 0.0;}
    else {AH_lwQ[i] = 0.5*per_ps/Gamma;}
    if(Gamma_C<1.0e-8) {AH_lwC[i] = 0.0;}
    else {AH_lwC[i] = 0.5*per_ps/Gamma_C;}
	}
	IO_Quant.close();
	IO_Class.close();
	//Storage for properties of phonons in current interation
	double* freq    = new double[3];
	double* ah_f    = new double[3];
	double* ah_f_C  = new double[3];
	double* ah_lw   = new double[3];
	double* ah_lw_C = new double[3];
	double*** mode  = new double**[3];
	mode[0] = new double*[3*UC_DOF];
	for(i=0;i<3;i++) {
		if(i<2) {mode[i+1] = mode[i] + UC_DOF;}
		mode[i][0] = new double [UC_DOF*2];
		for(nu[2]=0;nu[2]<UC_DOF-1;nu[2]++) {mode[i][nu[2]+1] = mode[i][nu[2]] + 2;}
	}


	//CALCULATE FREQUENCY SHIFT AND LINEWIDTH
	if(PE==0) {
	  if(wv_start==NULL) {
	    sprintf(filename, "%s%i.txt", "Shift_Width", iter);
	    S_W_Out.open(filename, ios::trunc);
	    sprintf(filename, "%s%i.txt", "Classical_Shift_Width", iter);
	    S_W_Out_C.open(filename, ios::trunc);
      S_W_Out   <<Sym_F<<" Out of "<<Symmetry->TotalDOF
        <<" QH, Third, Fourth Freq in rad/ps and Lifetimes in ps for a "<<UC[0];
      S_W_Out_C <<Sym_F<<" Out of "<<Symmetry->TotalDOF
        <<" QH, Third, Fourth Freq in rad/ps and Lifetimes in ps for a "<<UC[0];
      for(i=1;i<DIM;i++) {
        S_W_Out   <<"x"<<UC[i];
        S_W_Out_C <<"x"<<UC[i];
      }
      for(i=0;i<Material->n_mat;i++) {
        S_W_Out   <<" "<<Material->symb[i];
        S_W_Out_C <<" "<<Material->symb[i];
      }
      S_W_Out   <<" Lattice at "<<Temperature*Parameter->temp
        <<"K. (atoms/UC="<<Unit_Cell->natom<<")\n";
      S_W_Out_C <<" Lattice at "<<Temperature*Parameter->temp
        <<"K. (atoms/UC="<<Unit_Cell->natom<<")\n";
    }
    else {
	    Log  <<"  Continuing from previous unfinished calculation of iteration "<<iter<<".\n"<<endl;
	    cout <<"  Continuing from previous unfinished calculation of iteration "<<iter<<".\n"<<endl;
      sprintf(filename, "%s%i.txt", "Shift_Width", iter);
      S_W_Out.open(filename, ios::app);
      sprintf(filename, "%s%i.txt", "Classical_Shift_Width", iter);
      S_W_Out_C.open(filename, ios::app);
    }
    if( (!S_W_Out.is_open())||(!S_W_Out_C.is_open()) ) {
	    Log <<"Unable to open output files for anharmonic calculation."<<endl;
	    exit(0);
	  }
	}
	if(wv_start==NULL) {
    wv_start = new int[DIM];
    for(i=0;i<DIM;i++) {wv_start[i] = upper[i];}
    nu_start = 0;
	}
	per_ps = 1.0e-12*1.0546e-34/(16.0*double(Symmetry->TotalDOF/UC_DOF))
		/(Parameter->energy*Parameter->time*Parameter->time);
	mntr = 1;
	outer_F = PE;
	outer_E = PE;
	//Outer Loop
	for(int loop=PE;loop<Sym_F;loop+=nPE) {
	  nu[0] = loop%UC_DOF;
	  nu[1] = loop/UC_DOF;
	  for(i=0;i<DIM;i++) {
	    wv_N[0][i] = Symmetry->unique[nu[1]][i];
	    wv[0][i] = double(wv_N[0][i])/double(UC[i]);
    }
    //Read proper data from arrays
    inner_E = Symmetry->E_map[upper[0]-wv_N[0][0]][upper[1]-wv_N[0][1]][upper[2]-wv_N[0][2]] + nu[0];
    Delta3 = Delta4 = Gamma = 0.0;
    Delta3_C = Delta4_C = Gamma_C = 0.0;
    freq[0] = QH_f[loop];
    ah_f[0] = AH_fQ[loop];
    ah_f_C[0] = AH_fC[loop];
    ah_lw[0] = AH_lwQ[loop];
    ah_lw_C[0] = AH_lwC[loop];
    for(i=0;i<UC_DOF;i++) {
      mode[0][i][0] = E[inner_E][i][0];
      mode[0][i][1] = E[inner_E][i][1];
    }

/*    //Check if we are continuing from a partial calculation
    if(wv_N[0][0]>wv_start[0]) {Progress(mntr, loop, Sym_F, begin);continue;}
    else if(wv_N[0][0]==wv_start[0]) {
      if(wv_N[0][1]>wv_start[1]) {Progress(mntr, loop, Sym_F, begin);continue;}
      else if(wv_N[0][1]==wv_start[1]) {
        if(wv_N[0][2]>wv_start[2]) {Progress(mntr, loop, Sym_F, begin);continue;}
        else if(wv_N[0][2]==wv_start[2]) {
          if(nu[0]<nu_start) {Progress(mntr, loop, Sym_F, begin);continue;}
        }
      }
    }*/
    cout <<wv_N[0][0]<<" "<<wv_N[0][1]<<" "<<wv_N[0][2]<<" ("<<nu[0]<<")"<<endl;

    if(ah_f_C[0]<1.0e-8) {  //Skip bulk translation
      T = freq[0]*1.0e-12/Parameter->time;
      Write(loop, T, Delta3, Delta3_C, Delta4, Delta4_C, Gamma, Gamma_C, S_W_Out, S_W_Out_C);
      Progress(mntr, loop+nPE, begin);
      continue;
    }


    //Inner Loop
    for(wv_N[1][0]=upper[0];wv_N[1][0]>=lower[0];wv_N[1][0]--) {
    wv[1][0] = double(wv_N[1][0])/double(UC[0]);
    for(wv_N[1][1]=upper[1];wv_N[1][1]>=lower[1];wv_N[1][1]--) {
    wv[1][1] = double(wv_N[1][1])/double(UC[1]);
    for(wv_N[1][2]=upper[2];wv_N[1][2]>=lower[2];wv_N[1][2]--) {
    wv[1][2] = double(wv_N[1][2])/double(UC[2]);

      inner_F = Symmetry->F_map[upper[0]-wv_N[1][0]][upper[1]-wv_N[1][1]][upper[2]-wv_N[1][2]];
      inner_E = Symmetry->E_map[upper[0]-wv_N[1][0]][upper[1]-wv_N[1][1]][upper[2]-wv_N[1][2]];

      //Compute third wave vector
      wv_N[2][0] = (UC[0] - wv_N[0][0] - wv_N[1][0] - lower[0])%UC[0] + lower[0];
      wv_N[2][1] = (UC[1] - wv_N[0][1] - wv_N[1][1] - lower[1])%UC[1] + lower[1];
      wv_N[2][2] = (UC[2] - wv_N[0][2] - wv_N[1][2] - lower[2])%UC[2] + lower[2];
      wv[2][0] = double(wv_N[2][0])/double(UC[0]);
      wv[2][1] = double(wv_N[2][1])/double(UC[1]);
      wv[2][2] = double(wv_N[2][2])/double(UC[2]);
      line_F = Symmetry->F_map[upper[0]-wv_N[2][0]][upper[1]-wv_N[2][1]][upper[2]-wv_N[2][2]];
      line_E = Symmetry->E_map[upper[0]-wv_N[2][0]][upper[1]-wv_N[2][1]][upper[2]-wv_N[2][2]];


      //Identify symmetries (need factor4 even if frequency shift is not computed)
      factor4 = Symmetry->Sym_Negative(wv_N[1]);
      factor3 = Symmetry->Sym_3Phonon(wv_N[1], wv_N[2]);
      if(Symmetry->Sym_Negative(wv_N[2])>0.0) {sym = false;}
      else {sym = true;}

#if defined FREQ_SHIFT
      if(factor4>0.0) {Transform->TransformG4(mode[0], wv);}
#endif
      if(factor3>0.0) {Transform->TransformG3(mode[0], wv);}
      for(nu[1]=0;nu[1]<UC_DOF;nu[1]++) {

        freq[1] = QH_f[inner_F];
        ah_f[1] = AH_fQ[inner_F];
        ah_f_C[1] = AH_fC[inner_F];
        ah_lw[1] = AH_lwQ[inner_F];
        ah_lw_C[1] = AH_lwC[inner_F];
        for(i=0;i<UC_DOF;i++) {
          mode[1][i][0] = E[inner_E][i][0];
          if(factor4>0.0) {mode[1][i][1] = E[inner_E][i][1];}
          else {mode[1][i][1] = -E[inner_E][i][1];}
        }
        inner_F += 1;
        inner_E += 1;
        if(ah_f_C[1]<1.0e-8) {continue;}  //Skip bulk translation

        //occupation number
        if((DV=ah_f[1]*B_E)<230.0) {n1 = 1.0/(exp(DV)-1.0);}
        else {n1 = 0.0;}	//n1<1/(exp(230)-1)
        n1_C = 1.0/(ah_f_C[1]*B_E);		//Classical value

#if defined FREQ_SHIFT
        if(factor4>0.0) {
          //Fourth Order//
          T = factor4*Transform->TransformT4(mode[1]);
          Delta4 += T/ah_f[1]*(n1+n1+1.0);
          Delta4_C += T/ah_f_C[1]*(n1_C+n1_C);

          //Third  Order//
          if(abs(wv_N[2][0])+abs(wv_N[2][1])+abs(wv_N[2][2])==0) {
            for(nu[2]=0;nu[2]<UC_DOF;nu[2]++) {
              if((ah_f_C[2] = AH_fC[line_F+nu[2]])<1.0e-8) {continue;}
              freq[2]=QH_f[line_F+nu[2]];
              ah_f[2] = AH_fQ[line_F+nu[2]];
              ah_f_C[2] = AH_fC[line_F+nu[2]];
              for(i=0;i<UC_DOF;i++) {
                mode[2][i][0] = E[line_E+nu[2]][i][0];
                if(sym) {mode[2][i][1] = -E[line_E+nu[2]][i][1];}
                else {mode[2][i][1] = E[line_E+nu[2]][i][1];}
              }
              T = 4.0*Transform->Transform_3_0(mode, wv);
              Delta3 -= T/ah_f[1]/ah_f[2]*(n1+n1+1.0);
              Delta3_C -= T/ah_f_C[1]/ah_f_C[2]*(n1_C+n1_C);
            }
          }
        }
#endif

        if(factor3<0.0) {continue;}
        Transform->TransformH3(mode[1]);
        for(nu[2]=0;nu[2]<UC_DOF;nu[2]++) {
          freq[2] = QH_f[line_F];
          ah_f[2]   = AH_fQ[line_F];
          ah_f_C[2] = AH_fC[line_F];
          ah_lw[2]   = AH_lwQ[line_F];
          ah_lw_C[2] = AH_lwC[line_F];
          for(i=0;i<UC_DOF;i++) {
            mode[2][i][0] = E[line_E][i][0];
            if(sym) {mode[2][i][1] = -E[line_E][i][1];}
            else {mode[2][i][1] = E[line_E][i][1];}
          }
          line_F += 1;
          line_E += 1;			//Do not place after continue
          if(ah_f_C[2]<1.0e-8) {continue;}  //Skip bulk translation

          //occupation numbers
          if((DV=ah_f[2]*B_E)<230.0) {n2 = 1.0/(exp(DV)-1.0);}
          else {n2 = 0.0;}		//n2<1/(exp(230)-1)
          n2_C = 1.0/(ah_f_C[2]*B_E);		//Classical value

          //
          PV_DV(ah_f, n1+0.5, n2+0.5, ah_lw[0]+ah_lw[1]+ah_lw[2], PV, DV);
          PV_DV(ah_f_C, n1_C, n2_C, ah_lw_C[0]+ah_lw_C[1]+ah_lw_C[2], PV_C, DV_C);
          T = factor3*Transform->TransformT3(mode[2]);
          T_C = T/ah_f_C[1]/ah_f_C[2];
          T = T/ah_f[1]/ah_f[2];
          //DV =1./DV/(1.0e-12/Parameter->time);
          //Gamma += T*DV;
          Gamma += 1;
          //Gamma += DV/(1.0e-12/Parameter->time);
          //S_W_Out<<Gamma<<endl;
          Gamma_C += T_C*DV_C;
#if defined FREQ_SHIFT
          Delta3 += T*PV;
          Delta3_C += T_C*PV_C;
#endif

        }  //for(nu[2])
        line_F -= UC_DOF;
        line_E -= UC_DOF;

      }  //for(nu[1])
    }  //for(wv_N[1][2])
    }  //for(wv_N[1][1])
    }  //for(wv_N[1][0])

    //Output line shift and width for phonon mode
    T = per_ps/ah_f[0];
    //T = per_ps;
    Delta3 *= T;
    Delta4 *= 2.0*T;
    //Gamma = 1.0/(2.0*Gamma*T);
    T = per_ps/ah_f_C[0];
    Delta3_C *= T;
    Delta4_C *= 2.0*T;
    Gamma_C = 1.0/(2.0*Gamma_C*T);
    T = freq[0]*1.0e-12/Parameter->time;
    //S_W_Out<<per_ps<<endl;
    Write(loop, T, Delta3, Delta3_C, Delta4, Delta4_C, Gamma, Gamma_C, S_W_Out, S_W_Out_C);
    Progress(mntr, loop+nPE, begin);

	}  //for(loop)
#if defined PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);	//Ensure nodes do not get too far out of sync
#endif


	//CLOSE FILES, DEALLOCATE MEMORY AND RETURN
	if(PE==0) {
	  S_W_Out.close();
	  S_W_Out_C.close();
	}
	delete[] wv_N[0];  wv_N[0] =NULL;
	delete[] wv_N;     wv_N    =NULL;
	delete[] wv[0];    wv[0]   =NULL;
	delete[] wv;       wv      =NULL;
	delete[] nu;       nu      =NULL;
	delete[] freq;     freq    =NULL;
	delete[] ah_f;     ah_f    =NULL;
	delete[] ah_f_C;   ah_f_C  =NULL;
	delete[] ah_lw;    ah_lw   =NULL;
	delete[] ah_lw_C;  ah_lw_C =NULL;
	for(i=0;i<3;i++) {delete[] mode[i][0];  mode[i][0] = NULL;}
	delete[] mode[0];  mode[0] =NULL;
	delete[] mode;     mode    =NULL;
	delete[] wv_start; wv_start=NULL;
	return;
}


/*SUBROUTINE PV_DV (COMPUTES COMPLETE PRINCIPAL VALUE AND DELTA VALUES)*/
void ANHARMONIC::PV_DV(double *F, double n1, double n2, double e, double& PV, double& DV) {
	//DECLARE LOCAL VARIABLES
	double n, a, b, A, B, e2;
	e2 = e*e;


	//COMPUTE (1/x)_P = x/(x*x+e*e)) AND D(x) = e/(x*x+e*e))
	n = n1 + n2;
	a = F[0] - F[1] - F[2];
	A = 1.0/(a*a+e2);
#if defined FREQ_SHIFT
	b = F[0] + F[1] + F[2];
	B = 1.0/(b*b+e2);
	PV = n * (a*A - b*B);
#endif
	DV = n * e * (A/* - B*/);

	n = n1 - n2;
	a = F[0] + F[1] - F[2];
	A = 1.0/(a*a+e2);
	b = F[0] - F[1] + F[2];
	B = 1.0/(b*b+e2);
#if defined FREQ_SHIFT
	PV += n * (a*A - b*B);
#endif
	DV += n * e * (A - B);
	//S_W_Out<<DV<<endl;

	return;
}


/*SUBROUTINE Write (OUTPUTS DATA)*/
void ANHARMONIC::Write(int n, double F_qh, double D3Q, double D3C, double D4Q,
  double D4C, double GQ, double GC, ofstream &Quant, ofstream &Class) {
#if defined PARALLEL
	MPI_Status status;
  if(PE==0) {
    Quant <<F_qh<<'\t'<<D3Q<<'\t'<<D4Q<<'\t'<<GQ<<endl;
    Class <<F_qh<<'\t'<<D3C<<'\t'<<D4C<<'\t'<<GC<<endl;
    for(int i=1;i<nPE;i++) {
      if(n+i>=Sym_F) {break;}
      MPI_Recv(&F_qh, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&D3Q,  1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&D4Q,  1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
      MPI_Recv(&GQ,   1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
      MPI_Recv(&D3C,  1, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &status);
      MPI_Recv(&D4C,  1, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
      MPI_Recv(&GC,   1, MPI_DOUBLE, i, 7, MPI_COMM_WORLD, &status);
      Quant <<F_qh<<'\t'<<D3Q<<'\t'<<D4Q<<'\t'<<GQ<<endl;
      Class <<F_qh<<'\t'<<D3C<<'\t'<<D4C<<'\t'<<GC<<endl;
    }
  }
  else {
    MPI_Ssend(&F_qh, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Ssend(&D3Q,  1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Ssend(&D4Q,  1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    MPI_Ssend(&GQ,   1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    MPI_Ssend(&D3C,  1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
    MPI_Ssend(&D4C,  1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
    MPI_Ssend(&GC,   1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
  }
#else
  Quant <<F_qh<<'\t'<<D3Q<<'\t'<<D4Q<<'\t'<<GQ<<endl;
  Class <<F_qh<<'\t'<<D3C<<'\t'<<D4C<<'\t'<<GC<<endl;
#endif

  return;
}


/*SUBROUTINE Progress (OUTPUT PROGRESS TO SCREEN AND LOG FILE)*/
void ANHARMONIC::Progress(int &mntr, int loop, int begin) {
  if(PE!=0) {return;}
  if(Sym_F*mntr<=100*loop) {
    if(loop>Sym_F) {if(loop-nPE<Sym_F) {loop = Sym_F;}}
    Log <<"Anharmonic calculation is "<<int(double(loop)/Sym_F*100)
      <<"% complete after "<<double(time(0)-begin)/60.0<<" minutes.";
    Log <<"  "<<loop<<" of "<<Sym_F<<" waves calculated."<<endl;
    cout <<int(double(loop)/Sym_F*100)<<"% done after "
      <<double(time(0)-begin)/60.0<<" minutes.";
    cout <<"  "<<loop<<" of "<<Sym_F<<" waves calculated."<<endl;
    mntr += 1;
  }
  return;
}


/*
  ANHARMONIC Diagnostics: DETERMINES IF AND WHERE A PREVIOUS
CALCULATION STOPPED DURING THE CURRENT CONTINUED ITERATION.
*/
void ANHARMONIC::Diagnostics(void) {
  //DETERMINE IF FILES FOR ITERATION continue_ EXIST
  Log  <<"\nChecking for Incomplete Anharmonic Calculation for Iteration "<<continue_<<"..."<<endl;
  cout <<"\nChecking for Incomplete Anharmonic Calculation for Iteration "<<continue_<<"..."<<endl;
  char filename_Q[256];
  char filename_C[256];
  char header[1024];
  fstream SWQ_io, SWC_io;
  sprintf(filename_Q, "%s%i.txt", "Shift_Width", continue_);
  SWQ_io.open(filename_Q, ios::in);
  if(!SWQ_io.is_open()) {return;}
  sprintf(filename_C, "%s%i.txt", "Classical_Shift_Width", continue_);
  SWC_io.open(filename_C, ios::in);
  if(!SWC_io.is_open()) {return;}
  SWQ_io.getline(header, 1024);
  SWC_io.ignore(int(1e6), '\n');
  Log  <<"\nDetected files from previous calculation of iteration "<<continue_<<".\n";
  cout <<"\nDetected files from previous calculation of iteration "<<continue_<<".\n";


  //READ DATA AND COMPARE QH FREQUENCIES
  int nu, j, i = 0;
  int wv[3];
  double per_ps = 1.0e12*Parameter->time;	//(Non-)dimensionalizing factor (1/ps)
  double *AH_3Q = new double[Sym_F];
  double *AH_3C = new double[Sym_F];
  double *AH_4Q = new double[Sym_F];
  double *AH_4C = new double[Sym_F];
	for(wv[0]=upper[0];wv[0]>=lower[0];wv[0]--) {
  for(wv[1]=upper[1];wv[1]>=lower[1];wv[1]--) {
  for(wv[2]=upper[2];wv[2]>=lower[2];wv[2]--) {
    //Test Symmetry Operations
    if(Symmetry->Sym_Negative(wv)<0.0) {continue;}
    if(Symmetry->Sym_Freq(wv)) {continue;}

    for(nu=0;nu<UC_DOF;nu++) {
      //Read Anharmonic Data
      SWQ_io >>AH_fQ[i]>>AH_3Q[i]>>AH_4Q[i]>>AH_lwQ[i];
      SWC_io >>AH_fC[i]>>AH_3C[i]>>AH_4C[i]>>AH_lwC[i];
      if( (fabs(AH_fQ[i]*per_ps-QH_f[i])>1.0e-4)||
          (fabs(AH_fC[i]*per_ps-QH_f[i])>1.0e-4)||
          (SWQ_io.eof())||(SWC_io.eof()) ) {
        //Exit loops
        for(j=0;j<DIM;j++) {wv[j] = lower[j] - 1;}
        nu = UC_DOF;
        continue;
      }

      i += 1;
    }
  }
  }
	}
	SWQ_io.clear();
	SWQ_io.close();
  SWC_io.clear();
	SWC_io.close();


	//WRITE NEW SHIFT_WIDTH FILES
	if(i==Sym_F) {
	  Log  <<"  No anomolies detected.  Continuing with next iteration."<<endl;
    cout <<"  No anomolies detected.  Continuing with next iteration."<<endl;
	}
	else {
	  Log  <<"  Anomolies detected.  Continuing current iteration from";
    cout <<"  Anomolies detected.  Continuing current iteration from";
    continue_ -= 1;
    i -= 2;
    if(i<0) {i = 0;}
    wv_start = new int[3];
    nu_start = UC_DOF;
    for(j=0;j<DIM;j++) {wv_start[j] = lower[j];}
    ofstream SWQ_io2, SWC_io2;
    SWQ_io2.open(filename_Q, ios::trunc);
    SWC_io2.open(filename_C, ios::trunc);
    SWQ_io2 <<header<<endl;
    SWC_io2 <<header<<endl;
    j = 0;
    for(wv[0]=upper[0];wv[0]>=lower[0];wv[0]--) {
    for(wv[1]=upper[1];wv[1]>=lower[1];wv[1]--) {
    for(wv[2]=upper[2];wv[2]>=lower[2];wv[2]--) {
      //Test Symmetry Operations
      if(Symmetry->Sym_Negative(wv)<0.0) {continue;}
      if(Symmetry->Sym_Freq(wv)) {continue;}

      for(nu=0;nu<UC_DOF;nu++) {
        if(j==i) {
          for(nu_start=0;nu_start<DIM;nu_start++) {
            wv_start[nu_start] = wv[nu_start];
            wv[nu_start] = lower[nu_start] - 1;
          }
          nu_start = nu;
          nu = UC_DOF;
          continue;
        }
        Write(j, AH_fQ[j], AH_3Q[j], AH_3C[j], AH_4Q[j], AH_4C[j],
              AH_lwQ[j], AH_lwC[j], SWQ_io2, SWC_io2);
        j += 1;
      }
    }
    }
    }
    SWQ_io2.close();
    SWC_io2.close();
    for(j=0;j<DIM;j++) {
      Log  <<wv_start[j]<<" ";
      cout <<wv_start[j]<<" ";
    }
    Log  <<"("<<nu_start<<")."<<endl;
    cout <<"("<<nu_start<<")."<<endl;
	}


  //DEALLOCATE MEMORY AND RETURN
  delete[] AH_3Q;  AH_3Q=NULL;
  delete[] AH_3C;  AH_3C=NULL;
  delete[] AH_4Q;  AH_4Q=NULL;
  delete[] AH_4C;  AH_4C=NULL;

  return;
}


/*ANHARMONIC DESTRUCTOR: DEALLOCATES MEMORY*/
ANHARMONIC::~ANHARMONIC() {
  if(iterations==0) {return;}
  int i;
  //DEALLOCATE MEMORY
  delete[] QH_f;      QH_f     =NULL;
  delete[] AH_fQ;     AH_fQ    =NULL;
  delete[] AH_fC;     AH_fC    =NULL;
  delete[] AH_lwQ;    AH_lwQ   =NULL;
  delete[] AH_lwC;    AH_lwC   =NULL;
  for(i=0;i<Sym_E;i++) {delete[] E[i][0];  E[i][0]=NULL;}
  delete[] E[0];      E[0]     =NULL;
  delete[] E;         E        =NULL;
  return;
}
