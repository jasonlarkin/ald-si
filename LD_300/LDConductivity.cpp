/*						             Conductivity.cpp	              					*/
/*							              02/29/2008		              					*/
/*********************************************************************
*    This subroutine computes the conductivity of a material and the *
*  contribution as a function of anharmonic freqeucny given the		   *
*  anharmonic frequencies, lifetimes, and goroup velocities of a  	 *
*  lattice.  k_x = sum_branch&wave vector(C_v*v_g^2*tau).			       *
*********************************************************************/


/*DECLARE HEADERS*/
#include "LDCode.h"
#include <cmath>


/*DECLARE PREPROCESSOR DEFINITIONS*/


/*DECLARE GLOBAL VARIABLES*/


/*MAIN ROUTINE*/
int Conductivity(int iter, double Temperature) {
	//DELCARE LOCAL VARIABLES
	int *UC = Lattice->N;
	int i, j, jj, dim, mntr;	  //Counters
	int *wv_N=new int[DIM];			//Wave vector counters
	int *wv_IBZ=new int[2*DIM]; //Wave vector in the irreducible Brillioun zone
	double **v_=new double*[DIM];//
	int MaxWave = Symmetry->TotalDOF, MaxSym = Symmetry->Sym_F;	//Phonons at each wave vector, total number of phonons
	int natom = Unit_Cell->natom;					//Number of: atoms per unit cell, unit cells
	int *lower = Lattice->lower, *upper = Lattice->upper;	//Lower and upper bounds of 1st Brillioun zone
	int bins = 1000;			      //Number of bins to use for histogram
	double F, F_C;				      //Current frequencies
	double V, x, y, z, e, C_v, n, n_C;//Volume, temp, combined variables, and specific heat
	double factor;				      //Multiplicity factor for symmetry conditions
	double *k, *k_C;	          //Lattice thermal conductivity
	double F_max;				        //Maximum frequency or bin size
	double *DOS, *DOS_C;		    //Density of states
	double **k_DOS, **k_DOS_C;	//Histogram (Conductivity as a function of frequency)
  double h_bar = 1.05457e-34; //Plank's constant [J-s]
	double RMS_G = 0.0, RMS_D = 0.0;
	double RMS_G_C = 0.0, RMS_D_C = 0.0;
  char filename[256];


	//INITIALIZE FILES AND DATA
	sprintf(filename, "%s%i.xls", "Conductivity", iter);
	ofstream Output(filename);
	sprintf(filename, "%s%i.txt", "Shift_Width", iter);
	ifstream S_WQ(filename);//Input file for frequencies/lifetimes (10^12 rad/s)
	if(!S_WQ.is_open()) {Log <<"Could not open file '"<<filename<<"'"<<endl;exit(0);}
	sprintf(filename, "%s%i.txt", "Classical_Shift_Width", iter);
	ifstream S_WC(filename);//Classical values
	if(!S_WC.is_open()) {Log <<"Could not open file '"<<filename<<"'"<<endl;exit(0);}


	//ALLOCATE MEMORY AND READ DATA
	S_WQ.ignore(10000, '\n');		//Skip rest of header line
	S_WC.ignore(10000, '\n');		//Skip entire header line
	ifstream Velocity("Velocity.txt");		//Input file for velocities
	if(!Velocity.is_open()) {Log <<"Could not open file 'Velocity.txt'"<<endl;exit(0);}
	Velocity.ignore(10000, '\n');			//Skip header line
	double *QH_F = new double[MaxSym];
	double *freq = new double[MaxSym];
	double *freq_C = new double[MaxSym];
	double *tau = new double[MaxSym];
	double *tau_C = new double[MaxSym];
	double **vel = new double*[MaxSym];
	vel[0] = new double[MaxSym*DIM];
	F_max = 0;
	for(i=0;i<MaxSym;i++) {
	  if(i<MaxSym-1) {vel[i+1] = vel[i] + DIM;}
	  S_WQ >>QH_F[i]>>x>>e>>tau[i];
	  if((freq[i]=QH_F[i]+x+e)<1.0e-8) {freq[i] = QH_F[i];}
	  if(F_max<freq[i]) {F_max = freq[i];}
	  S_WC >>freq_C[i]>>x>>e>>tau_C[i];
	  if((freq_C[i]=QH_F[i]+x+e)<1.0e-8) {freq_C[i] = QH_F[i];}
	  if(F_max<freq_C[i]) {F_max = freq_C[i];}
	  for(j=0;j<DIM;j++) {Velocity >>vel[i][j];}
	  for(dim=0;dim<DIM;dim++) {
	    if(Lattice->BC[dim]!=SCATTERING) {continue;}
      if( (fabs(vel[i][dim])>1.0e-5)&&(QH_F[i]>1.0e-5) ) {		//Boundary scattering
        x = 0.0;
        for(j=0;j<DIM;j++) {x += Lattice->a[dim][j]*Lattice->a[dim][j];}
        x = sqrt(x)*Lattice->N[dim]*Parameter->length/2.0/fabs(vel[i][dim])*1.0e12;
        if(freq[i]>1.0e-5) {e = x/freq[i]*QH_F[i];}
        else {e = x;}
        tau[i] = (tau[i]*e)/(tau[i]+e);
        if(freq_C[i]>1.0e-5) {e = x/freq_C[i]*QH_F[i];}
        else {e = x;}
        tau_C[i] = (tau_C[i]*e)/(tau_C[i]+e);
      }
	  }
	}
	S_WQ.close();
	S_WC.close();
	Velocity.close();
	double** u = new double*[natom];      //Quantum displacement vectors
  u[0] = new double[natom*DIM];
  double** u_C = new double*[natom];    //Classical displacement vectors
  u_C[0] = new double[natom*DIM];
  for(j=0;j<natom;j++) {
    if(j<natom-1) {
      u[j+1] = u[j] + DIM;
      u_C[j+1] = u_C[j] + DIM;
	  }
    for(dim=0;dim<DIM;dim++) {
      u[j][dim] = 0.0;
      u_C[j][dim] = 0.0;
	  }
	}


	//INITIALIZE MORE VARIABLES
	V = Lattice->V;
	for(dim=0;dim<DIM;dim++) {V *= Lattice->N[dim]*Parameter->length;}
	const double B_E = 7.63850231/Temperature/Parameter->temp;//B_E = h_bar/(k_B*T) freq in rad/ps
	const double k_B_V = 1.3806e-23/V*1.0e-12;//Boltzmann's constant divided by volume*1.0e-12 (from tau)
	bins = MaxWave/(DIM*2000);   //6000 phonons per bin
	if(bins<50) {bins = 50;}	//Minimum of 50 bins
	F_max = F_max/double(bins-1); //Bin size
	k = new double[DIM];
	k_C = new double[DIM];
	double *mfp = new double[DIM];
	double *mfp_C = new double[DIM];
	DOS = new double [bins+1];
	DOS_C = new double [bins+1];
	k_DOS = new double*[DIM];
	k_DOS_C = new double*[DIM];
	k_DOS[0] = new double[(bins+1)*DIM];
	k_DOS_C[0] = new double[(bins+1)*DIM];
	for(i=0;i<=bins;i++) {DOS[i] = DOS_C[i] = 0.0;}
  for(j=0;j<DIM;j++) {
    if(j<DIM-1) {
      k_DOS[j+1] = k_DOS[j] + bins+1;
      k_DOS_C[j+1] = k_DOS_C[j] + bins+1;
    }
    mfp[j] = mfp_C[j] = k[j] = k_C[j] = 0.0;
    for(i=0;i<=bins;i++) {k_DOS[j][i] = k_DOS_C[j][i] = 0.0;}
  }


	//CALCULATE CONDUCTIVITY
	ifstream Evector("Eigenvector.txt");
	if(!Evector.is_open()) {Log <<"Could not open file 'Eigenvector.txt'"<<endl;exit(0);}
	Evector.ignore(10000, '\n');

	//Loop Over Wave Vectors
	i = 0;
	mntr = 1;
	for(wv_N[0]=upper[0];wv_N[0]>=lower[0];wv_N[0]--) {
  for(wv_N[1]=upper[1];wv_N[1]>=lower[1];wv_N[1]--) {
  for(wv_N[2]=upper[2];wv_N[2]>=lower[2];wv_N[2]--) {

    //Perform symmetry operations
    if((factor=Symmetry->Sym_Negative(wv_N))<0.0) {continue;}
    Symmetry->Sym_Freq(wv_N, &i, wv_IBZ);
    for(j=0;j<DIM;j++) {v_[j] = vel[i] + wv_IBZ[j];}

    for(j=0;j<UC_DOF;j++) {
      //Compute The Heat Capacity and Occupation Number
      F = freq[i];				//Assign frequecies to varaibles
      F_C = freq_C[i];
      if(F<1.0e-8) {n=0.0;C_v=k_B_V;}//Account for zero frequency
      else {
        x = B_E * F;				  //h_bar * freq / (k_B * T)
        e = exp(x);				    //exp[x]
        if(e>1.0e200) {e=1.0e200;}//Limit exponential to avoid overflow
        n = 1.0/(e-1.0);			//Quantum occupation number
        C_v = (k_B_V*x*x) * e*n*n;//Quantum heat capacity [x^2*k_B * e/(e-1)^2]
      }
      if(F_C<1.0e-8) {n_C=0.0;}//Account for zero frequency
      else {n_C = 1.0/(B_E * F_C);}	//Classical occupation number


      //Determine RMS data
      e = 1.0;
      if(QH_F[i]>=1.0e-8) {
        x = (F-QH_F[i])/QH_F[i];
        RMS_D += factor*x*x;
        x = 1.0/(2.0*tau[i]*QH_F[i]);
        RMS_G += factor*x*x;
        x = (F_C-QH_F[i])/QH_F[i];
        RMS_D_C += factor*x*x;
        x = 1.0/(2.0*tau_C[i]*QH_F[i]);
        RMS_G_C += factor*x*x;

        e = 1.0/(QH_F[i]*QH_F[i]);
      }


      //Accumulate The Conductivity
      mfp[0] += factor * (F*F*e) * C_v * fabs(*v_[0]);
      mfp[1] += factor * (F*F*e) * C_v * fabs(*v_[1]);
      mfp[2] += factor * (F*F*e) * C_v * fabs(*v_[2]);
      x = factor * (F*F*e) * C_v * (*v_[0]) * (*v_[0]) * tau[i];
      y = factor * (F*F*e) * C_v * (*v_[1]) * (*v_[1]) * tau[i];
      z = factor * (F*F*e) * C_v * (*v_[2]) * (*v_[2]) * tau[i];
      k[0] += x;
      k[1] += y;
      k[2] += z;
      jj = 0;
      while(jj <= bins) {
        if( F<=(F_max*(jj+1.0)) ) {
          k_DOS[0][jj] += x;
          k_DOS[1][jj] += y;
          k_DOS[2][jj] += z;
          DOS[jj] += factor;
          break;
        }
        jj += 1;
      }
      mfp_C[0] += factor * (F_C*F_C*e) * k_B_V * fabs(*v_[0]);
      mfp_C[1] += factor * (F_C*F_C*e) * k_B_V * fabs(*v_[1]);
      mfp_C[2] += factor * (F_C*F_C*e) * k_B_V * fabs(*v_[2]);
      x = factor * (F_C*F_C*e) * k_B_V * (*v_[0]) * (*v_[0]) * tau_C[i];
      y = factor * (F_C*F_C*e) * k_B_V * (*v_[1]) * (*v_[1]) * tau_C[i];
      z = factor * (F_C*F_C*e) * k_B_V * (*v_[2]) * (*v_[2]) * tau_C[i];
      k_C[0] += x;
      k_C[1] += y;
      k_C[2] += z;
      jj = 0;
      while(jj <= bins) {
        if( F_C<=(F_max*(jj+0.5)) ) {
          k_DOS_C[0][jj] += x;
          k_DOS_C[1][jj] += y;
          k_DOS_C[2][jj] += z;
          DOS_C[jj] += factor;
          break;
        }
        jj += 1;
      }


      //Accumulate The RMS Displacement
      if(F>=1.0e-8) {F   = factor*(2.0*n + 1.0)/F;}
      else {F = 1.0e100;}
      if(F_C>=1.0e-8) {F_C = factor*(2.0*n_C)/F_C;}
      else {F_C = 1.0e100;}
      for(jj=0;jj<natom;jj++) {
        for(dim=0;dim<DIM;dim++) {
          Evector >>e>>x;
          e = e*e + x*x;		//Square eigenvector
          if(F<1.0e50) {u[jj][dim] += e*F;}
          if(F_C<1.0e50) {u_C[jj][dim] += e*F_C;}
        }
      }


      //Increment counters
      i += 1;
      v_[0] += DIM;
      v_[1] += DIM;
      v_[2] += DIM;


      //Output Progress To Screen
      if(MaxSym*mntr<=10*i) {
        cout <<int(double(i)/MaxSym*10)*10<<"% done.";
        cout <<"  "<<i<<" of "<<MaxSym<<" waves calculated."<<endl;
        while(MaxSym*mntr<=10*i) {mntr += 1;}
      }

    }  //for(j)
  }  //for(wv_N[2])
  }  //for(wv_N[1])
	}  //for(wv_N[0])


	//SCALE DISPLACEMENTS [h_bar/(2*N*mass)]
	x = h_bar/(2.0*UC[0]*UC[1]*UC[2]*Parameter->mass)*1.0e-12*1.0e10*1.0e10;
	for(jj=0;jj<natom;jj++) {
	  for(dim=0;dim<DIM;dim++) {
	    u[jj][dim]   *= x;
	    u_C[jj][dim] *= x;
	  }
	}


	//PRINT RESULTS AND CLOSE FILES
	Evector.close();
	Output <<"Total Quantum Conductivity (W/m-K) = \t"<<k[0]<<'\t'<<k[1]<<'\t'<<k[2]<<endl;
	Output <<"Total Classical Conductivity (W/m-K) = \t"<<k_C[0]<<'\t'<<k_C[1]<<'\t'<<k_C[2]<<endl;
	Output <<"\nQuantum <(Delta/omega)^2>^1/2 = \t"<<sqrt(RMS_D/double(MaxWave))
		   <<"\nQuantum <(Gamma/omega)^2>^1/2 = \t"<<sqrt(RMS_G/double(MaxWave))<<endl;
	Output <<"\nClassical <(Delta/omega)^2>^1/2 = \t"<<sqrt(RMS_D_C/double(MaxWave))
		   <<"\nClassical <(Gamma/omega)^2>^1/2 = \t"<<sqrt(RMS_G_C/double(MaxWave))<<endl;
  Output <<"\nQuantum effective MFP (nm) = \t"<<k[0]/mfp[0]*1e-3<<'\t'<<k[1]/mfp[1]*1e-3<<'\t'<<k[2]/mfp[2]*1e-3<<endl;
  Output <<"Classical effective MFP (nm) = \t"<<k_C[0]/mfp_C[0]*1e-3<<'\t'<<k_C[1]/mfp_C[1]*1e-3<<'\t'<<k_C[2]/mfp_C[2]*1e-3<<endl;
	Output <<"\nQuantum RMS Displacements in Angstroms"<<endl;
	n = 0.0;
	for(jj=0;jj<natom;jj++) {
      Output <<"Atom "<<jj<<":"<<endl;
      x = 0.0;
      for(dim=0;dim<DIM;dim++) {
        Output <<"  "<<sqrt(u[jj][dim]);
        x += u[jj][dim];
	  }
      x = sqrt(x);
      Output <<"  :  "<<x<<endl;
      n += x;
	}
	Output <<"Average Quantum RMS Displacement = "<<n/double(natom)<<endl;
	Output <<"\nClassical RMS Displacements in Angstroms"<<endl;
	n = 0.0;
	for(jj=0;jj<natom;jj++) {
      Output <<"Atom "<<jj<<":"<<endl;
      x = 0.0;
      for(dim=0;dim<DIM;dim++) {
        Output <<"  "<<sqrt(u_C[jj][dim]);
        x += u_C[jj][dim];
	  }
      x = sqrt(x);
      Output <<"  :  "<<x<<endl;
      n += x;
	}
	Output <<"Average Classical RMS Displacement = "<<n/double(natom)<<endl;
	Output <<"\n\n\nDensity of States and Frequency Dependence of "
		<<"Thermal Conductivity\nQuantum\t\t\t\t\t\tClassical\n"
		<<"Anharm Freq (1/ps)\tDOS (area = 1e12)\tk-DOS (area = 1e12)\t\t"
		<<"Anharm Freq (1/ps)\tDOS (area = 1e12)\tk-DOS (area = 1e12)\n";
	DOS[0] = DOS[0] - 3; //Remove 3 translational modes (k = [0 0 0])
	DOS_C[0] = DOS_C[0] - 3;
//	factor = 1.0e-18*(1.0e-12*1.0e30/F_max/V);  //scale to a reasonable range (1e12 Phonons/m^3)
	factor = 1.0/(UC[0]*UC[1]*UC[2]*DIM*natom*F_max);  //scale to integrate to 10^12
	k[0] = 1.0/k[0]/F_max;			//scale to integrate to 10^12
	k[1] = 1.0/k[1]/F_max;			//scale to integrate to 10^12
	k[2] = 1.0/k[2]/F_max;			//scale to integrate to 10^12
	k_C[0] = 1.0/k_C[0]/F_max;	//scale to integrate to 10^12
	k_C[1] = 1.0/k_C[1]/F_max;	//scale to integrate to 10^12
	k_C[2] = 1.0/k_C[2]/F_max;	//scale to integrate to 10^12
	for(i=0;i<=bins;i++) {
		Output <<F_max*(i+0.5)<<'\t'<<DOS[i]*factor<<'\t'<<k_DOS[0][i]*k[0]<<'\t'<<k_DOS[1][i]*k[1]<<'\t'<<k_DOS[2][i]*k[2]<<"\t\t";
		Output <<F_max*(i+0.5)<<'\t'<<DOS_C[i]*factor<<'\t'<<k_DOS_C[0][i]*k_C[0]<<'\t'<<k_DOS_C[1][i]*k_C[1]<<'\t'<<k_DOS_C[2][i]*k_C[2]<<endl;
	}
	Output.close();


	//DEALLOCATE MEMORY AND RETURN
	delete[] wv_N;       wv_N      =NULL;
	delete[] wv_IBZ;     wv_IBZ    =NULL;
	delete[] v_;         v_        =NULL;
	delete[] QH_F;       QH_F      =NULL;
	delete[] freq;       freq      =NULL;
	delete[] freq_C;     freq_C    =NULL;
	delete[] tau;        tau       =NULL;
	delete[] tau_C;      tau_C     =NULL;
	delete[] vel[0];     vel[0]    =NULL;
	delete[] vel;        vel       =NULL;
	delete[] u[0];       u[0]      =NULL;
	delete[] u;          u         =NULL;
	delete[] u_C[0];     u_C[0]    =NULL;
	delete[] u_C;        u_C       =NULL;
	delete[] mfp;        mfp       =NULL;
	delete[] mfp_C;      mfp_C     =NULL;
	delete[] k;          k         =NULL;
	delete[] k_C;        k_C       =NULL;
	delete[] DOS;        DOS       =NULL;
	delete[] DOS_C;      DOS_C     =NULL;
	delete[] k_DOS[0];   k_DOS[0]  =NULL;
	delete[] k_DOS;      k_DOS     =NULL;
	delete[] k_DOS_C[0]; k_DOS_C[0]=NULL;
	delete[] k_DOS_C;    k_DOS_C   =NULL;

	return(0);
}
