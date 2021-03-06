//COMPILE: g++ anharmonic.cpp -O3 -g -Wall -fopenmp -L/usr/local/lib -llapacke -llapack -lf77blas -lcblas -latlas -lgfortran -lm -lc
//         -lgsl -lcblas -latlas


//=================STANDARD HEADER===================
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <cmath>
#include <math.h>
#include <stdio.h>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <stdlib.h>
#include <complex.h>
//#include <complex>

//=========================GSL HEADER=================
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_eigen.h>
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_math.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_cblas.h>

//====================LAPACK HEADER==================
#include <lapack/lapacke.h>
extern "C" {
  #include <atlas/cblas.h>
}

using namespace std;

//=========================================================================================================================================================
//==========================================================================================================================================================
//======================INPUT/OUTPUT STREAMS=====================
ifstream data_in("input.dat");
ifstream phi2_in("PHI2.dat");
ifstream phi3_in("PHI3.dat");
ifstream phi4_in("PHI4.dat");
ofstream data_out("data.dat");
ofstream irorder2_out("IRBorder2.dat");
ofstream irorder3_out("IRBorder3.dat");
ofstream irorder4_out("IRBorder4.dat");
ofstream fullpos_out("FULLpos.dat");
ofstream fulldisp4_out("FULLorder4.dat");
ofstream fulldisp3_out("FULLorder3.dat");
ofstream fulldisp2_out("FULLorder2.dat");
//==============================================================


//======================STRUCTURE TYPES=========================
  struct LATTICE{                 //carries all information about the lattice
    double *r, *m, *R, ax, ay, az,*pos,*pR, *klist, *IBZ_klist;             // m is the mass of atoms in central unit cell, ri s position in central unit cell, pos is position of all atoms, pR is the position of primitive unit cell of that atom
      int natom, ncell,Nx,Ny,Nz,nNx,nNy,nNz,type,isConv,nlist,kx,ky,kz,kx_scale,ky_scale,kz_scale,nNeighbour2,nNeighbour3,nNeighbour4;
     int *index, ninterac,IBZ_nlist, *IBZ_num_interac;         // 1st coloumn is the atom index and 2nd is the atom type
      int *ikinterac, *kmap, *IBZ_kmap;
  };

  struct DYNA{                 //stores force constant and other information related with force constant like frequency, etc.
    double *order2, *order3, *order4;                     // force constant correspoing to displament in iorder lists
    double complex *eigvec;                                // eigen vectors of each mode
    double *freq;                                       // frequency (sqrt(eigval)) of each mode
     int *iorder2, *iorder3, *iorder4,norder2,norder3,norder4;        //iorder records indices and norder the number of corresponding indices in displcaement list
    double *linewidth, *epsi;                                 // linewidth and epsilon for each mode
    double *velocity;                                       // group velocities
     int *info2, *info3, *info4;                  // information (atom type) in order2/oreder3/order4 list
    double *pos2, *pos3, *pos4, *m2, *m3, *m4;               //positions and mass of atoms in order2/order3/order4 list
  };
/*
  struct disp{                   // useful for creating and storing irreducible list of force constant, etc
    double *order2, *order3, *order4;
     long int *iorder2, *iorder3, *iorder4;
     int norder2,norder3,norder4;
  };
*/

struct K_INFO{
    double complex  *exp_kr2, *exp_kr3;
    double *nbe;   //occupation of mode 1,2,3,4,5,6....bose einstein
     int *vecmap31, *vecmap32, *vecmap33;
  };

//===============================================================

//=====================FUNCTION PROTOTYPES=======================
void sub_readInput(LATTICE *lattice);
void sub_klist(LATTICE *lattice);
void sub_kinteraction(LATTICE *lattice);
int sub_search(long int *sklist,int s, int e,long int tmp);
void sub_ksymmetry(LATTICE *lattice);
void sub_posIndex(LATTICE *lattice);
void sub_fulldisp(LATTICE *lattice, DYNA *dynamic);
void sub_posPrimitive(LATTICE *lattice);
void sub_readPhi(DYNA *dynamic);
void sub_info(LATTICE *lattice, DYNA *dynamic);
void sub_kscale(LATTICE *lattice);
void sub_eigen(LATTICE *lattice, DYNA *dynamic);
void sub_prelinewidth(DYNA *dynamic, LATTICE *lattice, K_INFO *kprop);
double sub_delta(double w1,double w2,double w3,double e,double n2,double n3);
void sub_epsi(DYNA *dynamic,int i, LATTICE *lattice);
void sub_linewidth(LATTICE *lattice, DYNA *dynamic, int iteration,K_INFO *kprop);

//===============================================================

//=============global variables=================================
const double pi = 3.14159265;
int numthread;
double temperature;
double hbar = (6.626068*(1e-34))/(2*pi);
double Kb = 1.3806503*(1e-23);
int FC_available;             // check whether to calculate displacement list or to do full conductivity calculation
const double unit_length = 5.43*(1e-10);
const double mass = 28.0855*(1.67e-27);
const double eV = 1.6*(1e-19);
const double non_time = 1e-12;
//=================================================================

//===========================================================================================================================================================
//===========================================================================================================================================================
int main(){
  clock_t start, end,start1,end1;
  double cpuTime;
  start = clock();
  start1=start;

  //non dimensionalising variables
  Kb = Kb/((mass*unit_length*unit_length)/(non_time*non_time));
  hbar = hbar/((mass*unit_length*unit_length)/(non_time));

  //============================================================================================

  //================ variable declaration=======================================================
  LATTICE *lattice;
  DYNA *dynamic;
  K_INFO *kprop;
  //disp *irrdisp;
  //=============================================================================================

  //=========READING INPUTS AND DEVELOPING STRUCTURE=============================================
  lattice = new LATTICE;
  sub_readInput(lattice);

  end = clock();
  cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
  start = clock();
  cout<<"Structure created in "<<cpuTime<<" seconds."<<endl;
  //================================================================================================

  //======================= DEVELOPING RECIPROCAL SPACE VECTORS=====================================
  if(FC_available==1){
    sub_klist(lattice);

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"Reciprocal vectors created in "<<cpuTime<<" seconds."<<endl;
  }
  //=================================================================================================

  //======================= CREATING IRREDUCIBLE BRILLIOUN ZONE==================================
  if(FC_available==1){
    sub_ksymmetry(lattice);

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"Irreducible Brillioun zone created in "<<cpuTime<<" seconds."<<endl;
  }
  //=================================================================================================

  //============== CHECKING ALLOWED K COMBINATIONS K1 + K2 = K3 +G  ====================================
  if(FC_available==1){
    sub_kinteraction(lattice);
    //cout<<lattice->ninterac<<endl;
    //cout<<lattice->IBZ_nlist*lattice->nlist<<endl;

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"Allowed K points checked in "<<cpuTime<<" seconds."<<endl;
  }
  //===================================================================================================

  //============ Generating and writing full force constant list and index list======================
  dynamic = new DYNA;
  sub_posIndex(lattice);              // writing position file in increasing distance from origin
  sub_fulldisp(lattice,dynamic);                 // full 2nd, 3rd and 4th order displacement list
  sub_posPrimitive(lattice);           // primitive cells for dynamic matrix


  end = clock();
  cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
  start = clock();
  cout<<"Full displacement list created in "<<cpuTime<<" seconds."<<endl;
  //====================================================================================================

  /*
  //======GENERATING IRREDUCIBLE NUMBER OF DISPLACEMENT LIST FOR FORCE CONSTANTS=====================


  end = clock();
  cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
  start = clock();
  cout<<"Identifying symmetries in the structure "<<cpuTime<<" seconds."<<endl;
  //=================================================================================================
  */

  //=========READING FORCE CONSTANTS FROM FILE=======================================================
  if(FC_available==1){
    sub_readPhi(dynamic);

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"FCs read in "<<cpuTime<<" seconds."<<endl;
  }
  //=================================================================================================

  //===============creating info list (to get rid of array of array)=================================
  if(FC_available==1){
    sub_info(lattice,dynamic);

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"Info list created in "<<cpuTime<<" seconds."<<endl;
  }
  //===================================================================================================

  //=========Scaling k values to real values=======================================================
  if(FC_available==1){
    sub_kscale(lattice);

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"K values scaled in "<<cpuTime<<" seconds."<<endl;
  }
  //=================================================================================================

  //========calculating eigen vectors and frequencies for all k points==============================
  if(FC_available==1){
    sub_eigen(lattice,dynamic);       //calculates eigen vectors and eigen values for all k points

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"Eigen vectors and values calculated in "<<cpuTime<<" seconds."<<endl;
  }
  //================================================================================================

  //===========preprocessing for linewidths====================================================
  if(FC_available==1){
    kprop = new K_INFO;
    sub_prelinewidth(dynamic,lattice,kprop);

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"Preprocessing for linewidths done in "<<cpuTime<<" seconds."<<endl;
  }
  //======================================================================================================


  //================solving for linewidths===================
  if(FC_available==1){
    sub_epsi(dynamic,0,lattice);
    for(int i=0;i<5;i++){
      cout<<i<<endl;
      sub_linewidth(lattice, dynamic,i,kprop);
      sub_epsi(dynamic,i+1,lattice);


      /*for( int k=0;k<lattice->IBZ_nlist;k++){
	for( int l=0;l<6;l++){
	  data_out<<dynamic->freq[6*lattice->IBZ_kmap[k]+l]<<'\t'<<(dynamic->linewidth[6*k+l])<<endl;
	}
      }*/


    }

    for( int k=0;k<lattice->IBZ_nlist;k++){
      for( int l=0;l<6;l++){
	data_out<<dynamic->freq[6*lattice->IBZ_kmap[k]+l]<<'\t'<<1./(2.*dynamic->linewidth[6*k+l])<<endl;
      }
    }

    end = clock();
    cpuTime= (end-start)/ (double) (CLOCKS_PER_SEC);
    start = clock();
    cout<<"linewidths for different interactions (egn 16) calculated in "<<cpuTime<<" seconds."<<endl;
  }
  //================================================================================================


  end = clock();
  end1=end;
  cpuTime= (end1-start1)/ (double) (CLOCKS_PER_SEC);
  cout<<"TOTAL TIME IS "<<cpuTime<<" SECONDS."<<endl;
}



//========== FUNCTIONS===========================================================================================================================

//----------calculating epsilon for delta function----------------------------------------------------------------------------
void sub_epsi(DYNA *dynamic,int iteration, LATTICE *lattice){
  double epsi0,prefac;
  epsi0 = 0.005;
  prefac = hbar/(16.*lattice->nlist);

  if(iteration==0) {
    dynamic->epsi = new double [lattice->IBZ_nlist*6];
    for( int k=0;k<lattice->IBZ_nlist;k++){
      for(int i=0;i<6;i++){
	dynamic->epsi[6*k+i] = epsi0*dynamic->freq[6*k+i];
      }
    }
  }else{
    for( int k=0;k<lattice->IBZ_nlist;k++){
      for(int i=0;i<6;i++){
	dynamic->epsi[6*k+i] = prefac*dynamic->linewidth[6*k+i];
	dynamic->linewidth[6*k+i] *= prefac;
      }
    }
  }
}
//---------------------------------------------------------------------------------------------------------------------------

//-----------delta function-------------------------------------------------------------------------------------------------------
double sub_delta(double w1,double w2,double w3,double e,double n2,double n3){
  double result,a,b,c,e2;
  a = w1-w2-w3;
  b = w1+w2-w3;
  c = w1-w2+w3;
  e2 = e*e;
  result = (e*(((n2+n3)*(1./(a*a+e2))) + ((n2-n3)*((1./(b*b+e2))-(1./(c*c+e2))))));
  return result;
}
//---------------------------------------------------------------------------------------------------------------------------------

//--------------calculating linewidths--------------------------------------------------------------------------------------------
void sub_linewidth(LATTICE *lattice, DYNA *dynamic, int iteration, K_INFO *kprop){


  if(iteration==0) dynamic->linewidth = new double [6*lattice->IBZ_nlist];
  for( int i=0;i<lattice->IBZ_nlist;i++){
    for(int j=0;j<6;j++){
      dynamic->linewidth[6*i+j] = 0.;
    }
  }
  int i,nu1,nu2,nu3,i1,i2,i3,j,i11,factor,k,count;
  double tmp,e,w1,w2,w3,n2,n3,e1,e2,e3,tmp1;
  double complex tmpc1,*f3, *f31, *f32,tmpc;
  f3 = new double complex [6*dynamic->norder3];
  f31 = new double complex [6*dynamic->norder3];
  f32 = new double complex [6*6*dynamic->norder3];


  count =0;

  for(i = 0;i<lattice->IBZ_nlist;i++){
    i11 = lattice->ikinterac[3*count+0];
    i1 = lattice->IBZ_kmap[i11];
    for(nu1=0;nu1<6;nu1++){
      for(k=0;k<dynamic->norder3;k++){
	f3[nu1*dynamic->norder3 + k] = dynamic->order3[k] * dynamic->eigvec[36*i1 + 6*(kprop->vecmap31[k]) + nu1];
      }
    }
    for(j =0;j<lattice->IBZ_num_interac[i];j++){

      i2 = lattice->ikinterac[3*count+1];
      i3 = lattice->ikinterac[3*count+2];
      if(i2==i3) factor = 1;
      else if(i2>i3)factor = 2;
      else factor=0;

      if(factor>0){
	// force const * exp(ik1.r1) * exp(ik2.r2)===================================================
	nu2=0;
	nu3=0;
	tmpc = kprop->exp_kr2[dynamic->norder3*i2 + nu3]  * kprop->exp_kr3[dynamic->norder3*i3 + nu3];
	do{
	  for(nu1=0;nu1<27;nu1++){
	    for(k=0;k<6;k++){
	      f31[k*dynamic->norder3+nu2] = f3[k*dynamic->norder3+nu2]*tmpc;
	    }
	    nu2++;
	  }
	  nu3 += 27;
	  tmpc = kprop->exp_kr2[dynamic->norder3*i2 + nu3]  * kprop->exp_kr3[dynamic->norder3*i3 + nu3];
	}while(nu2<dynamic->norder3);
	//===============================================================================================

	for(nu1=0;nu1<6;nu1++){
	  for(nu2=0;nu2<6;nu2++){
	    for(k=0;k<dynamic->norder3;k++){
	      f32[nu1*(dynamic->norder3*6) + nu2*dynamic->norder3 + k] = f31[nu1*dynamic->norder3 + k] * dynamic->eigvec[36*i2 + 6*(kprop->vecmap32[k]) + nu2];
	    }
	  }
	}

	for(nu1=5;nu1>=0;--nu1){
	  w1 = dynamic->freq[6*i1+nu1];
	  e1 = dynamic->epsi[6*i11+nu1];
	  for(nu2=5;nu2>=0;--nu2){
	    w2 = dynamic->freq[6*i2+nu2];
	    n2 = kprop->nbe[6*i2+nu2];
	    e2 = dynamic->epsi[6*lattice->kmap[i2]+nu2];
	    e3 = e1+e2;
	    for(nu3 =5 ;nu3 >= 0;nu3--){
	      w3 = dynamic->freq[6*i3+nu3];
	      n3 = kprop->nbe[6*i3+nu3];
	      e = dynamic->epsi[6*lattice->kmap[i3]+nu3] + e3;
	      if((w1==0)||(w2==0)||(w3==0)) tmp=0;
	      else tmp = sub_delta(w1,w2,w3,e,n2 + 0.5,n3 + 0.5);
	      if(tmp!=0){
		tmpc1=0;
		for(k=0;k<dynamic->norder3;k++){
		  tmpc1 += f32[nu1*(dynamic->norder3*6) + nu2*dynamic->norder3 + k]*dynamic->eigvec[36*i3 + 6*(kprop->vecmap33[k]) + nu3];
		}
		tmp1 = cabs(tmpc1);
		tmp *= tmp1*tmp1;
	      }
	      dynamic->linewidth[6*i11 + nu1] +=  (factor * tmp);
	    }
	  }
	}

      }

      count++;
    }
  }

}
//----------------------------------------------------------------------------------------------------------------------------------


//--------------storing all information of k points for efficient calculation of linewidth-----------------------------------------
void sub_prelinewidth(DYNA *dynamic, LATTICE *lattice, K_INFO *kprop){
  double complex tt;
  tt =1I;
  kprop->exp_kr2 = new complex double [(lattice->nlist)*dynamic->norder3];
  kprop->exp_kr3 = new complex double [(lattice->nlist)*dynamic->norder3];
  kprop->nbe= new double [6*lattice->nlist];

  for(int i=0;i<lattice->nlist;i++){
    for(int j=0;j<6;j++){
      kprop->nbe[6*i+j] = 1./((exp((hbar*dynamic->freq[6*i+j])/(Kb*temperature)))-1);
    }
  }

  for( int i=0;i<lattice->nlist;i++){
    for( int j=0;j<dynamic->norder3;j++){
      kprop->exp_kr2[(dynamic->norder3)*i+j] = cexp(tt*((( (lattice->klist[3*i+0]*( dynamic->pos3[9*j + 3*(1) + 0] -  (lattice->pR[3*(dynamic->iorder3[6*j+1]) + 0 ])))   +   (lattice->klist[3*i+1]*( dynamic->pos3[9*j + 3*(1) + 1] - (lattice->pR[3*(dynamic->iorder3[6*j+1]) + 1 ])))  +     (lattice->klist[3*i+2]*( dynamic->pos3[9*j + 3*(1) + 2] - (lattice->pR[3*(dynamic->iorder3[6*j+1]) + 2 ])))))));
      kprop->exp_kr3[(dynamic->norder3)*i+j] = cexp(tt*((( (lattice->klist[3*i+0]*( dynamic->pos3[9*j + 3*(2) + 0] - (lattice->pR[3*(dynamic->iorder3[6*j+2]) + 0 ])))   +   (lattice->klist[3*i+1]*( dynamic->pos3[9*j + 3*(2) + 1] - (lattice->pR[3*(dynamic->iorder3[6*j+2]) + 1 ])))  +     (lattice->klist[3*i+2]*( dynamic->pos3[9*j + 3*(2) + 2] - (lattice->pR[3*(dynamic->iorder3[6*j+2]) + 2 ])))))));
    }
  }

  double r[2][3];
  r[0][0] = 0; r[0][1] = 0; r[0][2] = 0;
  r[1][0] = 0.25; r[1][1] = 0.25; r[1][2] = 0.25;
  for(int i=0;i<lattice->nlist;i++){
    for(int j=0;j<2;j++){
      for(int d=0;d<3;d++){
	for(int nu =0;nu<6;nu++){
	  dynamic->eigvec[36*i + 6*(3*j+d) + nu] *= (((cexp(tt*( (lattice->klist[3*i+0]*r[j][0]) + (lattice->klist[3*i+1]*r[j][1]) + (lattice->klist[3*i+2]*r[j][2]) ))))/(sqrt(lattice->m[j]*dynamic->freq[6*i+nu])));
	}
      }
    }
  }


  kprop->vecmap31 = new  int [dynamic->norder3];
  kprop->vecmap32 = new  int [dynamic->norder3];
  kprop->vecmap33 = new  int [dynamic->norder3];
  for( int i=0;i<dynamic->norder3;i++){
    kprop->vecmap31[i] = 3*(lattice->index[2*(dynamic->iorder3[6*i+0])+1]) + dynamic->iorder3[6*i+3] - 1;
    kprop->vecmap32[i] = 3*(lattice->index[2*(dynamic->iorder3[6*i+1])+1]) + dynamic->iorder3[6*i+4] - 1;
    kprop->vecmap33[i] = 3*(lattice->index[2*(dynamic->iorder3[6*i+2])+1]) + dynamic->iorder3[6*i+5] - 1;
  }
}
//----------------------------------------------------------------------------------------------------------------------------------

//-------calculating eigen values and eigen vectors---------------------------------------------------------------------------
void sub_eigen(LATTICE *lattice, DYNA *dynamic){
  //specific to diamond structure only with primitive unit cell
  double kx,ky,kz,px,py,pz;
  int x,y,info,k;
   int i;
  double m1,m2;
  double complex tt;
  tt = 1I;

  dynamic->freq = new double [6*lattice->nlist];         //6 branches for diamond structure with primitive cell
  dynamic->eigvec = new double complex [36*lattice->nlist];
  double complex *tmpdynax, *tmpdynay, *tmpdynaz;
  tmpdynax = new double complex [36*lattice->nlist];
  tmpdynay = new double complex [36*lattice->nlist];
  tmpdynaz = new double complex [36*lattice->nlist];


  for(k=0;k<lattice->nlist;k++){
    for(i=0;i<36;i++){
      dynamic->eigvec[36*k+i]=0;
      tmpdynax[36*k+i] = 0;
      tmpdynay[36*k+i] = 0;
      tmpdynaz[36*k+i] = 0;
    }
  }

  // creating dynamic matrices for all k points
    for(i=0;i<dynamic->norder2;i++){
      x = (3*(lattice->index[2*dynamic->iorder2[4*i+0] + 1])+dynamic->iorder2[4*i+2]-1);
      y = (3*(lattice->index[2*dynamic->iorder2[4*i+1] + 1])+dynamic->iorder2[4*i+3]-1);
      px = (dynamic->pos2[6*i+3*1+0] -dynamic->pos2[6*i+3*0+0]);
      py = (dynamic->pos2[6*i+3*1+1]- dynamic->pos2[6*i+3*0+1]);
      pz = (dynamic->pos2[6*i+3*1+2]- dynamic->pos2[6*i+3*0+2]);
      m1 = dynamic->m2[2*i+0];
      m2 = dynamic->m2[2*i+1];
      m1 = 1./(sqrt(m1*m2));

#pragma omp parallel for private(k,kx,ky,kz) num_threads(numthread)
      for(k=0;k<lattice->nlist;k++){
	kx = lattice->klist[3*k+0]; ky = lattice->klist[3*k+1]; kz = lattice->klist[3*k+2];
	dynamic->eigvec[36*k + 6*x + y] += ((m1 * (dynamic->order2[i])))* (cexpf((tt)*(((kx*px)+(ky*py)+(kz*pz)))));
      }
    }


    // solving eigen value problem
#pragma omp parallel for private(k,info) num_threads(numthread)
    for(k=0;k<lattice->nlist;k++){
      info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U',6,&dynamic->eigvec[36*k],6,&dynamic->freq[6*k]);
      if(info!=0) cout<<"ERROR: failed to solve for eigen values and eigenvectors."<<endl;
    }

    // frequencies from eigen values
    for(k=0;k<lattice->nlist;k++){
      //data_out<<lattice->klist[3*k+0]<<'\t'<<lattice->klist[3*k+1]<<'\t'<<lattice->klist[3*k+2]<<'\n';
      for(i=0;i<6;i++){
	dynamic->freq[6*k+i] = (dynamic->freq[6*k+i]>0)?(sqrt(dynamic->freq[6*k+i])):(0);
	//data_out<<dynamic->freq[6*k+i]<<'\t';
      }
     //data_out<<" "<<endl;
    }


    // group velocities
    dynamic->velocity = new double [lattice->nlist*6*3];
    for(i=0;i<dynamic->norder2;i++){
      x = (3*(lattice->index[2*dynamic->iorder2[4*i+0] + 1])+dynamic->iorder2[4*i+2]-1);
      y=(3*(lattice->index[2*dynamic->iorder2[4*i+1] + 1])+dynamic->iorder2[4*i+3]-1);
      px = (dynamic->pos2[6*i+3*1+0] -dynamic->pos2[6*i+3*0+0]);
      py = (dynamic->pos2[6*i+3*1+1] -dynamic->pos2[6*i+3*0+1]);
      pz = (dynamic->pos2[6*i+3*1+2] -dynamic->pos2[6*i+3*0+2]);
      m1 = dynamic->m2[2*i+0];
      m2 = dynamic->m2[2*i+1];
      m1 = 1./(sqrt(m1*m2));
#pragma omp parallel for private(k,kx,ky,kz) num_threads(numthread)
      for(k=0;k<lattice->nlist;k++){
	kx = lattice->klist[3*k+0]; ky = lattice->klist[3*k+1]; kz = lattice->klist[3*k+2];
	tmpdynax[36*k + 6*x + y] += (tt*px)*((m1 * (dynamic->order2[i])))* (cexpf((tt)*(((kx*px)+(ky*py)+(kz*pz)))));
	tmpdynay[36*k + 6*x + y] += (tt*py)*((m1 * (dynamic->order2[i])))* (cexpf((tt)*(((kx*px)+(ky*py)+(kz*pz)))));
	tmpdynaz[36*k + 6*x + y] += (tt*pz)*((m1 * (dynamic->order2[i])))* (cexpf((tt)*(((kx*px)+(ky*py)+(kz*pz)))));
      }
    }

    double complex tmpy[6],result;
    double alpha = 1.0;
    double beta =0.0;
    for (k=0;k<lattice->nlist;k++){
      for(i=0;i<6;i++){
	cblas_zgemv(CblasRowMajor,CblasNoTrans,6,6,&alpha,&tmpdynax[36*k],6,&dynamic->eigvec[36*k+i],6,&beta,tmpy,1);
	cblas_zdotc_sub(6,&dynamic->eigvec[36*k+i],6,tmpy,1,&result);
	dynamic->velocity[18*k+3*i+0] = (1./(2.0*dynamic->freq[6*k+i]))*(creal(result));


	cblas_zgemv(CblasRowMajor,CblasNoTrans,6,6,&alpha,&tmpdynay[36*k],6,&dynamic->eigvec[36*k+i],6,&beta,tmpy,1);
	cblas_zdotc_sub(6,&dynamic->eigvec[36*k+i],6,tmpy,1,&result);
	dynamic->velocity[18*k+3*i+1] = (1./(2.0*dynamic->freq[6*k+i]))*(creal(result));


	cblas_zgemv(CblasRowMajor,CblasNoTrans,6,6,&alpha,&tmpdynaz[36*k],6,&dynamic->eigvec[36*k+i],6,&beta,tmpy,1);
	cblas_zdotc_sub(6,&dynamic->eigvec[36*k+i],6,tmpy,1,&result);
	dynamic->velocity[18*k+3*i+2] = (1./(2.0*dynamic->freq[6*k+i]))*(creal(result));

      }
    }

    /*for(k=0;k<lattice->nlist;k++){
      data_out<<lattice->klist[3*k+0]<<'\t'<<lattice->klist[3*k+1]<<'\t'<<lattice->klist[3*k+2]<<'\n';
      for(i=0;i<6;i++){
	for(int j=0;j<6;j++){
	  data_out<<cimag (dynamic->eigvec[36*k + 6*j + i])<<'\t';
	}data_out<<endl;
      }
     data_out<<" "<<endl;
    }*/

    delete [] tmpdynax;
    delete [] tmpdynay;
    delete [] tmpdynaz;
}
//------------------------------------------------------------------------------------------------------------------------------------------

//----------------scaling k values by 2*pi/ax--------------------------------------------------------------------------------------------------
void sub_kscale(LATTICE *lattice){
  for(int i=0;i<lattice->nlist;i++){
    lattice->klist[3*i+0] *= 2.*pi/(lattice->ax);
    lattice->klist[3*i+1] *= 2.*pi/(lattice->ay);
    lattice->klist[3*i+2] *= 2.*pi/(lattice->az);
  }
  for( int i=0;i<lattice->IBZ_nlist;i++){
    lattice->IBZ_klist[3*i+0] *= 2.*pi/(lattice->ax);
    lattice->IBZ_klist[3*i+1] *= 2.*pi/(lattice->ay);
    lattice->IBZ_klist[3*i+2] *= 2.*pi/(lattice->az);
  }
}
//---------------------------------------------------------------------------------------------------------------------------------------------

//------------full info of displacement list------------------------------------------------------------------------------------------
void sub_info(LATTICE *lattice, DYNA *dynamic){
  dynamic->pos2 = new double [6*dynamic->norder2];
  dynamic->pos3 = new double [9*dynamic->norder3];
  dynamic->pos4 = new double [12*dynamic->norder4];
  dynamic->m2 = new double [2*dynamic->norder2];
  dynamic->m3 = new double [3*dynamic->norder3];
  dynamic->m4 = new double [4*dynamic->norder4];

  for( int i=0;i<dynamic->norder2;i++){             // i is the index of disp list...j is the coloumn in the list...and then 3 directions for it
    for(int j=0;j<2;j++){
      dynamic->pos2[6*i + 3*j + 0] = (lattice->pos[3*(dynamic->iorder2[4*i+j])+0]);
      dynamic->pos2[6*i + 3*j + 1] = (lattice->pos[3*(dynamic->iorder2[4*i+j])+1]);
      dynamic->pos2[6*i + 3*j + 2] = (lattice->pos[3*(dynamic->iorder2[4*i+j])+2]);
      dynamic->m2[2*i+j] = lattice->m[lattice->index[2*(dynamic->iorder2[4*i+j])+1]];
    }
  }

  for( int i=0;i<dynamic->norder3;i++){
    for(int j=0;j<3;j++){
      dynamic->pos3[9*i + 3*j + 0] = (lattice->pos[3*((dynamic->iorder3[6*i+j]))+0]);
      dynamic->pos3[9*i + 3*j + 1] = (lattice->pos[3*((dynamic->iorder3[6*i+j]))+1]);
      dynamic->pos3[9*i + 3*j + 2] = (lattice->pos[3*((dynamic->iorder3[6*i+j]))+2]);
      dynamic->m3[3*i+j] = lattice->m[lattice->index[2*(dynamic->iorder3[6*i+j])+1]];
    }
  }

  for( int i=0;i<dynamic->norder4;i++){
    for(int j=0;j<4;j++){
      dynamic->pos4[12*i + 3*j + 0] = (lattice->pos[3*((dynamic->iorder4[8*i+j]))+0]);
      dynamic->pos4[12*i + 3*j + 1] = (lattice->pos[3*((dynamic->iorder4[8*i+j]))+1]);
      dynamic->pos4[12*i + 3*j + 2] = (lattice->pos[3*((dynamic->iorder4[8*i+j]))+2]);
      dynamic->m4[4*i+j] = lattice->m[lattice->index[2*(dynamic->iorder4[8*i+j])+1]];
    }
  }



  dynamic->info2 = new  int [2*dynamic->norder2]; //1st two coloumns are atom type (0 or 1 in diamond structure)
  dynamic->info3 = new  int [3*dynamic->norder3];
  dynamic->info4 = new  int [4*dynamic->norder4];

  for ( int i=0;i<dynamic->norder2;i++){
    for(int j=0;j<2;j++){
      dynamic->info2[2*i+j] = lattice->index[2*(dynamic->iorder2[4*i+j])+1];
    }
  }

  for ( int i=0;i<dynamic->norder3;i++){
    for(int j=0;j<3;j++){
      dynamic->info3[3*i+j] = lattice->index[2*(dynamic->iorder3[6*i+j])+1];
    }
  }
  for ( int i=0;i<dynamic->norder4;i++){
    for(int j=0;j<4;j++){
      dynamic->info4[4*i+j] = lattice->index[2*(dynamic->iorder4[8*i+j])+1];
    }
  }


}
//---------------------------------------------------------------------------------------------------------------------------------------

//------------Reading Force Constants from File-------------------------------------------------------------------------------------
void sub_readPhi(DYNA *dynamic){
  dynamic->order2 = new double[dynamic->norder2];
  dynamic->order3 = new double[dynamic->norder3];
  dynamic->order4 = new double[dynamic->norder4];
   int i,j,k;
  i=0; j=0;k=0;
  while((!phi2_in.eof())&&(i<dynamic->norder2)){
    phi2_in>>dynamic->order2[i];
    i++;
  }
  while((!phi3_in.eof())&&(i<dynamic->norder3)){
    phi3_in>>dynamic->order3[j];
    j++;
  }
  while((!phi4_in.eof())&&(i<dynamic->norder4)){
    phi4_in>>dynamic->order4[k];
    k++;
  }
  if((i<dynamic->norder2)||(j<dynamic->norder3)||(k<dynamic->norder4)) cout<<"ERROR: Insufficient force constants in a file"<<endl;

/*cout<<i<<'\t'<<dynamic->norder2<<endl;
  cout<<j<<'\t'<<dynamic->norder3<<endl;
  cout<<k<<'\t'<<dynamic->norder4<<endl;
  */


  for(i=0;i<dynamic->norder2;i++){
    dynamic->order2[i] *= (eV*(1e+20))/(mass/(non_time*non_time));
  }
  for(i=0;i<dynamic->norder3;i++){
    dynamic->order3[i] *= (((eV*(1e+30)))/(mass/(unit_length*non_time*non_time)));
  }
  for(i=0;i<dynamic->norder4;i++){
    dynamic->order4[i] *= ((eV*(1e+40)))/(mass/(unit_length*unit_length*non_time*non_time));
  }

}
//-------------------------------------------------------------------------------------------------------------------------------------

//------------------------GENERATING PRIMITIVE STRUCTURE FOR DYNAMIC MATRIX FROM CONVENTIONAL STRUCTURE----------------------------------------------
void sub_posPrimitive(LATTICE *lattice){
  int count;
  double epsi = 1e-6;
  double tmp;
  count =0;
  for(int i=0;i<lattice->natom;i++){                     // CHECKING NUMBER OF ATOMS
    for (int j=0;j<lattice->ncell;j++){
      count++;
    }
  }

  lattice->pR = new double [count*3];                //primitive cell locations for all atoms in supercells
  for(int i=0;i<count;i++){
    tmp = lattice->pos[3*i+0] - (0.5*lattice->ax)*lround((lattice->pos[3*i+0])/(0.5*lattice->ax));
    tmp = (tmp>0)?(tmp):(-tmp);
    if(tmp<epsi){
      lattice->index[2*i+1] = 0;
      lattice->pR[3*i+0] = lattice->pos[3*i+0]-0;
      lattice->pR[3*i+1] = lattice->pos[3*i+1]-0;
      lattice->pR[3*i+2] = lattice->pos[3*i+2]-0;
    }else{
      lattice->index[2*i+1] = 1;
      lattice->pR[3*i+0] = lattice->pos[3*i+0]-0.25*lattice->ax;
      lattice->pR[3*i+1] = lattice->pos[3*i+1]-0.25*lattice->ay;
      lattice->pR[3*i+2] = lattice->pos[3*i+2]-0.25*lattice->az;
    }
  }
}
//-----------------------------------------------------------------------------------------------------------------------------------------

//--------------Creating full displacement list for 2nd,3rd and 4th order-----------------------------------------------------------------------
void sub_fulldisp(LATTICE *lattice, DYNA *dynamic){
   int i0[8], index1[5000],neighbour,count1,count2,index2[5000];
  double cutoff;

  //----checking indices for all atoms in central cell-----------------------------------------------------------------
  for(int i=0;i<((lattice->natom)*(lattice->ncell));i++){
    if((lattice->pos[3*i+0]==lattice->r[3*0+0]) && (lattice->pos[3*i+1]==lattice->r[3*0+1])  && (lattice->pos[3*i+2]==lattice->r[3*0+2])) i0[0] = lattice->index[2*i+0];
    if((lattice->pos[3*i+0]==lattice->r[3*1+0]) && (lattice->pos[3*i+1]==lattice->r[3*1+1])  && (lattice->pos[3*i+2]==lattice->r[3*1+2])) i0[1] = lattice->index[2*i+0];
    if((lattice->pos[3*i+0]==lattice->r[3*2+0]) && (lattice->pos[3*i+1]==lattice->r[3*2+1])  && (lattice->pos[3*i+2]==lattice->r[3*2+2])) i0[2] = lattice->index[2*i+0];
    if((lattice->pos[3*i+0]==lattice->r[3*3+0]) && (lattice->pos[3*i+1]==lattice->r[3*3+1])  && (lattice->pos[3*i+2]==lattice->r[3*3+2])) i0[3] = lattice->index[2*i+0];
    if((lattice->pos[3*i+0]==lattice->r[3*4+0]) && (lattice->pos[3*i+1]==lattice->r[3*4+1])  && (lattice->pos[3*i+2]==lattice->r[3*4+2])) i0[4] = lattice->index[2*i+0];
    if((lattice->pos[3*i+0]==lattice->r[3*5+0]) && (lattice->pos[3*i+1]==lattice->r[3*5+1])  && (lattice->pos[3*i+2]==lattice->r[3*5+2])) i0[5] = lattice->index[2*i+0];
    if((lattice->pos[3*i+0]==lattice->r[3*6+0]) && (lattice->pos[3*i+1]==lattice->r[3*6+1])  && (lattice->pos[3*i+2]==lattice->r[3*6+2])) i0[6] = lattice->index[2*i+0];
    if((lattice->pos[3*i+0]==lattice->r[3*7+0]) && (lattice->pos[3*i+1]==lattice->r[3*7+1])  && (lattice->pos[3*i+2]==lattice->r[3*7+2])) i0[7] = lattice->index[2*i+0];
  }
  //for(int i=0;i<8;i++){
  //  cout<<lattice->pos[3*(i0[i])+0]<<'\t'<<lattice->pos[3*(i0[i])+1]<<'\t'<<lattice->pos[3*(i0[i])+2]<<endl;
  //}
  //------------------------------------------------------------------


  //======================================================================================================================================
  //-------------creating 4th order force constant list for 2 atoms of unit cell-----------------------------------------
  int tmp1,tmp2;
  neighbour = lattice->nNeighbour4;
  if(neighbour==1) cutoff=0.1876;
  if(neighbour==2) cutoff=0.5001;
  if(neighbour==3) cutoff=0.6876;
  if(neighbour==4) cutoff=1.0001;
  if(neighbour==5) cutoff=1.1876;
  count1 =0; count2=0;
  for(int i=0;i<((lattice->natom)*(lattice->ncell));i++){

    if(((lattice->pos[3*i0[0]+0]-lattice->pos[3*i+0])*(lattice->pos[3*i0[0]+0]-lattice->pos[3*i+0]) + (lattice->pos[3*i0[0]+1]-lattice->pos[3*i+1])*(lattice->pos[3*i0[0]+1]-lattice->pos[3*i+1]) + (lattice->pos[3*i0[0]+2]-lattice->pos[3*i+2])*(lattice->pos[3*i0[0]+2]-lattice->pos[3*i+2])) < cutoff){
      index1[count1] = lattice->index[2*i+0];
      count1++;
    }

    if(((lattice->pos[3*i0[1]+0]-lattice->pos[3*i+0])*(lattice->pos[3*i0[1]+0]-lattice->pos[3*i+0]) + (lattice->pos[3*i0[1]+1]-lattice->pos[3*i+1])*(lattice->pos[3*i0[1]+1]-lattice->pos[3*i+1]) + (lattice->pos[3*i0[1]+2]-lattice->pos[3*i+2])*(lattice->pos[3*i0[1]+2]-lattice->pos[3*i+2])) < cutoff){
      index2[count2] = lattice->index[2*i+0];
      count2++;
    }
  }


  tmp1=0; tmp2=0;
  for( int i1 =0; i1<count1;i1++){
    for( int i2=0;i2<count1;i2++){
      for( int i3=0; i3<count1;i3++){
	for(int a=1;a<4;a++){
	  for(int b=1;b<4;b++){
	    for(int c=1;c<4;c++){
	      for(int d=1;d<4;d++){
		fulldisp4_out<<i0[0]<<'\t'<<index1[i1]<<'\t'<<index1[i2]<<'\t'<<index1[i3]<<'\t'<<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<endl;
		tmp1++;
	      }
	    }
	  }
	}
      }
    }
  }
  for( int i1 =0; i1<count2;i1++){
    for( int i2=0;i2<count2;i2++){
      for( int i3=0; i3<count2;i3++){
	for(int a=1;a<4;a++){
	  for(int b=1;b<4;b++){
	    for(int c=1;c<4;c++){
	      for(int d=1;d<4;d++){
		fulldisp4_out<<i0[1]<<'\t'<<index2[i1]<<'\t'<<index2[i2]<<'\t'<<index2[i3]<<'\t'<<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<endl;
		tmp2++;
	      }
	    }
	  }
	}
      }
    }
  }

  //------writing same interactions in array------------------------------------------
  dynamic->iorder4 = new  int [8*(tmp1+tmp2)];
  dynamic->norder4 = tmp1+tmp2;
  tmp1=0;
  for( int i1 =0; i1<count1;i1++){
    for( int i2=0;i2<count1;i2++){
      for( int i3=0; i3<count1;i3++){
	for(int a=1;a<4;a++){
	  for(int b=1;b<4;b++){
	    for(int c=1;c<4;c++){
	      for(int d=1;d<4;d++){
		dynamic->iorder4[tmp1*8+0] = i0[0];
		dynamic->iorder4[tmp1*8+1] = index1[i1];
		dynamic->iorder4[tmp1*8+2] = index1[i2];
		dynamic->iorder4[tmp1*8+3] = index1[i3];
		dynamic->iorder4[tmp1*8+4] = a;
		dynamic->iorder4[tmp1*8+5] = b;
		dynamic->iorder4[tmp1*8+6] = c;
		dynamic->iorder4[tmp1*8+7] = d;
		tmp1++;
	      }
	    }
	  }
	}
      }
    }
  }
  for( int i1 =0; i1<count2;i1++){
    for( int i2=0;i2<count2;i2++){
      for( int i3=0; i3<count2;i3++){
	for(int a=1;a<4;a++){
	  for(int b=1;b<4;b++){
	    for(int c=1;c<4;c++){
	      for(int d=1;d<4;d++){
		dynamic->iorder4[tmp1*8+0] = i0[1];
		dynamic->iorder4[tmp1*8+1] = index2[i1];
		dynamic->iorder4[tmp1*8+2] = index2[i2];
		dynamic->iorder4[tmp1*8+3] = index2[i3];
		dynamic->iorder4[tmp1*8+4] = a;
		dynamic->iorder4[tmp1*8+5] = b;
		dynamic->iorder4[tmp1*8+6] = c;
		dynamic->iorder4[tmp1*8+7] = d;
		tmp1++;
	      }
	    }
	  }
	}
      }
    }
  }
  //--------------------------------------------------------------------------------------------------------------------
  //===================================================================================================================================

  //======================================================================================================================================
  //-------------creating 3rd order force constant list for 2 atoms of unit cell-----------------------------------------
  neighbour = lattice->nNeighbour3;
  if(neighbour==1) cutoff=0.1876;
  if(neighbour==2) cutoff=0.5001;
  if(neighbour==3) cutoff=0.6876;
  if(neighbour==4) cutoff=1.0001;
  if(neighbour==5) cutoff=1.1876;
  count1 =0; count2=0;
  for(int i=0;i<((lattice->natom)*(lattice->ncell));i++){

    if(((lattice->pos[3*i0[0]+0]-lattice->pos[3*i+0])*(lattice->pos[3*i0[0]+0]-lattice->pos[3*i+0]) + (lattice->pos[3*i0[0]+1]-lattice->pos[3*i+1])*(lattice->pos[3*i0[0]+1]-lattice->pos[3*i+1]) + (lattice->pos[3*i0[0]+2]-lattice->pos[3*i+2])*(lattice->pos[3*i0[0]+2]-lattice->pos[3*i+2])) < cutoff){
      index1[count1] = lattice->index[2*i+0];
      count1++;
    }

    if(((lattice->pos[3*i0[1]+0]-lattice->pos[3*i+0])*(lattice->pos[3*i0[1]+0]-lattice->pos[3*i+0]) + (lattice->pos[3*i0[1]+1]-lattice->pos[3*i+1])*(lattice->pos[3*i0[1]+1]-lattice->pos[3*i+1]) + (lattice->pos[3*i0[1]+2]-lattice->pos[3*i+2])*(lattice->pos[3*i0[1]+2]-lattice->pos[3*i+2])) < cutoff){
      index2[count2] = lattice->index[2*i+0];
      count2++;
    }
  }
  tmp1=0;tmp2=0;
  for( int i1 =0; i1<count1;i1++){
    for( int i2=0;i2<count1;i2++){
      for(int a=1;a<4;a++){
	for(int b=1;b<4;b++){
	  for(int c=1;c<4;c++){
	    fulldisp3_out<<i0[0]<<'\t'<<index1[i1]<<'\t'<<index1[i2]<<'\t'<<a<<'\t'<<b<<'\t'<<c<<endl;
	    tmp1++;
	  }
	}
      }
    }
  }
  for( int i1 =0; i1<count2;i1++){
    for( int i2=0;i2<count2;i2++){
	for(int a=1;a<4;a++){
	  for(int b=1;b<4;b++){
	    for(int c=1;c<4;c++){
		fulldisp3_out<<i0[1]<<'\t'<<index2[i1]<<'\t'<<index2[i2]<<'\t'<<a<<'\t'<<b<<'\t'<<c<<endl;
		tmp2++;
	      }
	    }
	  }
	}
  }

  // writing the same output in array
  dynamic->iorder3 = new  int [6*(tmp1+tmp2)];
  dynamic->norder3 = tmp1+tmp2;
  tmp1=0;
  for( int i1 =0; i1<count1;i1++){
    for( int i2=0;i2<count1;i2++){
      for(int a=1;a<4;a++){
	for(int b=1;b<4;b++){
	  for(int c=1;c<4;c++){
	    dynamic->iorder3[tmp1*6+0] = i0[0];
	    dynamic->iorder3[tmp1*6+1] = index1[i1];
	    dynamic->iorder3[tmp1*6+2] = index1[i2];
	    dynamic->iorder3[tmp1*6+3] = a;
	    dynamic->iorder3[tmp1*6+4] = b;
	    dynamic->iorder3[tmp1*6+5] = c;
	    tmp1++;
	  }
	}
      }
    }
  }
  for( int i1 =0; i1<count2;i1++){
    for( int i2=0;i2<count2;i2++){
	for(int a=1;a<4;a++){
	  for(int b=1;b<4;b++){
	    for(int c=1;c<4;c++){
	      dynamic->iorder3[tmp1*6+0] = i0[1];
	      dynamic->iorder3[tmp1*6+1] = index2[i1];
	      dynamic->iorder3[tmp1*6+2] = index2[i2];
	      dynamic->iorder3[tmp1*6+3] = a;
	      dynamic->iorder3[tmp1*6+4] = b;
	      dynamic->iorder3[tmp1*6+5] = c;
	      tmp1++;
	      }
	    }
	  }
	}
  }
  //--------------------------------------------------------------------------------------------------------------------
  //===================================================================================================================================

  //======================================================================================================================================
  //-------------creating 2nd order force constant list for 2 atoms of unit cell-----------------------------------------
  neighbour = lattice->nNeighbour2;
  if(neighbour==1) cutoff=0.1876;
  if(neighbour==2) cutoff=0.5001;
  if(neighbour==3) cutoff=0.6876;
  if(neighbour==4) cutoff=1.0001;
  if(neighbour==5) cutoff=1.1876;
  count1 =0; count2=0;
  for(int i=0;i<((lattice->natom)*(lattice->ncell));i++){

    if(((lattice->pos[3*i0[0]+0]-lattice->pos[3*i+0])*(lattice->pos[3*i0[0]+0]-lattice->pos[3*i+0]) + (lattice->pos[3*i0[0]+1]-lattice->pos[3*i+1])*(lattice->pos[3*i0[0]+1]-lattice->pos[3*i+1]) + (lattice->pos[3*i0[0]+2]-lattice->pos[3*i+2])*(lattice->pos[3*i0[0]+2]-lattice->pos[3*i+2])) < cutoff){
      index1[count1] = lattice->index[2*i+0];
      count1++;
    }

    if(((lattice->pos[3*i0[1]+0]-lattice->pos[3*i+0])*(lattice->pos[3*i0[1]+0]-lattice->pos[3*i+0]) + (lattice->pos[3*i0[1]+1]-lattice->pos[3*i+1])*(lattice->pos[3*i0[1]+1]-lattice->pos[3*i+1]) + (lattice->pos[3*i0[1]+2]-lattice->pos[3*i+2])*(lattice->pos[3*i0[1]+2]-lattice->pos[3*i+2])) < cutoff){
      index2[count2] = lattice->index[2*i+0];
      count2++;
    }
  }
  tmp1=0;tmp2=0;
  for( int i1 =0; i1<count1;i1++){
    for(int a=1;a<4;a++){
      for(int b=1;b<4;b++){
	fulldisp2_out<<i0[0]<<'\t'<<index1[i1]<<'\t'<<a<<'\t'<<b<<endl;
	tmp1++;
      }
    }
  }
  for( int i1 =0; i1<count2;i1++){
    for(int a=1;a<4;a++){
      for(int b=1;b<4;b++){
	fulldisp2_out<<i0[1]<<'\t'<<index2[i1]<<'\t'<<a<<'\t'<<b<<endl;
	tmp2++;
      }
    }
  }

  //writing the same thing in array
  dynamic->iorder2 = new  int [4*(tmp1+tmp2)];
  dynamic->norder2 = tmp1+tmp2;
  tmp1=0;
  for( int i1 =0; i1<count1;i1++){
    for(int a=1;a<4;a++){
      for(int b=1;b<4;b++){
	dynamic->iorder2[4*tmp1+0] = i0[0];
	dynamic->iorder2[4*tmp1+1] = index1[i1];
	dynamic->iorder2[4*tmp1+2] = a;
	dynamic->iorder2[4*tmp1+3] = b;
	tmp1++;
      }
    }
  }
  for( int i1 =0; i1<count2;i1++){
    for(int a=1;a<4;a++){
      for(int b=1;b<4;b++){
	dynamic->iorder2[4*tmp1+0] = i0[1];
	dynamic->iorder2[4*tmp1+1] = index2[i1];
	dynamic->iorder2[4*tmp1+2] = a;
	dynamic->iorder2[4*tmp1+3] = b;
	tmp1++;
      }
    }
  }
  //--------------------------------------------------------------------------------------------------------------------
  //===================================================================================================================================

}
//-----------------------------------------------------------------------------------------------------------------------------

//-----------------------writing position for full list sorted in order of distance from center--------------------------------------
void sub_posIndex(LATTICE *lattice){
  //SPECIFIC FOR DIAMOND STRUCTURE WITH SAME ATOMS ONLY

  int count;
  double *tmp;

  //=======creating position file================================================
  count =0;
  for(int i=0;i<lattice->natom;i++){
    for (int j=0;j<lattice->ncell;j++){
      count++;
    }
  }
  lattice->pos = new double [count*3];                       //positions of all atoms arranged as their distance from origin....
  lattice->index = new  int [2*count];              // 1st column as index
  tmp = new double [count];

  count=0;
  for(int i=0;i<lattice->natom;i++){
    for (int j=0;j<lattice->ncell;j++){
      lattice->pos[3*count+0] = lattice->R[3*j+0]+lattice->r[3*i+0];
      lattice->pos[3*count+1] = lattice->R[3*j+1]+lattice->r[3*i+1];
      lattice->pos[3*count+2] = lattice->R[3*j+2]+lattice->r[3*i+2];
      tmp[count] = lattice->pos[3*count+0]*lattice->pos[3*count+0] + lattice->pos[3*count+1]*lattice->pos[3*count+1] + lattice->pos[3*count+2]*lattice->pos[3*count+2];
      count++;
    }
  }
  //======================================================================================

  //sorting distance wise-----------------------------------------------------------------
  bool check;
  double tmpm,tmpx,tmpy,tmpz;
  do{
    check =true;
    for( int i=0;i<count-1;i++){
      if(tmp[i+1]<tmp[i]){
	tmpm = tmp[i]; tmpx = lattice->pos[3*i+0]; tmpy = lattice->pos[3*i+1]; tmpz = lattice->pos[3*i+2];
	tmp[i] = tmp[i+1]; lattice->pos[3*(i)+0] = lattice->pos[3*(i+1)+0]; lattice->pos[3*(i)+1] = lattice->pos[3*(i+1)+1];  lattice->pos[3*(i)+2] = lattice->pos[3*(i+1)+2];
	tmp[i+1] = tmpm; lattice->pos[3*(i+1)+0] = tmpx; lattice->pos[3*(i+1)+1] = tmpy; lattice->pos[3*(i+1)+2] = tmpz;
	check =false;
      }
    }
  }while(check==false);
  for( int i=0;i<count;i++){
    lattice->index[2*i+0] = i;
  }
  //------------------------------------------------------------------------------------

  // writing positions............
  for( int i=0;i<count;i++){
    fullpos_out<<lattice->index[2*i+0]<<'\t';
    fullpos_out<<lattice->pos[3*i+0]<<'\t';
    fullpos_out<<lattice->pos[3*i+1]<<'\t';
    fullpos_out<<lattice->pos[3*i+2]<<'\t';
    fullpos_out<<tmp[i]<<'\n';
  }

  delete [] tmp;
}
//-----------------------------------------------------------------------------------------------------------------------------

//-------------CREATING IBZ FROM FBZ-------------------------------------------------------------------------------------------------------
void sub_ksymmetry(LATTICE *lattice){
  lattice->kmap = new   int [lattice->nlist];

    int *track;
   int *magk;
  magk = new  int [lattice->nlist];
  track =  new   int [lattice->nlist];

  for(int i=0;i<lattice->nlist;i++){              // id and magnitude of k pints.....
    track[i] = i;
    magk[i] = lround((1e+5)*(lattice->klist[3*i+0]*lattice->klist[3*i+0] + lattice->klist[3*i+1]*lattice->klist[3*i+1] + lattice->klist[3*i+2]*lattice->klist[3*i+2]));
  }

  bool check;
   int fswap;
    int iswap;

  do{                  // arranging k in increasing distance from gamma point
    check = true;
    for(int i=0;i<lattice->nlist-1;i++){
      if(magk[i+1]<magk[i]){
	check =false;
	fswap = magk[i];
	iswap = track[i];
	magk[i] = magk[i+1];
	track[i] = track[i+1];
	magk[i+1] = fswap;
	track[i+1] = iswap;
      }
    }
  }while(check==false);


   int unique,val,count,start[10],indi;
  double ax,ay,az,bx,by,bz,sum[10],tmp,mul[10],tmpm;
  unique =0;
  val = 0;
  count=0;

  for(int i=0;i<lattice->nlist;i++){             // counting unique k points.......................
    if((magk[i]!=val)||(i==0)){
      count=0;
      unique++;
      val = magk[i];
      ax = (lattice->klist[3*track[i]+0])*lattice->kx;
      ay = (lattice->klist[3*track[i]+1])*lattice->ky;
      az = (lattice->klist[3*track[i]+2])*lattice->kz;
      ax = (ax>0)?(ax):(-ax);
      ay = (ay>0)?(ay):(-ay);
      az = (az>0)?(az):(-az);
      sum[count] = ax+ay+az;
      mul[count] = ax*ay*az;
      count++;
    }else{
      bx = (lattice->klist[3*track[i]+0])*lattice->kx;
      by = (lattice->klist[3*track[i]+1])*lattice->ky;
      bz = (lattice->klist[3*track[i]+2])*lattice->kz;
      bx = (bx>0)?(bx):(-bx);
      by = (by>0)?(by):(-by);
      bz = (bz>0)?(bz):(-bz);
      tmp = bx+by+bz;
      tmpm = bx*by*bz;
      //data_out<<ax<<endl;
      check =false;
      for( int j=0;j<count;j++){
	if( ((sum[j]-tmp < 0.01)&&(sum[j]-tmp > -(0.01))) &&  ((mul[j]-tmpm < 0.01) && (mul[j]-tmpm > -0.01)) )	  check = true;
      }
      if(check==false){
	unique++;
	sum[count] = tmp;
	mul[count] = tmpm;
	count++;
      }
    }
  }
  //cout<<lattice->nlist<<'\t'<<unique<<endl;

  lattice->IBZ_nlist = unique;                             // storing unique k points and creating mapping
  lattice->IBZ_klist = new double [unique*3];
  unique=0;
  val = 0;
  count=0;;

  for(int i=0;i<lattice->nlist;i++){
    if((magk[i]!=val)||(i==0)){
      count=0;
      lattice->IBZ_klist[3*unique+0] = lattice->klist[3*track[i]+0];
      lattice->IBZ_klist[3*unique+1] = lattice->klist[3*track[i]+1];
      lattice->IBZ_klist[3*unique+2] = lattice->klist[3*track[i]+2];
      lattice->kmap[track[i]] = unique;
      unique++;
      val = magk[i];
      ax = (lattice->klist[3*track[i]+0])*lattice->kx;
      ay = (lattice->klist[3*track[i]+1])*lattice->ky;
      az = (lattice->klist[3*track[i]+2])*lattice->kz;
      ax = (ax>0)?(ax):(-ax);
      ay = (ay>0)?(ay):(-ay);
      az = (az>0)?(az):(-az);
      sum[count] = ax+ay+az;
      mul[count] = ax*ay*az;
      start[count]=i;
      indi = i;
      count++;
    }else{
      bx = (lattice->klist[3*track[i]+0])*lattice->kx;
      by = (lattice->klist[3*track[i]+1])*lattice->ky;
      bz = (lattice->klist[3*track[i]+2])*lattice->kz;
      bx = (bx>0)?(bx):(-bx);
      by = (by>0)?(by):(-by);
      bz = (bz>0)?(bz):(-bz);
      tmp = bx+by+bz;
      tmpm = bx*by*bz;
      //data_out<<ax<<endl;
      check =false;
      for( int j=0;j<count;j++){
	if(((sum[j]-tmp < 0.01)&&(sum[j]-tmp > -(0.01))) &&  ((mul[j]-tmpm < 0.01) && (mul[j]-tmpm > -0.01))  )	{  check = true; indi = start[j];}
      }
      if(check==false){
	lattice->IBZ_klist[3*unique+0] = lattice->klist[3*track[i]+0];
	lattice->IBZ_klist[3*unique+1] = lattice->klist[3*track[i]+1];
	lattice->IBZ_klist[3*unique+2] = lattice->klist[3*track[i]+2];
	lattice->kmap[track[i]] = unique;
	unique++;
	sum[count] = tmp;
	mul[count] = tmpm;
	start[count] = i;
	count++;
      }else{
	lattice->kmap[track[i]] = lattice->kmap[track[indi]];
      }
    }
  }

  /*for(int i=0;i<lattice->nlist;i++){
    data_out<<lattice->klist[3*track[i]+0]<<" "<<lattice->klist[3*track[i]+1]<<" "<<lattice->klist[3*track[i]+2]<<" ";
    data_out<<lattice->IBZ_klist[3*lattice->kmap[track[i]] + 0]<<" "<<lattice->IBZ_klist[3*lattice->kmap[track[i]] + 1]<<" "<<lattice->IBZ_klist[3*lattice->kmap[track[i]] + 2]<<" "<<lattice->kmap[track[i]]<<endl;
  }*/


  lattice->IBZ_kmap= new   int [lattice->IBZ_nlist];
  int j;
  for( int i=0;i<lattice->IBZ_nlist;i++){
    check = false;
    j=0;
    do{
      if(((lattice->IBZ_klist[3*i+0]-lattice->klist[3*j+0] < 1e-6)&&(lattice->IBZ_klist[3*i+0]-lattice->klist[3*j+0] > -(1e-6)))){
	if(((lattice->IBZ_klist[3*i+1]-lattice->klist[3*j+1] < 1e-6)&&(lattice->IBZ_klist[3*i+1]-lattice->klist[3*j+1] > -(1e-6)))){
	  if(((lattice->IBZ_klist[3*i+2]-lattice->klist[3*j+2] < 1e-6)&&(lattice->IBZ_klist[3*i+2]-lattice->klist[3*j+2] > -(1e-6)))){
	    lattice->IBZ_kmap[i]=j;
	    check = true;
	  }
	}
      }
      j++;
    }while(check==false);
  }


  delete [] magk;
  delete [] track;
}
//----------------------------------------------------------------------------------------------------------------------------------------

//---------------------SEARCHING FOR GIVEN VALUE IN A LIST OF POINTS----------------------------------------------------------------------------
int sub_search(long int *sklist,int s, int e, long int tmp){
  // this subroutine does searching only on intergers....
  int out,m;
  if((s-e)>0) return -1;

  m = floor((s+e)/2.);
  if(sklist[m]==tmp) return m;
  //if((sklist[m]>=(tmp-(1e-2)))&&(sklist[m]<=(tmp+(1e-2)))) return m;
  if(sklist[m]>tmp) out = sub_search(sklist,s,m-1,tmp);
  if(sklist[m]<tmp) out = sub_search(sklist,m+1,e,tmp);

  return out;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------

//-------------------- CHECKING MOMENTUM CONSERVATION K1+K2+K3-G = 0-------------------------------------------------------------------
void sub_kinteraction(LATTICE *lattice){
  //.......................This Subrotuine is for fcc only with primitive unit cell................................
  lattice->ikinterac = new   int [lattice->IBZ_nlist*lattice->nlist*3];
  lattice->ninterac = 0;
  long int *sklist, *IBZ_sklist;
    int *track, *IBZ_track;
  sklist = new long int [lattice->nlist];
  track =  new   int [lattice->nlist];
  IBZ_sklist = new long int [lattice->IBZ_nlist];
  IBZ_track =  new   int [lattice->IBZ_nlist];

  //=====generating sklist=================================
  long int mx,my,mz;
  mx =1000; my = 1000000; mz =1000000000;
  for(int i=lattice->nlist-1; i>=0; --i){
    sklist[i] = lround(mx*lattice->klist[3*i+0]) + lround(my*lattice->klist[3*i+1]) + lround(mz*lattice->klist[3*i+2]);
    track[i] = i;
    //data_out<<sklist[i]<<'\n';
  }
  for(int i=lattice->IBZ_nlist-1; i>=0; --i){
    IBZ_sklist[i] = lround(mx*lattice->IBZ_klist[3*i+0]) + lround(my*lattice->IBZ_klist[3*i+1]) + lround(mz*lattice->IBZ_klist[3*i+2]);
    IBZ_track[i] = i;
    //data_out<<sklist[i]<<'\n';
  }
  //=========================================================

  //=====sorting the sklist in ascending order
  bool check;
  long int tmp;
  int t;
  do{
    check =true;
    for(int i=lattice->nlist-2;i>=0;--i){
      if(sklist[i+1]<sklist[i]){
	tmp = sklist[i+1];
	t = track[i+1];
	sklist[i+1] = sklist[i];
	track[i+1] = track[i];
	sklist[i] = tmp;
	track[i] = t;
	check =false;
      }
    }
  }while(check==false);
  do{
    check =true;
    for(int i=lattice->IBZ_nlist-2;i>=0;--i){
      if(IBZ_sklist[i+1]<IBZ_sklist[i]){
	tmp = IBZ_sklist[i+1];
	t = IBZ_track[i+1];
	IBZ_sklist[i+1] = IBZ_sklist[i];
	IBZ_track[i+1] = IBZ_track[i];
	IBZ_sklist[i] = tmp;
	IBZ_track[i] = t;
	check =false;
      }
    }
  }while(check==false);
  //============================================

  //=======STARTING CHECKING MOMENTUM===========
  long int tmpx,tmpy,tmpz,tmp1,tmp2,tmp3;
  long int *dtmp;
  int res,s,e,i1,i2,ginterac,gsize,count;
  ginterac = 2;
  i2=0;
  gsize= round(pow((((ginterac)*2)+1),3));
  dtmp = new long int [gsize];
  for (int r1=-ginterac;r1<=ginterac;++r1){                             // G vectors to check with.....everything is to be scaled by 2*pi/ax
    for(int r2=-ginterac;r2<=ginterac;++r2){
      for(int r3=-ginterac;r3<=ginterac;++r3){
	tmpx = lround(mx*((-1)*r1 + (1)*r2 + (1)*r3));
	tmpy = lround(my*((1)*r1 + (-1)*r2 + (1)*r3));
	tmpz = lround(mz*((1)*r1 + (1)*r2 + (-1)*r3));
	tmp = tmpx+tmpy+tmpz;                                          //desired value from combination of sklist
	dtmp[i2] = tmp;
	i2++;
      }
    }
  }



	for(i1 = lattice->IBZ_nlist-1;i1>=0;--i1){
	  tmp1 = IBZ_sklist[i1];
	  for(i2 = lattice->nlist-1;i2>=0;--i2){
	    check =false;
	    count=0;
	      do{
		tmp2 = sklist[i2];
		tmp3 = -tmp1 - tmp2 + dtmp[count];
		s = 0; e = lattice->nlist-1;
		res = sub_search(sklist,s,e,tmp3);
		if(res>=0){
		    lattice->ikinterac[3*lattice->ninterac+0] = IBZ_track[i1]; lattice->ikinterac[3*lattice->ninterac+1] = track[i2]; lattice->ikinterac[3*lattice->ninterac+2] = track[res];
		    lattice->ninterac++;
		    check = true;
		}
	      count++;
	    }while((count<gsize)&&(check==false));
	  }
	}
	//================================
  //============================================


  lattice->IBZ_num_interac = new   int [lattice->IBZ_nlist];
  int ttt1,ttt2,ttt3,ttt4;
  ttt1 = lattice->ikinterac[3*0+0];
  ttt3=1;
  ttt4=0;
  for( int i=1;i<lattice->ninterac;i++){
    ttt2 = lattice->ikinterac[3*i+0];
    if(ttt1==ttt2){
      ttt3++;
    }else{
      lattice->IBZ_num_interac[ttt4] = ttt3;
      ttt3=1;
      ttt1=ttt2;
      ttt4++;
    }
  }lattice->IBZ_num_interac[ttt4] = ttt3;


  delete [] sklist;
  delete [] track;
  delete [] IBZ_sklist;
  delete [] IBZ_track;
  delete [] dtmp;
}
//------------------------------------------------------------------------------------------------------------------------------------------------

//---------------------------------- CONSTRUCTING RECIPROCAL STRUCTURE-----------------------------------------------------------------------------
void sub_klist(LATTICE *lattice){
  //this subroutine is for fcc with primitive unit cell only.....
  double b1[3],b2[3],b3[3];
  lattice->klist = new double [lattice->nlist*3];


  if((lattice->type==2)&&(lattice->isConv==0)){                                        // primitive fcc cell
    b1[0] = -1.; b1[1] = 1.; b1[2] = 1.;
    b2[0] = 1.; b2[1] = -1.; b2[2] = 1.;
    b3[0] = 1.; b3[1] = 1.; b3[2] = -1.;
    int i=0;
    for(int i1 = -round(lattice->kx/2)+1; i1<round(lattice->kx/2)+1;i1++){
      for(int i2 = -round(lattice->ky/2)+1; i2<round(lattice->ky/2)+1;i2++){
        for(int i3 = -round(lattice->kz/2)+1; i3<round(lattice->kz/2)+1;i3++){
          lattice->klist[3*i+0] = ((i1*b1[0])/(lattice->kx)) + ((i2*b2[0])/(lattice->ky)) +((b3[0]*i3)/(lattice->kz));
          lattice->klist[3*i+1] = ((i1*b1[1])/(lattice->kx)) + ((i2*b2[1])/(lattice->ky)) +((b3[1]*i3)/(lattice->kz));
          lattice->klist[3*i+2] = ((i1*b1[2])/(lattice->kx)) + ((i2*b2[2])/(lattice->ky)) +((b3[2]*i3)/(lattice->kz));
          i++;
        }
      }
    }

    //   checking if point is in first brillioun zone
    double p [14][3];
    p[0][0] = 1.; p[0][1]=1.; p[0][2]=1.;
    p[1][0] = 1.; p[1][1]=1.; p[1][2]=-1.;
    p[2][0] = 1.; p[2][1]=-1.; p[2][2]=1.;
    p[3][0] = 1.; p[3][1]=-1.; p[3][2]=-1.;
    p[4][0] = -1.; p[4][1]=1.; p[4][2]=1.;
    p[5][0] = -1.; p[5][1]=1.; p[5][2]=-1.;
    p[6][0] = -1.; p[6][1]=-1.; p[6][2]=1.;
    p[7][0] = -1.; p[7][1]=-1.; p[7][2]=-1.;
    p[8][0] =2*1.; p[8][1] = 0; p[8][2]=0;
    p[9][0] =-2*1.; p[9][1] = 0; p[9][2]=0;
    p[10][0]=0; p[10][1] = 2*1.; p[10][2]=0;
    p[11][0] =0; p[11][1] = -2*1.; p[11][2]=0;
    p[12][0] =0; p[12][1] = 0; p[12][2]=2*1.;
    p[13][0] =0; p[13][1] = 0; p[13][2]=-2*1.;

    bool check;
    double tmp1,tmp2;
    for(int i1 = 0; i1<lattice->nlist;i1++){
      do{
        check = true;
        for(int i2=0;i2<14;i2=i2+2){
          tmp1 = pow((lattice->klist[3*i1+0] - p[i2][0]),2) + pow((lattice->klist[3*i1+1] - p[i2][1]),2) + pow((lattice->klist[3*i1+2] - p[i2][2]),2);
          tmp2 = pow((lattice->klist[3*i1+0] ),2) + pow((lattice->klist[3*i1+1]),2) + pow((lattice->klist[3*i1+2] ),2);
          if(tmp2>tmp1){
            check = false;
            lattice->klist[3*i1+0] -=p[i2][0];
            lattice->klist[3*i1+1] -=p[i2][1];
            lattice->klist[3*i1+2] -=p[i2][2];
          }
        }
        for(int i2=1;i2<14;i2=i2+2){
          tmp1 = pow((lattice->klist[3*i1+0] - p[i2][0]),2) + pow((lattice->klist[3*i1+1] - p[i2][1]),2) + pow((lattice->klist[3*i1+2] - p[i2][2]),2);
          tmp2 = pow((lattice->klist[3*i1+0] ),2) + pow((lattice->klist[3*i1+1]),2) + pow((lattice->klist[3*i1+2] ),2);
          if(tmp2>tmp1){
            check = false;
            lattice->klist[3*i1+0] -=p[i2][0];
            lattice->klist[3*i1+1] -=p[i2][1];
            lattice->klist[3*i1+2] -=p[i2][2];
          }
        }
      }while(check==false);
    }
  }
}
//---------------------------------------------------------------------------------------------------------------------------------------------------

//---------------------- constructing real lattice----------------------------------------------------------------------------------------------
void sub_readInput(LATTICE *lattice){
  //this is a general subroutine....valid for all types and formats.....

  std::string line;
    int type,isConv;

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> lattice->ax >> lattice->ay >> lattice->az;                      // unit cell size in x,y,z direction

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> type;                                                          // lattice type

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> isConv;                                                          // unit cell  type for reciprocal space

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> lattice->natom;                                               // # atoms in unit cell

  lattice->r = new double [lattice->natom*3];
  lattice->m = new double [lattice->natom];


  // lattice vectors for fcc cell
  double a1[3],a2[3],a3[3];


  if((lattice->natom==2)&&(type==2)){
    std::getline(data_in,line);
    std::getline(data_in,line);
    std::getline(data_in,line);
    data_in>> lattice->r[3*0+0] >> lattice->r[3*0+1] >> lattice->r[3*0+2] >> lattice->m[0];
    std::getline(data_in,line);
    data_in>> lattice->r[3*1+0] >> lattice->r[3*1+1] >> lattice->r[3*1+2] >> lattice->m[1];

    // lattice vectors for fcc cell
    a1[0] = 0.5; a1[1] = 0.5; a1[2] = 0.0;
    a2[0] = 0.0; a2[1] = 0.5; a2[2] = 0.5;
    a3[0] = 0.5; a3[1] = 0.0; a3[2] = 0.5;

    for(int i=0;i<lattice->natom;i++){
      lattice->r[3*i+0] *= lattice->ax;
      lattice->r[3*i+1] *= lattice->ay;
      lattice->r[3*i+2] *= lattice->az;
    }
  }

  if((lattice->natom==8)&&(type==2)){
    std::getline(data_in,line);
    std::getline(data_in,line);
    std::getline(data_in,line);
    data_in>> lattice->r[3*0+0] >> lattice->r[3*0+1] >> lattice->r[3*0+2] >> lattice->m[0];
    std::getline(data_in,line);
    data_in>> lattice->r[3*1+0] >> lattice->r[3*1+1] >> lattice->r[3*1+2] >> lattice->m[1];
    std::getline(data_in,line);
    data_in>> lattice->r[3*2+0] >> lattice->r[3*2+1] >> lattice->r[3*2+2] >> lattice->m[2];
    std::getline(data_in,line);
    data_in>> lattice->r[3*3+0] >> lattice->r[3*3+1] >> lattice->r[3*3+2] >> lattice->m[3];
    std::getline(data_in,line);
    data_in>> lattice->r[3*4+0] >> lattice->r[3*4+1] >> lattice->r[3*4+2] >> lattice->m[4];
    std::getline(data_in,line);
    data_in>> lattice->r[3*5+0] >> lattice->r[3*5+1] >> lattice->r[3*5+2] >> lattice->m[5];
    std::getline(data_in,line);
    data_in>> lattice->r[3*6+0] >> lattice->r[3*6+1] >> lattice->r[3*6+2] >> lattice->m[6];
    std::getline(data_in,line);
    data_in>> lattice->r[3*7+0] >> lattice->r[3*7+1] >> lattice->r[3*7+2] >> lattice->m[7];

    // lattice vectors for fcc cell
    a1[0] = 1.0; a1[1] = 0.0; a1[2] = 0.0;
    a2[0] = 0.0; a2[1] = 1.0; a2[2] = 0.0;
    a3[0] = 0.0; a3[1] = 0.0; a3[2] = 1.0;

    for(int i=0;i<lattice->natom;i++){
      lattice->r[3*i+0] *= lattice->ax;
      lattice->r[3*i+1] *= lattice->ay;
      lattice->r[3*i+2] *= lattice->az;
    }
  }
  //cout<<lattice->r[3*7+0]<<'\t'<<lattice->r[3*7+1]<<'\t'<<lattice->r[3*7+2]<<endl;


  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> lattice->Nx>>lattice->Ny>>lattice->Nz>>lattice->nNx>>lattice->nNy>>lattice->nNz;                                                            // # unit cells in each direction
  //lattice->Nx = (lattice->Nx/2); lattice->Ny = (lattice->Ny/2); lattice->Nz = (lattice->Nz/2);


  lattice->ncell = (((lattice->Nx+lattice->nNx+1)*(lattice->Ny+lattice->nNy+1)*(lattice->Nz+lattice->nNz+1)));                                        // 8 because 2*2*2 as lattice->Nx =lattice->Nx_required/2 and 1 to include central cell as well

  lattice->R = new double [(lattice->ncell)*3];
  lattice->R[3*0+0] = 0.0;
  lattice->R[3*0+1] = 0.0;
  lattice->R[3*0+1] = 0.0;                                               // central cell
  int i=1;                                                               // generating cloud around central cell.

  for(int i1 = -lattice->nNx; i1<lattice->Nx+1;i1++){
    for(int i2 = -lattice->nNy; i2<lattice->Ny+1;i2++){
      for(int i3 = -lattice->nNz; i3<lattice->Nz+1;i3++){
	if(!((i1==0)&&(i2==0)&&(i3==0))){
	  lattice->R[3*i+0] = lattice->ax*(i1*a1[0] + i2*a2[0] + i3*a3[0]);
	  lattice->R[3*i+1] = lattice->ay*(i1*a1[1] + i2*a2[1] + i3*a3[1]);
	  lattice->R[3*i+2] = lattice->az*(i1*a1[2] + i2*a2[2] + i3*a3[2]);
	  i++;
	}
      }
    }
  }


  lattice->Nx = lattice->Nx + lattice->nNx + 1 ;
  lattice->Ny = lattice->Ny + lattice->nNy + 1 ;
  lattice->Nz = lattice->Nz + lattice->nNz + 1 ;
  lattice->type = type;
  lattice->isConv = isConv;

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>>lattice->kx>>lattice->ky>>lattice->kz;

  lattice->nlist = (lattice->kx*lattice->ky*lattice->kz);
  lattice->kx_scale = 2.*pi/(lattice->ax);
  lattice->ky_scale = 2.*pi/(lattice->ay);
  lattice->kz_scale = 2.*pi/(lattice->az);

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> numthread;
  //cout<<lattice->ncell<<endl;

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> lattice->nNeighbour2>> lattice->nNeighbour3>> lattice->nNeighbour4;
  //cout<<lattice->nNeighbour<<endl;

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> temperature;

  std::getline(data_in,line);
  std::getline(data_in,line);
  std::getline(data_in,line);
  data_in>> FC_available;   // checking whether to run full code or only FC part...
}
//--------------------------------------------------------------------------------------------------------------------------------------------------
//=================================================================================================================================================

