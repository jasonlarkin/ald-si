/*								            LDGuess.cpp		                        */
/*                            07/02/2007                            */
/*********************************************************************
*    This program provides a guess for the frequency shift and       *
*  linewidth to use in an anharmonic lattice dynamics calculation.   *                            *
*********************************************************************/

#include "LDCode.h"
#include "Anharmonic.h"
#include <cmath>


void ANHARMONIC::LDGuess(void) {
  //CREATE FILES WITH INITIAL GUESSES
  Log  <<"\nBuilding Files of Initial Guesses..."<<endl;
  cout <<"\nBuilding Files of Initial Guesses..."<<endl;
  Log  <<"  For quantum system from ";
  cout <<"  For quantum system from ";
  ofstream SW("Shift_Width0.txt");
  if(Q_Interpolate.size()==0) {Guess(SW);}
  else {Interpolate(Q_Interpolate, SW);}
  SW.close();
  Log  <<"  For classical system from ";
  cout <<"  For classical system from ";
  ofstream CSW("Classical_Shift_Width0.txt");
  if(C_Interpolate.size()==0) {Guess(CSW);}
  else {Interpolate(C_Interpolate, CSW);}
  CSW.close();

  return;
}


void ANHARMONIC::Guess(ofstream &Out) {
  Out  <<"Initial Guess from Harmonic Frequencies"<<endl;
  Log  <<"harmonic frequencies."<<endl;
  cout <<"harmonic frequencies."<<endl;
  double qh;
  ifstream Freq("Frequency.txt");
  Freq.ignore(10000, '\n');
  while( (Freq.peek()=='\n')||(Freq.peek()=='\t')||(Freq.peek()==' ') ) {Freq.ignore(1, '@');}
  while(!Freq.eof()) {
    Freq >>qh;
    if(qh>=1.0e-5) {Out <<qh<<"\t0.0\t"<<qh*fs_guess<<'\t'<<0.5/(qh*lw_guess)<<endl;}
    else {Out <<qh<<"\t0.0\t"<<qh*fs_guess<<'\t'<<0.0<<endl;}
    while( (Freq.peek()=='\n')||(Freq.peek()=='\t')||(Freq.peek()==' ') ) {Freq.ignore(1, '@');}
  }

  Freq.close();
  return;
}


void ANHARMONIC::Interpolate(string In_str, ofstream &Out) {
  Out  <<"Interpolation"<<endl;
  Log  <<"interpolating file '"<<In_str<<"'"<<endl;
  cout <<"interpolating file '"<<In_str<<"'"<<endl;
  //DECLARE LOCAL VARIABLES
  int rows, r, i, j;
  int top, bot;
  double qh, fs, lw=0.0;
  double top_qh, bot_qh;
  ifstream In(In_str.c_str());
  if(!In.is_open()) {Log <<"Could not open interpolation file '"<<In_str<<"'."<<endl;exit(0);}


	//COUNT ROWS OF INPUT DATA
	Log  <<"   * Counting lines..."<<endl;
	cout <<"   * Counting lines..."<<endl;
	In.ignore(int(1e6), '\n');
	rows = -1;
	while(++rows<1e9) {
		In.ignore(int(1e6), '\n');
		if(In.peek()==EOF) {break;}
	}
	if(rows>=1e9) {Log <<"Interpolation file exceeds size limit."<<endl;exit(0);}


	//ALLOCATE DYNAMIC MEMORY
	double *QH = new double[rows];
	double *FS = new double[rows];
	double *LW = new double[rows];


	//RETURN TO BEGINNING OF INPUT FILE, READ DATA, AND ORDER
	Log <<"   * Reading data..."<<endl;
	cout <<"   * Reading data..."<<endl;
	In.clear();
	In.seekg(0, ios::beg);
	In.ignore(int(1e6), '\n');
	for(r=0;r<rows;r++) {
	  In >>QH[r]>>FS[r]>>fs>>LW[r];
	  FS[r] += fs;
	  if(LW[r]>1.0e-5) {LW[r] = 1.0/LW[r];}
	  else {LW[r] = 0.0;}
	  for(i=0;i<r;i++) {
	    if(QH[r]<QH[i]) {
	      qh = QH[r];
	      fs = FS[r];
	      lw = LW[r];
	      for(j=r;j>i;j--) {
	        QH[j] = QH[j-1];
	        FS[j] = FS[j-1];
	        LW[j] = LW[j-1];
        }
	      QH[i] = qh;
	      FS[i] = fs;
	      LW[i] = lw;
	      break;
	    }
	  }
  }
	In.close();


	//REDUCE DATA BY ELIMINATING MULTIPLE FREQUENCIES
	r = rows;
	while(--r>0) {
	  i = r-1;
	  fs = FS[r];
	  lw = LW[r];
	  while(fabs(QH[r]-QH[i])<1.0e-5) {
	    fs += FS[i];
	    lw += LW[i];
	    i -= 1;
	    if(i<0) {break;}
	  }
	  if((top=r-i)<=1) {continue;}
	  i += 1;
	  FS[i] = fs/double(top);
	  LW[i] = lw/double(top);
	  j = r+1;
	  r = i;
	  while((i+=1)<rows) {
	    QH[i] = QH[j];
	    FS[i] = FS[j];
	    LW[i] = LW[j];
	    j += 1;
	  }
	  rows = rows - top + 1;
	}


	//INTERPOLATE FOR NEW VALUES
	Log <<"   * Interpolating values..."<<endl;
	cout <<"   * Interpolating values..."<<endl;
	ifstream Freq("Frequency.txt");
	Freq.ignore(int(1e6), '\n');
  while( (Freq.peek()=='\n')||(Freq.peek()=='\t')||(Freq.peek()==' ') ) {Freq.ignore(1, '@');}
  double last = -100.0;
  while(!Freq.eof()) {
    Freq >>qh;
    if(fabs(qh-last)>1.0e-5) {
      //Bisection
      top = rows;
      bot = -1;
      if(qh>QH[rows-1]) {bot=rows-1;}
      else if(qh<QH[0]) {top=0;}
      while(top-bot>1) {
        i = (top+bot) >> 1;   //Divide by 2
        if(qh==QH[i]) {
          fs = FS[i];
          if(qh<1.0e-5) {lw = 0.0;}
          else {lw = 1.0/LW[i];}
          Out <<qh<<"\t0.0\t"<<fs<<'\t'<<lw<<endl;
          last = qh;
          while( (Freq.peek()=='\n')||(Freq.peek()=='\t')||(Freq.peek()==' ') ) {Freq.ignore(1, '@');}
          if(Freq.eof()) {return;}
          Freq >>qh;
          while(fabs(qh-last)<1.0e-5) {
            Out <<qh<<"\t0.0\t"<<fs<<'\t'<<lw<<endl;
            while( (Freq.peek()=='\n')||(Freq.peek()=='\t')||(Freq.peek()==' ') ) {Freq.ignore(1, '@');}
            if(Freq.eof()) {return;}
            Freq >>qh;
          }
          top = rows;
          bot = -1;
          if(qh>QH[rows-1]) {bot=rows-1;}
          else if(qh<QH[0]) {top=0;}
        }
        else if(qh>QH[i]) {bot = i;}
        else {top = i;}
      }
      if(top>=rows) {fs=FS[bot];  lw=LW[bot];}
      else if(bot<0) {fs=FS[0];  lw=LW[0];}
      else {
        bot_qh = qh-QH[bot];
        top_qh = QH[top] - QH[bot];
        fs = bot_qh*(FS[top]-FS[bot])/top_qh + FS[bot];
        lw = bot_qh*(LW[top]-LW[bot])/top_qh + LW[bot];
      }
      if(qh<1.0e-5) {lw = 0.0;}
      else {lw = 1.0/lw;}
    }
    Out <<qh<<"\t0.0\t"<<fs<<'\t'<<lw<<endl;
    last = qh;
    while( (Freq.peek()=='\n')||(Freq.peek()=='\t')||(Freq.peek()==' ') ) {Freq.ignore(1, '@');}
  }

  Freq.close();


  return;
}
