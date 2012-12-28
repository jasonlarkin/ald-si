#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

ofstream structure_out("structure.dat");
ifstream supercell_in;
ifstream displacement_in;

int main(){
   double delta = 0.01;
   int precision;
   precision = 15;

   char line[512], word[512];
   stringstream iss;
   
   double a1[3], a2[3], a3[3];
   double non_mass, non_length;
   double *pos;
   int *atomid;
   int natom;
   double x1,y1,z1;

   //......reading structure from file...............................
   supercell_in.open("./supercell.dat",ifstream::in);
   if(!supercell_in) cout<<"Error in opening input file: supercell.dat "<<endl;
   
   supercell_in.getline(line,512,'\n');        //comment line
   supercell_in.getline(line,512,'\n');        //comment line
   iss.clear();
   iss << line;
   iss >> non_length >> non_mass;
   
   supercell_in.getline(line,512,'\n');        //comment line
   iss.clear();
   iss << line;
   iss >> a1[0] >> a1[1] >>a1[2];
   
   supercell_in.getline(line,512,'\n');        //comment line
   iss.clear();
   iss << line;
   iss >> a2[0] >> a2[1] >>a2[2];
   
   supercell_in.getline(line,512,'\n');        //comment line
   iss.clear();
   iss << line;
   iss >> a3[0] >> a3[1] >>a3[2];
   supercell_in.getline(line,512,'\n');        //comment line
   
   natom = 0;
   pos = new double [10000*4];
   atomid = new int [10000];
   while(supercell_in.getline(line,512,'\n')){
      iss.clear();
      iss << line;
      iss >> atomid[natom] >> pos[4*natom+0] >> pos[4*natom+1] >> pos[4*natom+2] >>pos[4*natom+3];
      natom++;
      if(natom==10000) cout<<"Please increase array size for atomid and pos"<<endl;
   }


   for(int i=0;i<natom;i++){
      x1 = (pos[4*i+0]*a1[0] + pos[4*i+1]*a2[0]  + pos[4*i+2]*a3[0])*(non_length/(1e-10)); 
      y1 = (pos[4*i+0]*a1[1] + pos[4*i+1]*a2[1]  + pos[4*i+2]*a3[1])*(non_length/(1e-10)); 
      z1 = (pos[4*i+0]*a1[2] + pos[4*i+1]*a2[2]  + pos[4*i+2]*a3[2])*(non_length/(1e-10));
      pos[4*i+0] = x1;
      pos[4*i+1] = y1;
      pos[4*i+2] = z1;
   }
   for(int i=0;i<3;i++){
      a1[i] *= (non_length/(1e-10));
      a2[i] *= (non_length/(1e-10));
      a3[i] *= (non_length/(1e-10));
   }
   non_mass = non_mass*(1e+3);
   //............structure read..........................................


   //.........reading displacements.....................................
   int ndisp, *dispid;
   int ndisp2, ndisp3, ndisp4;
   int *displacement;
   int tmp[6], count;

   ndisp = 0;
   ndisp2 = 0;
   ndisp3 = 0;
   ndisp4 = 0;
   displacement = new int [8*1000];
   dispid = new int [1000];

   displacement_in.open("./displacement.dat", ifstream::in);
   if(!displacement_in) cout<<"Error in opening displacement file: displacement.dat"<<endl;

   displacement_in.getline(line,512,'\n');
   displacement_in.getline(line,512,'\n');

   while(displacement_in.getline(line,512,'\n')){  
      iss.clear();
      iss << line;
      iss.getline(word,512,'\t');
      count = 0;
      while(iss >> tmp[count]){
         count++;
      }
      if(count==2){
         dispid[ndisp] = ndisp;
         displacement[8*ndisp+0] = tmp[0];
         displacement[8*ndisp+4] = tmp[1];
         ndisp2++;
      }
      if(count==4){
         dispid[ndisp] = ndisp;
         displacement[8*ndisp+0] = tmp[0];
         displacement[8*ndisp+1] = tmp[1];
         displacement[8*ndisp+4] = tmp[2];
         displacement[8*ndisp+5] = tmp[3];
         ndisp3++;
      }
      if(count==6){
         dispid[ndisp] = ndisp;
         displacement[8*ndisp+0] = tmp[0];
         displacement[8*ndisp+1] = tmp[1];
         displacement[8*ndisp+2] = tmp[2];
         displacement[8*ndisp+4] = tmp[3];
         displacement[8*ndisp+5] = tmp[4];
         displacement[8*ndisp+6] = tmp[5];
         ndisp4++;
      }

      ndisp++;
      if(ndisp==1000) cout<<"Please increase array size for displacement"<<endl;
   }

   //..........displacements read.......................................


   //..........now creating structure files...........................
   string filename;
   int d1, d2, d3, d4;
   int id1, id2, id3, id4;
   int s1,s2,s3,s4;
   
   for(int i =0; i<ndisp2;i++){
      stringstream steam;
      steam<<"disp_"<<i<<".dat";
      filename = steam.str();
            
      for(int j=0;j<natom;j++){
         if(atomid[j]==displacement[8*i+0]){
            id1 = j;
            d1 = displacement[8*i+4];
            d1 = (d1>0)?(d1):(-d1);
            s1 = displacement[8*i+4]/d1;
         }
      }

      d1--;
      pos[4*id1 + d1] += delta*s1;
     
      //...now writing into file.....................
      ofstream output;
      output.open(filename.c_str());
      output<<"Comment"<<'\n';
      output<<" "<<'\n';
      output<<'\n';
      output<<natom<<" atoms"<<'\n';
      output<<"1 atom types"<<'\n';
      output<<'\n';
      output<<"0 "<<setprecision(precision)<<a1[0]<<" xlo xhi"<<'\n';
      output<<"0 "<<setprecision(precision)<<a2[1]<<" ylo yhi"<<'\n';
      output<<"0 "<<setprecision(precision)<<a3[2]<<" zlo zhi"<<'\n';
      output<<'\n';
      output<<"Masses"<<'\n';
      output<<'\n';
      output<<"1 "<<28.0855<<'\n';
      output<<'\n';
      output<<"Atoms"<<'\n';
      output<<'\n';
      for(int l=0;l<natom;l++){
         output<<l+1<<" 1 "<<((int)(round(1e+5*pos[4*l+0])))/(1e+5)<<'\t'<<setprecision(precision)<<pos[4*l+1]<<'\t'<<setprecision(precision)<<pos[4*l+2]<<'\n';
      }


      output.close();
      //..............................................

      pos[4*id1 + d1] -= delta*s1;

   }

   for(int i = ndisp2; i< ndisp2 + ndisp3;i++){
      stringstream steam;
      steam<<"disp_"<<i<<".dat";
      filename = steam.str();
            
      for(int j=0;j<natom;j++){
         if(atomid[j]==displacement[8*i+0]){
            id1 = j;
            d1 = displacement[8*i+4];
            d1 = (d1>0)?(d1):(-d1);
            s1 = displacement[8*i+4]/d1;
         }
      }
      for(int j=0;j<natom;j++){
         if(atomid[j]==displacement[8*i+1]){
            id2 = j;
            d2 = displacement[8*i+5];
            d2 = (d2>0)?(d2):(-d2);
            s2 = displacement[8*i+5]/d2;
         }
      }

      d1--;
      d2--;
      pos[4*id1 + d1] += delta*s1;
      pos[4*id2 + d2] += delta*s2;
     
      //...now writing into file.....................
      ofstream output;
      output.open(filename.c_str());
      output<<"Comment"<<'\n';
      output<<" "<<'\n';
      output<<'\n';
      output<<natom<<" atoms"<<'\n';
      output<<"1 atom types"<<'\n';
      output<<'\n';
      output<<"0 "<<setprecision(precision)<<a1[0]<<" xlo xhi"<<'\n';
      output<<"0 "<<setprecision(precision)<<a2[1]<<" ylo yhi"<<'\n';
      output<<"0 "<<setprecision(precision)<<a3[2]<<" zlo zhi"<<'\n';
      output<<'\n';
      output<<"Masses"<<'\n';
      output<<'\n';
      output<<"1 "<<non_mass<<'\n';
      output<<'\n';
      output<<"Atoms"<<'\n';
      output<<'\n';
      for(int l=0;l<natom;l++){
         output<<l+1<<" 1 "<<((int)(round(1e+5*pos[4*l+0])))/(1e+5)<<'\t'<<setprecision(precision)<<pos[4*l+1]<<'\t'<<setprecision(precision)<<pos[4*l+2]<<'\n';
      }


      output.close();
      //..............................................

      pos[4*id1 + d1] -= delta*s1;
      pos[4*id2 + d2] -= delta*s2;
   }

   for(int i = ndisp2 + ndisp3; i< ndisp2 + ndisp3 + ndisp4;i++){
      stringstream steam;
      steam<<"disp_"<<i<<".dat";
      filename = steam.str();
            
      for(int j=0;j<natom;j++){
         if(atomid[j]==displacement[8*i+0]){
            id1 = j;
            d1 = displacement[8*i+4];
            d1 = (d1>0)?(d1):(-d1);
            s1 = displacement[8*i+4]/d1;
         }
      }
      for(int j=0;j<natom;j++){
         if(atomid[j]==displacement[8*i+1]){
            id2 = j;
            d2 = displacement[8*i+5];
            d2 = (d2>0)?(d2):(-d2);
            s2 = displacement[8*i+5]/d2;
         }
      }
      for(int j=0;j<natom;j++){
         if(atomid[j]==displacement[8*i+2]){
            id3 = j;
            d3 = displacement[8*i+6];
            d3 = (d3>0)?(d3):(-d3);
            s3 = displacement[8*i+6]/d3;
         }
      }

      d1--;
      d2--;
      d3--;
      pos[4*id1 + d1] += delta*s1;
      pos[4*id2 + d2] += delta*s2;
      pos[4*id3 + d3] += delta*s3;
     
      //...now writing into file.....................
      ofstream output;
      output.open(filename.c_str());
      output<<"Comment"<<'\n';
      output<<" "<<'\n';
      output<<'\n';
      output<<natom<<" atoms"<<'\n';
      output<<"1 atom types"<<'\n';
      output<<'\n';
      output<<"0 "<<setprecision(precision)<<a1[0]<<" xlo xhi"<<'\n';
      output<<"0 "<<setprecision(precision)<<a2[1]<<" ylo yhi"<<'\n';
      output<<"0 "<<setprecision(precision)<<a3[2]<<" zlo zhi"<<'\n';
      output<<'\n';
      output<<"Masses"<<'\n';
      output<<'\n';
      output<<"1 "<<non_mass<<'\n';
      output<<'\n';
      output<<"Atoms"<<'\n';
      output<<'\n';
      for(int l=0;l<natom;l++){
         output<<l+1<<" 1 "<<((int)(round(1e+5*pos[4*l+0])))/(1e+5)<<'\t'<<setprecision(precision)<<pos[4*l+1]<<'\t'<<setprecision(precision)<<pos[4*l+2]<<'\n';
      }


      output.close();
      //..............................................

      pos[4*id1 + d1] -= delta*s1;
      pos[4*id2 + d2] -= delta*s2;
      pos[4*id3 + d3] -= delta*s3;
   }
   //..................................................................



}
