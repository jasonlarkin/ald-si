/*                       BuildNeighborList.cpp                      */
/*                            12/02/2008                            */
/*********************************************************************
*    Source code file for functions in the NEIGHBOR class.           *
*********************************************************************/

/*DEFINE HEADERS*/
#include "LDCode.h"
#include <cmath>


/*DECLARE PREPROCESSOR COMMANDS*/
//#define ZERO 1.0e-5


/*PRIVATE SUBROUTINE Two_Body_List*/
N_LIST *BuildNeighborList_R(int &b0, int Mat, double cutoff) {
	//DECLARE LOCAL VARIABLES
	int n_neig;                 //Counts number of neighbors
	int b1;                     //Atomic counter
	int d0, d1;                 //Directional counters
	int natom=Unit_Cell->natom; //Copy number of atoms in unit cell
	int *n_UC;                  //Copy of unit cells
	int *L1;                    //Unit cell counter
	double *r_L1, *r_01;        //Separation distances
	double *X0, *X1;            //Atomic positions
	double R2;                  //Square of separation distance/generic double
	N_LIST *start, *end;        //Linked lists


	//ALLOCATE MEMORY AND INITIALIZE VARIABLES
	//Allocate Memory
	n_UC = new int[DIM];
	L1 = new int[DIM];
	r_L1 = new double[DIM];
	r_01 = new double[DIM];
	X0 = new double[DIM];
	X1 = new double[DIM];
	//Find Appropriate Values for Unit Cell Range
	for(d0=0;d0<DIM;d0++) {
	  R2 = 0.0;
	  for(d1=0;d1<DIM;d1++) {R2 += (Lattice->a[d0][d1])*(Lattice->a[d0][d1]);}
	  R2 = sqrt(R2);
	  n_UC[d0] = int(ceil(cutoff/R2)+1.1);
	  if(Lattice->BC[d0]==FREE) {n_UC[d0] = 0;}
	}
	//Compute Position of b0
	for(d0=0;d0<DIM;d0++) {X0[d0] = Unit_Cell->X[b0][d0];}
	//
	cutoff *= cutoff;           //Need squared value
	n_neig = 0;
  start = new N_LIST;
  end = start;


	//COMPUTE NEIGHBOR LIST, LIST POINTER, AND DATA LIST
  for(L1[0]=n_UC[0];L1[0]>=-n_UC[0];L1[0]--) {//Loop over unit cells
    for(L1[1]=n_UC[1];L1[1]>=-n_UC[1];L1[1]--) {
      for(L1[2]=n_UC[2];L1[2]>=-n_UC[2];L1[2]--) {
        //Compute Vector to Lattice Position
        for(d0=0;d0<DIM;d0++) {
          r_L1[d0] = 0.0;
          for(d1=0;d1<DIM;d1++) {r_L1[d0] += L1[d1]*Lattice->a[d1][d0];}
        }

        for(b1=0;b1<natom;b1++) {//Loop over atoms in unit cell
          //Check for Duplicate, Material, and Distance
          if((L1[0]==0)&&(L1[1]==0)&&(L1[2]==0)&&(b1==b0)) {continue;}
          if(Unit_Cell->mat[b1]!=Mat) {continue;}
          R2 = 0.0;
          for(d0=0;d0<DIM;d0++) {
            X1[d0] = r_L1[d0]+Unit_Cell->X[b1][d0];
            r_01[d0] = X1[d0] - X0[d0];	//Reversed for direct use in Dynamical Matrix
            R2 += r_01[d0] * r_01[d0];
          }

          if(cutoff<R2) {continue;}

          //Restart If n_UC Is Too Small
          for(d0=0;d0<DIM;d0++) {
            if( (abs(L1[d0])==-n_UC[d0])&&(Lattice->BC[d0]!=FREE) ) {
              Log  <<d0<<" INCREASING LATTICE UC_DOF"<<endl;
              cout <<d0<<" INCREASING LATTICE UC_DOF"<<endl;
              n_UC[d0] += 1;
              L1[0] = n_UC[0] + 1;
              L1[1] = n_UC[1] + 1;
              L1[2] = n_UC[2] + 1;
              b1 = natom;
              n_neig = 0;
              delete start;
              start = new N_LIST;
              end = start;
              continue;
            }
          }

          //Store Info in Linked List
          n_neig += 1;
          end->b  = b1;
          for(d0=0;d0<DIM;d0++) {
            end->l[d0] = L1[d0];
            end->x[d0] = X1[d0];
          }
          end->next = new N_LIST;
          end = end->next;

        }  //for(b1...)
      }  //for(l[2]...)
    } //for(l[1]...)
  } //for(l[0]...)


	//DEALLOCATE MEMORY
	delete[] L1;    L1   = NULL;
	delete[] n_UC;  n_UC = NULL;
	delete[] r_L1;  r_L1 = NULL;
	delete[] r_01;  r_01 = NULL;
	delete[] X0;  X0 = NULL;
	delete[] X1;  X1 = NULL;


	//ASSIGN NUMBER OF NEIGHBORS TO b0 AND RETURN LINKED LIST OF NEIGHBORS
	b0 = n_neig;
	return(start);

}


/*
SUBROUTINE ManageNeighborList: SHIFTS UNIT CELL INDICIES BY THOSE
  OF THE CENTER ATOM AND REMOVES THE ASSOCIATED ATOM (IF ANY) GIVEN
  IN A PREVIOUS LIST.  IF NO ATOMS NEED REMOVED ENTER THE CENTER ATOM
  TWICE.
*/
void ManageNeighborList(int &n_neig, N_LIST *n_list, int b_center,
                        int *l_center, int b_remove, int *l_remove) {
  //DECLARE LOCAL VARIABLES
  int i, j, n, l_total;       //Counters
  N_LIST *n_list_ptr;         //Pointer to current entry in n_list


  //CORRECT HIGHER ORDER NEIGHBOR LISTS AND REMOVE DUPLICATES
  n_list_ptr = n_list;
  for(n=0;n<n_neig;n++) {
    l_total = 0;
    //Correct unit cells
    for(i=0;i<DIM;i++) {
      n_list_ptr->l[i] += l_center[i];
      for(j=0;j<DIM;j++) {
        n_list_ptr->x[i] += double(l_center[j])*Lattice->a[j][i];
      }
      if(n_list_ptr->l[i]==l_remove[i]) {l_total+=1;}
    }
    //Check for duplicate (from previous list) and remove if needed
    if( (n_list_ptr->b==b_remove)&&(l_total==DIM) ) {
      if(n_list_ptr->next!=NULL) {
        n_list_ptr->b = n_list_ptr->next->b;
        for(i=0;i<DIM;i++) {
          n_list_ptr->l[i] = n_list_ptr->next->l[i];
          n_list_ptr->x[i] = n_list_ptr->next->x[i];
        }
        N_LIST *temp = n_list_ptr->next->next;
        n_list_ptr->next->next = NULL;
        delete n_list_ptr->next;
        n_list_ptr->next = temp;
      }
      else {delete n_list_ptr;  n_list_ptr=NULL;}
      n_neig -= 1;
      n -= 1;
      continue;
    }
    //Increment neighbor list
    n_list_ptr = n_list_ptr->next;
  }

	return;
}
