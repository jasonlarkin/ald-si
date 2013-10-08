/*                         Neighbor_Force.h                         */
/*                            06/06/2009                            */
/*********************************************************************
*    Header file for the Neighbor_Force template.                    *
*********************************************************************/

#ifndef NEIGHBOR_FORCE_H_INCLUDED
#define NEIGHBOR_FORCE_H_INCLUDED

/*DEFINE HEADERS*/
#include "LDCode.h"


template <class T>
/*Neighbor_Force: COMPUTES NEIGHBOR LISTS AND FORCE CONSTANTS*/
T **Neighbor_Force(T **F_C, POTENTIAL **Potential) {
  //DECLARE LOCAL VARIABLES
  int b0;                     //Counter
  int n_mat = Material->n_mat;//Number of materials
  int M;                      //Material of atom b0
  POTENTIAL *P_ptr = NULL;    //Potential pointer


  //ALLOCATE MEMORY
  F_C = new T*[Unit_Cell->natom];
  for(b0=0;b0<Unit_Cell->natom;b0++) {F_C[b0] = new T(b0);}


  //LOOP OVER ALL ATOMS IN UNIT CELL COMPUTING FORCE CONSTANTS
  for(b0=PE;b0<Unit_Cell->natom;b0+=nPE) {
    M = Unit_Cell->mat[b0];
    P_ptr = Potential[M];
    while(P_ptr!=NULL) {
      P_ptr->ForceConstants(b0, F_C[b0]);
      P_ptr = P_ptr->next;
    }
  }


  //SCALE FORCE CONSTANTS AND UPDATE
  for(b0=0;b0<Unit_Cell->natom;b0++) {
    if((b0%nPE)==PE) {F_C[b0]->Send();}
    else {F_C[b0]->Recv(b0%nPE);}
  }

  return(F_C);
}

#endif // NEIGHBOR_FORCE_H_INCLUDED
