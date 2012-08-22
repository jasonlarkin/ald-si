#ifndef LIST_H_INCLUDED
#define LIST_H_INCLUDED


/*LIST CLASS USED AS LINKED LIST IN NEIGHBOR LIST*/
class LIST {
public:
  LIST(void) {
    l = new int[DIM];
    r = new double[DIM];
    next = NULL;
    return;
  }
  ~LIST(void) {
    delete[] l;  l = NULL;
    delete[] r;  r = NULL;
    return;
  }
	int n;
	int *l;
	double *r;
	LIST *next;
};


#endif // LIST_H_INCLUDED
