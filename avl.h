#include "structure.h"

// creates a new Avl
Avl avl_create();

// add a new map in the Avl.
// if a Map with the same signature already exist in the Avl, nothing is added and M is freed
void avl_insert(Avl *A, Map *M);

// print an Avl
void avl_print(Avl A);

// Free the avl
void avl_free(struct AvlNode *A);

// Free the abr
void abr_free(struct abr *S);
