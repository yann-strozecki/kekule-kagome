#include "structure.h"
#include "signature.h"
#include "index.h"
#include "avl.h"
#include "chrono.h"
// ######
// 1. Print of Avl (for debug)
// ######

void avl_node_print(struct AvlNode *A, int prof) {
int i;
if (A) {
	avl_node_print(A->l,prof+1);
	for (i=0 ; i<prof ; i++) printf("   ");  
	printf("%x.%x.%x.%x.%x\n",A->M->signature[4]%0x10,A->M->signature[5]%0x10,A->M->signature[6]%0x1000,A->M->signature[7]%0x1000,A->M->signature[8]%0x1000);
	avl_node_print(A->r,prof+1);
	}
}

void avl_print(Avl A) {
printf("\n### AVL ###\n");
printf("%ld nodes\n",A.nbnodes);
avl_node_print(A.root,0);
printf("###########\n");
}

// ######
// 2. Avl
// ######

// ## 2.0 Free a Map


// ## 2.1 Macros used to deal with height and equilibrium

#define sup(x,y)( (x)>(y)?(x):(y) )
#define lh(A) ( ((A)->l)?(A)->l->h:0 )   // Height of the left subtree
#define rh(A) ( ((A)->r)?(A)->r->h:0 )   // Height of the right subtree
#define height(A) ( 1+sup(lh(A),rh(A)) ) // Height of a node
#define eq(A) ( (A)?(lh(A)-rh(A)):0 )    // Equilibrium of a node

// ## 2.2 Creates a new Avl
Avl avl_create() {
Avl A;
A.nbnodes = 0;
A.nbcalls = 0;
A.lower_cutindex  = -1.0; // None yet
A.limit_cutindex  = -1.0; // None yet
A.higher_cutindex = -1.0; // None yet
A.root = NULL; 
A.S = NULL;
return A;
}

// ## 2.2 Insertion of a new node

// ## 2.2.1 Balance the tree after insertion (see section 2.2.4 for the body)
struct AvlNode* balance(struct AvlNode *A);

// ## 2.2.2 Duplicate a map before insertion in the avl
Map *duplicate(Map *M) {
Map *newM;
newM = malloc(sizeof(Map)); if (!newM) {fprintf(stderr,"cannot malloc new map in avlc.c\n"); exit(61); }
*newM = *M; // Copy all the fields
// debug qst 
// static int MMM=0;
// printf("%3d newM->signature = %p %p\n",MMM++,newM->signature, M->signature);
M->signature = NULL; // The signature is usefull only in the avl
newM->vertexarray = malloc(newM->vertexnumber*sizeof(struct VertexMap));
if (!newM->vertexarray) {fprintf(stderr,"cannot malloc new vertexarray in avlc.c\n"); exit(62); }
int v;
for (v=0 ; v<newM->vertexnumber ; v++) {
	newM->vertexarray[v] = M->vertexarray[v];
 	newM->vertexarray[v].edges = malloc(newM->vertexarray[v].degree*sizeof(struct Edge));
 	if (!newM->vertexarray[v].edges) {fprintf(stderr,"cannot malloc new edges in avl.c\n"); exit(63); }
 	int e;
	for (e=0 ; e<newM->vertexarray[v].degree ; e++) newM->vertexarray[v].edges[e] = M->vertexarray[v].edges[e];
	}
newM->vertexautomorphic = NULL;
newM->faces.distrib = NULL;
return newM;
}

// ## 2.2.3 Insertion of the Map M (if not present) in the AvlNode A
//          the argument inserted will be 1 if M is inserted and 0 is M is found to already exist and is not inserted
struct AvlNode* insert(struct AvlNode *A, Map *M, int *inserted) {
int delta;
if (!A) {	// new node to insert
	A = malloc(sizeof(struct AvlNode)); if (!A) { fprintf(stderr,"cannot malloc AvlNode\n"); exit(64); }
	A->h = 1; A->M = duplicate(M); A->nb_iso = 1; A->l = A->r = NULL;
	*inserted = 1;
	return A;
	}
delta = signature_compare(M,A->M);
if (delta == 0) { (A->nb_iso)++; (A->M)->average_complexity += M->average_complexity; *inserted = 0; } // already present nothing to insert
if (delta <  0) { A->l = insert(A->l,M,inserted); } // try to insert in left subtree
if (delta >  0) { A->r = insert(A->r,M,inserted); } // try to insert in right subtree
if (*inserted)  { A = balance(A); } // balance the tree after insertion
return A;
}

// ## 2.2.4 Free a node
void avl_node_free(struct AvlNode *A) {
int v;
if (A->M->vertexarray) {
	for (v=0 ; v<A->M->vertexnumber ; v++) free(A->M->vertexarray[v].edges);
	free(A->M->vertexarray);
	}
// qst debug
// fprintf(stderr,"free signature in avl %p\n",A->M->signature); 
if (A->M->signature) free(A->M->signature);
if (A->M->vertexautomorphic) free(A->M->vertexautomorphic);
if (A->M->faces.distrib) free(A->M->faces.distrib);
free(A->M);
free(A);
}

void avl_free(struct AvlNode *A) {
if (A) {
	avl_free (A->l);
	avl_free (A->r);
	avl_node_free(A);
	}
}

void abr_free(struct abr *S) {
if (S) {
	abr_free (S->l);
	abr_free (S->r);
	if (S->signature) free(S->signature);
	free(S);
	}
}

// ## 2.2.5 Remove the limit map (the one to the left) in the AVL A
struct AvlNode *avl_remove(struct AvlNode *A, float *wci) {
if (A->l==NULL) { // limit map to remove found
	struct AvlNode *B = A->r;
	while (B) { *wci = B->M->cutindex; B = B->l; } // if A->r is NULL *wci has the value computed in the father
	B = A->r;
	avl_node_free(A);
	return B;
	}
else {
	*wci = A->M->cutindex;
	A->l = avl_remove(A->l,wci);
	A->h = height(A);
	A = balance(A);
	return A;
	}
}


struct abr *is_in_abr(struct abr *S, Map *M, int *inserted) {
if (S==NULL) { 
	S = malloc(sizeof(struct abr)); if (!S) exit(87);
	S->l = S->r = NULL;
	S->signature = M->signature;
	*inserted=1;
	return S;
	}
int delta = 0;
int i = 0;
while ((i<M->edgenumber) && (delta==0)) {
	if (S->signature[i] < M->signature[i]) delta=-1;
	if (S->signature[i] > M->signature[i]) delta=+1;
	i++;
	}
if (delta==0) *inserted=0;
if (delta<0) S->l = is_in_abr(S->l,M,inserted);
if (delta>0) S->r = is_in_abr(S->r,M,inserted);
return S;
}



// ## 2.2.6 Insertion of the Map M in the Avl A
void avl_insert(Avl *A, Map *M) {
(A->nbcalls)++;
int inserted = 0;
M->cutindex  = 0;
M->signature = NULL;
//Uncomment to test cut index speed
/*	double t_begin,t_end;
	t_begin = what_time_is_it();
	M->cutindex  = sparsity_poly(M);
	t_end = what_time_is_it();
	printf("time : %f \n \n \n \n",t_end-t_begin);
	exit(0);
*/
	// Limitation of the number of nodes in the AVL
if (SELECTIONCUTINDEX) {
	signature_compute(M);
	A->S = is_in_abr(A->S,M,&inserted);
	if (!inserted) {free(M->signature); M->signature=NULL; return;}
	}
if (SELECTIONCUTINDEX) M->cutindex = sparsity_poly(M);
	// For the first one
if (A->lower_cutindex == -1.0)  A->lower_cutindex = A->limit_cutindex = A->higher_cutindex = M->cutindex;
if (M->cutindex < A->lower_cutindex)  A->lower_cutindex  = M->cutindex;
if (M->cutindex > A->higher_cutindex) A->higher_cutindex = M->cutindex;
if ((A->nbnodes <= NBMAXNODESINAVL) || (M->cutindex > A->limit_cutindex) || !SELECTIONCUTINDEX) {
	if (!SELECTIONCUTINDEX) signature_compute(M);
	inserted = 0;
	A->root = insert(A->root,M,&inserted); // inserted returns 1 if M is inserted and 0 if M is already present
	A->nbnodes += inserted;
	if (inserted) {
		if (A->limit_cutindex > M->cutindex) A->limit_cutindex = M->cutindex;
		if ((A->nbnodes >= NBMAXNODESINAVL+1) && SELECTIONCUTINDEX) { A->root = avl_remove(A->root,&(A->limit_cutindex)); (A->nbnodes)--; }
		}
	else {free(M->signature); M->signature=NULL;}
    }
}

// ## 2.2.7 Functions to balance the tree
struct AvlNode* rotate_right(struct AvlNode *b) {
/*
    b             a
   / \           / \
  a   Z   ==>   X   b
 / \               / \
X   Y             Y   Z
*/
struct AvlNode *a=b->l;
if (!a) return b; // nothing to rotate
b->l = a->r;
a->r = b;
b->h = height(b); // Height of b may have changed
a->h = height(a);
return a;
}

struct AvlNode* rotate_left(struct AvlNode *b) {
/*
    b             c
   / \           / \
  X   c   ==>   b   Z
     / \       / \
    Y   Z     X   Y
*/
struct AvlNode *c=b->r;
if (!c) return b; // nothing to rotate
b->r = c->l;
c->l = b;
b->h = height(b); // Height of b may have changed
c->h = height(c);
return c;
}

struct AvlNode* balance(struct AvlNode *A) {
/*
         a (+2)
        / \     => right rotation on a
(+1,0) b   c

         a (+2)
        / \     => left rotation on b AND right rotation on a
  (-1) b   c
Symetrically if equilibrium of a is -2
*/
int eqA;
eqA = eq(A);
if (eqA == +2) {
	if (eq(A->l)==-1) A->l = rotate_left(A->l);
	A = rotate_right(A);
	}
if (eqA == -2) {
	if (eq(A->r)==+1) A->r = rotate_right(A->r);
	A = rotate_left(A);
	}
return A;
}

