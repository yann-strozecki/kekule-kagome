#include "structure.h"


void print_signature(Map *M);

// Computes the minimal signature of M
void signature_compute(Map *M);

// Computes the signature starting at vertex_deb, edge_deb
void signature_vertex_edge(Map *M, int vertex_deb, int edge_deb);


// if signature(M1) < signature(M2) returns a negative value
// if signature(M1) = signature(M2) returns 0
// if signature(M1) < signature(M2) returns a positive value
int signature_compare(Map *M1, Map *M2);
