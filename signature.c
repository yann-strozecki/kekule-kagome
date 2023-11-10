#include <string.h>
#include "structure.h"
#include "signature.h"
#include "global.h"

// ######
// 1. Print a signature
// ######



void print_vertexsignature(unsigned int sig){
  int id,label,index;
  label = sig & 0xff; sig = sig>>8;
  id = sig & 0xff; sig = sig>>8;
  index = sig & 0xff;
  printf("id of the vertex %d, label of the edge %d, index of the next vertex %d \n",id,label,index);
}


void print_signature(Map *M) {
int i;
 printf("%d edges\n Début de la signature \n",M->edgenumber);
 for (i=0 ; i<M->edgenumber ; i++) {
   print_vertexsignature(M->signature[i]);
 }
 printf("Fin de la signature \n \n");
}
 
// ######
// 2. Computes the signature starting from the vertex and the edge given in args.
// ######


void signature_vertex_edge(Map *M, int vertex_deb, int edge_deb) {
  // Rank is an array that gives for each vertex its rank in the signature, -1 if the vertex has no rank yet
  if( M->vertexarray[vertex_deb].edges[edge_deb].vertexindex == vertex_deb) return; //do not compute a signature beginning by a an edge going to itself
  memcpy(rank,rank_init,sizeof(int)*mapsize);//initialize the rank array
  rank[vertex_deb] = 0;
  int rankindex = 1;
  int stack_length = 1;
  int lower_signature = 0;
  int S_index,nextvertex, val;
  edge_stack[0] = edge_deb;
  int currentvertex = vertex_deb, currentedge = edge_deb;
  for(S_index = 0; S_index < M->edgenumber;S_index++){ //stop when the DFS have explored all edges
    nextvertex = M->vertexarray[currentvertex].edges[currentedge].vertexindex;
    val = vertices[M->vertexarray[currentvertex].type].edges[currentedge]
      +  (vertices[M->vertexarray[currentvertex].type].id<<8); //this is a redondant information but can be useful to cut the computation earlier
    //val contains the information on the label of the current edge, and the 
    //rank of nextvertex (added later). 
    if(rank[nextvertex]==-1){//index a new found vertex 
      rank[nextvertex] = rankindex; 
      rankindex++;
      //           Move to the next vertex              //
      currentedge = M->vertexarray[currentvertex].edges[currentedge].edgeindex;
      currentvertex = nextvertex;
      edge_stack[stack_length] = currentedge;
      stack_length++;
    }
    else{
      if(currentedge == edge_stack[stack_length-1]){ //go back in the DFS
	currentedge = M->vertexarray[currentvertex].edges[currentedge].edgeindex;
	currentvertex = nextvertex;
	stack_length--;
      }
    }
    val += rank[nextvertex]<<16;//update the val using the rank of the next vertex
    //M->average_complexity++;
    if (!lower_signature && (val > M->signature[S_index])) {return;} // Signature is higher ==> stop
    lower_signature = lower_signature | (val < M->signature[S_index]); // Signature is lower ==> new better signature found, cannot be stopped
    M->signature[S_index] = val;
    currentedge++;//move to the next edge
    if(currentedge == M->vertexarray[currentvertex].degree) currentedge = 0;//modulo degree without %
  }
}


// ######
// 3. Computes the minimal signature of a map
// ######
void signature_compute(Map *M) {
  int vertex, edge;
  // Computes the number of edges (in fact 2 times the number of edges)
  M->edgenumber=0; 
  for (vertex = 0 ; vertex<M->vertexnumber; vertex++) M->edgenumber += M->vertexarray[vertex].degree;//could be computed along the creation of the backbone
  // Intialize the signature
  M->signature = calloc(M->edgenumber,sizeof(unsigned int));
  M->signature[0] = 0x7fffffff; 	// The left most bit is always 0 in a signature. 
  // This initialization garanties that the first computed signature will be lower
  // Computes the signature for each starting position (vertex,edge)
  //M->average_complexity =0;
  for (vertex=0 ; vertex < M->vertexnumber ; vertex++){
    for (edge=0 ; edge < M->vertexarray[vertex].degree ; edge++){
      if(vertices[M->vertexarray[vertex].type].id == 0) signature_vertex_edge(M,vertex,edge);
    }
  }
}



// ######
// 4. Signature comparison
// ######
// if signature(M1) < signature(M2) returns a negative value
// if signature(M1) = signature(M2) returns 0
// if signature(M1) < signature(M2) returns a positive value
int signature_compare(Map *M1, Map *M2) {
    if (M1->edgenumber != M2->edgenumber) {fprintf(stderr,"Strange: signature_compare of 2 maps not of the same size\n"); exit(3);}//it is possible to have different number of edges for the same number of vertices, but it does not happen in our examples
  if(SELECTIONCUTINDEX){
    float delta = M1->cutindex - M2->cutindex; //we may compile that only if the cut index is computed for every map
    if (delta < -1e-6) return -1; 
    if (1e-6 < delta)  return +1; 
  }
  int i;
  for (i=0 ; i<M1->edgenumber ; i++) {
    if (M1->signature[i] < M2->signature[i]) return -1;
    if (M1->signature[i] > M2->signature[i]) return +1;
  }
  return 0;
}
