#include <stdio.h>
#include <stdlib.h>


#ifndef ___STRUCTURE_H
#define ___STRUCTURE_H

#define MAGICNUMBER 256 //used as a basis to encode informations about the free labels in a map into a single integer
#define BACKBONETYPE 1 // 0 for generating tree, 1 for generating path, 2 for generating cycle
#define FOLDSTEP 1 // 0 for generating backbone only, 1 to get the folding also
#define ALMOSTFOLDABLE 1 //0 for generating all paths, 1 for generating the almost foldable only
#define NBMAXNODESINAVL 1000000 // Maximum number of nodes stored in the AVL
#define NBFILEMAX 2000 // Maximum number of files written in the Results directory
#define SELECTIONCUTINDEX 0// Use the cutindex as the first criteria to compare signatures, option to remove
#define METAMOTIF 0 //do not use metamotifs with trees, it is not useful and not implemented
//a define to control the output ?
//most of the int should be char to make it more memory efficient (cache efficient ?)


typedef struct pair {
	int first;
	int second;
} pair;

typedef struct Label {
	int *list;// array of the value of the labels in the .mot (only the positive ones)
	int size; //number of different labels
} Label; //Label and compatibleVertices could be unified



typedef struct Vertex {  //An array of all possible vertices is built at the beginning of the algorithm 
	int letter; // Integer that gives the rank of the letter in the alphabet 9 for I, 10 for J, 22 for V, ...
	int degree; // Number of edges
	int id; // Position in the .mot file 
	int *edges; // Array containing the label of the edges
	int labelvalue;
} Vertex;


typedef struct CompatibleVertices {
	int *list; //array of the index of vertices which can be connected to some letter by their first edge 
	int vertexnumber; //size of the list
} CompatibleVertices;


// The outline structure //
typedef struct FreeEdge {  //a free edge in the graph
	int label; //number which labels the edge
	int vertexindex;  //index of the corresponding vertex in the vertexarray of the  Map
	int edgeindex; //number of the edge in the edges array  
} FreeEdge; //to merge with Edge

// The Map Structure//
typedef struct Edge {
	int vertexindex;//vertex connected to the edge
	int edgeindex;//position of this edge in the neighborood of the end vertex
} Edge; 

typedef struct VertexMap {
  int degree;//degree of the vertex
  int type;//position in the list of vertices of type Vertex -> should be id
  Edge *edges;//array of the neighbour given by their index in VertexMap 
} VertexMap;//we should store the id instead of the type since it is what we really want

struct Faces {
  int nb;       // The number of faces
  int max;      // The size of the largest face size
  float mean;   // The mean of the faces size
  float index;  // (max - second greatest size)/mean 
  unsigned long long int *distrib; // distrib[i] is the number of faces of size i (distrib[0] is always 0)
};

typedef struct Map{
  int vertexnumber; //number of vertices in the Map
  int edgenumber; // number of edges in the Map, it is also the size of the signature
  VertexMap *vertexarray;//Array of the vertices with their neighbours
  unsigned int *signature;
  int average_complexity;
  float cutindex; // Indice de coupe;
  int *vertexautomorphic; // Gives for each vertex the smallest indices among all the automorphic vertices (gives itself if no automorphism)	
  struct Faces faces; // La distribution des faces
} Map;	

//The dual Map Structure used to compute the cut index

typedef struct DualMap{
  int vertexnumber;
  struct DualVertex *vertexarray;
} DualMap;

struct DualVertex{
  int degree;
  struct DualEdge *edges;
};

struct DualEdge{
  int vertexindex;
  int edgeindex;
  int region;
  int weigth;
};

// The Avl structure
struct abr {  // The abr in the avl to store the signatures
	unsigned int *signature;
	struct abr *l, *r;
};

struct AvlNode { // One node of the avl
	int h; // height of the tree
	Map *M; // the data of the node
	int nb_iso; // number of isomorphism of the same map obtained during construction 
	struct AvlNode *l, *r; // left and right trees
};

typedef struct Avl{
	unsigned long int nbnodes; // number of elements in the avl   ##
	unsigned long long int nbcalls; // number of call to the insert function ##
	float lower_cutindex; // lower cutindex
	float limit_cutindex; // limit cutindex stored in the avl
	float higher_cutindex; // higher cutindex
	struct AvlNode *root; // pointer to the root of the avl
	struct abr *S;
} Avl;


//All variables defined here globally correspond to unique instance of their type. They are used to store immutable informations all along the algorithm.



#endif
