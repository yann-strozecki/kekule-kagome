#include "structure.h"
#include <math.h>
#include "fold.h"
#include "global.h"

//#################### Print functions ###########################

void printcompatible(int size){
  int i,j;
  for(i=0;i<size;i++)
    {
      printf("Label %d:",i);
      for(j=0; j< connection[i].vertexnumber;j++)
	{
	  printf("Vertex %d ",connection[i].list[j]);
	}
      printf("\n");
    }
}

void printvertices(int size, Vertex *vert){
int i,j;
printf("%d vertices \n", size);
for(i=0; i<size;i++) {
	printf("Degree: %d, Id:  %d, Value: %d Edges: ",vert[i].degree,vert[i].id,vert[i].labelvalue);
	for(j=0; j < vert[i].degree; j++) {
		printf("%d ",vert[i].edges[j]);
		}
	printf("\n");
	}
}

void printbinary(int code, int size)
{
  int i;
  for(i=0;i<size;i++)
    {
      printf("%d", code%2);
      code /= 2;
  }
}

void printindex(int code, int codesize, int base){
  int i,temp;
  printf("Index %d (",code);
  for(i=0;i<codesize;i++)
    {
      temp = code%base - base/2;
      printf("%d ", temp);
      code /= base;
    }
  printf(")");
}

void printmatrix(unsigned int ***mat, int dim1, int dim2,int dim3,int magicnumber)
{
  int i,j,k;
  for(i=dim1-1; i>=0;i--)
    {
      printf("\n Size %d\n",dim1 - i);
      for(j=0;j<dim2;j++){
	printf("\n Complementary label %d\n",j);
	for(k=0;k<dim3;k++)
	  {
	    if(mat[i][j][k]){
	      printindex(k,dim2/2,magicnumber);
	      printf("Set ");
	      printbinary(mat[i][j][k],dim2);
	      printf("; ");
	    }
	  }
      }
      printf("\n"); printf("\n");
    }
}



void printmatrixtree(unsigned int **mat, int dim1, int dim2,int magicnumber,int codesize){
  int i,j;
  for(i=0; i<dim1;i++)
    {
      printf("\n Complementary size %d\n",i);
      for(j=0;j<dim2;j++){
	if( mat[i][j])	printindex(j,codesize,magicnumber);
      }
      printf("\n");
    }
}

/************************************** Input functions **********************************/

//Create all vertices given in the .mot, return the number of vertices
int read_input(char* inputname, int *maxdegree){
  FILE* inputfile = NULL;
  int verticesnumber, degree, i, j;
  int ret;
  int period; //is not used, should be removed from the file format
  inputfile = fopen(inputname, "r");
  ret = fscanf(inputfile, "%d %d", &verticesnumber, &period);
  vertices = malloc(sizeof(Vertex)*verticesnumber);
  *maxdegree = 0;
  for(i=0; i < verticesnumber;i++) {
    ret = fscanf(inputfile, "%d", &(vertices[i].letter));
    ret = fscanf(inputfile, "%d", &degree);
    if (degree > *maxdegree) *maxdegree = degree;
    vertices[i].degree = degree;
    vertices[i].id = i;
    vertices[i].edges = malloc(sizeof(int)*degree);
    for(j=0; j < degree; j++) {
      ret = fscanf(inputfile,"%d", &(vertices[i].edges[j]));
    } 
  }
  fclose(inputfile);
  if (ret<0) exit(11);
  return verticesnumber;
} 

int read_concat(char* inputname){
  int i,j,k;
  int motifnumber;
  FILE *inputfile = fopen(inputname, "r");
  int ret = fscanf(inputfile, "%d", &motifnumber);
  translation = malloc(motifnumber*sizeof(Map));
  for(i = 0; i< motifnumber; i++){
    ret = fscanf(inputfile, "%d", &(translation[i].vertexnumber)); 
    translation[i].vertexarray = malloc(sizeof(VertexMap)*translation[i].vertexnumber);
    for(j = 0; j < translation[i].vertexnumber;j++){
      ret = fscanf(inputfile, "%d %d", &(translation[i].vertexarray[j].type), &(translation[i].vertexarray[j].degree) );
      translation[i].vertexarray[j].edges = malloc(sizeof(Edge)*translation[i].vertexarray[j].degree);
      for(k=0;k < translation[i].vertexarray[j].degree;k++){
	ret = fscanf(inputfile, "%d %d", &(translation[i].vertexarray[j].edges[k].vertexindex), &(translation[i].vertexarray[j].edges[k].edgeindex) );
      }
    }
  }
  fclose(inputfile);
  if (ret<0) exit(12);
  return motifnumber;
}

/************************************************ Construction of the datastructures ****************************/

Vertex* create_metavertices(int motifnumber){//
  int i,j,k,outlinesize;
  outlinearray = malloc(motifnumber*sizeof(FreeEdge*));
  Vertex *alternatevertices = malloc(motifnumber*sizeof(Vertex));
  //normalize the label of the motifs in alternatevertices
  for(i=0;i<motifnumber;i++){ //for each metamotif compute the outline and create a metavertex
    //store a map beteen the edges of the meta vertex and of the map which represents the elements of the matavertex
    alternatevertices[i].labelvalue = 0;
    for(j=0; j < translation[i].vertexnumber;j++){alternatevertices[i].labelvalue += vertices[translation[i].vertexarray[j].type].labelvalue; }
    outlinearray[i] = malloc(sizeof(FreeEdge)*translation[i].vertexnumber*4); //here the degree is bounded by 6, could be changed in other applications 
    outlinesize = computeoutline(&(translation[i]), outlinearray[i],vertices);  
    alternatevertices[i].id = i;
    alternatevertices[i].degree = outlinesize;
    alternatevertices[i].edges = malloc(sizeof(int)*outlinesize);
    for(j=0;j<outlinesize;j++){
      k=0;
      while(label.list[k] != outlinearray[i][j].label) k++;
      outlinearray[i][j].label = k;
      alternatevertices[i].edges[j] = k;
    }
  }
  return alternatevertices;
}



/////////////Fonctions to create datastructures with important precomputed informations to speed up path enumeration////

//Index the labels by their order of apparition in vertices and store that in label. 
//The negative labels are stored after the positive ones in the same order
//Return the number of different labels.
//It also compute the fields labelvalue
void normalize_labels(int size, Vertex *vert ){ 
  int i,j,k,labelnumber=0,test;
  label.list = malloc(sizeof(int)*10); //we have at most 10 different labels (can be changed)  
  for(i=0; i < size; i++) {//create label.list
    for(j=0; j < vert[i].degree; j++) {
      test = 1;
      for(k=0; k< labelnumber;k++) {
	if(label.list[k] == abs(vert[i].edges[j])) {
	  test = 0;
	  break;
	}
      }
      if(test) {
	label.list[labelnumber] = abs(vert[i].edges[j]);
	labelnumber++;
      }
    }
  }
  for(i=0;i <labelnumber;i++){label.list[labelnumber+i] = -label.list[i];}
  // change the labels to their indices in label.list and compute labelvalue
  label.size = 2*labelnumber;
  for(i=0; i < size; i++) {
    vert[i].labelvalue=0;
    for(j=0; j < vert[i].degree; j++) {
      k=0;
      while(label.list[k] != vert[i].edges[j]) k++;
      vert[i].edges[j] = k;
      if(k < labelnumber) {
	vert[i].labelvalue += (int)pow(MAGICNUMBER,k);//MAGIC number is large enough to encode without collisions the number of each label
      }
      else {
	vert[i].labelvalue += -(int) pow(MAGICNUMBER,k - labelnumber);
      }
    }
  }
}



int non_isomorph(int *list, int size, int shift) { //return false if the list of edges shifted by shift is isomorph to the list shifted by less
int i,j, test;
for(i=0; i < shift; i++) {
	test = 1;
	for(j=0; j<size; j++) {
		if ( list[(j+i) % size] != list[(i + shift) % size]) {
			test = 0;
			break;
			}
		}
	if(test) return 0;
	}
return 1;
}

 //Add to the list of vertices all vertices with a rotation on their edges (and do not add isomorphic copies of a vertex), return the number of vertices
////when the vertices are meta vertices add also the relevant maps and outlines to translation and outlinearray 
Vertex *create_rotated_vertices(int size, int *newsize, Vertex *vert){
  int i,j,k,maxverticesnumber=0,verticesnumber=size;
  for(i=0;i<size;i++)
    {
      maxverticesnumber+= vert[i].degree;
    }
  vert = realloc(vert, maxverticesnumber*sizeof(Vertex));
  if(METAMOTIF){
    outlinearray = realloc(outlinearray, maxverticesnumber*sizeof(FreeEdge*));
    translation =  realloc(translation, maxverticesnumber*sizeof(Map));
  }
  for(i=0;i<size;i++)//for each vertex add a rotated version if it is not already there 
    {
      for(j=1; j < vert[i].degree;j++)
	{
	  if(non_isomorph(vert[i].edges, vert[i].degree,j))
	    {
	      vert[verticesnumber] =  vert[i];
	      vert[verticesnumber].edges = malloc(sizeof(int)*vert[i].degree);
	      if(METAMOTIF){
		outlinearray[verticesnumber] = malloc(sizeof(FreeEdge)*vert[i].degree);
		translation[verticesnumber] = translation[i];
	      }
		for(k=0; k < vert[i].degree; k++)
		{
		  vert[verticesnumber].edges[k] =  vert[i].edges[(j+k) % vert[i].degree];
		  if(METAMOTIF) outlinearray[verticesnumber][k] =  outlinearray[i][(j+k) % vert[i].degree];
		}
	      verticesnumber++;
	    }
	}
    }
  *newsize = verticesnumber;
  return vert;
}

//Create an helper array with for each label the vertices which can be connected by their first edge to this label
void create_concatenation_helper(int size, Vertex * vert)
{ 
  int i,j;
  int opposite = label.size/2;
  connection = malloc(sizeof(CompatibleVertices)*label.size);
  for(i=0;i<label.size;i++)
    {
      connection[i].list = malloc(sizeof(int)*size);
      connection[i].vertexnumber = 0;
      for(j=0;j<size;j++)
	{
	  if( abs(vert[j].edges[0] - i) ==  opposite) {
	    connection[i].list[connection[i].vertexnumber] = j;
	    connection[i].vertexnumber++;
	  }
	}
    }
}


//Create a multidimensional table which contains the information about the almost foldable trees


void almost_foldable_tree(int vertexnumber, Vertex *vert)
{
  int i,j,k;
  int alphabetsize = label.size/2;
  int maxsize = (int) pow(MAGICNUMBER, alphabetsize);
  almostfoldabletree = malloc(sizeof(unsigned int*)*mapsize);
  for(i=0;i<mapsize;i++)//initialise mat to 0
    {
      almostfoldabletree[i] = calloc(maxsize,sizeof(unsigned int));
    }
  for(i=0;i<vertexnumber;i++) almostfoldabletree[mapsize-1][vert[i].labelvalue + shift] = 1;
  for(i=1;i<mapsize;i++) {
    for(j=0;j<maxsize;j++) 
      {
	if(almostfoldabletree[mapsize -i][j]){
	  for(k=0;k<vertexnumber;k++) almostfoldabletree[mapsize-i-1][vert[k].labelvalue + j] = 1;
	}
      }
  }
  //printmatrixtree(mat,size,maxsize,MAGICNUMBER,alphabetsize);
}

//Create a multidimensional table which contains the information about the almost foldable paths

void almost_foldable_path(int vertexnumber, Vertex *vert)
{
  int i,j,k,index,value=0,l;// magic number is defined in structure.h, make it larger when labelnumber is large
  pair **add_after_a_label = malloc(sizeof(pair*)*label.size); //will be used to add all possible vertices after an edge of given label
  int **metavertexsize =  malloc(sizeof(int*)*label.size);
  for(i=0;i<label.size;i++) {
    add_after_a_label[i] = malloc(sizeof(pair)*vertexnumber);
    metavertexsize[i] = malloc(sizeof(int)*vertexnumber);
  }
  int *sizes = calloc(label.size,sizeof(int));// sizes[i] is the size of add_after_a_label[i]
  int alphabetsize = label.size/2;
  int maxsize = (int) pow(MAGICNUMBER, alphabetsize);
  almostfoldablepath = malloc(sizeof(unsigned int**)*mapsize);
  for(i=0;i<mapsize;i++)//initialise mat to 0
    {
      almostfoldablepath[i] = malloc(sizeof(unsigned int*)*label.size);
      for(j=0;j<label.size;j++)
		{
		almostfoldablepath[i][j] = calloc(maxsize,sizeof(unsigned int));
		}
    }
  //Compute add_after_a_label and initialize almostfoldablepath[0]
  for(i=0;i<vertexnumber;i++)
    {
      value = vert[i].edges[0] + (vert[i].edges[0] >= alphabetsize ? -alphabetsize : alphabetsize); //label of the first edge (we index by its complement to facilitate the use of almostfoldablepath in concat)
      for(j=1;j<vert[i].degree;j++)//j=0 is used for the beginning of the path
	{
	  if(METAMOTIF){almostfoldablepath[mapsize-translation[i].vertexnumber][value][vert[i].labelvalue + shift] |= 1<<vert[i].edges[j];}
	  else {almostfoldablepath[mapsize-1][value][vert[i].labelvalue + shift] |= 1<<vert[i].edges[j];}//add a free edge to the possible end of the path of size 1
	}
      add_after_a_label[value][sizes[value]].first = vert[i].labelvalue;//first is the labelvalue
      if(METAMOTIF) add_after_a_label[value][sizes[value]].second = almostfoldablepath[mapsize-translation[i].vertexnumber][value][vert[i].labelvalue + shift];
      else add_after_a_label[value][sizes[value]].second = almostfoldablepath[mapsize-1][value][vert[i].labelvalue + shift];//second is the set of labels of ending edges we have just computed
      if (METAMOTIF) metavertexsize[value][sizes[value]] = translation[i].vertexnumber;
      sizes[value]++;
      //printf("add labelvalue %d label %d size %d vertexnumber %d\n",vert[i].labelvalue,almostfoldablepath[mapsize-translation[i].vertexnumber][value][vert[i].labelvalue + shift],translation[i].vertexnumber,i);
    }
  //Loop to build almostfoldablepath from size i by adding a (meta)vertex
  for(i = mapsize-1;i>0;i--)
    {
      for(j=0;j<label.size;j++)//for all first label of the path
	{
	  for(index=0;index < maxsize;index++)//for all values of Cg (vector of number of free element in the alphabet)
	    {
	      if(almostfoldablepath[i][j][index])//if there is a compatible path
		{
		  for(k=0;k<label.size;k++)//for all possible last edge of the path add the possible end motifs
		    {
		      if(1<<k & almostfoldablepath[i][j][index])//if the label k is a valid end label
			for(l=0;l < sizes[k]; l++)
			  {
			    if(METAMOTIF && (i >= metavertexsize[k][l])) almostfoldablepath[i - metavertexsize[k][l]][j][index + add_after_a_label[k][l].first] |=  add_after_a_label[k][l].second;
			    else almostfoldablepath[i - 1][j][index + add_after_a_label[k][l].first] |=  add_after_a_label[k][l].second;
			  }
		    }
		}
	    }
	}
    }
  //printmatrix(almostfoldablepath,mapsize,label.size,maxsize,MAGICNUMBER);
  for(i=0;i<label.size;i++){free(add_after_a_label[i]); free(metavertexsize[i]);}
  free(add_after_a_label);
  free(metavertexsize);
  free(sizes);
}
