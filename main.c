#include <math.h>
#include <unistd.h>
#include "structure.h"
#include "input.h"
#include "concat.h"
#include "avl.h"
#include "signature.h"
#include "index.h"
#include "chrono.h"
#include "output.h"


Vertex *vertices;// Vertices given in the mot and their rotation
Map *translation; //store the metamotifs as maps
FreeEdge **outlinearray;//free edge of the metamotifs
Label label; //Label of the edges
CompatibleVertices *connection;//helper structure to do the concatenation fast
unsigned int *** almostfoldablepath;//structure to cut in the enumeration of paths
unsigned int ** almostfoldabletree;//structure to cut in the enumeration of trees
Vertex *concatvertices;  
int mapsize,maxdegree,maxoutlinesize,shift;
int *previous,*mat,*foldable_with,**fold_matrix,*call_stack;
pair *to_fold,*folded;
pair *face ,**edge, *tree,*queue;
DualMap *dual; 
int **edgebis;
int *rank, *rank_init, *edge_stack;





void prettyprint(unsigned long long int  x) {
  int t,g,m,k,u;
  t = x/1e12; x = x-t*1e12;
  g = x/1e9;  x = x-g*1e9;
  m = x/1e6;  x = x-m*1e6;
  k = x/1e3;  u = x-k*1e3;
  if (t) printf("%5d,%03d,%03d,%03d,%03d",t,g,m,k,u);
  else if (g) printf("      %3d,%03d,%03d,%03d",g,m,k,u);
  else if (m) printf("          %3d,%03d,%03d",m,k,u);
  else if (k) printf("              %3d,%03d",k,u);
  else        printf("                  %3d",u);
}


int compute_average(struct AvlNode *root)
{
  if (root == NULL) return 0;
  return  (((root->M)->average_complexity / (root->nb_iso)) + compute_average(root->l) +  compute_average(root->r));
}

int main (int argc, char *argv[])
{
  int i,j;
  printf("## START ##\n");
  chrono();
  if (argc < 3)  { fprintf(stderr,"usage: %s mot_file size\n",argv[0]); exit(1); }
  mapsize = atoi(argv[2]);

  printf("1. Read of the input file: %s, Size = %d. \n Creation of the data structures.\n",argv[1],mapsize);

  Avl avl = avl_create();
  int vertexnumber, vertexconcatnumber;
  vertexnumber = read_input (argv[1], &maxdegree);
  //printvertices(vertexnumber,vertices);
  normalize_labels(vertexnumber, vertices);
  //for(i=0;i<label.size;i++) printf("label %d : %d ",i,label.list[i]);
  //printvertices(vertexnumber,vertices);
  concatvertices = vertices;
  if (METAMOTIF){
    vertexnumber = read_concat(argv[3],translation);
    concatvertices = create_metavertices(vertexnumber); //we use the metavertices for the concat part
    //printvertices(vertexnumber,concatvertices);
  }
  concatvertices = create_rotated_vertices(vertexnumber,&vertexconcatnumber,concatvertices);
  //for(i=0;i<vertexconcatnumber;i++){printmap(translation[i]);}
  if(!METAMOTIF) vertices = concatvertices;
  //hack while I have not deal with a better gestion of the type /id
  //printvertices(vertexconcatnumber,concatvertices);
  //for(i=0;i<vertexconcatnumber;i++){
  //  for(j=0;j<concatvertices[i].degree;j++){
  //    printf("%d ", outlinearray[i][j].label);
  //  }
  //  printf("\n");
  //}
  create_concatenation_helper(vertexconcatnumber, concatvertices);
  //printcompatible(label.size);
  shift = 0;
  for(i=0; i<label.size/2;i++){shift+= (int) (pow(MAGICNUMBER,i+1)/2) ;}//used to create and  read the almost_foldable structures (base to encode vectors of small dimension)
  if(ALMOSTFOLDABLE){
    if(BACKBONETYPE){almost_foldable_path(vertexconcatnumber, concatvertices);}
    else{almost_foldable_tree(vertexconcatnumber, concatvertices);}
  }
  
  /********* structures useful for the concat and the fold, memory allocation once and for all here **********/

  //Some of the memory allocation are used for several variables which are not used at the same time
  
  maxoutlinesize = (maxdegree-1)*mapsize;
  Map map = initializemap();//create the map 
  FreeEdge *outline = malloc(maxoutlinesize*sizeof(FreeEdge));
  if (!outline) { fprintf(stderr,"cannot map.outline in concat.c\n"); exit(41); }
  previous = malloc(maxoutlinesize*sizeof(int));// the same memory is used for several arrays
  mat = malloc(maxoutlinesize*maxoutlinesize*sizeof(int));
  foldable_with = previous;
  fold_matrix = malloc((maxoutlinesize-1)*sizeof(int*));//compute the array for the dynamic programming
  if (!fold_matrix) { fprintf(stderr,"cannot malloc res in fold.c\n"); exit(51); }
  for(i=0; i < maxoutlinesize-1;i++) 
    {
      fold_matrix[i]=malloc(sizeof(int)*(maxoutlinesize/2+1));
      if (!fold_matrix[i]) { fprintf(stderr,"cannot malloc res[%d] in fold.c\n",i); exit(52);}
    }
  to_fold = malloc(maxoutlinesize/2*sizeof(pair));
  if (!to_fold) { fprintf(stderr,"cannot malloc to_fold in fold.c\n"); exit(53); }
  folded = malloc(maxoutlinesize/2*sizeof(pair));
  if (!folded) { fprintf(stderr,"cannot malloc folded in fold.c\n"); exit(54); }
  call_stack = previous;

  /************************** Allocation of the structures used to compute the cut-index *****************************************************/
  face = malloc(sizeof(pair)*mapsize*2);//at most twice more faces than vertices (reused)
  edge = malloc(sizeof(pair*)*mapsize*2);//used to mark the edges of the graph (should remove the 2)
  for(i=0;i<2*mapsize;i++){edge[i] = malloc(sizeof(pair)*mapsize*2);}
  dual = malloc(sizeof(DualMap));
  dual->vertexarray = malloc(sizeof(struct DualVertex)*mapsize*2);
  for(i=0;i<2*mapsize;i++){dual->vertexarray[i].edges = malloc(sizeof(struct DualEdge)*mapsize*6);}//in a planar graph at most 3 times more edges than vertices and we have the same number of edges as M
  edgebis = malloc(sizeof(int*)*mapsize*6);//used to mark the edges of the dual graph
  for(i=0;i<2*mapsize;i++){edgebis[i] = malloc(sizeof(int)*mapsize*6);}//same number of edges as M 
  tree = face;
  queue = malloc(sizeof(pair)*mapsize*2*64);//64 is a bound on  weigthmax, should be changed when we deal with large graphs, two dimension vector stored into a one dimensional vector

  /**************************** Allocation of the structures used to compute the signature ******************/
  
  rank = malloc(mapsize*sizeof(int)); 
  rank_init = malloc(mapsize*sizeof(int));
  for(i=0;i<mapsize;i++) rank_init[i] = -1;
  edge_stack = malloc(mapsize*sizeof(int));

  /************************************** The generation procedure *************************************************/
  chrono();

  printf("2. Build of the maps\n");

  unsigned long long int backbonenumber;
  
  if(BACKBONETYPE) backbonenumber = generate_paths(vertexnumber, &avl, map, outline);
  else backbonenumber =  generate_trees(vertexnumber, &avl, map, outline);
  
  printf("\r\033[2K   Nb backbones: "); prettyprint(backbonenumber); printf("\n");
  printf("  Nb folded maps       : "); prettyprint(avl.nbcalls); printf("\n");
  printf("  Nb unique maps         : "); prettyprint(avl.nbnodes); printf("\n");
  /* printf("  Cut Index lower      : "); printf("%25.3lf\n",avl.lower_cutindex);
  printf("  Cut Index limit      : "); printf("%25.3lf\n",avl.limit_cutindex);
  printf("  Cut Index higher     : "); printf("%25.3lf\n",avl.higher_cutindex);*/
  double tps = chrono();


/*************************** Compute the indices for each generated graph ****************************/ 


  printf("3. Computation of the indices\n");
  compute_all_indices(avl);
  double tps_indices = chrono();
   /*************************** Creation of the chemdraw and  html files ****************************/ 

  printf("4. Write of the results\n");
  char LongDirName[512];
  DirToWrite(argv[1], mapsize, LongDirName);
  Output(argv[1],LongDirName,mapsize,avl,vertices);
  print_all_indices(LongDirName,avl);
  //printf("Complexité moyenne de la signature sur les maps de taille %d :%d  \n",mapsize, compute_average(avl.root)/avl.nbnodes);
  chrono();
  print_global_indices(LongDirName,backbonenumber,mapsize,avl,vertices,tps,tps_indices);

  /*********************************** Free of all allocated memory ****************************/ 
  printf("5. Free of the memory\n");
  
  abr_free(avl.S);
  avl_free(avl.root);
  //Immutable structure
  if(ALMOSTFOLDABLE){
    if(BACKBONETYPE){
      for (i=0 ; i < mapsize ; i++) {
	for (j=0 ; j<label.size ; j++) free(almostfoldablepath[i][j]);
	free(almostfoldablepath[i]);
      }
      free(almostfoldablepath);
    }
    else{
      for (i=0 ; i < mapsize ; i++) {free(almostfoldabletree[i]);}
      free(almostfoldabletree);
    }
  }
  for (i=0 ; i < vertexconcatnumber ; i++) free(concatvertices[i].edges);
  free(concatvertices);
  if(METAMOTIF) free(vertices);
  for (i=0 ; i<label.size ; i++) free(connection[i].list);
  free(connection);
  free(label.list);

  // mutable structures used for the map generation or index computation
  free(outline);
  for (i=0 ; i < map.vertexnumber ; i++) {free(map.vertexarray[i].edges);}
  free(map.vertexarray);//free the map
  free(previous);
  free(mat);
  free(to_fold);
  free(folded);
  for (i=0 ; i<maxoutlinesize-1 ; i++) {free(fold_matrix[i]);}
  free(fold_matrix);
  free(face);
  for(i=0;i<2*mapsize;i++){free(edge[i]);free(dual->vertexarray[i].edges);free(edgebis[i]);}
  free(edgebis);
  free(dual->vertexarray);
  free(edge);
  free(dual);
  free(queue);
  free(rank);
  free(rank_init);
  free(edge_stack);
  
  chrono();
  //print_global_indices(LongDirName,backbonenumber,mapsize,avl,vertices,tps,tps_indices);
  printf("## END ##\n");
  return 0;
}
