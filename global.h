extern Vertex *vertices;// Vertices given in the mot and their rotation
extern Map *translation; //store the metamotifs as maps
extern FreeEdge **outlinearray;//free edge of the metamotifs
extern Label label; //Label of the edges
extern CompatibleVertices *connection;//helper structure to do the concatenation fast
extern unsigned int *** almostfoldablepath;//structure to cut in the enumeration of paths
extern unsigned int ** almostfoldabletree;//structure to cut in the enumeration of trees
extern Vertex *concatvertices;	
extern int mapsize,maxdegree,maxoutlinesize,shift;
extern int *previous,*mat,*foldable_with,**fold_matrix,*call_stack;
extern pair *to_fold,*folded;
extern pair *face ,**edge, *tree,*queue;
extern DualMap *dual; 
extern int **edgebis;
extern int *rank, *rank_init, *edge_stack;