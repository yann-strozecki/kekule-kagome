#include "structure.h"
#include "input.h"
#include "fold.h"
#include "avl.h"
#include "chrono.h"
#include "global.h"


#define CURRENTLABEL concatvertices[map.vertexarray[currentsize].type].edges[currentedge]



int smallest_representation(Map *map)
{//the id should be in the map instead of the type: then simplify the condition
  int i;
  for(i=0; i < map->vertexnumber/2 ; i++)
    {
      if(concatvertices[map->vertexarray[i].type].id  > concatvertices[map->vertexarray[map->vertexnumber-i-1].type].id)
	{
	  return 0;
	}
      if(concatvertices[map->vertexarray[i].type].id  < concatvertices[map->vertexarray[map->vertexnumber-i-1].type].id)
	{
	  return 1;
	}
    }
  return 1;
}

int smallest_rep(int *metastack){
  int i;
  for(i=0;i<mapsize/2;i++){
    if(metastack[i] > metastack[mapsize -i -1]) return 0;
    if(metastack[i] < metastack[mapsize -i -1]) return 1;
  }
  return 1;
}

Map initializemap(){//create a new map with value -1 on the edges
	int i,j;
	Map map;
	map.vertexnumber = mapsize;
	map.edgenumber = 0;
	map.vertexarray = malloc(sizeof(struct VertexMap)*mapsize);
	if (!map.vertexarray) { fprintf(stderr,"cannot malloc map.vertexarray in concat.c\n"); exit(42); }
	for (i=0 ; i<mapsize ; i++) {
		map.vertexarray[i].edges = malloc(sizeof(struct Edge)*maxdegree);
		if (!map.vertexarray[i].edges) { fprintf(stderr,"cannot malloc map.vertexarray[%d].edges in concat.c\n",i); exit(43); }
		for(j=0;j<maxdegree;j++) {
			map.vertexarray[i].edges[j].vertexindex = -1;
			map.vertexarray[i].edges[j].edgeindex = -1;
			}
		map.signature = NULL;
		}
	return map;
}


void prettyprint2(FILE *F, unsigned long long int x, char *s) {
	int t,g,m,k,u;
	t = x/1e12; x = x-t*1e12;
	g = x/1e9;  x = x-g*1e9;
	m = x/1e6;  x = x-m*1e6;
	k = x/1e3;  u = x-k*1e3;
	if (t) fprintf(F,"%d,%03d,%03d,%03d,%03d",t,g,m,k,u);
	else if (g) fprintf(F,"%d,%03d,%03d,%03d",g,m,k,u);
	else if (m) fprintf(F,"%d,%03d,%03d",m,k,u);
	else if (k) fprintf(F,"%d,%03d",k,u);
	else        fprintf(F,"%3d",u);
	fprintf(F," %s",s);
}

void affichage_dynamique(unsigned long long int numberofbackbone, Avl *avl) {
	static unsigned int lemodulo=0xff; 
	static double start_time;
	static double last_time;
	double current_time = what_time_is_it();
	if (!numberofbackbone) { start_time = last_time = current_time; }
	if (!(numberofbackbone & lemodulo)) {
		if (current_time-last_time > 1.0) lemodulo = (lemodulo>>1);
		if (current_time-last_time < 0.5) lemodulo = (lemodulo<<1)+1;
		last_time = current_time;
		printf("\r\033[2K");
		printf("   %.1lf sec   ",current_time-start_time);
		prettyprint2(stdout,numberofbackbone,"quasi-uniq backbones   ");
		prettyprint2(stdout,avl->nbcalls,"folded maps   ");
		prettyprint2(stdout,avl->nbnodes,"uniq maps");
		fprintf(stdout," [%.3f %.3f %.3f] cutindex\r",avl->lower_cutindex,avl->limit_cutindex,avl->higher_cutindex);
		fflush(stdout); 
		}
}

//Function to add a metamotif to the Map we are building

void add_metamotif(Map *m, int motifindex, int shift, int vertex, int edge){//connect meta in m at the shift position through edge
  
    int i,j;
    for(i = 0; i < translation[motifindex].vertexnumber; i++){
      m->vertexarray[i+shift].degree = translation[motifindex].vertexarray[i].degree;
      m->vertexarray[i+shift].type = translation[motifindex].vertexarray[i].type;
      for(j=0; j < translation[motifindex].vertexarray[i].degree; j++){
	if(translation[motifindex].vertexarray[i].edges[j].vertexindex != -1){
	  m->vertexarray[i+shift].edges[j].vertexindex = translation[motifindex].vertexarray[i].edges[j].vertexindex + shift;
	  m->vertexarray[i+shift].edges[j].edgeindex = translation[motifindex].vertexarray[i].edges[j].edgeindex;
	}
	else{m->vertexarray[i+shift].edges[j].vertexindex = -1;}
      }
    }
    if(shift != 0){//otherwise, it is the first vertex and no edge has to be changed
    int attachedvertex = outlinearray[motifindex][0].vertexindex + shift;
    int attachededge = outlinearray[motifindex][0].edgeindex;
    m->vertexarray[vertex].edges[edge].vertexindex = attachedvertex;
    m->vertexarray[vertex].edges[edge].edgeindex = attachededge;
    m->vertexarray[attachedvertex].edges[attachededge].vertexindex = vertex; 
    m->vertexarray[attachedvertex].edges[attachededge].edgeindex = edge;
    }
}


//Main function which create every possible path and then for each call folding_enumeration to create all possible folding
unsigned long long int generate_paths(int verticesnumber, Avl *avl, Map map, FreeEdge *outline)
{
  static unsigned long long int backbonenumber = 0;
  int i,j, currentsize, currentedge, labelvalues,currentlabel;
  int *connectionstack = malloc(mapsize*sizeof(int));
  if (!connectionstack) { fprintf(stderr,"cannot malloc connectionstack in concat.c\n"); exit(41); }
  int* metastack;
  if(METAMOTIF){//code path when we concatenate metamotif 
    int vertex,edge,currenttype,newtype;
    int * edgestack = malloc(mapsize*sizeof(int));
    metastack = malloc(mapsize*sizeof(int));//save the type of meta motif in the path already built
    //The concatenation begin here
    int bound = (BACKBONETYPE==2) ? 1 : verticesnumber; //optimization for cycle begin by one type of vertex only: we assume that all motifs will be used in the maps 
    for(i=0; i<bound; i++)
      {//optimization : go over the first vertices only (not the rotated ones)
	//printf("nouveau premier vertex \n");
	for(j=0; j < concatvertices[i].degree; j++)
	  { 
	    if (non_isomorph(concatvertices[i].edges,concatvertices[i].degree,j))//optimization: do not consider two edges which are equivalent in their neighborood for the first vertex
	      {
		add_metamotif(&map,i,0,0,0);
		currentsize = translation[i].vertexnumber -1;//size -1 of the path
		metastack[currentsize]= i;
		connectionstack[currentsize] = 0;
		labelvalues = shift - concatvertices[i].labelvalue;
		currentedge = j; //the edge we want to connect the next vertex to
		while(currentsize >=0)
		  {
		    currenttype = metastack[currentsize];
		    currentlabel = outlinearray[currenttype][currentedge].label;
		    //printf("size % d, type %d, currentedge %d \n",currentsize, currenttype,currentedge);
		    //printf(" %d \n",!almostfoldablepath[currentsize+1][currentlabel][labelvalues]);
		    while(currentedge < concatvertices[currenttype].degree && (connectionstack[currentsize] == connection[currentlabel].vertexnumber || (ALMOSTFOLDABLE && !almostfoldablepath[currentsize+1][currentlabel][labelvalues]) ))
		      //move to the next free edge in the current vertex if the previous has been used or cannot be used to generate an almost foldable path
		      {
			connectionstack[currentsize]=0;//reinitialize the possible connections
			currentedge++;//use the next edge
			if(currentedge < concatvertices[currenttype].degree) currentlabel = outlinearray[currenttype][currentedge].label;
			if(currentsize == 0){currentedge =  map.vertexarray[currentsize].degree;}//when dealing with the first vertex skip to the last edge so that the program exit the while loop
		      } 
		    if(currentedge < concatvertices[currenttype].degree)//test whether there is an edge left to add some motif 
		      {
			newtype = connection[currentlabel].list[connectionstack[currentsize]];
			//printf("newtype %d currentedge %d currentlabel %d, stack %d \n",newtype,currentedge, currentlabel,connectionstack[currentsize]);
			if(currentsize + translation[newtype].vertexnumber < mapsize){
			  //the vertex and edge used to do the connection
			  edgestack[currentsize]= currentedge;
			  vertex = outlinearray[currenttype][currentedge].vertexindex + currentsize + 1 - translation[currenttype].vertexnumber;
			  edge = outlinearray[currenttype][currentedge].edgeindex;
			  labelvalues +=  -concatvertices[newtype].labelvalue;//add the value of the labels of the new vertex
			  add_metamotif(&map,newtype,currentsize+1,vertex,edge);
			  currentsize += translation[newtype].vertexnumber;
			  metastack[currentsize] = newtype;
			  currentedge = 1; //the edge 0 is connected to the previous vertex therefore we consider the next one
			  connectionstack[currentsize]=0;
			}
			else{ //skip to the next meta motif if it makes the path larger than mapsize
			  connectionstack[currentsize]++;
			}
		      }
		    if(currentsize == mapsize-1 && (!ALMOSTFOLDABLE || (labelvalues == shift)) && smallest_rep(metastack))// 
		      //smallest representation devrait etre faite sur metastack
		      //produce a complete backbone when the path is of the right size and the path is almost foldable
		      {
			//printmap(map);
			affichage_dynamique(backbonenumber, avl);
			backbonenumber++;//should increment this value before testing smallest representation and almost foldable
			if(FOLDSTEP) folding_enumeration(map, avl, outline);//i and j are here to help produce cycles from paths
		      }
		    currenttype = metastack[currentsize];
		    if(currentsize == mapsize-1 || currentedge == concatvertices[currenttype].degree)//remove the current vertex if the path is finished or if we have connected all possible concatvertices to the current vertex
		      {
			labelvalues +=  concatvertices[currenttype].labelvalue;//subtract the labelvalue of the deleted vertex
			currentsize -= translation[currenttype].vertexnumber;
			if(currentsize >= 0)//remove an edge to the new current motif
			  {
			    currentedge = edgestack[currentsize];
			    currenttype = metastack[currentsize];
			    vertex = outlinearray[currenttype][currentedge].vertexindex + currentsize + 1 - translation[currenttype].vertexnumber;
			    edge = outlinearray[currenttype][currentedge].edgeindex;
			    map.vertexarray[vertex].edges[edge].vertexindex=-1;//remove the edge to the previous current vertex 
			    connectionstack[currentsize]++; //use the next possible vertex in the future connection
			  }
		      }
		  }
	      } 
	  }
      }
  }
  else{//code path when we concatenate single motifs
    //The concatenation begin here
    int bound = (BACKBONETYPE==2) ? 1 : verticesnumber; //optimization for cycle begin by one type of vertex only: we assume that all motifs will be used in the maps 
    for(i=0; i<bound; i++)
      {//optimization : go over the first vertices only (not the rotated ones)
	for(j=0; j < concatvertices[i].degree; j++)
	  { 
	    if (non_isomorph(concatvertices[i].edges,concatvertices[i].degree,j))//optimization: do not consider two edges which are equivalent in their neighborood for the first vertex
	      {
		map.vertexarray[0].type = i;   //initialize the first vertex of the path
		map.vertexarray[0].degree = concatvertices[i].degree;
		currentsize=0; //size -1 of the current path
		connectionstack[0] = 0;
		labelvalues = shift - concatvertices[i].labelvalue;//on devrait intégrer le shift ici
		currentedge = j; //the edge we want to connect the next vertex to
		while(currentsize >=0)
		  {
		    while(currentedge < map.vertexarray[currentsize].degree && (connectionstack[currentsize] == connection[CURRENTLABEL].vertexnumber || (ALMOSTFOLDABLE && !almostfoldablepath[currentsize+1][CURRENTLABEL][labelvalues]) ))//move to the next free edge in the current vertex if the previous has been used or cannot be used to generate an almost foldable path
		      //should add a condition with the size in the almostfoldable part
		      {
			connectionstack[currentsize]=0;//reinitialize the possible connections
			currentedge++;//use the next edge
			if(currentsize == 0){currentedge =  map.vertexarray[currentsize].degree;}//when dealing with the first vertex skip to the last edge so that the program exit the while loop
		      } 
		    if(currentedge < map.vertexarray[currentsize].degree)//test whether there is an edge left to add some motif
		      {
		
			//create the new vertex in the map and the corresponding edge 
			map.vertexarray[currentsize].edges[currentedge].vertexindex = currentsize+1;//create the new edge in the current vertex,  
			map.vertexarray[currentsize].edges[currentedge].edgeindex=0;  
			map.vertexarray[currentsize+1].type = connection[CURRENTLABEL].list[connectionstack[currentsize]];
			currentsize++;//the map has now one more vertex
			labelvalues +=  -concatvertices[map.vertexarray[currentsize].type].labelvalue;//add the value of the labels of the new vertex
			map.vertexarray[currentsize].degree = concatvertices[map.vertexarray[currentsize].type].degree;
			map.vertexarray[currentsize].edges[0].vertexindex = currentsize -1;//could be put in the map once and for all
			map.vertexarray[currentsize].edges[0].edgeindex = currentedge;//create the edge in the new vertex
			currentedge = 1; //the edge 0 is connected to the previous vertex therefore we consider the next one
			connectionstack[currentsize]=0;
		      }
		    if(currentsize == mapsize-1 && (!ALMOSTFOLDABLE || (labelvalues==shift)) && smallest_representation(&map) ) //produce a complete backbone when the path is of the right size and the path is almost foldable
		      {
			affichage_dynamique(backbonenumber, avl);
			backbonenumber++;//should increment this value before testing smallest representation and almost foldable
			if(FOLDSTEP) folding_enumeration(map, avl, outline);
		      }
		    if(currentsize == mapsize-1 || currentedge == map.vertexarray[currentsize].degree)//remove the current vertex if the path is finished or if we have connected all possible concatvertices to the current vertex
		      {
			labelvalues +=  concatvertices[map.vertexarray[currentsize].type].labelvalue;//subtract the labelvalue of the deleted vertex
			currentedge = map.vertexarray[currentsize].edges[0].edgeindex;
			currentsize--;
			if(currentsize >= 0)
			  {
			    map.vertexarray[currentsize].edges[currentedge].vertexindex=-1;//remove the edge to the previous current vertex 
			    connectionstack[currentsize]++; //use the next possible vertex in the future connection
			  }
		      }
		  }
	      } 
	  }
      }
  }
  free(connectionstack);
  return backbonenumber;
}

 
 



//Function which create every possible tree and then for each call folding_enumeration to create all possible folding

unsigned long long int generate_trees(int verticesnumber, Avl *avl, Map map, FreeEdge *outline)
{
  int size=1, free_edge_number=0,i,j,currentvertex,currentedge,currentlabel,newtype;
  static unsigned long long int backbonenumber = 0;
  pair *edgestack = malloc(maxdegree*mapsize*sizeof(pair));//first field vertex, second field edge
  pair **enumeration_position = malloc(mapsize*sizeof(pair*));//what is the next connection and in which context does it holds 
  for(i=0;i<mapsize;i++) enumeration_position[i] = malloc(maxdegree*sizeof(pair));
  
  map.vertexarray[0].type = 0;   //initialize the first vertex of the tree as vertex of type 0
  map.vertexarray[0].degree = concatvertices[0].degree;
  int labelvalues = concatvertices[0].labelvalue;
  for(i=0;i<map.vertexarray[0].degree;i++) //use free_edge_number instead of i ?
    {
      edgestack[free_edge_number].first = 0;
      edgestack[free_edge_number].second = i;
      enumeration_position[0][i].first=0;
      free_edge_number++;
    }//begin by the first type of vertex only 
  while(1){
    //printf("size %d free_edge_number %d\n",size,free_edge_number);
    //for(i=0;i<free_edge_number;i++) printf("(%d,%d)  ", edgestack[i].first, edgestack[i].second);
    //printf("\n");
    if(size < mapsize && free_edge_number > 0 && (!ALMOSTFOLDABLE || almostfoldabletree[size][shift - labelvalues]))
      {//we add a motif to the first free edge in the tree
	free_edge_number--;//the current edge will not be free anymore
	currentvertex = edgestack[free_edge_number].first;
	currentedge = edgestack[free_edge_number].second;
	map.vertexarray[currentvertex].edges[currentedge].vertexindex = size;
	map.vertexarray[currentvertex].edges[currentedge].edgeindex=0;
	map.vertexarray[size].edges[0].vertexindex= currentvertex;
	map.vertexarray[size].edges[0].edgeindex= currentedge;
	currentlabel = concatvertices[map.vertexarray[currentvertex].type].edges[currentedge];
	newtype = connection[currentlabel].list[enumeration_position[currentvertex][currentedge].first];
	map.vertexarray[size].type = newtype;
	labelvalues += concatvertices[newtype].labelvalue;
	map.vertexarray[size].degree = concatvertices[newtype].degree;
	//printf("currentvertex %d currentedge %d currentlabel %d type %d\n",currentvertex,currentedge,currentlabel,newtype);
	for(i=1;i<concatvertices[newtype].degree;i++) {
	  enumeration_position[size][i].first = 0;
	  edgestack[free_edge_number].first = size;
	  edgestack[free_edge_number].second = i;
	  free_edge_number ++;
	}
	size++;
      }
    else
      {//remove the last vertex
	//remove the edge of the vertex
	//printf("remove \n");
	size--;
	if(size == 0) break;
	currentvertex = map.vertexarray[size].edges[0].vertexindex;
	currentedge = map.vertexarray[size].edges[0].edgeindex;
	//printf("currentvertex % d currentedge %d\n",currentvertex,currentedge);
	enumeration_position[currentvertex][currentedge].first++;
	if(free_edge_number == 0){//case with no free edge usable : we must add back the frozen edges
	  //printf("add back frozen edges \n");
	  for(i=0;i<size;i++){
	    for(j=0;j<map.vertexarray[i].degree;j++){
	      if(map.vertexarray[i].edges[j].vertexindex == -1 && size < enumeration_position[i][j].second){//if an edge is used in a context we have just changed, free it back
		edgestack[free_edge_number].first = i;
		edgestack[free_edge_number].second = j;
		free_edge_number ++;
		enumeration_position[i][j].first = 0;
		enumeration_position[i][j].second = 0;
	      }
	    }
	  }
	}
	else{free_edge_number += 1 - concatvertices[map.vertexarray[size].type].degree;}//remove the unused free edges
	if(size == mapsize-1 ){//case with a full size map
	  affichage_dynamique(backbonenumber, avl);
	  backbonenumber++;//should increment this value before testing smallest representation and almost foldable
	  if(FOLDSTEP && !(ALMOSTFOLDABLE && labelvalues))
	    {
	      folding_enumeration(map, avl, outline);
	    }
	  //printmap(map);
	}//remove the free edges of the motif 
	labelvalues -= concatvertices[map.vertexarray[size].type].labelvalue;
	currentlabel =  concatvertices[map.vertexarray[currentvertex].type].edges[currentedge];
	map.vertexarray[currentvertex].edges[currentedge].vertexindex =-1;
	//printf("currentvertex %d currentedge %d pos %d degree %d currentlabel %d \n",currentvertex,currentedge,enumeration_position[currentvertex][currentedge], connection[currentlabel].vertexnumber,currentlabel);
	if(enumeration_position[currentvertex][currentedge].first < connection[currentlabel].vertexnumber) //if there is still possible connections add back the edge we have freed by removing the last motif
	  {
	    //printf("add edge back \n");
	    edgestack[free_edge_number].first = currentvertex;
	    edgestack[free_edge_number].second = currentedge;
	    free_edge_number++;
	  }
	else{enumeration_position[currentvertex][currentedge].second = size;}
      }
  }
  free(edgestack);
  for(i=0;i<mapsize;i++) free(enumeration_position[i]);
  free(enumeration_position);
  return backbonenumber;
}



