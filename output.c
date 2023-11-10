#include "output.h"
#include "colordot.h"
#include "chemdraw.h"
#include "string.h"
#include "global.h"

void Create_Dir(char *DirName)
{
  char Name[strlen(DirName)+10];
  
  mkdir("Results", 0755);
  sprintf(Name,"Results/%s",DirName);
  mkdir(Name,0755);
}


void GetName(char *InputFile, char *DirName)
{
  int i,j;

  i = 0;
  j = 0;
  while (InputFile[i] != '.')
    {
      if (InputFile[i] == '/') j = 0;
      else {
		DirName[j] = InputFile[i];
		j++;
      }
      i++;
    }
  DirName[j] = '-';j++; 
  sprintf(&DirName[j],"%d",BACKBONETYPE); j++;
  sprintf(&DirName[j],"%d",FOLDSTEP); j++;
  sprintf(&DirName[j],"%d",ALMOSTFOLDABLE); j++;
  DirName[j] = '\0';  
}

int EdgeArity(Map *M, int n, int k)
{
	int i;
	int nb = 0;
	for(i=0;i< M->vertexarray[n].degree; i++)
	{
		if (M->vertexarray[n].edges[i].vertexindex+1 == k)
			nb ++;
	}
	return nb;
}


void WriteFileCDXML(char *FileName, Map *M, Vertex *vertices)
{
  FILE *file;
  int i, coordX, coordY, temp, j;
  int compteur = 1;
  int last,nb;

  file = fopen(FileName,"w");
  fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
  fprintf(file, "<!DOCTYPE CDXML SYSTEM \"http://www.camsoft.com/xml/cdxml.dtd\" >\n");
  fprintf(file, "<CDXML\n");
  fprintf(file, "CreationProgram=\"ChemDraw 11.0\">\n");
  fprintf(file, "<page BoundingBox=\"0 0 600 600\" Header=\"%s\">\n", FileName);
  fprintf(file, "<fragment>\n");

  coordX = 100;
  coordY = 210;
  
  for ( i = 0; i < M->vertexnumber; i++ ) {
    fprintf(file, "<n\n");
    fprintf(file, "id= \"%d\"\n", i+1);
    fprintf(file, "p=\"%d %d\"\n", coordX, coordY);  
    fprintf(file, "NodeType = \"Element\"\n");
    fprintf(file, "Element = \"%d\"\n", chemdraw[vertices[M->vertexarray[i].type].letter-1]);
    fprintf(file, " LabelDisplay = \"Auto\">\n");
    fprintf(file, "</n>\n");
    if ( compteur == 1 ) {
      temp = coordX ;
      coordX = coordY;
      coordY = temp;
    }
    else if ( compteur == 2 ) {
      temp = coordX ;
      coordX = coordY;
      coordY = temp;
      coordY += 30;
      compteur = 0;
    }
      
    compteur++;
  }
  for ( i = 0; i < M->vertexnumber; i++ ) 
    {
		last = -1;
      for ( j = 0; j < M->vertexarray[i].degree; j++ )
	{ 
	  if (M->vertexarray[i].edges[j].vertexindex > i && (M->vertexarray[i].edges[j].vertexindex+1) != last)
	    {
	      fprintf(file, "<b\n");
	      fprintf(file, "B=\"%d\"\n", i+1);
	      
	      fprintf(file, "E=\"%d\"\n", M->vertexarray[i].edges[j].vertexindex+1);
	      nb = EdgeArity(M,i,M->vertexarray[i].edges[j].vertexindex+1);
		  if (nb > 1)
		  {
			fprintf(file, "Order=\"%d\"\n", nb);
		  fprintf(file, "DoublePosition=\"0\"\n");
		  }
	      fprintf(file, "/>\n");
	      last = M->vertexarray[i].edges[j].vertexindex+1;
	    }
	}
    }
  fprintf(file, "</fragment>\n");
  fprintf(file, "</page> \n");
  fprintf(file, "</CDXML> \n");
  fclose(file);
}

void WriteFileDOT(char *FileName, Map *M, Vertex *vertices)
{
  FILE *f;
  int i,j;
  char origine[64];
  char dest[64];
  int the_type;
  char the_letter;
  
  f = fopen(FileName,"w");
  fprintf(f,"graph G {\n");
  fprintf(f,"bgcolor=\"transparent\"");
  fprintf(f, "    node [shape=circle, height=0.6, width=0.6, color=lightblue2, style=filled];\n");
  fprintf(f, "    edge [color=blue];\n");
  fprintf(f, "\n");
  for(i=0;i < M->vertexnumber;i++)
  {
	the_type = M->vertexarray[i].type;
	the_letter = 'A' - 1 + vertices[the_type].letter;
    fprintf(f,"\" %d \" [",i);
    fprintf(f,"color=\"%s\", style=filled,",color[M->vertexautomorphic[i]]);
    fprintf(f,"label=\"%c%d\"];\n",the_letter,vertices[the_type].id);
  } 
  for(i=0;i < M->vertexnumber;i++)
    {
      sprintf(origine,"\" %d \"",i);
      for (j=0;j < M->vertexarray[i].degree;j++)
	{
	  if (M->vertexarray[i].edges[j].vertexindex > i)
	    {
	      sprintf(dest,"\" %d \"",M->vertexarray[i].edges[j].vertexindex);
	      fprintf(f,"%s -- %s;\n",origine,dest);
	    }
	}
    }
  fprintf(f, "}\n");
  
  fclose(f);
}

void WriteCDXML(char* LongDirName, char *DirName,struct AvlNode *a, int MapSize, Vertex *vertices, short lim, float high)
{
  char FileName[512];
  static int Nb = 1;
  static int Nb_Write = 0;
  if (a != NULL)
    {
      sprintf(FileName,"%s/%s_SAT_%02d_%06d.cdxml",LongDirName,DirName,MapSize,Nb);
      if ((lim == 0) || (lim ==1 && a->M->cutindex == high && Nb_Write < NBFILEMAX)) {
			WriteFileCDXML(FileName,a->M,vertices);
			Nb_Write ++;
	  }
      Nb++;
      WriteCDXML(LongDirName, DirName,a->l, MapSize,vertices,lim,high);
      WriteCDXML(LongDirName,DirName,a->r,MapSize,vertices,lim,high);
    }
}

void WriteDOT(char *LongDirName, char* DirName, struct AvlNode *a, int MapSize, Vertex *vertices, short lim, float high)
{
  char FileName[512];
  static int Nb = 1;
  static int Nb_Write = 0;

  if (a != NULL)
    {
      sprintf(FileName,"%s/%s_SAT_%02d_%06d.dot",LongDirName,DirName,MapSize,Nb);
    if ((lim == 0) || (lim ==1 && a->M->cutindex == high && Nb_Write < NBFILEMAX)) {
		  WriteFileDOT(FileName,a->M,vertices);
		Nb_Write ++;
	  }
      Nb++;
      WriteDOT(LongDirName,DirName, a->l, MapSize,vertices,lim,high);
      WriteDOT(LongDirName,DirName,a->r,MapSize,vertices,lim,high);
    }
}

int nb_size_faces(struct Faces F)
{
  int i;
  int cpt=0;

  for(i=1;i<=F.max;i++)
    {
      if (F.distrib[i]) cpt++;
    }
  return cpt;
}

void WriteIsoClass(FILE *f, Map *M, Vertex *vertices)
{
  int i,j,cpt;
  int the_type,the_id;
  int type[M->vertexnumber];
  char the_letter;
  int nb[M->vertexnumber];
  //int n;
  
  //n = 0;
  for(i=0;i<M->vertexnumber;i++) type[i] = 0;
  for(i=0;i<M->vertexnumber;i++)
  {
    the_type = M->vertexarray[i].type;
    the_id = vertices[M->vertexarray[i].type].id;
    if (type[the_id] == 0)
	{
		type[the_id] = 1;
		the_letter = 'A' - 1 + vertices[the_type].letter;
		cpt = 0;
		for (j=0;j<M->vertexnumber;j++) nb[j] = 0;
		for (j=0;j<M->vertexnumber;j++) 
	    { 
			if (vertices[M->vertexarray[j].type].id ==  the_id && nb[M->vertexautomorphic[j]] == 0) 
			{
				nb[M->vertexautomorphic[j]] = 1; 
				cpt ++;
			}
	    }
		//n++;
		fprintf(f,"%c%d : %d\n",the_letter, the_id, cpt);
	}
  }
  //fprintf(f,"-- %2d ",n);
}

void WriteComposition(FILE *f, Map *M, Vertex *vertices)
{
  int the_id,the_type,i,j,cpt;
  int type[M->vertexnumber];
  char the_letter;


  for(i=0;i<M->vertexnumber;i++) type[i] = 0;
  
  for(i=0;i<M->vertexnumber;i++)
    {
      the_type = M->vertexarray[i].type;
      the_id = vertices[M->vertexarray[i].type].id;
      the_letter = 'A' - 1 + vertices[the_type].letter;
      if (type[the_id] == 0)
	{
	  type[the_id] = 1;
	  cpt = 0;
	  for (j=0;j<M->vertexnumber;j++) if (vertices[M->vertexarray[j].type].id == the_id) cpt ++;
	  fprintf(f,"%c%d : %d\n",the_letter, the_id, cpt);
	}
    }
}

// A RAJOUTER : IMAGE CLIQUABLE pour devenir plus grande. Generation automatique des images. Différents indices.
void WriteIND(char *LongDirName,char *MotName, struct AvlNode *a, int MapSize, Vertex *vertices, short lim, float high)
{
	char FileName[512];
	static int Nb = 1;
	int i;
	FILE *f;
	static int Nb_Write = 0;

	sprintf(FileName,"%s/%s_SAT_%02d_%06d.ind",LongDirName,MotName,MapSize,Nb);
	
  if (a != NULL)
    {
		if ((lim == 0) || (lim ==1 && a->M->cutindex == high && Nb_Write < NBFILEMAX)) {
			f = fopen(FileName,"w");
			fprintf(f,"%s/SAT_%02d/%s_SAT_%02d_%06d.cdxml\n",MotName,MapSize,MotName,MapSize,Nb);
			fprintf(f,"%s/SAT_%02d/%s_SAT_%02d_%06d.dot.svg\n", MotName,MapSize,MotName,MapSize,Nb);
			fprintf(f,"%s/SAT_%02d/%s_SAT_%02d_%06d.dot\n",MotName,MapSize,MotName,MapSize,Nb);
			fprintf(f,"%5.2f \n",a->M->cutindex);
			WriteIsoClass(f,a->M,vertices);
			WriteComposition(f,a->M,vertices);
			fprintf(f,"%3d \n",a->M->faces.nb);
			fprintf(f,"%3d \n",nb_size_faces(a->M->faces));
			for(i=0;i<=a->M->faces.max;i++)
			{
				if (a->M->faces.distrib[i] == 0) fprintf(f,". ");
				else fprintf(f,"%llu ",a->M->faces.distrib[i]);
			}
			fclose(f);
			Nb_Write++;
	  }
	  Nb = Nb + 1;
      WriteIND(LongDirName,MotName,a->l,MapSize,vertices,lim,high);
      WriteIND(LongDirName,MotName,a->r,MapSize,vertices,lim,high);
    }
}



void print_global_indices(char *LongDirName, unsigned long long int backbonenumber, int MapSize, Avl avl, Vertex *vertices, double tps,double tps_indices)
{ 
	FILE *f;
	char FileName[512];

	sprintf(FileName,"%s/indices.txt",LongDirName);
	f = fopen(FileName,"w");
	fprintf(f,"%d\n",BACKBONETYPE);
	fprintf(f,"%d\n",ALMOSTFOLDABLE);
	fprintf(f,"%lf\n",tps); // temps de generation + éventuellement folding et isomorphisme
	fprintf(f,"%llu\n",backbonenumber);
	fprintf(f,"%llu\n",avl.nbcalls); 
 	fprintf(f,"%lu\n",avl.nbnodes);
	fprintf(f,"%lf\n",tps_indices); // temps jusqu'au calcul d'indices	
	fprintf(f,"%25.3lf\n",avl.lower_cutindex);
	fprintf(f,"%25.3lf\n",avl.limit_cutindex);
	fprintf(f,"%25.3lf\n",avl.higher_cutindex);
	fclose(f);
}
void DirToWrite(char *InputFile, int MapSize, char *LongDirName)
{
	char DirName[256];
	char MotName[256];
	char cmd[512];
	
	GetName(InputFile,MotName);
	Create_Dir(MotName);
	sprintf(LongDirName,"Results/%s",MotName);
	sprintf(cmd,"cp %s %s",InputFile,LongDirName);
	if (system(cmd) == -1) {
		printf("%s",cmd);
		exit(56);
	}
	strcat(LongDirName,"/");
	sprintf(DirName,"SAT_%02d",MapSize);
	strcat(LongDirName,DirName);
	mkdir(LongDirName,0755);
}

float HigherCutIndex (struct AvlNode *a)
{
	float cuitl,cuitr,max;
	if (a != NULL)
	{
		cuitl = HigherCutIndex(a->l);
		cuitr = HigherCutIndex(a->r);
		max = (cuitr<cuitl?cuitl:cuitr);
		if (max > a->M->cutindex) return max;
		else return a->M->cutindex;
	}
	return 0.0;
}

void Output(char *InputFile, char *LongDirName, int MapSize, Avl avl, Vertex *vertices)
{
  char MotName[256];
    
  // Recuperer le dir name à l'aide du nom du fichier
 // fprintf(stderr,"Récupération du dirname \n");
  GetName(InputFile,MotName);
  avl.higher_cutindex = HigherCutIndex(avl.root);
  //fprintf(stderr,"%f\n",avl.higher_cutindex);
  //exit(1);
  if (avl.nbnodes < NBFILEMAX) {
    fprintf(stderr,"\t writing CDXML \n");
    WriteCDXML(LongDirName,MotName,avl.root,MapSize,vertices,0,avl.higher_cutindex);
    fprintf(stderr,"\t writing DOT \n");
    WriteDOT(LongDirName,MotName,avl.root,MapSize,vertices,0,avl.higher_cutindex);
    fprintf(stderr,"\t writing IND \n");
    WriteIND(LongDirName,MotName,avl.root,MapSize,vertices,0,avl.higher_cutindex);
  }
  else 
  {
	fprintf(stderr,"\t writing CDXML \n");
    WriteCDXML(LongDirName,MotName,avl.root,MapSize,vertices,1,avl.higher_cutindex);
    fprintf(stderr,"\t writing DOT \n");
    WriteDOT(LongDirName,MotName,avl.root,MapSize,vertices,1,avl.higher_cutindex);
    fprintf(stderr,"\t writing IND \n");
    WriteIND(LongDirName,MotName,avl.root,MapSize,vertices,1,avl.higher_cutindex);
  }
}
