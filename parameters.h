//all options of the code, used to select different modes of the program 
//or to execute only part of the pipeline 


#define MAGICNUMBER 256 //used as a basis to encode informations about the free labels in a map into a single integer
#define BACKBONETYPE 1 // 0 for generating tree, 1 for generating path, 2 for generating cycle
#define FOLDSTEP 1 // 0 for generating backbone only, 1 to get the folding also
#define ALMOSTFOLDABLE 1 //0 for generating all paths, 1 for generating the almost foldable only
#define NBMAXNODESINAVL 1000000 // Maximum number of nodes stored in the AVL
#define NBFILEMAX 2000 // Maximum number of files written in the Results directory
#define SELECTIONCUTINDEX 0// Use the cutindex as the first criteria to compare signatures, option to remove
#define METAMOTIF 0 //do not use metamotifs with trees, it is not useful and not implemented
