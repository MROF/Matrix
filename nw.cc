#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define MAX(x,y) ((x) > (y) ? (x) : (y))

//==============================================
//=================== Main =====================

int main(int argc, char * argv[] )
{
  char *a;
  char *b;

  if (argc != 3)
   {
     fprintf (stderr, "syntax: [str1] [str2]\n");
     return (EXIT_FAILURE);
   }

  b = argv[1];
  a = argv[2];

// b = strdup("AAAC");
// a = strdup("AGC");

  int match=1;
  int mismatch=-1;
  int gap=-2;

  int i, j;
  int score;

  int alen= strlen(a)+1;
  int blen= strlen(b)+1;
  
  int matrix[alen][blen];

  for(i=0;i<alen;i++)
    for(j=0;j<blen;j++)
      matrix[i][j]=0;

//-------------------Initializing The Matrix-------------------
  matrix[0][0]=0;

  for(i=1;i<alen;i++)
    matrix[i][0] = matrix[i-1][0]+gap;

  for(j=1;j<blen;j++)
    matrix[0][j] = matrix[0][j-1]+gap;

//---------------------Filling The Matrix----------------------
  for(i=1;i<alen;i++)
    for(j=1;j<blen;j++)
    {
      if (a[i-1]==b[j-1]) score=match;  
      else score=mismatch;

      matrix[i][j]= MAX( matrix[i-1][j-1]+score,  MAX(matrix[i-1][j]+gap,  matrix[i][j-1]+gap));
    }

//---------------------printing the matrix---------------------
  for(i=0;i<alen;i++)
  {
    for(j=0;j<blen;j++)
      printf("%4d", matrix[i][j]);
    printf("\n");
  }

 return(0) ;
}
