#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define MAX(x,y) ((x) > (y) ? (x) : (y))

int main(int argc, char * argv[] )
{
  char * a;
  char *b;

  if (argc !=3) 
  {
    fprintf(stderr, "Syntax: [str1] [str2]\n");
    return(EXIT_FAILURE);
  }

  a = argv[1];
  b = argv[2];

//  b = strdup("ATTCA");
//  a = strdup("TCA");

  int match=1;
  int mismatch=-1;
  int d=-3;
  int e=-1;

  int i, j;
  int score;
  int inf=-99999;

  int alen= strlen(a)+1;
  int blen= strlen(b)+1;
  
  int mm[alen][blen];
  int lx[alen][blen];
  int ly[alen][blen];

//------------------Initializing The Matricies-------------------

  mm[0][0]=0;

  for(i=1;i<alen;i++) 
  { 
    mm[i][0] = d+(i*e);
    ly[i][0] = inf;
    lx[i][0] = d+(i*e);
  }
		
  for(j=1;j<blen;j++)
  {
    mm[0][j] = d+(j*e);
    lx[0][j]= inf;
    ly[0][j]= d+(j*e);
  }
    
//---------------------Filling The Matricies----------------------

  for(i=1;i<alen;i++)
    for(j=1;j<blen;j++)
    {
      if (a[i-1]==b[j-1]) score=match;  
      else score=mismatch;

      lx[i][j]= e+ MAX( mm[i-1][j]+d, lx[i-1][j]);     //check the previous rows
		
      ly[i][j]= e+ MAX( mm[i][j-1]+d, ly[i][j-1]);    //check the previous column   
  
      mm[i][j]= MAX( mm[i-1][j-1] + score, MAX ( lx[i][j], ly[i][j] ));
    }
//---------------------printing the matricies---------------------

  printf("__________________________________\n");
  for(i=0;i<alen;i++)   
  {
    for(j=0;j<blen;j++)
      printf("%4d", mm[i][j]);
    printf("\n");
  }

 return(0) ;
}
