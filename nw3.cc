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

  int alen= strlen(a)+1;
  int blen= strlen(b)+1;
  
  int mm[alen][blen];
  int lxy[alen][blen];

//--------------Initializing And Filling The Matricies------------

  mm[0][0]=0;

  for(i=1;i<alen;i++)
  {
    mm[i][0] = d+(i*e);
    lxy[i][0] = d+(i*e);

    for(j=1;j<blen;j++)
    {
      mm[0][j] = d+(j*e);
      lxy[0][j]= d+(j*e);

      if (a[i-1]==b[j-1]) score=match;  
      else score=mismatch;

      lxy[i][j]= e+ MAX( mm[i-1][j]+d, MAX( lxy[i-1][j], MAX( mm[i][j-1]+d, lxy[i][j-1]) ) );    
      mm[i][j]= MAX( mm[i-1][j-1] + score, lxy[i][j]);
    }
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
