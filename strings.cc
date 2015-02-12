/*
gcc -Wall strings.cc -o strings

./strings > s.txt

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main (void)
{
  int elements=6000;			//desired elements
  int size=4;    			//DNA size
  char text[]="ACGT";			//DNA characters
//  char text[4]={'A','C','G','T'};	//DNA characters
  
  char *matrix;
  posix_memalign ((void **) &matrix, 16, 30 * sizeof(char) );
  int r, j=0;
  srand ( time(0) );

  printf("./nw-sse3 -a ");

  for ( int i = 0; i < elements; ++i )
  {  
    matrix[i]= text[ rand() % size ];
    printf("%c", matrix[i] );  
  }

  printf(" -b ");

  for ( int i = 0; i < elements; ++i )
  {  
    r= rand() % size;

    if (r==0)
    {     
      matrix[i]= text[ rand() % size ];
      j++;
    }

    printf("%c", matrix[i] );  
  }

  printf(" -m 0 -s 1 -g 1 -e 1");

  return 0;
}
