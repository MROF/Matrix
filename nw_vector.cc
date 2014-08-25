//  ./nw-avx -a AAAA -b ADERRTTYU -m 0 -s -1 -g -3 -e -1

#include <x86intrin.h>
#include <emmintrin.h>
#include <immintrin.h>

#ifdef __AVX
#include <x86intrin.h>
#include <emmintrin.h>
#include <immintrin.h>
#endif
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef __AVX
#define BYTE_ALIGNMENT 32
#else
#define BYTE_ALIGNMENT 16
#endif

int inverse(float arr[],int x,int y); // function to reverse the array

void print(char ch, float *arr)
{
  printf("%c  ",ch);
  for(int k=0; k<8; k++)
    printf("%4.0f ", arr[k] );
  printf("\n");
}

int main(int argc, char * argv[] )
{
  printf("\n--------------------------------------\n");

  char *a;
  char *b;
  int match;  
  int mismatch;
  int gapopen;
  int gapextend;

  int option=0;
  int padding=BYTE_ALIGNMENT / 4;   //padding = 8 or 4 
  int qq[padding];      //gapopen and extention vectors
  int rr[padding];
  int score[padding];
  int i, j;

  char c;
  #if (defined(__AVX))
   __m256 h, t1, x;
   x= _mm256_setzero_ps();
  #endif

  if (argc == 1)
   {
     fprintf(stderr, " syntax: %s arguments\n", argv[0]);
     printf(" use the following parameters -a : -b : -m : -s : -g : -e :\n");
     return (0);
   }

  while ((option = getopt(argc, argv,"a:b:m:s:g:e:")) != -1)
  {
    switch (option) 
    {
      case 'a' : a = optarg;		  break;
      case 'b' : b = optarg;		  break;
      case 'm' : match = atoi(optarg); 	  break;
      case 's' : mismatch = atoi(optarg); break;
      case 'g' : gapopen= atoi(optarg);   break;
      case 'e' : gapextend= atoi(optarg); break;
      default  : printf(" use the following parameters -a : -b : -m : -s : -g : -e :\n");
		 exit(EXIT_FAILURE);
    }
  }

  printf ("Match:%4d   Mismatch:%d\n", match, mismatch);
  printf ("Gapopen:%2d   gapextend:%d\n\n", gapopen, gapextend);

  for (i=0; i<padding;i++)   //fill the gap vectors will penalties
  {
    qq[i]=gapopen;
    rr[i]=gapextend;
    score[i]=0;
  }

  int n= strlen(a);
  int m= strlen(b);
  int mm=m;

  if (m%padding) // not 4 or 8 greater common divisor    
  {
    mm = m + (padding - (m % padding)); //new size of the array
    for (j=m; j<mm;j++)
      b[j]='Z';      //padding letters
  }

  float *hh;
  float *f;
  posix_memalign ((void **)&hh, BYTE_ALIGNMENT, n * mm * sizeof(float));
  posix_memalign ((void **)&f, BYTE_ALIGNMENT,  n * mm * sizeof(float));

  for(i=0;i<n*mm;i++)
    hh[i]=f[i]=-99;

  hh[0]=0;

  for(i=1; i<mm; i++)
  {
    hh[n*i]= f[n*i]= gapopen+(i*gapextend);
  
    for (j=1; j<n;j++)
    {
      c=a[j-1];
      hh[j]= f[j]= gapopen+(j*gapextend);
    }

  }

  int y= mm/padding;
  float *output;

  posix_memalign ((void **)&output, BYTE_ALIGNMENT, padding * sizeof(float));

  inverse(hh,mm,n); // inverse the matrix to work with vector approach 

  for (j=0; j<n;j++)
  {
    for (i=0;i<y;i++)
    { 
      h= _mm256_load_ps( &hh [(j*padding*2)+(padding*i)]) ; 


      _mm256_store_ps(output, h);       //|
      output[0]=output[7];              //| shiftr7 x0000000
      for(int k=1;k<padding;k++)        //|    
        output[k]=0;                    //|
/**/  print('t', output);               //|
      t1= _mm256_load_ps(&output [0]) ; //| 
 
      _mm256_store_ps(output, h);       //|
/**/  print('h', output);               //|
      for(int k=padding; k>0; k--)      //|shiftl 0xxxxxxx
        output[k]=output[k-1];          //|
      output[0]=0;                      //|
      h= _mm256_load_ps( &output[0]);   //|

      h= _mm256_or_ps(h, x);
/**/  print('+', output);

      _mm256_store_ps(output, t1);     //copy t1 int x 
      x= _mm256_load_ps(&output [0]);
/**/  _mm256_store_ps(output, x);
      print('x', output);
   


      printf("\n");
    }
    printf("\n");
  }

 free(hh);
 free(f);
 free(output);
 return(0) ;
}

//============ Reversing the array =========
int inverse(float *arr,int x,int y)
{
  int l = 0;

  float *tmp = NULL;      // temp matrix
  posix_memalign ((void **)&tmp, BYTE_ALIGNMENT, x * y * sizeof(float)); 
  for (int j=0;j<y;j++)
    for (int i=0;i<x;i++)	
    {
      tmp[l]=arr[(y*i)+j];
      l++;
    }

  for (int i=0;i<x*y;i++)
   {
    arr[i]=tmp[i];
   }

  free(tmp);
  return (0);
}
