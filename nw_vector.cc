//  ./nw-avx -a ABAA -b ADERRTTYU -m 0 -s -1 -g -3 -e -1

#ifdef __AVX
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
  float gapopen= 0;
  float gapextend= 0;

  int option=0;
  int padding=BYTE_ALIGNMENT / 4;   //padding = 8 or 4 
  char c;
  char cc;

  #if (defined(__AVX))
    __m256 h;
    __m256 e;
    __m256 t1;
    __m256 x;
    __m256 f;
    __m256 score_vector;
    __m256 qq;    // gap open vector
    __m256 rr;    // gap extention vector

    x= _mm256_setzero_ps();             //initialize vectors 
    score_vector= _mm256_setzero_ps();
    qq= _mm256_set1_ps (gapopen);
    rr= _mm256_set1_ps (gapextend);
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
  printf ("Gapopen:%2.0f  gapextend:%2.0f\n\n", gapopen, gapextend);

  int n= strlen(a);
  int m= strlen(b);
  int mm=m;

  if (m%padding) // not 4 or 8 greater common divisor    
  {
    mm = m + (padding - (m % padding)); //new size of the array
    for (int j=m; j<mm;j++)
      b[j]='Z';      //padding letters
  }

  float *scorematrix; // matrix filled with match and mismatchi costs
  float *hh;
  float *ee;
  float *output;

  posix_memalign ((void **) &hh, BYTE_ALIGNMENT, n * mm * sizeof(float));
  posix_memalign ((void **) &ee, BYTE_ALIGNMENT, n * mm * sizeof(float));
  posix_memalign ((void **) &scorematrix, BYTE_ALIGNMENT,  n * mm * sizeof(float));
  posix_memalign ((void **)&output, BYTE_ALIGNMENT, padding * sizeof(float));

  int l=0;                        // Filling the score matrix
  for (int j=0; j<n; j++)
  {
    c=a[j];
    for (int i=0; i<mm; i++)
    {
      cc= b[i];
     
      if (c==cc)
        scorematrix[l] = match;
      else
        scorematrix[l] = mismatch;

      l++;      
    }
  }

  for(int i=0;i<n*mm;i++)
    hh[i]=ee[i]=-99;

//  hh[0]=0;       // if on then next i and j should start from 1

  for(int i=0; i<mm; i++)   //initializing the HH and EE matricies
  {
    hh[n*i]= ee[n*i]= gapopen+(i*gapextend);
  
    for (int j=0; j<n;j++)
      hh[j]= ee[j]= gapopen+(j*gapextend);
  }

  int y= mm/padding;

  inverse(hh,mm,n); // inverse the matrix to work with vector approach 

  for (int j=0; j<n;j++)
  {
    c=a[j];

    for (int i=0;i<y;i++)
    { 
      h = _mm256_load_ps( &hh [(j*padding*2)+(padding*i)]); 
      e = _mm256_load_ps( &ee [(j*padding*2)+(padding*i)]); 

      _mm256_store_ps(output, h);       //|
      output[0]=output[7];              //| shiftr7 x0000000
      for(int k=1;k<padding;k++)        //|    
        output[k]=0;                    //|
      print('t', output);               //|                                   //X
      t1= _mm256_load_ps(&output [0]) ; //| 
 
      _mm256_store_ps(output, h);       //|
      print('h', output);               //|                                   //x
      for(int k=padding; k>0; k--)      //|shiftl 0xxxxxxx
        output[k]=output[k-1];          //|
      output[0]=0;                      //|
      h= _mm256_load_ps( &output[0]);   //|

      h= _mm256_or_ps(h, x);
      print('+', output);                                                     //x

      _mm256_store_ps(output, t1);      //copy t1 into x 
      x= _mm256_load_ps(&output [0]);
      _mm256_store_ps(output, x);                                             //x
      print('x', output);                                                     //x

      score_vector= _mm256_load_ps( &scorematrix [(j*padding*2)+(padding*i)]) ; 
      h= _mm256_add_ps(h, score_vector);

      h= _mm256_max_ps (h, e);




      printf("\n");
    }
    printf("_______________\n");
  }

 free(hh);
 free(ee);
 free(output);
 free(scorematrix);

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
