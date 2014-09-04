//  ./nw-avx -a ATAGAAGTAG -b TCAGTCAGTCAGAGC -m 0 -s 1 -g 1 -e 1

#include <emmintrin.h>
#include <immintrin.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define BYTE_ALIGNMENT 32

//int inverse(float arr[],int x,int y); // function to reverse the array
void print(__m256 vector, int padding);
__m256 shiftr7(__m256 vector, int padding);
__m256 shiftl1(__m256 vector, int padding);

//================ Main ===================
//=========================================

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
  int n;
  int m;
  int mm;

//  __m256 x;
  __m256 h;
  __m256 e;
  __m256 t1;
//    __m256 f;
  __m256 score;
  __m256 qq;    // gap open vector
  __m256 rr;    // gap extention vector
  __m256 n_vector;

  float *score_matrix; // matrix filled with match and mismatchi costs
  float *HH;
  float *EE;
  float *output;

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

  qq= _mm256_set1_ps (gapopen);
  rr= _mm256_set1_ps (gapextend);

  n= strlen(a); // database / columns
  m= strlen(b); // query / rows  
  mm=m;

  if (m%padding) // not 4 or 8 greater common divisor    
  {
    mm = m + (padding - (m % padding)); //new size of the array
    for (int j=m; j<mm;j++)
      a[j]='Z';      //padding letters
  }

  posix_memalign ((void **) &HH, BYTE_ALIGNMENT, mm * sizeof(float));
  posix_memalign ((void **) &EE, BYTE_ALIGNMENT, mm * sizeof(float));
  posix_memalign ((void **) &score_matrix, BYTE_ALIGNMENT,  n * mm * sizeof(float));
  posix_memalign ((void **)&output, BYTE_ALIGNMENT, padding * sizeof(float));

//-----------------------------------------------------------
  int l=0;                        // Filling the score matrix

  for (int i=0; i<n; i++)
  {
    cc= a[i];
    for (int j=0; j<mm; j++)
    {
      c=b[j];
      score_matrix[l] = (c==cc) ? match : mismatch;
      l++;      
    }
  }

//-----------------------------------------------------------

  int y= mm/padding;
  int flag=0; 

  for(int i=0; i<mm; i++)   //initializing the HH and EE matricies
  {
    HH[i]= 1 * gapopen + ( (i+1) * gapextend );
    EE[i]= 2 * gapopen + ( (i+2) * gapextend );
  }

//----------------------------------------------------------
  for (int j=0; j<n;j++)
  {
    if (j==0)
      h= _mm256_setzero_ps();
    else
      h= _mm256_set_ps( 0, 0, 0, 0, 0, 0, 0, (gapopen+ j*gapextend) );

    for (int i=0;i<y;i++)
    { 
      score = _mm256_load_ps(&score_matrix[flag]); //pull the scores 
      flag += padding;

      n_vector= _mm256_load_ps( &HH[ i*padding ] );

      n_vector= shiftl1( n_vector, padding );

//       t1= shiftr7( n_vector, padding );
       
      h= _mm256_or_ps( h, n_vector );
      e= _mm256_load_ps( &EE[ i*padding ] );
      h= _mm256_add_ps (h, score);

      h= _mm256_min_ps (h, e);
      _mm256_store_ps(HH +(padding*i), h );

      print (h, padding);

      h= _mm256_add_ps (h,qq);   // h= h+(q+r)
      h= _mm256_add_ps (h,rr);  

      e= _mm256_add_ps (e,rr);   //e= e+r
      e= _mm256_min_ps (h, e);   //e= min(h,e)
      _mm256_store_ps( EE+(padding*i), e );

//       h= _mm256_set_ps(0,0,0,0,0,0,0,0);
//       h= _mm256_or_ps(h,t1);
      printf("\t");
    }
    printf("\n");
  }

 free(score_matrix);
 free(HH);
 free(EE);
 free(output);
 
 return(0) ;
}

/*
      _mm256_store_ps(output, t1);      //copy t1 into x 
      x= _mm256_load_ps(&output [0]);
      _mm256_store_ps(output, x);                                             //x
      print('x', output);                                                     //x

*/

//==========================================
//----------- shiftl 0xxxxxxx---------------

__m256 shiftl1(__m256 vector, int padding)
{
  float * temp;
  posix_memalign ((void **) &temp, BYTE_ALIGNMENT, padding * sizeof(float));
 
  _mm256_store_ps( temp, vector );

  for(int k=padding; k>0; k--)
    temp[k]=temp[k-1];

  temp[0]=0;

  vector= _mm256_load_ps( &temp[0]);
  free( temp);

  return vector;
}

//------------ shiftr7 x0000000-------------

__m256 shiftr7(__m256 vector, int padding)
{
  float *temp;
  posix_memalign ((void **) &temp, BYTE_ALIGNMENT, padding * sizeof(float));

  _mm256_store_ps(temp, vector);

  temp[0]= temp[7];
 
  for(int k=1;k<padding;k++) 
    temp[k]= 0;
  vector= _mm256_load_ps( &temp[0] ); 
  free (temp);

 return vector;
}

//-----------------------------------------
void print ( __m256 vector, int padding )
{
  float * temp;

  posix_memalign ((void **) &temp, BYTE_ALIGNMENT, padding * sizeof(float));
 
  _mm256_store_ps( temp, vector );

  for(int k=0; k<padding; k++)
    printf( "%2.0f ", temp[k] );

  free( temp );
}
