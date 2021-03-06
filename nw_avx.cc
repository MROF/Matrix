//////////////////////////////////////////////////////////////////////////
//  Needleman/Wunsch using (32-bit) floating-point elements		//
//									//
//  Example   ./nw-avx -a ATAGAAGTAG -b TCAGTCAG -m 0 -s 1 -g 1 -e 1	//
//  where the variables refer to: 					//
//  a: database (horisontal)						//
//  b: query (Vertical)							//
//  m: match								//
//  s: mismatch								//
//  g: gapopen								//
//  e: gapextend							//
//////////////////////////////////////////////////////////////////////////

#include <x86intrin.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#define BYTE_ALIGNMENT 32

int padding = BYTE_ALIGNMENT / 4;    //new vecotr size = multiples of padding
int check2(__m256 vector1, __m256 vector2);
void print(__m256 vector);
void convert(char * s);
__m256 shiftr7(__m256 vector);
__m256 shiftl1(__m256 vector);

char map_nt[256] =
  {
    // A=0, C=1, G=2, T=3

    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };
//========================= Main ========================
//=======================================================

int main(int argc, char * argv[] )
{
  printf("====================================================================\n");
  struct timeval start_t, end_t;

  char *a;
  char *b;
  int match;  
  int mismatch;
  float q= 0;   // gap open
  float r= 0; // gap extend 
  int option=0;
  int n;
  int m;
  int mm;
  int y;
//  int flag;
  float *HH;
  float *EE;
  int score_matrix[4][4];

  __m256 x;
  __m256 E;
  __m256 T1;
  __m256 score;
  __m256 qq;    // gap open vector
  __m256 rr;    // gap extention vector
  __m256 qr;    // q+r
  __m256 H;
  __m256 H0;
  __m256 F0;
  __m256 F;
  __m256 T2;
  __m256 init;
  __m256 padding_mul_r;      // 8*rr

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
      case 'g' : q= atoi(optarg);   	  break;
      case 'e' : r= atoi(optarg); 	  break;
      default  : printf(" use the following parameters -a : -b : -m : -s : -g : -e :\n");
		 exit(EXIT_FAILURE);
    }
  }

  printf ("Match:%4d   Mismatch:%d\n", match, mismatch);
  printf ("Gapopen:%2.0f  gapextend:%2.0f\n\n", q, r);

  qq= _mm256_set1_ps (q);

  rr= _mm256_set1_ps (r);
  qr= _mm256_add_ps(qq, rr);
  padding_mul_r= _mm256_mul_ps( _mm256_set_ps(8, 8, 8, 8, 8, 8, 8, 8), rr);
 
  n= strlen(a);
  m= strlen(b);
  mm=m;

//  mm = ((m - 1) | (BYTE_ALIGNMENT - 1)) + 1;		//the new zize of the vector - multiples of 16
  mm = ((m - 1) | ( padding - 1)) + 1;       		//the new zize of the vector - multiples of 8

  posix_memalign ((void **) &HH, BYTE_ALIGNMENT, mm * sizeof(float));
  posix_memalign ((void **) &EE, BYTE_ALIGNMENT, mm * sizeof(float));

//---------------- Filling the score matrix ----------------

  convert(a);
  convert(b);	

  for(int i=0; i<4; i++)		
    for(int j=0; j<4; j++)
      score_matrix[i][j] = (i==j) ? match : mismatch;
//----------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@ start counting the time @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  gettimeofday(&start_t, NULL);
for(int counter=0; counter<1; counter++)
{
  y= mm/padding;

//================== initializing HH and EE ==================
  init= _mm256_set_ps( 8, 7, 6, 5, 4, 3, 2, 1);       // HH= q+ ( (i+1) * r )
  init= _mm256_mul_ps( init, rr);

  init= _mm256_add_ps( init, qq);
  _mm256_store_ps(HH, init );
  _mm256_store_ps(EE, _mm256_add_ps( init, qq) );     // EE= HH+ QQ

  for(int i=1; i<y; i++)
  {
    init = _mm256_load_ps(&HH[ (i-1) * padding]); 
    init = _mm256_add_ps( init, padding_mul_r );    // HH[i-1] + 8*rr
    _mm256_store_ps(HH + (padding*i), init );

    init = _mm256_load_ps(&EE[ (i-1) * padding]);
    init = _mm256_add_ps( init, padding_mul_r );    // EE[i-1] + 8*rr
    _mm256_store_ps(EE + (padding*i), init );
  }

//=============================================================

  for (int j=0; j<n;j++)
  {  
    x= ( (j==0) ? _mm256_setzero_ps() : _mm256_set_ps( 0, 0, 0, 0, 0, 0, 0, (q+ j*r) ) );
    H0= _mm256_set_ps( 0, 0, 0, 0, 0, 0, 0,  q+ (j+1)*r  );
    F0= _mm256_set_ps( 0, 0, 0, 0, 0, 0, 0, (2*q) + (j+1)*r );
    F= _mm256_setzero_ps();

//    flag=0;
    for (int i=0;i<y;i++)
    { 
/*
      score = _mm256_set_ps( score_matrix[ int(a[j]) ][ int(b[flag+7]) ], 
 			     score_matrix[ int(a[j]) ][ int(b[flag+6]) ], 
			     score_matrix[ int(a[j]) ][ int(b[flag+5]) ], 
			     score_matrix[ int(a[j]) ][ int(b[flag+4]) ], 
			     score_matrix[ int(a[j]) ][ int(b[flag+3]) ], 
			     score_matrix[ int(a[j]) ][ int(b[flag+2]) ], 
			     score_matrix[ int(a[j]) ][ int(b[flag+1]) ], 
			     score_matrix[ int(a[j]) ][ int(b[flag+0]) ] );
      flag += padding;
*/
      score = _mm256_set_ps( score_matrix[ int(a[j]) ][ int(b[padding*i +7]) ],
 			     score_matrix[ int(a[j]) ][ int(b[padding*i +6]) ], 
			     score_matrix[ int(a[j]) ][ int(b[padding*i +5]) ], 
			     score_matrix[ int(a[j]) ][ int(b[padding*i +4]) ], 
			     score_matrix[ int(a[j]) ][ int(b[padding*i +3]) ], 
			     score_matrix[ int(a[j]) ][ int(b[padding*i +2]) ], 
			     score_matrix[ int(a[j]) ][ int(b[padding*i +1]) ], 
			     score_matrix[ int(a[j]) ][ int(b[padding*i +0]) ] );

      H= _mm256_load_ps( &HH[ i*padding ] );
      E= _mm256_load_ps( &EE[ i*padding ] );

      T1= shiftr7( H );

      H= _mm256_or_ps( shiftl1(H), x );
      x= T1;                                     //pading for the H and E part

      H= _mm256_add_ps (H, score);
      H= _mm256_min_ps ( H, E);                  // min( H+score, E)

//     print (H);                               //    open to print just H,E 
//    _mm256_store_ps(HH +(padding*i), H );     //    open to print just H,E 

//*********************************************

      T2= _mm256_add_ps( _mm256_or_ps(shiftl1(H), H0), qr);    // the H+qr vector to be compared with F 
      if (check2(H, T2) or check2(H, F) )
      {
        F=T2;  
        do
  	{
	  F=_mm256_min_ps(F,T2);
          T2= _mm256_add_ps( _mm256_or_ps(shiftl1(F), F0), rr);
	} while ( check2(F, T2) );

        H= _mm256_min_ps(H,F);
      }

      else 
        F=_mm256_set_ps( (2*q) + (j+1)*r, 0, 0, 0, 0, 0, 0, 0);

      H0=shiftr7(H);
      F0=shiftr7(F);

      print (H);
      _mm256_store_ps( HH +(padding*i), H );

      E= _mm256_min_ps ( _mm256_add_ps (H,qr), _mm256_add_ps (E,rr) ); //E=min(H+q+r ,E+r)
      _mm256_store_ps( EE+(padding*i), E );
    }
    printf("\n");
  }
}
  gettimeofday(&end_t, NULL);
  printf("Total time taken by CPU: %ld \n", ((end_t.tv_sec - start_t.tv_sec)* 1000000 + end_t.tv_usec - start_t.tv_usec));

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  free(HH);
  free(EE);
  return(0) ;
}

//==========================================
//---------------- check2 -------------------
int check2(__m256 vector1, __m256 vector2)
{
  __m256 vcmp = _mm256_cmp_ps(vector2, vector1, _CMP_LT_OQ);
  int cmp = _mm256_movemask_ps(vcmp);
  return ((cmp>0) ? 1 : 0) ;
}

//----------- shiftl 0xxxxxxx---------------
__m256 shiftl1(__m256 vector)
{
  __m256 u = _mm256_permute_ps(vector, 0x93);
  __m256 v = _mm256_permute2f128_ps(u, u, 0x08);
  __m256 ans  = _mm256_blend_ps(u, v, 0x11);
  return ans;
}

//------------ shiftr7 x0000000-------------
__m256 shiftr7(__m256 vector)
{
  __m256 u = _mm256_permute_ps(vector, 0x27);
  __m256 v = _mm256_permute2f128_ps(u, u, 0x81);
  __m256 ans  = _mm256_blend_ps(_mm256_setzero_ps(), v, 0x01);
  return ans;
}

//--------------- print --------------------
void print ( __m256 vector)
{
  float * temp;
  posix_memalign ((void **) &temp, BYTE_ALIGNMENT, padding * sizeof(float));
 
  _mm256_store_ps( temp, vector );

  for(int k=0; k<padding; k++)
    printf( "%2.0f ", temp[k] );

  free( temp );
}

//--------------- convert ------------------
void convert(char * s)
{
  char c;
  char m;
  char * p = s;
  int i = 0;

  while ((c = *p++))
   {
     if ((m = map_nt[(int)c]) >= 0)
       *(s + i++) = m;
     else
     {        
       printf("Illegal character in sequence.\n");
       exit(EXIT_FAILURE);
     }
   }
}

