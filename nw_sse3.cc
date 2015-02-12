//////////////////////////////////////////////////////////////////////////
//  Needleman-Wunsch using (16-bit) integer elements			//
//									//
//  Example   ./nw-sse3 -a ATAGAAGTAG -b TCAGTCAG -m 0 -s 1 -g 1 -e 1	//
//  where the variables refer to: 					//
//  a: database (horisontal)						//
//  b: query (Vertical)							//
//  m: match								//
//  s: mismatch								//
//  g: gapopen								//
//  e: gapextend							//
//////////////////////////////////////////////////////////////////////////

#include<emmintrin.h>
#include <x86intrin.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#define BYTE_ALIGNMENT 16

int padding = BYTE_ALIGNMENT / 2;   //new vecotr size = multiples of padding
int check2(__m128i vector1, __m128i vector2);
void print(__m128i vector);
void convert(char * s);
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
  short int q= 0;   // gap open
  short int r= 0; // gap extend 
  int option=0;
  int n;
  int m;
  int mm;
  int y;
  int score_matrix[4][4];    //DNA letters 
  int flag;

  int16_t *s_matrix; // matrix filled with match and mismatchi costs
  __m128i *HH;
  __m128i *EE;
  __m128i x;
  __m128i E;
  __m128i T1;
  __m128i score;
  __m128i qq;    // gap open vector
  __m128i rr;    // gap extention vector
  __m128i qr;    // q+r
  __m128i H;
  __m128i H0;
  __m128i F0;
  __m128i F;
  __m128i T2;
  __m128i init;
  __m128i padding_mul_r;      // 8*rr
//  int mycounter=0;

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
  printf ("Gapopen:%2d  gapextend:%2d\n\n", q, r);

  qq= _mm_set1_epi16 (q);
  rr= _mm_set1_epi16 (r);
  qr= _mm_add_epi16(qq, rr);
  
  padding_mul_r= _mm_mullo_epi16( _mm_set1_epi16(8), rr);

  n= strlen(a);
  m= strlen(b);
  mm=m;

//  mm = ((m - 1) | (BYTE_ALIGNMENT - 1)) + 1;		//the new zize of the vector - multiples of 16
  mm = ((m - 1) | ( padding - 1)) + 1;       		//the new zize of the vector - multiples of 8

  posix_memalign ((void **) &HH, BYTE_ALIGNMENT, mm * sizeof(__m128i) );
  posix_memalign ((void **) &EE, BYTE_ALIGNMENT, mm * sizeof(__m128i) );
  posix_memalign ((void **) &s_matrix, BYTE_ALIGNMENT, 4 * m * sizeof(int16_t) );

//---------------- Filling the score matrix ----------------

  convert(a);
  convert(b);	

  for(int i=0; i<4; i++)		
    for(int j=0; j<4; j++)
      score_matrix[i][j] = (i==j) ? match : mismatch;

  int l=0;
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<m; j++)
    {
      s_matrix[l] = score_matrix[ i ][ (int)b[j] ];
      l++;
    }
  }

//----------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@ start counting the time @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  gettimeofday(&start_t, NULL);
for(int counter=0; counter<10; counter++)
{
  y= mm/padding;
//================== initializing HH and EE ==================
  init= _mm_set_epi16( 8, 7, 6, 5, 4, 3, 2, 1);       // HH= q+ ( (i+1) * r )

  init= _mm_mullo_epi16(init, rr);
  init= _mm_add_epi16( init, qq);
 
  _mm_store_si128(  HH, init );
  _mm_store_si128(  EE, _mm_add_epi16( init, qq) );     // EE= HH+ QQ

  for(int i=1; i<y; i++)
  {
    init = _mm_load_si128(  &HH[ (i-1) * padding]);
    init = _mm_add_epi16( init, padding_mul_r );    // HH[i-1] + 8*rr
    _mm_store_si128(  HH + (padding*i), init );

    _mm_store_si128( EE + (padding*i), _mm_add_epi16( init, qq)  );	// EE= HH+ QQ
  }

//=============================================================
  for (int j=0; j<n;j++)
  {  
    x= ( (j==0) ? _mm_setzero_si128() : _mm_set_epi16( 0, 0, 0, 0, 0, 0, 0, (q+ j*r) ) );
    H0= _mm_set_epi16( 0, 0, 0, 0, 0, 0, 0,  q+ (j+1)*r  );
    F0= _mm_set_epi16( 0, 0, 0, 0, 0, 0, 0, (2*q) + (j+1)*r );
    F= _mm_setzero_si128();
    flag=0;

    for (int i=0;i<y;i++)
    { 
      score = _mm_load_si128( (__m128i*)&s_matrix[ flag + (y*padding*a[j]) ]); //pull the scores

      H= _mm_load_si128(  &HH[ flag ] );
      E= _mm_load_si128(  &EE[ flag ] );

      T1= _mm_srli_si128(H,14);   //shiftr7(H)

      H= _mm_or_si128( _mm_slli_si128(H,2), x );  //shiftl1(H)
      x= T1;                                      //pading for the H and E part
      H= _mm_add_epi16 (H, score);
      H= _mm_min_epi16 ( H, E);                  // min( H+score, E)

//     print (H);                               //    open to print just H,E 
//    _mm_store_si128( HH +(padding*i), H );     //    open to print just H,E 

//*********************************************

      T2= _mm_add_epi16( _mm_or_si128( _mm_slli_si128(H,2), H0), qr);    // the H+qr vector to be compared with F 

//      if (check2(H, T2) or check2(H, F) )
//      {
        F=T2;  
        do
  	{
	  F=_mm_min_epi16(F,T2);
          T2= _mm_add_epi16( _mm_or_si128( _mm_slli_si128(F,2), F0), rr);
	} while ( check2(F, T2) );

        H= _mm_min_epi16(H,F);
/*      }

      else 
        F=_mm_set_epi16( (2*q) + (j+1)*r, 0, 0, 0, 0, 0, 0, 0);  
*/    
      H0=_mm_srli_si128(H,14);   //shiftr7(H)
      F0=_mm_srli_si128(F,14);   //shiftr7(F)

//      print (H);

      _mm_store_si128( HH +(flag), H );
      E= _mm_min_epi16 ( _mm_add_epi16 (H,qr), _mm_add_epi16 (E,rr) ); //E=min(H+q+r ,E+r)
      _mm_store_si128( EE+(flag), E );
      flag += padding;
    }
//    printf("\n");
  }

}
  gettimeofday(&end_t, NULL);
  printf("Total time taken by CPU: %ld \n", ( (end_t.tv_sec - start_t.tv_sec)*1000000 + end_t.tv_usec - start_t.tv_usec) /10);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  free(HH);
  free(EE);
  return(0) ;
}

//==========================================
//---------------- check2 -------------------
int check2(__m128i vector1, __m128i vector2)
{
  __m128i vcmp =  _mm_cmplt_epi16(vector2, vector1);
  int cmp = _mm_movemask_epi8(vcmp);
  return ((cmp>0) ? 1 : 0) ;
}

//--------------- print --------------------
void print (__m128i vector)
{
  int16_t *val = (int16_t*) &vector;
  printf("%2i %2i %2i %2i %2i %2i %2i %2i ", val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
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

