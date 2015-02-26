//////////////////////////////////////////////////////////////////////////
//  Needleman-Wunsch using (16-bit) integer elements			//
//									//
//  Example   ./nw-avx2 -a ATAGAAGTAG -b TCAGTCAG -m 0 -s 1 -g 1 -e 1	//
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

int padding = BYTE_ALIGNMENT / 2;   //new vecotr size = multiples of padding
int check2(__m256i vector1, __m256i vector2);
void print(__m256i vector);
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
  short int r= 0;   // gap extend 
  int option=0;
  int n;
  int m;
  int mm;
  int y;
  int score_matrix[4][4];
  int flag;
  int padding_times;

  int16_t *s_matrix; // matrix filled with match and mismatchi costs
  __m256i *HH;   
  __m256i *EE;
  __m256i x;
  __m256i E;
  __m256i T1;
  __m256i score;
  __m256i qq;    // gap open vector
  __m256i rr;    // gap extention vector
  __m256i qr;    // q+r
  __m256i H;
  __m256i H0;
  __m256i F0;
  __m256i F;
  __m256i T2;
  __m256i init;
  __m256i padding_mul_r;      // 8*rr

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

//  printf ("Match:%4d   Mismatch:%d\n", match, mismatch);
//  printf ("Gapopen:%2d  gapextend:%2d\n\n", q, r);

  qq= _mm256_set1_epi16 (q);
  rr= _mm256_set1_epi16 (r);
  qr= _mm256_add_epi16(qq, rr);

  padding_mul_r= _mm256_mullo_epi16( _mm256_set1_epi16(padding), rr);
 
  n= strlen(a);
  m= strlen(b);  

//  mm = ((m - 1) | (BYTE_ALIGNMENT - 1)) + 1;
  mm = ((m - 1) | ( padding - 1)) + 1;       //the new zize of the vector

  posix_memalign ( (void **) &HH, BYTE_ALIGNMENT, mm * sizeof(__m256i) );
  posix_memalign ( (void **) &EE, BYTE_ALIGNMENT, mm * sizeof(__m256i) );
  posix_memalign ((void **) &s_matrix, BYTE_ALIGNMENT, 4 *mm * sizeof(int16_t) );

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
for(int counter=0; counter<20; counter++)
{
  y= mm/padding;
//================== initializing HH and EE ==================
  init= _mm256_set_epi16( 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);       // HH= q+ ( (i+1) * r )
  init= _mm256_mullo_epi16( init, rr);
  init= _mm256_add_epi16( init, qq);

  _mm256_store_si256( HH, init );
  _mm256_store_si256( EE, _mm256_add_epi16( init, qq) );     // EE= HH+ QQ

  padding_times=padding;  //(i*padding)
  for(int i=1; i<y; i++)
  {
    init = _mm256_load_si256( &HH[ padding_times - padding]); 
    init = _mm256_add_epi16( init, padding_mul_r );    // HH[i-1] + 8*rr
    _mm256_store_si256( HH + (padding_times), init );
    _mm256_store_si256( EE + (padding_times), _mm256_add_epi16( init, qq ) );  // EE= HH+ QQ

    padding_times += padding;
  }
//=============================================================
  for (int j=0; j<n;j++)
  {  
    int addmul= q+ (j+1)*r;
    x= ( (j==0) ? _mm256_setzero_si256() : _mm256_set_epi16( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (q+ j*r) ) );
    H0= _mm256_set_epi16( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  addmul  );
    F0= _mm256_set_epi16( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q+addmul );
    F=  _mm256_setzero_si256();
    int index=mm*a[j];
    flag=0;
   
    for (int i=0;i<y;i++)
    { 
      score = _mm256_load_si256( (__m256i*)&s_matrix[ flag+index ]); //pull the scores
      H= _mm256_load_si256( &HH[ flag ] );
      E= _mm256_load_si256( &EE[ flag ] );

      T1= _mm256_alignr_epi8(_mm256_permute2x128_si256(H, H, _MM_SHUFFLE(2, 0, 0, 1)), H, 30);   //shiftr7(H)

      H= _mm256_alignr_epi8(H, _mm256_permute2x128_si256(H, H, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);   //shiftl1(H);
      H= _mm256_or_si256( H, x );	

      x= T1;                                     // pading for the H and E part
      H= _mm256_add_epi16 (H, score);
      H= _mm256_min_epi16 ( H, E);                  // min( H+score, E)
                    
//     print (H);                                   // open to print just H,E 
//    _mm256_store_si256( HH +(flag), H );     // open to print just H,E 

//*********************************************

      T2= _mm256_alignr_epi8(H, _mm256_permute2x128_si256(H, H, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);
      T2= _mm256_add_epi16( _mm256_or_si256( T2, H0), qr);    // the H+qr vector to be compared with F 

//      if (check2(H, T2) or check2(H, F) )
//      {
        F=T2;  
        do
  	{
	  F=_mm256_min_epi16(F,T2);
          T2= _mm256_alignr_epi8(F, _mm256_permute2x128_si256(F, F, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);   //shiftl1(F);
          T2= _mm256_add_epi16( _mm256_or_si256( T2, F0), rr);

	} while ( check2(F, T2) );

        H= _mm256_min_epi16(H,F);
/*      }

      else 
      {
        F= _mm256_set_epi16( q + addmul, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      }
*/
      H0= _mm256_alignr_epi8(_mm256_permute2x128_si256(H, H, _MM_SHUFFLE(2, 0, 0, 1)), H, 30);   //shiftr7(H)
      F0= _mm256_alignr_epi8(_mm256_permute2x128_si256(F, F, _MM_SHUFFLE(2, 0, 0, 1)), F, 30);   //shiftr7(H)

//      print (H);

      _mm256_store_si256( HH +(flag), H );
      E= _mm256_min_epi16 ( _mm256_add_epi16 (H,qr), _mm256_add_epi16 (E,rr) ); //E=min(H+q+r ,E+r)
      _mm256_store_si256( EE+(flag), E );
      flag += padding;
    }
//    printf("\n");
  }
}

  gettimeofday(&end_t, NULL);
  printf("Total time taken by CPU: %ld \n", ( (end_t.tv_sec - start_t.tv_sec)*1000000 + end_t.tv_usec - start_t.tv_usec) /20);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  free(HH);
  free(EE);
  free(s_matrix);
  return(0) ;
}

//==========================================
//---------------- check2 -------------------
int check2(__m256i vector1, __m256i vector2)
{

  __m256i vcmp = _mm256_cmpgt_epi16(vector1, vector2);
  int cmp = _mm256_movemask_epi8(vcmp);
  return cmp !=0;
}

//--------------- print ------------------
void print (__m256i vector)
{
  int16_t *val = (int16_t*) &vector;
  printf("%2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i %2i ", val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7],
  val[8], val[9], val[10], val[11], val[12], val[13], val[14], val[15]);
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

