#include <x86intrin.h>
#if (defined(__SSE) || defined(__AVX))
#include <emmintrin.h>
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

  float *a={0};
  float *b={0};

  float *c= {0};
  float *tmp = NULL;      //used for matrix invertion
  float * a4 = NULL;        // replacement arrays
  float * b4 = NULL;

  int i, j, k, count, count_max, count2, next=0;

#if (defined(__SSE))
__m128 va1, vb1, vc1;

#elif (defined(__AVX))
  __m256 va2, vb2, vc2;
#endif

  int print(float * c, int x, int y);     // function for printing the arrays
  int num_digits(int temp);		  // function to calculate the spasec 
  int inverse(float b[],int x,int y);	  // function to reverse the second array
  float * fill_array(char *file, int elms, int alignment);
  
  void naive(int n, int m, int p);
  void sse(int n, int m, int p, int mm);
  void avx(int n, int m, int p, int mm);

//==============================================
//=================== Main =====================

int main(int argc, char *argv[])
{
  int
    n,  //Number of rows of the first matrix
    m,  //Number of columns of the first matrix / number of rows of the second matrix
    p,  //Number of columns of the second matrix
#if (defined(__SSE) || defined(__AVX))
    size_a,
    size_b,
    padding,
    mm,
#endif
    option = 0;
  char * file_a, * file_b;

  if (argc == 1)
   {
     fprintf(stderr, " syntax: %s arguments\n", argv[0]);
     fprintf(stderr, "    Arguments:\n");
     return (0);
   }

  while ((option = getopt(argc, argv,"n:m:p:a:b:")) != -1)
  {
    switch (option) 
    {
      case 'n' : n = atoi(optarg);	break;
      case 'm' : m = atoi(optarg);	break;
      case 'p' : p = atoi(optarg); 	break;
      case 'a' : file_a = optarg; 	               
                 break;

      case 'b' : file_b= optarg; 	               
                 break;

      default  : printf("use the following parameters -n : -m : -p : -a a.txt -b b.txt\n");
		 exit(EXIT_FAILURE);
    }
  }
  
  a = fill_array(file_a, n * m, BYTE_ALIGNMENT);
  b = fill_array(file_b, m * p, BYTE_ALIGNMENT);
  if (!a)
   {
     printf ("Error while reading files!\n");
     exit(0);
   }
#if (defined(__SSE) || defined(__AVX))
  padding = BYTE_ALIGNMENT / 4;
  size_a = n * m;
  size_b = p * m;
  mm = m;
  if (m % padding)
   {
    size_a = n * (m + (BYTE_ALIGNMENT / 4 - (m % padding)));
    size_b = p * (m + (BYTE_ALIGNMENT / 4 - (m % padding)));
    mm     = m + (BYTE_ALIGNMENT / 4 - (m % padding));
   }
  posix_memalign ((void **)&a4, BYTE_ALIGNMENT, size_a * sizeof(float));

  padding = BYTE_ALIGNMENT / 4;
  posix_memalign ((void **)&b4,  BYTE_ALIGNMENT, size_b * sizeof(float)); 
  posix_memalign ((void **)&tmp, BYTE_ALIGNMENT, size_b * sizeof(float)); 
  for (i = 0; i < size_a; ++i)
    a4[i] = 0;
  for (i = 0; i < size_b; ++i)
    b4[i] = 0;
#endif

  posix_memalign((void **)&c, BYTE_ALIGNMENT, n * p * sizeof(float));
  for (i = 0; i < n * p; ++ i) c[i] = 0;

  printf("\n======================================================\n");

  print (a, n, m);
  printf("\nx\n");
  print (b, m, p);
  printf("\n=\n");

#if (!defined(__SSE) && !defined(__AVX))
  naive(n, m, p);
#elif (defined(__SSE))
  sse(n,m,p,mm);
#elif (defined(__AVX))
  avx(n,m,p,mm);
#endif

#if (defined(__SSE) || defined(__AVX))
  free(a4);
  free(b4);
  free(tmp);
#endif
  free(a);
  free(b);
  free(c);

  return (0);
}

//=================== Fill array ========================
float * fill_array(char *file_name, int elms, int alignment)
{

  FILE* fp;
  size_t filesize;
  char * buf;
  char * eptr;
  char ch;
  int nComma = 0;
  float value;
  int i = 0;
  float * arr = NULL;
  
  fp = fopen(file_name, "r");
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  rewind(fp);
   
  buf = (char *) malloc ((filesize + 1)* sizeof (char));
  fread (buf, sizeof(char), filesize, fp);
  buf[filesize] = 0;

  //arr = (float *) calloc (elms, sizeof (float));
  if (posix_memalign ((void **)&arr, alignment, elms * sizeof (float)))
   {
     fprintf (stderr, "Error\n");
     exit (0);
   }

  eptr = buf;
  while ((ch = *eptr++))
   if (ch == ',')  ++nComma;

  if (nComma + 1 != elms) return (NULL);
  
  eptr = buf;
  while (*eptr)
   {
     if ((*eptr >= '0' && *eptr <= '9') || (*eptr == '-')) 
        sscanf(eptr, "%f", &value);
     else
       assert(0);

     arr[i++] = value;
     if (*eptr == '-') eptr++;
     while (*eptr && ((*eptr >= '0' && *eptr <= '9') || (*eptr == '.'))) ++eptr;
     while (*eptr && ( (*eptr != '-') && (!(*eptr >= '0' && *eptr <= '9')) ))
          eptr++;
   }

  free(buf);

  fclose(fp);
  return (arr);
}
//============== printing function =============

int print(float * c, int x, int y)		
{	
  int max = 0; 

  for(i=0;i<x*y;i++)   //Counting the max number to count its length later
  {
    if(c[i]>max) 
      max=c[i];
  }

  count_max=num_digits(max);	//the lenght of the max number

  printf("\n");
  for (i=0;i<x;i++)   
  {	
    for (j=0;j<y;j++)
    {
      count2=num_digits(c[(i*y)+j]);   // The length of the current number 

      if (j>0 and j<=y-1)       
        printf(",  ");

      if (c[(i*y)+j]<0)   
	count2++;

//      int space=count_max-count2+2;
 
//      for (k=0;k<space;k++)  // printing spaces
//        printf(" ");

      printf("%5.3f",c[(i*y)+j]);
    }
    printf("\n");
  }
  return (0);
}

//============== counting the spaces ============== 
int num_digits(int temp)
{
  count=1;
  temp /= 10;
  while(temp > 0)
  {
    count++;
    temp = temp/10;
  }

  return (count); 
}

#if (defined(__AVX) || defined(__SSE))
//============ Reversing the second array =========
int inverse(float * b,int x,int y)
{
  int l = 0;

  for (j=0;j<y;j++)
    for (i=0;i<x;i++)	
    {
      tmp[l]=b[(y*i)+j];
      l++;
    }
  for (i=0;i<x*y;i++)
   {
    b[i]=tmp[i];
   }
  return (0);
}
#endif

//============== Using non-intrinsic =============
void naive(int n, int m, int p) 
{
  for (k=0;k<n;k++) 
  {
    for (i=0;i<p;i++)
     { 
      for (j=0;j<m;j++)
       {
        c[(k*p)+i]= c[(k*p)+i] + ( a[(k*m)+j] * b[(j*p)+i] );
       }
     } 
   }

   print (c, n, p);
}

#if (defined(__SSE))
//============== Using SSE Instruction =============
void sse(int n, int m, int p, int mm)
{
  float *output;

  int next = 0;

  posix_memalign ((void **)&output, BYTE_ALIGNMENT, (BYTE_ALIGNMENT / 4) * sizeof(float));

  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      a4[i+(j*mm)]=a[i+(j*m)];	

  for (i=0;i<m;i++)
    for (j=0;j<p;j++)
      b4[(i*p)+j]=b[(i*p)+j];

  inverse(b4,mm,p);

  for (j=0;j<mm/4;j++)   //to count the parts of the victor
  {
    next=0;
    for (k=0;k<n;k++)
    { 
      va1= _mm_load_ps(&a4[(mm*k)+(j*4)]);

      for (i=0;i<p;i++)    
      { 
        vb1= _mm_load_ps(&b4[(mm*i)+(j*4)]); 
        vc1 = _mm_mul_ps(va1,vb1);
        vc1 = _mm_hadd_ps(vc1,vc1);
        vc1 = _mm_hadd_ps(vc1,vc1);
        _mm_store_ps(output, vc1);

        c[next]=c[next]+output[0];

        next++;
      }
    }
  }

  free(output);
  print (c, n, p);
}
#endif

#if (defined(__AVX))
//============== Using AVX Instruction =============
void avx(int n, int m, int p, int mm)
{
  int next = 0;
  float * output;

  posix_memalign ((void **)&output, BYTE_ALIGNMENT, (BYTE_ALIGNMENT / 4) * sizeof(float));

  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      a4[i+(j*mm)]=a[i+(j*m)];	

  for (i=0;i<m;i++)
    for (j=0;j<p;j++)
      b4[(i*p)+j]=b[(i*p)+j];

  inverse(b4,mm,p);

  for (j=0;j<mm/8;j++)   //to count the parts of the victor
  {
    next=0;
    for (k=0;k<n;k++)
    {
      va2= _mm256_load_ps(&a4[(mm*k)+(j*8)]);    //to load the vector from the first array
      for (i=0;i<p;i++)
      { 
        vb2= _mm256_load_ps(&b4[(mm*i)+(j*8)]);  //to load the vector from the second array
        vc2 = _mm256_mul_ps(va2	,vb2);
        vc2 = _mm256_hadd_ps(vc2,vc2);
        vc2 = _mm256_hadd_ps(vc2,vc2);
        _mm256_store_ps(output, vc2);

        c[next]=c[next]+output[0]+output[7];

        next++;
      }
    }
  }
  
  free(output);
  print (c, n, p);
}
#endif
