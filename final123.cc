#include <x86intrin.h>          										         
#include <emmintrin.h>
#include <unistd.h>
#include <getopt.h> 
#include <stdlib.h>
#include <stdio.h>

/*
n=2
m=4;
p=2;
a= {0, 2, -1, 5, 7, 6, 1, 3};
b= {2, 1, 5, 3, 1, 4, 1, 2};

n=3
m=4
p=2
a= {0, 2, -1, 5, 7, 6, 1, -2.22, 3.54, 4, 7, 1};
b= {2, 1, 5, 3, 1, 4, 1, 2};

n=3;
m=2;
p=4;

a= {0, 2, -1, 5, 7, 6};
b= {2, 1, 5, 3, 1, 4, 1, 2};

n=7
m=2;
p=5;
a= {0, 2, -1, 5, 7, 6, 8, 3, 4, 8, 6, 7, -1, 0};
b= {2, 1, 5, 3, 1, 4, 1, 2, 4, -2};

n=3;
m=9;
p=2;
a= {0, 2, -1, 5, 7, 6, 8, 3, 4, 8, 6, 7, -1, 0, 5, 6, 1, 3, 4, 2, 0, 7, 9, 0, 5, 7, 1}
b= {2, 1, 5.5, 3, 1, 4, 1, 2, 4, -2, 2, 1, 5, 3, 1, 4, 4, -2};
*/

  int n;  //Number of rows of the first matrix
  int m;  //Number of columns of the first matrix / number of rows of the second matrix
  int p;  //Number of columns of the second matrix

  float *a={0};
  float *b={0};

  float *c= {0};
  float *tmp= {0};      //used for matrix invertion
  float *output= {0};
  float *a4={0};        // replacement arrays
  float *b4={0};

  float *c2= {0};
  float *tmp2= {0};      //used for matrix invertion
  float *output2= {0};
  float *a42={0};        // replacement arrays
  float *b42={0};

  int option = 0, i, j, k, l=0, count, count_max, count2, next=0, mm=1;
  float max=0;

  __m128 va1, vb1, vc1;
  __m256 va2, vb2, vc2;

  int print(float c[], int x, int y);     // function for printing the arrays
  int num_digits(int temp);		  // function to calculate the spasec 
  int inverse(float b[],int x,int y);	  // function to reverse the second array
  int clear(float arr[],int x,int y);
  int fill_array(float arr[], char *file, int x, int y);
  
  int naive();
  int sse();
  int avx();

  char * file_name;

//==============================================
//=================== Main =====================

int main(int argc, char *argv[])
{
  posix_memalign((void **)&a, 128, 256 * sizeof(float));						//memory allignment
  posix_memalign((void **)&b, 128, 256 * sizeof(float));


  posix_memalign((void **)&a4, 128, 128 * sizeof(float));
  posix_memalign((void **)&b4, 128, 128 * sizeof(float));
  posix_memalign((void **)&c, 128, 128 * sizeof(float));
  posix_memalign((void **)&tmp, 128, 128 * sizeof(float));
  posix_memalign((void **)&output, 128, 128 * sizeof(float));

  posix_memalign((void **)&a42, 128, 128 * sizeof(float));
  posix_memalign((void **)&b42, 128, 128 * sizeof(float));
  posix_memalign((void **)&c2, 128, 128 * sizeof(float));
  posix_memalign((void **)&tmp2, 128, 128 * sizeof(float));
  posix_memalign((void **)&output2, 128, 128 * sizeof(float));

  while ((option = getopt(argc, argv,"n:m:p:a:b:")) != -1)
  {
    switch (option) 
    {
      case 'n' : n = atoi(optarg);	break;
      case 'm' : m = atoi(optarg);	break;
      case 'p' : p = atoi(optarg); 	break;
      case 'a' : file_name = optarg; 	               
                 fill_array(a,file_name,n ,m);
                 break;

      case 'b' : file_name= optarg; 	               
                 fill_array(b,file_name,m, p);
                 break;

      default  : printf("use the following parameters -n : -m : -p : -a a.txt -b b.txt\n");
		 exit(EXIT_FAILURE);
    }
  }

  printf("\n======================================================\n");
  print (a, n, m);
  printf("\nx\n");
  print (b, m, p);
  printf("\n=\n");

  naive();
  sse();
  avx();

  return (0);
}

//=================== Fill array ========================

int fill_array(float arr[], char *file_name, int x, int y)
{
  FILE* myfile;
  myfile = fopen(file_name, "r");
  int index=0;
  char character;
  float my_arr;

  if (myfile != NULL)
  {
    while (index<(x*y))
    {
      character = fgetc(myfile);	
      if (character!=',')
      { 
        fscanf(myfile, "%f", &arr[index]);
 	index++;
      }
    }
  }

  fclose(myfile);

  return (0);
}

//============== printing function =============

int print(float c[], int x, int y)		
{	
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
        printf(", ");

      if (c[(i*y)+j]<0)   
	count2++;

      int space=count_max-count2+2;
 
      for (k=0;k<space;k++)  // printing spaces
        printf(" ");

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

//============ Reversing the second array =========
int inverse(float b[],int x,int y)
{
  for (j=0;j<y;j++)
    for (i=0;i<x;i++)	
    {
      tmp[l]=b[(y*i)+j];
      l++;
    }

  for (i=0;i<x*y;i++)
    b[i]=tmp[i];
  return (0);
}

//============== Clear =============

int clear(float arr[],int x,int y)
{
  for (j=0;j<x*y;j++)
    arr[j]=0;

  l=0,next=0, mm=1;
  max=0;

  return (0);
}

//============== Using non-intrinsic =============

int naive() 
{
  for (k=0;k<n;k++) 
    for (i=0;i<p;i++)     
      for (j=0;j<m;j++)
        c[(k*p)+i]= c[(k*p)+i] + ( a[(k*m)+j] * b[(j*p)+i] );

  printf("\n========= naive ========= \n");
  print (c, n, p);
}

//==================================================
//============== Using SSE Instruction =============
int sse()
{
  clear(c,n,p);
  
  if (m%4!=0)   //filling the (<4) arrays
  {
    mm=(m/4)*4+4; // the new size that should be used

    for (i=0;i<m;i++)
      for (j=0;j<n;j++)
        a4[i+(j*mm)]=a[i+(j*m)];	

    for (i=0;i<m;i++)
      for (j=0;j<p;j++)
        b4[(i*p)+j]=b[(i*p)+j];

    inverse(b4,mm,p);
  }

  if (m%4==0) inverse(b,m,p);
  
  for (j=0;j<mm/4;j++)   //to count the parts of the victor
  {
    next=0;
    for (k=0;k<n;k++)
    { 
      if (m%4==0) va1= _mm_load_ps(&a[(mm*k)*(j*4)]);   //to load the vector from the first array 
      if (m%4!=0) va1= _mm_load_ps(&a4[(mm*k)+(j*4)]);

      for (i=0;i<p;i++)    
      { 
        if (m%4==0) vb1= _mm_load_ps(&b[(mm*i)+(j*4)]);  //to load the vector from the second array
        if (m%4!=0) vb1= _mm_load_ps(&b4[(mm*i)+(j*4)]);

        vc1 = _mm_mul_ps(va1,vb1);
        vc1 = _mm_hadd_ps(vc1,vc1);
        vc1 = _mm_hadd_ps(vc1,vc1);

        _mm_store_ps(output, vc1);

        c[next]=c[next]+output[0];//+output[3];

        next++;
      }
    }
  }

  printf("\n========== SSE ========== \n");
  print (c, n, p);
}

//==================================================
//============== Using AVX Instruction =============
int avx()
{
  clear(c,n,p);

  if (m%8!=0)   //filling the (<8) arrays
  {
    mm=(m/8)*8+8; // the new size that should be used
  
    for (j=0;j<n;j++)
      for (i=0;i<m;i++)
        a42[i+(j*mm)]=a[i+(j*m)];	

    for (i=0;i<m;i++)
      for (j=0;j<p;j++)
        b42[(i*p)+j]=b[(i*p)+j];

    inverse(b42,mm,p);
  }

  if (m%8==0)  inverse(b,m,p);  

  for (j=0;j<mm/8;j++)   //to count the parts of the victor
  {
    next=0;
    for (k=0;k<n;k++)
    {
      if (m%8==0) va2= _mm256_load_ps(&a[(mm*k)+(j*8)]);
      if (m%8!=0) va2= _mm256_load_ps(&a42[(mm*k)+(j*8)]);    //to load the vector from the first array
      for (i=0;i<p;i++)
      { 
        if (m%8==0) vb2= _mm256_load_ps(&b[(mm*i)+(j*8)]);   //to load the vector from the second array
        if (m%8!=0) vb2= _mm256_load_ps(&b42[(mm*i)+(j*8)]);  //to load the vector from the second array
        vc2 = _mm256_mul_ps(va2	,vb2);
        vc2 = _mm256_hadd_ps(vc2,vc2);
        vc2 = _mm256_hadd_ps(vc2,vc2);
        _mm256_store_ps(output2, vc2);

        c2[next]=c2[next]+output2[0]+output2[7];

        next++;
      }
    }
  }

  printf("\n========== AVX ========== \n");
  print (c2, n, p);
  printf("\n");
}
