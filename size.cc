#include <x86intrin.h>          										         
#include <emmintrin.h>
#include <unistd.h>
#include <getopt.h> 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


int main(int argc, char *argv[])
{
  FILE* myfile;
  myfile = fopen("a.txt", "r");
  int loc1;
  int size;
//-----------------------------------------
  struct stat st;
  stat("a.txt", &st);
  size = st.st_size;
  std::cout<<"size "<<size<<std::endl; 

  stat("b.txt", &st);
  size = st.st_size;
  std::cout<<"size "<<size<<std::endl;

//-----------------------------------------

  fseek( myfile, 0, SEEK_END ); 
  loc1= ftell(myfile);   
  fseek( myfile, 0, SEEK_SET );
  std::cout<<"size "<< loc1<<std::endl;

  myfile = fopen("b.txt", "r");
  fseek( myfile, 0, SEEK_END ); 
  loc1= ftell(myfile);   
  fseek( myfile, 0, SEEK_SET );
  std::cout<<"size "<< loc1<<std::endl;

  fclose(myfile);
  return(0);
}

