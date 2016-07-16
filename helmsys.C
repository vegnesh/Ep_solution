#include "helmsys.h"
#include <stdio.h>
#include <string.h>

void getfilename(char* filename)
{
  char files[] = "coordinates_custom.dat";
  strcpy(filename,files);
}

void getcoords(double coords[2][2])
{
  coords[0][0] = 0.0  ;
  coords[0][1] = 0.5  ;
  coords[1][0] =-1.0 ;
  coords[1][1] = 1.0 ;

}
void getdomain(double coords[2][2])
{
  coords[0][0] = 0.0  ;
  coords[0][1] = 0.08  ;
  coords[1][0] = 0.0  ;
  coords[1][1] = 0.486 ;

}


void getnnodes(int &n)
{
  n = 15;
}
