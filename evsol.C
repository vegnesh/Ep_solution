#include <iostream>
#include "evsol.h"
#include <math.h>
#include <fstream>
#include <cstring>
// Contains the function definitons required for the solution of Ev

using namespace std;
 
double PI=atan(1.0)*4.0; // The Value of pi





/*-----------------NOTE -------------------------------
We use Simpsons Rule for Integration. For Simpson Rule, the total number of points must be odd.
*/

//PI = atan(1.0)*4.0;

/*-----------------------------------------------------------------------------------------*/
	// This function returns the integrand of the elliptic integrals
/*-----------------------------------------------------------------------------------------*/

double ellipint(double k, double psi)
{

  return sqrt(1-pow(k*sin(psi),2.0));

}

/*-----------------------------------------------------------------------------------------*/
//			Gets the 'k' Value required for elliptic Integral			
/*-----------------------------------------------------------------------------------------*/

double getk(double Rc, double Zc, double r, double z)
{

      return sqrt(4*Rc*r/(pow(Rc+r,2.0)+pow(z-Zc,2.0)));
}


/*-----------------------------------------------------------------------------------------*/
//		Returns the value of G(k)
/*-----------------------------------------------------------------------------------------*/
double getG(double Rc, double Zc, double r, double z)
{

	double k = getk(Rc,Zc,r,z);
        if(k>-1e-6 && k<1e-6)
            return 0.0;
        else 
	return getG2(k);
}

double getG2(double k)
{
return ((2-k*k)*getKint(k) - 2*getEint(k))/k;
}

/*-----------------------------------------------------------------------------------------*/

// --- Following function evaluates the K(k) elliptic integral --------------------------------

/*-----------------------------------------------------------------------------------------*/

double getKint(double k)
{
	
	int npoints = 10001; // NUmber of points used for quadrature

	double minval=0.0, maxval = PI/2.0;
	double dx = (maxval-minval)/(npoints-1);


	double i=minval;
	int index = 0;
	double sum;

	while(i<maxval) // Did not put for equal checking conditiion because of double values being used.
	{
	  if (index ==0)
	  sum = 1.0/ellipint(k,i);
	  else if(index % 2 == 0)
	  {
	  sum += 2.0/ ellipint(k,i);	    
	  }
	  else if(index %2 !=0)
	  sum += 4.0/ ellipint(k,i);

	  i+= dx;
	  ++index;

	}

	sum += 1.0/ellipint(k,i); // Account for the last number
	sum *= dx/3.0;

	return sum;
}

/*-----------------------------------------------------------------------------------------*/

// ---------------- Following Function evaluates the E elliptic Integral -------------------

/*-----------------------------------------------------------------------------------------*/
double getEint(double k)
{
	
	int npoints = 10001; // NUmber of points used for quadrature

	double minval=0.0, maxval = PI/2.0;
	double dx = (maxval-minval)/(npoints-1);

	double i=minval;
	int index = 0;
	double sum;

	while(i<maxval) // Did not put for equal checking conditiion because of double values being used.
	{
	  if (index ==0)
	  sum = 1.0*ellipint(k,i);
	  else if(index % 2 == 0)
	  {
	  sum += 2.0* ellipint(k,i);	    
	  }
	  else if(index %2 !=0)
	  sum += 4.0* ellipint(k,i);

	  i+= dx;
	  ++index;

	}

	sum += 1.0*ellipint(k,i); // Account for the last number
	sum *= dx/3.0;

	return sum;
}


/*------------------------------------------i-----------------------------------------------*/
//Give the value of the EV solution for the given coordinates considering all the points 
//where the current is flowing through
/*------------------------------------------i-----------------------------------------------*/

void getdatafromfile (char *filename, double c[15][2], int &nval)
{
/* 	 Data Variables	    */

    int noc = 0;
    bool flag = false;
    double coords[15][2];
   nval = 0;
/*	 Reading from a file to get data	*/
   if(!flag)
{
    flag = !flag;

/* File Part */
char line[100];
int pos=0,flag=0;
  ifstream myfile (filename);
    if (myfile.is_open())
      {
          while ( myfile )
      { 
          pos = myfile.tellg();     
          myfile.getline(line,100);
          if (line[0]!='#' )
          if (strlen(line)!=0)
         {
           myfile.seekg(pos);
           myfile>>coords[noc][0]>>coords[noc][1];
	   ++noc;   
	   cout<<"\n R="<<coords[noc-1][0]<<" Z="<<coords[noc-1][1]<<"\n";
         }

    }

       myfile.close();
      for(int i=0;i<noc;i++)
{
coords[i][0]=coords[i][0]/1000.0;
coords[i][1]=coords[i][1]/1000.0;
}    
      
           

}

	     else cout << "Unable to open file";
    nval = noc;
    for(int j=0;j<nval;j++)
{
        c[j][0] = coords[j][0];
       c[j][1] = coords[j][1];
}             
/* File Part Ends */
}

}




double solev (double r,double z,char *filename)
{
    

    if (r==0.0)
    return 0.0;  

   double result = 0.0;
   double coords[15][2];
   int ncoords;
   double k = 0.0;
   getdatafromfile (filename,coords,ncoords);

   for(int i = 0;i<ncoords; i++)
{
    //cout<<"\n R="<<coords[i][0]<<" Z="<<coords[i][1]<<"\n";
    k = getk(coords[i][0],coords[i][1],r,z);
    result =result + sqrt(coords[i][0]/r)*getG2(k);
    //std::cout<<"\nk:"<<k<<"\n";
}
   
return result;

}


