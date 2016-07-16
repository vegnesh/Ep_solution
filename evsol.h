#ifndef EVSOL_H_INCLUDED
#define EVSOL_H_INCLUDED

/* -------Gives the function Definitions for getting the Ev solution------------------------

*/

extern double PI;

double getk(double Rc, double Zc, double r, double z);

double getG(double Rc, double Zc, double r, double z);

double getG2(double k);


void getdatafromfile (char *filename, double c[15][2], int &nval);

double getKint(double k);

double getEint(double k);

double ellipint(double k,double psi);

double solev(double r,double z,char * filename);
#endif
