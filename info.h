#ifndef GUARD_info_h
#define GUARD_info_h

#include<vector>

//data
class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;}
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j)
   double *y; // ith y is *(y+i) or y[i]
};

//prior and mcmc
class pinfo
{
public:
   pinfo() {pbd=1.0;pb=.5;alpha=.95;mybeta=.5;tau=1.0;sigma=1.0;}
//mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
//prior info
   double alpha;
   double mybeta;
   double tau; //mu ~ N(0,tau^2)
   std::vector<double> *mg;  //grid for mu.
   std::vector<double> *pm;  //prior for mu on grid.
//sigma
   double sigma;
};

//sufficient statistics for 1 node
class sinfo
{
public:
   sinfo() {n=0;sy=0.0;sy2=0.0;}
   size_t n;
   double sy;
   double sy2;
};

#endif
