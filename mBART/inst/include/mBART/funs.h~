#include <cmath>
#include <iostream>
#include "tree.h"
#include "info.h"
#include "rn.h" 

using std::cout;
using std::endl;

//pi and log(2*pi)
//extern "C" {
#include <R.h>
#include <Rmath.h>
//};
#ifndef PI
#define PI 3.1415926535897931
#endif
#define LTPI 1.83787706640934536

//--------------------------------------------------
//normal density
double pn(
   double x,    //variate
   double m,    //mean
   double v     //variance
);
//--------------------------------------------------
//draw from a discrete distribution
int rdisc(
   double *p,   //vector of probabilities
   rn& gen     //random number generator
);
//--------------------------------------------------
//evaluate tree tr on grid xi, write to os
void grm(tree& tr, xinfo& xi, std::ostream& os);
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else a/(1+d)^b
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi);
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv);
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
//--------------------------------------------------
//log of the integreted likelihood
double lil(size_t n, double sy, double sy2, double sigma, double tau);
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv);
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, double* fv);
//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv);
//--------------------------------------------------
// draw all the bottom node mu's
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, rn& gen);
//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi);
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
//--------------------------------------------------
//(ml,mu) will be the interval for the mu of the bottom node n, to keep monotonicity
void conint(double& ml, double& mu, tree::tree_p n, tree& t, xinfo& xi);
//--------------------------------------------------
// draw all the bottom node mu's with the monotonicity contstraint
void drmuc(tree& t, xinfo& xi, dinfo& di, pinfo& pi, rn& gen);
