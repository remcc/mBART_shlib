#ifndef RRN_H
#define RRN_H

#ifdef Rcpp_hpp
using R::rchisq;
#else
extern "C" {
#include <R.h>
#include <Rmath.h>
};
#endif

#include "rn.h"

class rrn: public rn
{
public:
//constructor
   rrn():df(1) {}
//virtual
   virtual ~rrn() {}
   virtual double normal() {return norm_rand();}
   virtual double uniform() { return unif_rand();}
   virtual double chi_square() {return rchisq((double)df);}
   virtual double exp() {return exp_rand();}
   virtual void set_df(int df) {this->df=df;}
//get,set
   int get_df() {return df;}
private:
   int df;
};

#endif
