#ifndef GUARD_funs_h
#define GUARD_funs_h

/*
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
*/

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

/*
#include <cmath>
#include "funs.h"
#include <map>
#include <limits>
*/

//--------------------------------------------------
// normal density N(x, mean, variance)
double pn(double x, double m, double v)
{
   double dif = x-m;
   return exp(-.5*dif*dif/v)/sqrt(2*PI*v);
}
//--------------------------------------------------
// draw from discrete distributin given by p, return index
int rdisc(double *p, rn& gen)
{

   double sum;
   double u = gen.uniform();

    int i=0;
    sum=p[0];
    while(sum<u) {
       i += 1;
       sum += p[i];
    }
    return i;
}
//--------------------------------------------------
//evalute tree tr on grid given by xi and write to os
void grm(tree& tr, xinfo& xi, std::ostream& os) 
{
   size_t p = xi.size();
   if(p!=2) {
      cout << "error in grm, p !=2\n";
      return;
   }
   size_t n1 = xi[0].size();
   size_t n2 = xi[1].size();
   tree::tree_cp bp; //pointer to bottom node
   double *x = new double[2];
   for(size_t i=0;i!=n1;i++) {
      for(size_t j=0;j!=n2;j++) {
         x[0] = xi[0][i]; 
         x[1] = xi[1][j]; 
         bp = tr.bn(x,xi);
         os << x[0] << " " << x[1] << " " << bp->getm() << " " << bp->nid() << endl;
      }
   }
   delete[] x;
}
//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
   int L,U;
   bool v_found = false; //have you found a variable you can split on
   size_t v=0;
   while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) v_found=true;
      v++;
   }
   return v_found;
}
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots)
{
   double pb;  //prob of birth to be returned
   tree::npv bnv; //all the bottom nodes
   t.getbots(bnv);
   for(size_t i=0;i!=bnv.size();i++) 
      if(cansplit(bnv[i],xi)) goodbots.push_back(bnv[i]);
   if(goodbots.size()==0) { //are there any bottom nodes you can split on?
      pb=0.0;
   } else { 
      if(t.treesize()==1) pb=1.0; //is there just one node?
      else pb=pi.pb;
   }
   return pb;
}
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
   int L,U;
   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi)
{
   if(cansplit(n,xi)) {
      return pi.alpha/pow(1.0+n->depth(),pi.mybeta);
   } else {
      return 0.0;
   }
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y

   bnv.clear();
   x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   sv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      y=di.y[i];

      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      ++(sv[ni].n);
      sv[ni].sy += y;
      sv[ni].sy2 += y*y;
   }
}
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         y = di.y[i];
         if(xx[v] < xi[v][c]) {
               sl.n++;
               sl.sy += y;
               sl.sy2 += y*y;
          } else {
               sr.n++;
               sr.sy += y;
               sr.sy2 += y*y;
          }
      }
   }
}
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      if(bn==nl) {
         y = di.y[i];
         sl.n++;
         sl.sy += y;
         sl.sy2 += y*y;
      }
      if(bn==nr) {
         y = di.y[i];
         sr.n++;
         sr.sy += y;
         sr.sy2 += y*y;
      }
   }
}
//--------------------------------------------------
//log of the integrated likelihood
double lil(size_t n, double sy, double sy2, double sigma, double tau)
{
   double yb,yb2,S,sig2,d;
   double sum, rv;

   yb = sy/n;
   yb2 = yb*yb;
   S = sy2 - (n*yb2);
   sig2 = sigma*sigma;
   d = n*tau*tau + sig2;
   sum = S/sig2 + (n*yb2)/d;
   rv = -(n*LTPI/2.0) - (n-1)*log(sigma) -log(d)/2.0;
   rv = rv -sum/2.0;
   return rv;
}
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv)
{
   double *xx;
   tree::tree_cp bn;
   fv.resize(di.n);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      fv[i] = bn->getm();
   }
}
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, double* fv)
{
   double *xx;
   tree::tree_cp bn;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      fv[i] = bn->getm();
   }
}
//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv)
{
   double *xx;
   tree::tree_cp bn;
   pv.resize(di.n);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      pv[i] = bn->nid();
   }
}
//--------------------------------------------------
// draw all the bottom node mu's
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, rn& gen)
{
   tree::npv bnv;
   std::vector<sinfo> sv;
   allsuff(t,xi,di,bnv,sv);

   double a = 1.0/(pi.tau * pi.tau);
   double sig2 = pi.sigma * pi.sigma;
   double b,ybar;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {
      b = sv[i].n/sig2;
      ybar = sv[i].sy/sv[i].n;
      bnv[i]->setm(b*ybar/(a+b) + gen.normal()/sqrt(a+b));
   }
}
//--------------------------------------------------
// draw all the bottom node mu's with the monotonicity contstraint
void drmuc(tree& t, xinfo& xi, dinfo& di, pinfo& pi, rn& gen)
{
   double bigval = std::numeric_limits<double>::infinity();

   tree::npv bnv; //bottom nodes
   std::vector<sinfo> sv; //sufficient stats for bottom nodes
   allsuff(t,xi,di,bnv,sv); //fill bnv and sv

   size_t nm = (pi.mg)->size(); //size of unconstrained mu grid
   double mla=(*pi.mg)[0]; //unrestricted grid for mu in [mla,mua]
   double mua=(*pi.mg)[nm-1];
   double ml; //restricted mu vals will be [ml,mu] in each bottom node 
   double mu;

   std::vector<double> pv;  //probs for mu vals
   std::vector<double> muv; //mu vals satisfying the constraint
   std::vector<double> priv; //prior at unrestricted mu vals
   double tmu; //current mu value
   size_t cnt; //count how many good mu values there are
   double temp;
   double lMax;
   size_t ii;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {
      ml=mla;mu=mua;
      conint(ml,mu,bnv[i],t,xi);
      cnt=0;lMax=-bigval;
      pv.clear(); muv.clear(); priv.clear();
      for(size_t j=0;j!=nm;j++) {
         tmu=(*pi.mg)[j];
         if((tmu >= ml) && (tmu<=mu)) {
            cnt++;
            temp = -.5*(-2.0*tmu*(sv[i].sy) + (sv[i].n)*tmu*tmu)/(pi.sigma*pi.sigma);
            if(temp>lMax) lMax=temp;
            pv.push_back(temp);
            muv.push_back(tmu);
            priv.push_back((*pi.pm)[j]);
         }
      }
      double il = 0.0;
      for(size_t j=0;j!=pv.size();j++) {
         pv[j] = priv[j]*exp(pv[j]-lMax);
         il += pv[j];
      }
      for(size_t j=0;j!=pv.size();j++) pv[j] /= il;
      ii = rdisc(&pv[0],gen);
      bnv[i]->setm(muv[ii]);
   }
}
//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
   cout << "xinfo: \n";
   for(size_t v=0;v!=xi.size();v++) {
      cout << "v: " << v << ", ncuts: " << xi[v].size() << endl;
      for(size_t j=0;j!=xi[v].size();j++) cout << "j,xi[v][j]: " << j << ", " << xi[v][j] << endl;
   }
   cout << "\n\n";
}
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
   double xinc;
   double ival = std::numeric_limits<double>::infinity();

   //compute min and max for each x
   std::vector<double> minx(p,ival);
   std::vector<double> maxx(p,-ival);
   double xx;
   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }
   }
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      xinc = (maxx[i]-minx[i])/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
   }
}
//--------------------------------------------------
//(ml,mu) will be the interval for the mu of the bottom node n, to keep monotonicity
// loop over bottom nodes:
//     if a node is below mu must be bigger than its
//     if a node is above mu must be smaller than its
void conint(double& ml, double& mu, tree::tree_p n, tree& t, xinfo& xi)
{
   double nmu; //mu for a bottom node
   char s;
   tree::npv bn;  //bottom nodes
   t.getbots(bn);
   for(size_t i=0;i!=bn.size();i++) {
      nmu = bn[i]->getm();
      if(!(bn[i]==n)) {
         s = n->nhb(bn[i],xi);
         if((s=='b') && (nmu>ml)) ml = nmu;
         if((s=='a') && (nmu<mu)) mu = nmu;
      }
   }
}

#endif
