#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

#include <fstream>
#include <vector>
#include <limits>

#include <Rcpp.h>

#include "info.h"
#include "tree.h"
#include "funs.h"
#include "bd.h"
#include "rrn.h"

using std::cout;
using std::endl;



RcppExport SEXP cmonbart(
   SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,
   SEXP _ixp,
   SEXP _itau,
   SEXP _inu,
   SEXP _ilambda,
   SEXP _ibase,
   SEXP _ipower,
   SEXP _ind,
   SEXP _iburn,
   SEXP _im,
   SEXP _imgsize,
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery
)
{
   Rprintf("*****Into main of monotonic bart\n");
   //-----------------------------------------------------------
   //random number generation
   GetRNGstate();
   //Rcpp::RNGScope scope;
   //rcpprn gen;
   rrn gen;


   //--------------------------------------------------
   //process args
   Rcpp::NumericMatrix xm(_ix);
   double *x = &xm[0];
   Rcpp::NumericVector yv(_iy);
   double *y = &yv[0];
   Rcpp::NumericMatrix xpm(_ixp);


   size_t p = xm.nrow();
   size_t n = xm.ncol();
   size_t np = xpm.ncol();
   double *xp;
   if(np)  xp = &xpm[0];

   double tau = Rcpp::as<double>(_itau);
   double nu = Rcpp::as<double>(_inu);
   double lambda = Rcpp::as<double>(_ilambda);
   double alpha = Rcpp::as<double>(_ibase);
   double mybeta = Rcpp::as<double>(_ipower);
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   size_t m = Rcpp::as<int>(_im);
   size_t nm = Rcpp::as<int>(_imgsize);

   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);

   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;


   //--------------------------------------------------
   /*
   //write out args
   cout << "p: " << p  << endl;
   cout << "n: " << n  << endl;
   cout << "np: " << np  << endl;
   cout << "y, first, last: "  << y[0] << ", " << y[n-1] << endl;
   cout << "nrows xp: " << xpm.nrow()  << endl;
   cout << "ncols xp: " << xpm.ncol()  << endl;
   cout << "first row x: " << x[0] << ", " << x[1] << endl;
   cout << "last x: " << x[(n-1)*p] << ", " << x[(n-1)*p+1] << endl;
   if(np) {
      cout << "first row xp: " << xp[0] << ", " << xp[1] << endl;
      cout << "last xp: " << xp[(np-1)*p] << ", " << xp[(np-1)*p+1] << endl;
   } else {
      cout << "no test observations\n";
   }
   Rprintf("tau: %lf\n",tau);
   Rprintf("nu: %lf\n",nu);
   Rprintf("lambda: %lf\n",lambda);
   Rprintf("tree prior base: %lf\n",alpha);
   Rprintf("tree prior power: %lf\n",mybeta);
   Rprintf("burn (nskip): %ld\n",burn);
   Rprintf("nd (ndpost): %ld\n",nd);
   Rprintf("m (ntree): %ld\n",m);
   Rprintf("nm (mu grid size): %ld\n",nm);
   */

   Rprintf("**********************\n");
   Rprintf("n: %ld\n",n);
   Rprintf("p: %ld\n",p);
   Rprintf("first and last y: %lf, %lf\n",y[0],y[n-1]);
   Rprintf("first row: %lf, %lf\n",x[0],x[p-1]);
   Rprintf("second row: %lf, %lf\n",x[p],x[2*p-1]);
   Rprintf("last row: %lf, %lf\n",x[(n-1)*p],x[n*p-1]);
   if(np) {
      Rprintf("np: %d\n",np);
      Rprintf("first row xp: %lf, %lf\n",xp[0],xp[p-1]);
      Rprintf("second row xp: %lf, %lf\n",xp[p],xp[2*p-1]);
      Rprintf("last row xp : %lf, %lf\n",xp[(np-1)*p],xp[np*p-1]);
   } else {
      Rprintf("no test observations\n");
   }
   Rprintf("tau: %lf\n",tau);
   Rprintf("nu: %lf\n",nu);
   Rprintf("lambda: %lf\n",lambda);
   //Rprintf("sigest: %lf\n",sigest);
   Rprintf("tree prior base: %lf\n",alpha);
   Rprintf("tree prior power: %lf\n",mybeta);
   Rprintf("burn (nskip): %ld\n",burn);
   Rprintf("nd (ndpost): %ld\n",nd);
   Rprintf("m (ntree): %ld\n",m);
   Rprintf("nm (mu grid size): %ld\n",nm);
   Rprintf("*****nkeeptrain,nkeeptest,nkeeptestme, nkeeptreedraws: %d, %d, %d, %d\n",
               nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
   Rprintf("*****printevery: %d\n",printevery);
   Rprintf("*****skiptr,skipte,skipteme,skiptreedraws: %d,%d,%d,%d\n",
               skiptr,skipte,skipteme,skiptreedraws);
   Rprintf("**********************\n");



   //--------------------------------------------------
   //--------------------------------------------------
   // main code
   //--------------------------------------------------
   //process train data
   double bigval = std::numeric_limits<double>::infinity();
   double miny = bigval; //use range of y to calibrate prior 
   double maxy = -bigval;
   sinfo allys;
   for(size_t i=0;i<n;i++) {
      if(y[i]<miny) miny=y[i];
      if(y[i]>maxy) maxy=y[i];
      allys.sy += y[i]; // sum of y
      allys.sy2 += y[i]*y[i]; // sum of y^2
   }
   allys.n = n;
   double ybar = allys.sy/n; //sample mean
   double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation
   cout << "ybar,shat: " << ybar << ", " << shat <<  endl;
 

   //--------------------------------------------------
   //process test data
   dinfo dip; //data information for prediction
   dip.n=np; dip.p=p; dip.x = &xp[0]; dip.y=0;
   Rprintf("dip.n: %ld\n",dip.n);


   //--------------------------------------------------
   // xinfo
   xinfo xi;
   size_t nc=100; //100 equally spaced cutpoints from min to max.
   makexinfo(p,n,&x[0],xi,nc);
   Rprintf("x1 cuts: %lf ... %lf\n",xi[0][0],xi[0][nc-1]);
   if(p>1) {
      Rprintf("xp cuts: %lf ... %lf\n",xi[p-1][0],xi[p-1][nc-1]);
   }

   //--------------------------------------------------
   //trees
   std::vector<tree> t(m);
   for(size_t i=0;i<m;i++) t[i].setm(ybar/m); 


   //--------------------------------------------------
   //prior and mcmc
   pinfo pi;
   pi.pbd=1.0; //prob of birth/death move
   pi.pb=.5; //prob of birth given  birth/death

   pi.alpha=alpha; //prior prob a bot node splits is alpha/(1+d)^beta
   pi.mybeta=mybeta; 
   pi.tau=tau;
   pi.sigma=shat;

   //***** discrete prior for constained model
   std::vector<double> mg(nm,0.0);  //grid for mu.
   double pridel=3*pi.tau;
   for(size_t i=0;i!=mg.size();i++) mg[i] = -pridel + 2*pridel*(i+1)/(nm+1);
   std::vector<double> pm(nm,0.0);  //prior for mu.

   double sum=0.0;
   for(size_t i=0;i!=mg.size();i++) {
      pm[i] = pn(mg[i],0.0,pi.tau*pi.tau);
      sum += pm[i];
   }
   for(size_t i=0;i!=mg.size();i++)  pm[i] /= sum;
   pi.mg = &mg;
   pi.pm = &pm;

   //--------------------------------------------------
   //dinfo
   double* allfit = new double[n]; //sum of fit of all trees
   for(size_t i=0;i<n;i++) allfit[i]=ybar;
   double* r = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
   double* ftemp = new double[n]; //fit of current tree
   dinfo di;
   di.n=n; di.p=p; di.x = &x[0]; di.y=r; //the y will be the residual

   //--------------------------------------------------
   //storage for ouput
   //in sample fit

   //out of sample fit
   double* ppredmean=0; //posterior mean for prediction
   double* fpredtemp=0; //temporary fit vector to compute prediction
   if(dip.n) {
      ppredmean = new double[dip.n];
      fpredtemp = new double[dip.n];
      for(size_t i=0;i<dip.n;i++) ppredmean[i]=0.0;
   }
   //for sigma draw
   double rss;  //residual sum of squares
   double restemp; //a residual

   //--------------------------------------------------
   //return data structures using Rcpp
   //draws
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
   //means
   Rcpp::NumericVector trmean(n); //train
   for(int i=0;i<n;i++) trmean[i]=0.0;
   Rcpp::NumericVector temean(np);
   for(int i=0;i<np;i++) temean[i]=0.0;
   //trees
   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   //--------------------------------------------------
   //mcmc
   cout << "\nMCMC:\n";
   time_t tp;
   int time1 = time(&tp);
   gen.set_df(n+nu);
   size_t trcnt=0;
   size_t tecnt=0;
   size_t temecnt=0;
   size_t treedrawscnt=0;
   bool keeptest,keeptestme,keeptreedraw;

   for(size_t i=0;i<(nd+burn);i++) {

      if(i%printevery==0) cout << "i: " << i << ", out of " << nd+burn << endl;
      //draw trees
      for(size_t j=0;j<m;j++) {
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) {
            allfit[k] = allfit[k]-ftemp[k];
            r[k] = y[k]-allfit[k];
         }
         bdc(t[j],xi,di,pi,gen);
         drmuc(t[j],xi,di,pi,gen);
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
      }
      //draw sigma
      rss=0.0;
      for(size_t k=0;k<n;k++) {restemp=y[k]-allfit[k]; rss += restemp*restemp;}
      pi.sigma = sqrt((nu*lambda + rss)/gen.chi_square());
      sdraw[i]=pi.sigma;
      if(i>=burn) {
         for(size_t k=0;k<n;k++) trmean[k]+=allfit[k];  

         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            for(size_t k=0;k<n;k++) trdraw(trcnt,k)=allfit[k];
            trcnt+=1;
         }

         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
         if(keeptest || keeptestme) {
            for(size_t j=0;j<dip.n;j++) ppredmean[j]=0.0;
            for(size_t j=0;j<m;j++) {
               fit(t[j],xi,dip,fpredtemp);
               for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
            }
         }
         if(keeptest) {
            for(size_t k=0;k<np;k++) tedraw(tecnt,k)=ppredmean[k];
            tecnt+=1;
         }
         if(keeptestme) {
            for(size_t k=0;k<np;k++) temean[k]+=ppredmean[k];
            temecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
            for(size_t jj=0;jj<m;jj++) treess << t[jj];
            treedrawscnt +=1;
         }
      }
   }
   int time2 = time(&tp);
   cout << "time for loop: " << time2-time1 << endl;
   Rprintf("check counts\n");
   Rprintf("trcnt,tecnt,temecnt,treedrawscnt: %d,%d,%d, %d\n",trcnt,tecnt,temecnt,treedrawscnt);

   for(size_t k=0;k<n;k++) trmean[k]/=nd;
   for(size_t k=0;k<np;k++) temean[k]/=temecnt;

   //--------------------------------------------------
   PutRNGstate();

   //draws of f(x)
   Rcpp::List ret;
   ret["sigma"]=sdraw;
   ret["yhat.train"]=trdraw;
   ret["yhat.test"]=tedraw;
   ret["yhat.train.mean"]=trmean;
   ret["yhat.test.mean"]=temean;

   //trees
   Rcpp::List treesL;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }
   treesL["cutpoints"] = xiret;
   ret["treedraws"] = treesL;

   return ret;
}
