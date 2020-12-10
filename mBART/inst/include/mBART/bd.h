#ifndef GUARD_bd_h
#define GUARD_bd_h

/*
#include "rn.h"
#include "info.h"
#include "tree.h"
*/

bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, rn& gen);
bool bdc(tree& x, xinfo& xi, dinfo& di, pinfo& pi, rn& gen);

/*
#include <iostream>
#include <limits>

#include "info.h"
#include "tree.h"
#include "bd.h"
#include "funs.h"

using std::cout;
using std::endl;
*/

/*
notation: (as in old code): going from state x to state y (eg, incoming tree is x).

note: rather than have x and making a tree y
we just figure out what we need from x, the drawn bottom node,the drawn (v,c).
note sure what the right thing to do is.
Could make y (using a birth) and figure stuff out from y.
That is how the old code works.
*/

/*
unconstrained is bd.
constrained is bdc.
*/

/*
note: the difference between coninteg1 and coninteg
is just numeric stability (I think).
*/


double coninteg1(sinfo& sl, sinfo& sr, xinfo& xi, pinfo& pi, tree x, tree::tree_p& nx, std::vector<double>& pv, std::vector<double>& muv)
{
   pv.clear();
   muv.clear();

   double syall = sl.sy + sr.sy;
   size_t nall = sl.n + sr.n;
   //cout << "syall,nall: " << syall << ", " << nall << endl;

   size_t nm = (pi.mg)->size();
   //cout << "nm: " << nm << endl;

   double mla=(*pi.mg)[0];
   double mua=(*pi.mg)[nm-1];

   //cout << "mla,mua: " << mla << ", " << mua << endl;

   double ml=mla;
   double mu=mua;
   conint(ml,mu,nx,x,xi);
   //cout << "ml,mu: " << ml << ", " << mu << endl;

   double tmu; //temp mu
   double iltc = 0.0;
   size_t cnt=0;
   double temp;
   double sumpri=0.0;
   double tpri;
   for(size_t i=0;i!=nm;i++) {
      tmu = (*pi.mg)[i];
      if((tmu >= ml) && (tmu<=mu)) {
         /*
         cnt++;
         temp = (*pi.pm)[i];
         temp *= temp*exp(-.5*(-2.0*tmu*syall + nall*tmu*tmu)/(pi.sigma*pi.sigma));
         iltc += temp;
         */
         cnt++;
         tpri = (*pi.pm)[i];
         sumpri += tpri;
         temp = tpri*exp(-.5*(-2.0*tmu*syall + nall*tmu*tmu)/(pi.sigma*pi.sigma));
         iltc += temp;
         pv.push_back(temp);
         muv.push_back(tmu);
      }
   }
   //cout << "in coninteg1, iltc: " << iltc << endl;
   for(size_t i=0;i!=pv.size();i++) {
      pv[i] /= iltc;
   //   cout << "in coninteg1, i, pv: " << i << ", " << pv[i] << ", " << muv[i] << endl;
   }
   //cout << "\ncnt in coninteg1: " << cnt << endl;
   return iltc;
}
double coninteg1(sinfo& sl, sinfo& sr, xinfo& xi, pinfo& pi, tree x, tree::tree_p& nx, std::vector<double>& pv, std::vector<double>& muv, double& lcfac)
{
   //cout << "into new coninteg1\n";

   double bigval = std::numeric_limits<double>::infinity();
   pv.clear();
   muv.clear();

   double syall = sl.sy + sr.sy;
   size_t nall = sl.n + sr.n;
   //cout << "syall,nall: " << syall << ", " << nall << endl;

   size_t nm = (pi.mg)->size();
   //cout << "nm: " << nm << endl;

   double mla=(*pi.mg)[0];
   double mua=(*pi.mg)[nm-1];

   //cout << "mla,mua: " << mla << ", " << mua << endl;

   double ml=mla;
   double mu=mua;
   conint(ml,mu,nx,x,xi);
   //cout << "ml,mu: " << ml << ", " << mu << endl;

   double tmu; //temp mu
   double iltc = 0.0;
   size_t cnt=0;
   double temp;
   double tpri;
   double lMax = -bigval;
   std::vector<double> priv;
   for(size_t i=0;i!=nm;i++) {
      tmu = (*pi.mg)[i];
      if((tmu >= ml) && (tmu<=mu)) {
         cnt++;
         tpri = (*pi.pm)[i];
         //temp = tpri*exp(-.5*(-2.0*tmu*syall + nall*tmu*tmu)/(pi.sigma*pi.sigma));
         temp = -.5*(-2.0*tmu*syall + nall*tmu*tmu)/(pi.sigma*pi.sigma);
         if(temp>lMax) lMax=temp;
         pv.push_back(temp);
         muv.push_back(tmu);
         priv.push_back(tpri);
      }
   }

   //cout << "lMax: " << lMax << endl;
   //for(size_t i=0;i!=pv.size();i++) {
    //  cout << "in coninteg1, i, pv, priv,muv: " << i << ", " << pv[i] << ", " << priv[i] << ", "  << muv[i] << endl;
   //}

   iltc=0.0;
   for(size_t i=0;i!=pv.size();i++) {
      pv[i] = priv[i]*exp(pv[i]-lMax);
      iltc += pv[i];
   }
   for(size_t i=0;i!=pv.size();i++) pv[i] /= iltc;
   lcfac=lMax;

   //cout << "lMax: " << lMax << endl;
   //for(size_t i=0;i!=pv.size();i++) {
    //  cout << "in coninteg1, i, pv, priv,muv: " << i << ", " << pv[i] << ", " << priv[i] << ", "  << muv[i] << endl;
   //}

   return iltc;
}


double coninteg2(sinfo& sl, sinfo& sr, xinfo& xi, pinfo& pi, tree y, size_t nxid, std::vector<double>& pv, std::vector<double>& mulv, std::vector<double>& murv)
//nxid is the nog node
{
   pv.clear();
   mulv.clear();
   murv.clear();

   size_t nm = (pi.mg)->size();
   double mla=(*pi.mg)[0];
   double mua=(*pi.mg)[nm-1];
   double ml;
   double mu;

   //get pointers to the nog node and its children
   tree::tree_p ny= y.getptr(nxid);
   tree::tree_p nyl = ny->getl();
   tree::tree_p nyr = ny->getr();

   double tmul,tmur;
   double illrc=0.0;
   int cnt=0;
   bool lg,rg;  //are right and left mu's good
   double tpril,tprir;
   double sumpr=0.0;
   double temp;
   for(size_t i=0;i!=nm;i++) { //i indexes mul
      tmul = (*pi.mg)[i];
      tpril = (*pi.pm)[i];
      for(size_t j=i;j<nm;j++) { //j indexes mur
         tmur = (*pi.mg)[j];
         lg=false; rg=false;
         //fix right, check left
         nyr->setm(tmur);
         ml=mla; mu=mua;
         conint(ml,mu,nyl,y,xi);
         if((tmul >= ml) && (tmul<=mu)) lg=true; 
         //fix left, check right
         nyl->setm(tmul);
         ml=mla; mu=mua;
         conint(ml,mu,nyr,y,xi);
         if((tmur >= ml) && (tmur<=mu)) rg=true; 
         if(lg && rg) {
            cnt++;
            tprir = (*pi.pm)[j];
            temp = tpril*tprir;
            sumpr += temp;
            temp *= exp(-.5*(-2.0*tmul*sl.sy + sl.n*tmul*tmul)/(pi.sigma*pi.sigma)); 
            temp *= exp(-.5*(-2.0*tmur*sr.sy + sr.n*tmur*tmur)/(pi.sigma*pi.sigma)); 
            illrc += temp;
            pv.push_back(temp);
            mulv.push_back(tmul);
            murv.push_back(tmur);
         }
      }
   }
   for(size_t i=0;i!=pv.size();i++) {
      pv[i] /= illrc;
      //cout << "i,pv,mul,mur: " << i << ", " << pv[i] << ", " << mulv[i] << ", " << murv[i] << endl;
   }
   //cout << "cnt,sumpr: " << cnt << ", " << sumpr << endl;
   //cout << "illrc: " << illrc << endl;
   return illrc;
}
double coninteg2(sinfo& sl, sinfo& sr, xinfo& xi, pinfo& pi, tree y, size_t nxid, std::vector<double>& pv, std::vector<double>& mulv, std::vector<double>& murv, double& lcfac)
//nxid is the nog node
{
   double bigval = std::numeric_limits<double>::infinity();

   pv.clear();
   mulv.clear();
   murv.clear();

   size_t nm = (pi.mg)->size();
   double mla=(*pi.mg)[0];
   double mua=(*pi.mg)[nm-1];
   double ml;
   double mu;

   //get pointers to the nog node and its children
   tree::tree_p ny= y.getptr(nxid);
   tree::tree_p nyl = ny->getl();
   tree::tree_p nyr = ny->getr();

   double tmul,tmur;
   double illrc=0.0;
   int cnt=0;
   bool lg,rg;  //are right and left mu's good
   double tpril,tprir;
   double sumpr=0.0;
   double temp;
   double lMax = -bigval;
   std::vector<double> priv;
   for(size_t i=0;i!=nm;i++) { //i indexes mul
      tmul = (*pi.mg)[i];
      tpril = (*pi.pm)[i];
      for(size_t j=i;j<nm;j++) { //j indexes mur
         tmur = (*pi.mg)[j];
         lg=false; rg=false;
         //fix right, check left
         nyr->setm(tmur);
         ml=mla; mu=mua;
         conint(ml,mu,nyl,y,xi);
         if((tmul >= ml) && (tmul<=mu)) lg=true; 
         //fix left, check right
         nyl->setm(tmul);
         ml=mla; mu=mua;
         conint(ml,mu,nyr,y,xi);
         if((tmur >= ml) && (tmur<=mu)) rg=true; 
         if(lg && rg) {
            cnt++;
            tprir = (*pi.pm)[j];
            temp = tpril*tprir;
            sumpr += temp;
            priv.push_back(temp);
            temp = -.5*(-2.0*tmul*sl.sy + sl.n*tmul*tmul)/(pi.sigma*pi.sigma); 
            temp += -.5*(-2.0*tmur*sr.sy + sr.n*tmur*tmur)/(pi.sigma*pi.sigma); 
            if(temp>lMax) lMax=temp;
            pv.push_back(temp);
            mulv.push_back(tmul);
            murv.push_back(tmur);
         }
      }
   }
   illrc=0.0;
   for(size_t i=0;i!=pv.size();i++) {
      pv[i] = priv[i]*exp(pv[i]-lMax);
      illrc += pv[i];
   }
   for(size_t i=0;i!=pv.size();i++) pv[i] /= illrc;
   lcfac=lMax;
   return illrc;
}

bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, rn& gen)
{
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal

      //draw bottom node, choose node index ni from list in goodbots 
      size_t ni = floor(gen.uniform()*goodbots.size()); 
      tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

      //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      size_t v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.mybeta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_cp nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
           Pnogy = 1.0/(nnogs+1.0);
         }
      }
  
      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuff(x,nx,v,c,xi,di,sl,sr);

      //--------------------------------------------------
      //compute alpha

      double alpha=0.0,alpha1=0.0,alpha2=0.0;
      double lill=0.0,lilr=0.0,lilt=0.0;
      if((sl.n>=5) && (sr.n>=5)) { //cloodge?
         lill = lil(sl.n,sl.sy,sl.sy2,pi.sigma,pi.tau);
         lilr = lil(sr.n,sr.sy,sr.sy2,pi.sigma,pi.tau);
         lilt = lil(sl.n+sr.n,sl.sy+sr.sy,sl.sy2+sr.sy2,pi.sigma,pi.tau);
   
         alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
         alpha2 = alpha1*exp(lill+lilr-lilt);
         alpha = std::min(1.0,alpha2);
      } else {
         alpha=0.0;
      }
      
      /*
      cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
      cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
      cout << "L,U: " << L << "," << U << endl;
      cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
         PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
         ", " << Pnogy << "," << Pbotx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop
      double a,b,s2,yb;
      double mul,mur; //means for new bottom nodes, left and right
      if(gen.uniform() < alpha) {
         //draw mul, mean for left node
         a= 1.0/(pi.tau*pi.tau); //a = 1/tau^2
         s2 = pi.sigma*pi.sigma; // sigma^2
         //left mean
         yb = sl.sy/sl.n;
         b = sl.n/s2; // b=n/sigma^2
         mul = b*yb/(a+b) + gen.normal()/sqrt(a+b);
         //draw mul, mean for left node
         yb = sr.sy/sr.n;
         b = sr.n/s2; // b=n/sigma^2
         mur = b*yb/(a+b) + gen.normal()/sqrt(a+b);
         //do birth
         x.birthp(nx,v,c,mul,mur);
         return true;
      } else {
         return false;
      }
   } else {
      //--------------------------------------------------
      //draw proposal

      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size()); 
      tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.mybeta);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,pi);
      double PGrx = pgrow(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
      if(!(nx->p)) { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one 
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuff(x,nx->getl(),nx->getr(),xi,di,sl,sr);

      //--------------------------------------------------
      //compute alpha

      double lill = lil(sl.n,sl.sy,sl.sy2,pi.sigma,pi.tau);
      double lilr = lil(sr.n,sr.sy,sr.sy2,pi.sigma,pi.tau);
      double lilt = lil(sl.n+sr.n,sl.sy+sr.sy,sl.sy2+sr.sy2,pi.sigma,pi.tau);

      double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
      double alpha2 = alpha1*exp(lilt - lill - lilr);
      double alpha = std::min(1.0,alpha2);

      /*
      cout << "death prop: " << nx->nid() << endl;
      cout << "nognds.size(), ni, nx: " << nognds.size() << ", " << ni << ", " << nx << endl;
      cout << "depth of nog node: " << dny << endl;
      cout << "PGny: " << PGny << endl;
      cout << "PGlx: " << PGlx << endl;
      cout << "PGrx: " << PGrx << endl;
      cout << "PBy: " << PBy << endl;
      cout << "Pboty: " << Pboty << endl;
      cout << "PDx: " << PDx << endl;
      cout << "Pnogx: " << Pnogx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "sigma: " << pi.sigma << endl;
      cout << "tau: " << pi.tau << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop
      double a,b,s2,yb;
      double mu;
      size_t n;
      if(gen.uniform()<alpha) {
         //draw mu for nog (which will be bot)
         n = sl.n + sr.n;
         a= 1.0/(pi.tau*pi.tau); //a = 1/tau^2
         s2 = pi.sigma*pi.sigma; // sigma^2
         yb = (sl.sy+sr.sy)/n;
         b = n/s2; // b=n/sigma^2
         mu = b*yb/(a+b) + gen.normal()/sqrt(a+b);
         //do death
         x.deathp(nx,mu);
         return true;
      } else {
         return false;
      }
   }
}

bool bdc(tree& x, xinfo& xi, dinfo& di, pinfo& pi, rn& gen)
{
   std::vector<double> pv;
   std::vector<double> mulv;
   std::vector<double> murv;
   std::vector<double> muv;
   double lfac1,lfac2;

   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death
      //cout << "\n**BIRTH\n";

      //--------------------------------------------------
      //draw proposal

      //draw bottom node, choose node index ni from list in goodbots 
      size_t ni = floor(gen.uniform()*goodbots.size()); 
      tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

      //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      size_t v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.mybeta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_cp nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
           Pnogy = 1.0/(nnogs+1.0);
         }
      }
      //cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
      //cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
      //cout << "L,U: " << L << "," << U << endl;
      //cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
      //   PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
      //   ", " << Pnogy << "," << Pbotx << endl;
  
      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuff(x,nx,v,c,xi,di,sl,sr);

      //cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      //cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;


      double alpha;
      double alpha1,alpha2;
      if((sl.n>=5) && (sr.n>=5)) { //cludge?

         //double iltc = coninteg1(sl,sr,xi,pi,x,nx,pv,muv);
         //cout << "first iltc: " << iltc << endl;
         double iltc = coninteg1(sl,sr,xi,pi,x,nx,pv,muv,lfac1);
         //cout << "second iltc: " << iltc*exp(lfac1) << endl;
         //cout << "in bdc birth, lfac1: " << lfac1 << endl;
         //iltc *= exp(lfac1);

         size_t nxid = nx->nid();
         tree y=x;
         y.birth(nxid,v,c,0.0,0.0);
         //double illrc = coninteg2(sl,sr,xi,pi,y,nxid,pv,mulv,murv);
         //cout << "illrc first: " << illrc << endl;
         //cout << "size pv 2: " << pv.size() << endl;
         double illrc = coninteg2(sl,sr,xi,pi,y,nxid,pv,mulv,murv,lfac2);
         //cout << "illrc second: " << illrc*exp(lfac2) << endl;
         //illrc *= exp(lfac2);
         
         alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
         //alpha2 = alpha1*illrc/iltc;
         alpha2 = alpha1*(illrc/iltc)*exp(lfac2-lfac1);
         alpha = std::min(1.0,alpha2);
      } else {
         alpha=0.0;
      }
      //cout << "alpha: " << alpha << endl;

      //--------------------------------------------------
      //finally ready to try metrop
      if(gen.uniform() < alpha) {
         size_t i = rdisc(&pv[0],gen);
         //cout << "i, mul, mur: " << i << ", " << mulv[i] << ", " << murv[i] << endl;
         x.birthp(nx,v,c,mulv[i],murv[i]);
         return true;
      } else {
         return false;
      }
   } else {
      //cout << "\n**DEATH\n";
      //--------------------------------------------------
      //draw proposal

      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size()); 
      tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.mybeta);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,pi);
      double PGrx = pgrow(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
      if(!(nx->p)) { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one 
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuff(x,nx->getl(),nx->getr(),xi,di,sl,sr);

      //cout << "nog node: " << nx << endl;
      //cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      //cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;

      //--------------------------------------------------
      //compute alpha

      size_t nidx = nx->nid();
      //double illrc = coninteg2(sl,sr,xi,pi,x,nidx,pv,mulv,murv);
      double illrc = coninteg2(sl,sr,xi,pi,x,nidx,pv,mulv,murv,lfac2);
      //illrc *= exp(lfac2);
      //cout << "size pv 2: " << pv.size() << endl;
      //cout << "death illrc: " << illrc << endl;

      tree y = x;
      y.death(nidx,0.0);
      tree::tree_p ny = y.getptr(nidx);
      //y.pr();
      //cout << "ny: " << ny;
      //double iltc = coninteg1(sl,sr,xi,pi,y,ny,pv,muv);
        double iltc = coninteg1(sl,sr,xi,pi,x,ny,pv,muv,lfac1);
         //cout << "second iltc: " << iltc*exp(lfac1) << endl;
         //cout << "in bdc birth, lfac1: " << lfac1 << endl;
         //iltc = iltc*exp(lfac1);
      //cout << "\ndeath iltc: " << iltc << endl;
      //cout << "size pv 1: " << pv.size() << endl;

      //double lill = lil(sl.n,sl.sy,sl.sy2,pi.sigma,pi.tau);
      //double lilr = lil(sr.n,sr.sy,sr.sy2,pi.sigma,pi.tau);
      //double lilt = lil(sl.n+sr.n,sl.sy+sr.sy,sl.sy2+sr.sy2,pi.sigma,pi.tau);

      double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
      //double alpha2 = alpha1*exp(lilt - lill - lilr);
      //double alpha2 = alpha1*iltc/illrc;
      double alpha2 = alpha1*(iltc/illrc)*exp(lfac1-lfac2);
      double alpha = std::min(1.0,alpha2);
      //cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;

      /*
      cout << "death prop: " << nx->nid() << endl;
      cout << "nognds.size(), ni, nx: " << nognds.size() << ", " << ni << ", " << nx << endl;
      cout << "depth of nog node: " << dny << endl;
      cout << "PGny: " << PGny << endl;
      cout << "PGlx: " << PGlx << endl;
      cout << "PGrx: " << PGrx << endl;
      cout << "PBy: " << PBy << endl;
      cout << "Pboty: " << Pboty << endl;
      cout << "PDx: " << PDx << endl;
      cout << "Pnogx: " << Pnogx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "sigma: " << pi.sigma << endl;
      cout << "tau: " << pi.tau << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop
      size_t ii; 
      if(gen.uniform()<alpha) {
         ii = rdisc(&pv[0],gen);
         //cout << "accept death, ii, muv.size(): " << ii << ", " << pv.size() << endl;
         //do death
         x.deathp(nx,muv[ii]);
         return true;
      } else {
         return false;
      }
   }
}

#endif
