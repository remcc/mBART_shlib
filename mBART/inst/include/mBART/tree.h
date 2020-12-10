#ifndef GUARD_tree_h
#define GUARD_tree_h

/*
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <vector>

#include "info.h"
#include "rn.h"
*/

/*
Class is a bit confusing because sometimes you
want to think about an instance as a node and
sometimes as the whole tree including the node
and all its children.
Fundamentally, we have a tree not a node.
*/

/*
three ways to access a node:
(i) node id: is the integer assigned by the node numbering system 
     assuming they are all there
(ii) node ind: is the index into the array of 
   node(tree) pointers returned by getnodes or getbots or getnogs
   which means you go left to right across the bottom of the tree
(iii) by its pointer (should only "give out" const pointers)
*/

//info contained in a node, used by input operator
struct node_info {
   std::size_t id; //node id
   std::size_t v;  //variable
   std::size_t c;  //cut point
   double m;       //mu
};

//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
// left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules

class tree {
public:
   //------------------------------
   //friends
   friend std::istream& operator>>(std::istream&, tree&);
   friend bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, rn& gen);
   friend bool bdc(tree& x, xinfo& xi, dinfo& di, pinfo& pi, rn& gen);

   //------------------------------
   //typedefs
   typedef tree* tree_p;
   typedef const tree* tree_cp;
   typedef std::vector<tree_p> npv; //Node Pointer Vector
   typedef std::vector<tree_cp> cnpv; //const Node Pointer Vector

   //------------------------------
   //tree constructors, destructors
   tree();
   tree(const tree&);
   tree(double);
   ~tree() {tonull();}

   //------------------------------
   //operators
   tree& operator=(const tree&);

   //------------------------------
   //node access
   //you are freely allowed to change mu, v, c
   //set----------
   void setm(double mu) {this->mu=mu;}
   void setv(size_t v) {this->v = v;}
   void setc(size_t c) {this->c = c;}
   //get----------
   double getm() const {return mu;} 
   size_t getv() const {return v;}
   size_t getc() const {return c;}
   tree_p getp() const {return p;}  //should this be tree_cp? 
   tree_p getl() const {return l;}
   tree_p getr() const {return r;}

   //------------------------------
   //tree functions
   //stuff about the tree----------
   size_t treesize() const; //number of nodes in tree
   size_t nnogs() const;    //number of nog nodes (no grandchildren nodes)
   size_t nbots() const;    //number of bottom nodes
   void pr(bool pc=true) const; //to screen, pc is "print children"
   //birth death using nid----------
   bool birth(size_t nid,size_t v, size_t c, double ml, double mr); 
   bool death(size_t nid, double mu); 
   //vectors of node pointers----------
   void getbots(npv& bv);         //get bottom nodes
   void getnogs(npv& nv);         //get nog nodes (no granchildren)
   void getnodes(npv& v);         //get vector of all nodes
   void getnodes(cnpv& v) const;  //get all nodes
   //find node from x and region for var----------
   tree_cp bn(double *x,xinfo& xi);  //find bottom node for x
   void rg(size_t v, int* L, int* U); //find region [L,U] for var v.
   //------------------------------
   //node functions
   size_t depth() const; //depth of a node
   size_t nid() const;   //node id
   char ntype() const;   //t:top;b:bottom;n:nog;i:interior, carefull a t can be bot
   tree_p getptr(size_t nid); //get node pointer from node id, 0 if not there.
   bool isnog() const;

   //------------------------------
   //constraint functions
   //'a': neighbor above, 'b': neighbor below, 'd': disjoint, 'x': nothing
   char nhb(tree_p n, xinfo& xi); 

private:
   //------------------------------
   //parameter for node
   double mu; 
   //------------------------------
   //rule: left if x[v] < xinfo[v][c]
   size_t v; 
   size_t c; 
   //------------------------------
   //tree structure
   tree_p p; //parent
   tree_p l; //left child
   tree_p r; //right child
   //------------------------------
   //utiity functions
   void cp(tree_p n,  tree_cp o); //copy tree
   void tonull(); //like a "clear", null tree has just one node
   void birthp(tree_p np,size_t v, size_t c, double ml, double mr); 
   void deathp(tree_p nb, double mu); //kill children of nog node nb 
};
std::istream& operator>>(std::istream&, tree&);
std::ostream& operator<<(std::ostream&, const tree&);

/*
#include <string>
#include <vector>
#include <map>
#include "tree.h"

using std::string;
using std::cout;
using std::endl;
*/

//--------------------------------------------------
// constructors
tree::tree(): mu(0.0),v(0),c(0),p(0),l(0),r(0) {}
tree::tree(double m): mu(m),v(0),c(0),p(0),l(0),r(0) {}
tree::tree(const tree& n): mu(0.0),v(0),c(0),p(0),l(0),r(0) { cp(this,&n); }
//--------------------------------------------------
//operators
tree& tree::operator=(const tree& rhs)
{
   if(&rhs != this) {
      tonull(); //kill left hand side (this)
      cp(this,&rhs); //copy right hand side to left hand side
   }
   return *this;
}
//--------------------------------------------------
//public functions
// find bottom node pointer given x
//--------------------
tree::tree_cp tree::bn(double *x,xinfo& xi)
{
   if(l==0) return this; //bottom node
   if(x[v] < xi[v][c]) {
      return l->bn(x,xi);
   } else {
      return r->bn(x,xi);
   }
}
//--------------------
//find region for a given variable
void tree::rg(size_t v, int* L, int* U)
{ 
   if(p==0)  { //no parent
      return;
   }
   if(p->v == v) { //does my parent use v?
      if(this == p->l) { //am I left or right child
         if((int)(p->c) <= (*U)) *U = (p->c)-1; 
      } else {
         if((int)(p->c) >= *L) *L = (p->c)+1; 
      }
   } 
   p->rg(v,L,U);
}
//--------------------
//tree size
size_t tree::treesize() const
{
   if(l==0) return 1;  //if bottom node, tree size is 1
   else return (1+l->treesize()+r->treesize());
}
//--------------------
size_t tree::nnogs() const
{
   if(!l) return 0; //bottom node
   if(l->l || r->l) { //not a nog
      return (l->nnogs() + r->nnogs()); 
   } else { //is a nog
      return 1;
   }
}
//--------------------
size_t tree::nbots() const
{
   if(l==0) { //if a bottom node
      return 1;
   } else {
      return l->nbots() + r->nbots();
   }
}
//--------------------
//depth of node
size_t tree::depth() const
{
   if(!p) return 0; //no parents
   else return (1+p->depth());
}
//--------------------
// node id
size_t tree::nid() const
//recursion up the tree
{
   if(!p) return 1; //if you don't have a parent, you are the top
   if(this==p->l) return 2*(p->nid()); //if you are a left child
   else return 2*(p->nid())+1; //else you are a right child
}
//--------------------
//node type
char tree::ntype() const
{
   //t:top, b:bottom, n:no grandchildren, i:internal
   if(!p) return 't';
   if(!l) return 'b';
   if(!(l->l) && !(r->l)) return 'n';
   return 'i';
}
//--------------------
//get bottom nodes
//recursion down the tree
void tree::getbots(npv& bv) 
{
   if(l) { //have children
      l->getbots(bv);
      r->getbots(bv);
   } else {
      bv.push_back(this);
   }
}
//--------------------
//get nog nodes
//recursion down the tree
void tree::getnogs(npv& nv) 
{
   if(l) { //have children
      if((l->l) || (r->l)) {  //have grandchildren
         if(l->l) l->getnogs(nv);
         if(r->l) r->getnogs(nv);
      } else {
         nv.push_back(this);
      }
   }
}
//--------------------
//get all nodes
//recursion down the tree
void tree::getnodes(npv& v) 
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
void tree::getnodes(cnpv& v)  const
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
//--------------------
//add children to  bot node nid
bool tree::birth(size_t nid,size_t v, size_t c, double ml, double mr)
{
   tree_p np = getptr(nid);
   if(np==0) {
      cout << "error in birth: bottom node not found\n";
      return false; //did not find note with that nid
   }
   if(np->l) {
      cout << "error in birth: found node has children\n";
      return false; //node is not a bottom node
   }

   //add children to bottom node np 
   tree_p l = new tree;
   l->mu=ml;
   tree_p r = new tree;
   r->mu=mr;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;

   return true;
}
//--------------------
//is the node a nog node
bool tree::isnog() const
{
   bool isnog=true;
   if(l) {
      if(l->l || r->l) isnog=false; //one of the children has children.
   } else {
      isnog=false; //no children
   }
   return isnog;
}
//--------------------
//kill children of  nog node nid
bool tree::death(size_t nid, double mu) 
{
   tree_p nb = getptr(nid);
   if(nb==0) {
      cout << "error in death, nid invalid\n";
      return false;
   }
   if(nb->isnog()) {
      delete nb->l;
      delete nb->r;
      nb->l=0;
      nb->r=0;
      nb->v=0;
      nb->c=0;
      nb->mu=mu;
      return true;
   } else {
      cout << "error in death, node is not a nog node\n";
      return false;
   }
}
//--------------------
//add children to bot node *np
void tree::birthp(tree_p np,size_t v, size_t c, double ml, double mr)
{
   tree_p l = new tree;
   l->mu=ml;
   tree_p r = new tree;
   r->mu=mr;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;
}
//--------------------
//kill children of  nog node *nb
void tree::deathp(tree_p nb, double mu) 
{
   delete nb->l;
   delete nb->r;
   nb->l=0;
   nb->r=0;
   nb->v=0;
   nb->c=0;
   nb->mu=mu;
}
//--------------------
//print out tree(pc=true) or node(pc=false) information
//uses recursion down
void tree::pr(bool pc) const
{
   size_t d = depth();
   size_t id = nid();

   size_t pid; 
   if(!p) pid=0; //parent of top node
   else pid = p->nid();

   string pad(2*d,' ');
   string sp(", ");
   if(pc && (ntype()=='t')) 
      cout << "tree size: " << treesize() << endl;
   //cout << pad << "(id,parent): " << id << sp << pid;
   cout << pad << "id: " << id;
   cout << sp << "(v,c): " << v << sp << c;
   cout << sp << "mu: " << mu; 
   cout << sp << "type: " << ntype();
   cout << sp << "depth: " << depth();
   cout << sp << "pointer: " << this << endl;

   if(pc) {
      if(l) {
         l->pr(pc);
         r->pr(pc);
      }
   }
}
//--------------------------------------------------
//private functions
//--------------------
//copy tree o to tree n
void tree::cp(tree_p n, tree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
{
   if(n->l) {
      cout << "cp:error node has children\n";
      return;
   }

   n->mu = o->mu;
   n->v = o->v;
   n->c = o->c;

   if(o->l) { //if o has children
      n->l = new tree;
      (n->l)->p = n;
      cp(n->l,o->l);
      n->r = new tree;
      (n->r)->p = n;
      cp(n->r,o->r);
   }
}
//--------------------
//cut back to one node
void tree::tonull()
{
   size_t ts = treesize();
   while(ts>1) { //if false ts=1
      npv nv;
      getnogs(nv);
      for(size_t i=0;i<nv.size();i++) {
         delete nv[i]->l;
         delete nv[i]->r;
         nv[i]->l=0;
         nv[i]->r=0;
      }
      ts = treesize(); 
   }
   mu=0.0;
   v=0;c=0;
   p=0;l=0;r=0;
}
//--------------------------------------------------
//functions
//--------------------
//input operator
std::ostream& operator<<(std::ostream& os, const tree& t)
{
   tree::cnpv nds;
   t.getnodes(nds);
   os << nds.size() << endl;
   for(size_t i=0;i<nds.size();i++) {
      os << nds[i]->nid() << " ";
      os << nds[i]->getv() << " ";
      os << nds[i]->getc() << " ";
      os << nds[i]->getm() << endl;
   }
   return os;
}
//--------------------
//ouput operator
std::istream& operator>>(std::istream& is, tree& t)
{
   size_t tid,pid; //tid: id of current node, pid: parent's id
   std::map<size_t,tree::tree_p> pts;  //pointers to nodes indexed by node id
   size_t nn; //number of nodes

   t.tonull(); // obliterate old tree (if there)

   //read number of nodes----------
   is >> nn;
   if(!is) {
      //cout << ">> error: unable to read number of nodes" << endl;
      return is;
   }

   //read in vector of node information----------
   std::vector<node_info> nv(nn);
   for(size_t i=0;i!=nn;i++) {
      is >> nv[i].id >> nv[i].v >> nv[i].c >> nv[i].m;
      if(!is) {
         //cout << ">> error: unable to read node info, on node  " << i+1 << endl;
         return is;
      }
   }

   //first node has to be the top one
   pts[1] = &t; //careful! this is not the first pts, it is pointer of id 1.
   t.setv(nv[0].v); t.setc(nv[0].c); t.setm(nv[0].m);
   t.p=0;

   //now loop through the rest of the nodes knowing parent is already there.
   for(size_t i=1;i!=nv.size();i++) {
      tree::tree_p np = new tree;
      np->v = nv[i].v; np->c=nv[i].c; np->mu=nv[i].m;
      tid = nv[i].id;
      pts[tid] = np;
      pid = tid/2;
      // set pointers
      if(tid % 2 == 0) { //left child has even id
         pts[pid]->l = np;
      } else {
         pts[pid]->r = np;
      }
      np->p = pts[pid];
   }
   return is;
}
//--------------------
// get node pointer from its nid
tree::tree_p tree::getptr(size_t nid)
{
   if(this->nid() == nid) return this; //found it
   if(l==0) return 0; //no children, did not find it
   tree_p lp = l->getptr(nid);
   if(lp) return lp; //found on left
   tree_p rp = r->getptr(nid); 
   if(rp) return rp; //found on right
   return 0; //never found it
}
//--------------------
// how does node n relate to this
//'a': neighbor above, 'b': neighbor below, 'd': disjoint, 'x': nothing
char tree::nhb(tree_p n, xinfo& xi)
//note:
//  each node is a region of the form \cap [Lv,Uv] v=0,1,..,(p-1)
{
   if(this==n) return 's'; //s means self
   size_t p = xi.size(); //need to loop over the p variables

   int mL,mU,oL,oU; //my range, other range
   short nind=0;
   //for each var=v, see if regions is disjoint (then done), or above or below
   // is it true that any n could only be a (say below neighor) on one v?
   //   is so could return 'b' or 'a' right away instead of using nind.
   //  consider a box below a box in three 3 then you could have it.
   //  could return as soon as you find below anyway
   for(size_t v=0;v<p;v++) { 
      mL=0; oL=0;
      mU = oU = xi[v].size()-1;
      rg(v,&mL,&mU);
      n->rg(v,&oL,&oU);
      //a neighbor will be 2 away since [L,U] are the usable ones, 
         //for example: 'b': [oL,oU] used cut point=C, [mL,mU] 
         //would have oU,C,mL in sequence
      if(oU < mL-2 || oL > mU+2 ) return 'd';
      if(oU == mL-2) nind = 1;
      if(oL == mU+2) nind = 2;
   }
   if(nind==1) return 'b';
   if(nind==2) return 'a';
   return 'x';
}

#endif
