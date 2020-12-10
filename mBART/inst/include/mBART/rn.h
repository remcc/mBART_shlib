#ifndef GUARD_rn
#define GUARD_rn

//pure virtual base class for random numbers
class rn
{
public:
   virtual double normal() = 0; //standard normal
   virtual double uniform() = 0; //uniform(0,1)
   virtual double chi_square() = 0; //chi-square
   virtual void set_df(int df) = 0; //set df for chi-square
   virtual ~rn() {}
};

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
