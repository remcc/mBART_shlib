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

#endif
