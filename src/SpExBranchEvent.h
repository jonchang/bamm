#ifndef SP_EX_BRANCH_EVENT_H
#define SP_EX_BRANCH_EVENT_H

#include "BranchEvent.h"

class Tree;
class Node;
class MbRandom;


class SpExBranchEvent : public BranchEvent
{

private:

    double _lamInit;   // Initial speciation rate at event
    double _lamShift;  // magnitude & direction of speciation shift
    double _muInit;    // Initial Mu rate at event
    double _muShift;   // magnitude & direction of mu shift

public:

    // constructors, depending on whether you want trait rate or lambda/mu
    SpExBranchEvent(double speciation, double lamshift, double extinction,
        double mushift, Node* x, Tree* tp, MbRandom*rp, double map);

    void   setLamInit(double x);
    double getLamInit();

    void   setMuInit(double x);
    double getMuInit();

    void   setLamShift(double x);
    double getLamShift();

    void   setMuShift(double x);
    double getMuShift();

};


inline void SpExBranchEvent::setLamInit(double x)
{
    _lamInit = x;
}


inline double SpExBranchEvent::getLamInit()
{
    return _lamInit;
}


inline void SpExBranchEvent::setMuInit(double x)
{
    _muInit = x;
}


inline double SpExBranchEvent::getMuInit()
{
    return _muInit;
}


inline void SpExBranchEvent::setLamShift(double x)
{
    _lamShift = x;
}


inline double SpExBranchEvent::getLamShift()
{
    return _lamShift;
}


inline void SpExBranchEvent::setMuShift(double x)
{
    _muShift = x;
}


inline double SpExBranchEvent::getMuShift()
{
    return _muShift;
}


#endif
