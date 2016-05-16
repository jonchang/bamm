#ifndef TRAIT_BRANCH_EVENT_H
#define TRAIT_BRANCH_EVENT_H

#include "BranchEvent.h"

class Tree;
class Node;
class Random;


class TraitBranchEvent : public BranchEvent
{

private:

    double _betaInit;  // initial beta value.
    double _betaShift; // temporal shift parameter of trait evolution rate.
    bool _isTimeVariable;
    bool _isJump;
    double _jump;
public:

    // constructors, depending on whether you want trait rate or lambda/mu
    TraitBranchEvent(double beta, double shift, bool isTimeVariable,
            bool isJump, double jump, Node* x, Tree* tp, Random& random, double map);
    virtual ~TraitBranchEvent() {};

    void   setBetaInit(double x);
    double getBetaInit();

    void   setBetaShift(double x);
    double getBetaShift();

    void setTimeVariable(bool isTimeVariable);
    bool isTimeVariable();

    // Jump-model parameters
    bool isJump();
    void setJump(double x);
    double getJump();
    
    //virtual void setIsEventValidForNode(bool x);
    virtual bool getIsEventValidForNode();
};


inline void TraitBranchEvent::setBetaInit(double x)
{
    _betaInit = x;
}


inline double TraitBranchEvent::getBetaInit()
{
    return _betaInit;
}


inline void TraitBranchEvent::setBetaShift(double x)
{
    _betaShift = x;
}


inline double TraitBranchEvent::getBetaShift()
{
    return _betaShift;
}


inline void TraitBranchEvent::setTimeVariable(bool isTimeVariable)
{
    _isTimeVariable = isTimeVariable;
}


inline bool TraitBranchEvent::isTimeVariable()
{
    return _isTimeVariable;
}

inline bool TraitBranchEvent::isJump()
{
    return _isJump;
}

inline double TraitBranchEvent::getJump()
{
    return _jump;
}

inline void TraitBranchEvent::setJump(double x)
{
    _jump = x;
}

#endif
