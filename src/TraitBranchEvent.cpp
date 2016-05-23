#include "TraitBranchEvent.h"

#include "global_macros.h"


class Random;
class Node;
class Tree;


TraitBranchEvent::TraitBranchEvent(double beta, double shift,
    bool isTimeVariable, bool isJump, double jump, Node* x, Tree* tp, Random& random, double map) :
        BranchEvent(x, tp, random, map),
        _betaInit(beta), _betaShift(shift), _isTimeVariable(isTimeVariable),
            _isJump(isJump), _jump(jump)
{
}


bool TraitBranchEvent::getIsEventValidForNode()
{
    return (!isJump());
}
