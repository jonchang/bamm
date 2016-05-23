#include "SpExBranchEvent.h"

#include "global_macros.h"


class Node;
class Tree;
class Random;


SpExBranchEvent::SpExBranchEvent(double speciation, double lamshift,
        double extinction, double mushift, bool isTimeVariable,
        Node* x, Tree* tp, Random& random, double map) :
    BranchEvent(x, tp, random, map), _lamInit(speciation), _lamShift(lamshift),
        _muInit(extinction), _muShift(mushift), _isTimeVariable(isTimeVariable)
{
}


// This is implemented solely to use this function to avoid mapping
// phenotypic jumps to nodes in the phenotypic jump model.
bool SpExBranchEvent::getIsEventValidForNode()
{
    return true;
}
