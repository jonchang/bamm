#ifndef SP_EX_BRANCH_EVENT_H
#define SP_EX_BRANCH_EVENT_H


// Forward declarations:
class Tree;
class Node;
class MbRandom;


/*
class SpExBranchEvent contains:
    (1) the node associated with the event
    (2) the map position of the event
    (3) parameters associated with the event, e.g., lambda

The initial event will always be the root node with a map position of 0.
This event is immutable.
*/

class SpExBranchEvent
{

public:

    class PtrCompare
    {
    public:
        bool operator()(const SpExBranchEvent* m1,
                        const SpExBranchEvent* m2) const {
            return *m1 < *m2;
        }
    };

private:

    double mapTime;
    Node* nodeptr;
    Tree* treePtr;
	MbRandom* ranPtr;

    // Keep values of the old pointer and old maptime associated
    // with the event for FAST reference if a LOCAL proposal is rejected.

    Node* oldNodePtr;
    double oldMapTime;

    /******************/
    // New event parameters, March 25, 2012
    double _absTime;   // real time, measured with t = 0 at root.
    double _lamInit;   // Initial speciation rate at event
    double _muInit;    // Initial Mu rate at event
    double _lamShift;  // magnitude & direction of speciation shift
    double _muShift;   // magnitude & direction of mu shift

    /*****************/
    // New parameters June 12 2012
    // allow rjMCMC to move between time-varying and time-constant partitions.
    bool _isEventTimeVariable;

public:

    // constructors, depending on whether you want trait rate or lambda/mu
    SpExBranchEvent(double speciation, double lamshift, double extinction,
        double mushift, Node* x, Tree* tp, MbRandom*rp, double map);

    ~SpExBranchEvent();

    void   setMapTime(double x);
    double getMapTime();

    void  setEventNode(Node* x);
    Node* getEventNode();

    void   setAbsoluteTime(double x);
    double getAbsoluteTime();

    void   setLamInit(double x);
    double getLamInit();

    void   setMuInit(double x);
    double getMuInit();

    void   setLamShift(double x);
    double getLamShift();

    void   setMuShift(double x);
    double getMuShift();

    void incrementMapPosition(double ink);
    void moveEventLocal(double stepsize);
    void moveEventGlobal();
    void setEventByMapPosition(double x);

    // Functions to set and manipulate OLD events:
    void  setOldEventNode(Node* x);
    Node* getOldEventNode();

    void   setOldMapTime(double x);
    double getOldMapTime();

    // Revert to old map position using oldPtr and oldMapTime
    // this only works if you have changed the nodeptr and maptime
    // relative to the values of oldNodePtr and oldMapTime
    void revertOldMapPosition();

    // Overloading comparision operator:
    bool operator<(const SpExBranchEvent& a) const;

    // For time-varying rjMCMC:
    void setIsEventTimeVariable(bool x);
    bool getIsEventTimeVariable();
};


inline void SpExBranchEvent::setMapTime(double x)
{
    mapTime = x;
}


inline double SpExBranchEvent::getMapTime()
{
    return mapTime;
}


inline void SpExBranchEvent::setEventNode(Node* x)
{
    nodeptr = x;
}


inline Node* SpExBranchEvent::getEventNode()
{
    return nodeptr;
}


inline void SpExBranchEvent::setAbsoluteTime(double x)
{
    _absTime = x;
}


inline double SpExBranchEvent::getAbsoluteTime()
{
    return _absTime;
}


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


inline void SpExBranchEvent::setOldEventNode(Node* x)
{
    oldNodePtr = x;
}


inline Node* SpExBranchEvent::getOldEventNode()
{
    return oldNodePtr;
}


inline void SpExBranchEvent::setOldMapTime(double x)
{
    oldMapTime = x;
}


inline double SpExBranchEvent::getOldMapTime()
{
    return oldMapTime;
}


inline bool SpExBranchEvent::getIsEventTimeVariable()
{
    return _isEventTimeVariable;
}


inline void SpExBranchEvent::setIsEventTimeVariable(bool x)
{
    _isEventTimeVariable = x;
}


#endif
