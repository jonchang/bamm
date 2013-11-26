#include <iostream>
#include <stdlib.h>

#include "SpExBranchEvent.h"
#include "Node.h"
#include "Tree.h"
#include "MbRandom.h"


/*
    Note that SpExBranchEvent inherits shift parameter,
    thus addition of event need not change likelihood of tree.
*/


SpExBranchEvent::SpExBranchEvent(double speciation, double lamshift,
    double extinction, double mushift, Node* x, Tree* tp,
    MbRandom* rp, double map)
{
    // Speciation and extinction from new event assumed constant in time...

    //std::cout << "Event ctor: map " << map << std::endl;


    _lamInit = speciation;
    _lamShift = lamshift;
    _muInit = extinction;
    _muShift = mushift;

    if (tp->getRoot() == x) {
        _absTime = 0.0;
    } else {
        _absTime = tp->getAbsoluteTimeFromMapTime(map);
    }

    nodeptr = x;
    mapTime = map;
    treePtr = tp;
	ranPtr = rp;

    oldMapTime = map;
    oldNodePtr = x;

    _isEventTimeVariable = false;

}


SpExBranchEvent::~SpExBranchEvent(void)
{


}

/*
 bool SpExBranchEvent::operator<(const BranchEvent& a) const
 Modified on 9.17.2012.
 if maptime is > than a.maptime
 THEN event must have occurred earlier than a.maptime event
 THus operator < would be TRUE. for a->maptime > b->maptime,
                        then a < b evaluates to TRUE

 */

bool SpExBranchEvent::operator<(const SpExBranchEvent& a) const
{

    if ( mapTime > a.mapTime ) {
        return true;
    } else {
        return false;
    }
}


/* this (along with incrementMapPosition) implements a uniform proposal mechanism
    with reflecting boundaries at all major tree boundaries: at the root, the proposal
    is bounced back onto the parent branch.
 At the tips OR at invalid nodes (defined in tree constructor)


 This mechanism of local move involves actually moving the position of the event
 If the proposed move is accepted, this is fine. If the proposed move is rejected,
    then the branchEvent will need to be reset to its original value.
 IN practice, this should be done by:
    oldMapPosition = getMapTime();
    moveEventLocal();
    calculate likelihood etc
    if (accept)
        event has already been updated.
    if (reject)
        setEventByMapPosition(oldMapPosition)
            // should undo it.

*/


void SpExBranchEvent::moveEventLocal(double stepsize)
{

    oldNodePtr = nodeptr;
    oldMapTime = mapTime;

    incrementMapPosition(stepsize);

    // shouldn't actually need to update model attributes.
    // if node position shifts, it just brings its current attributes
    // along, e.g., lambda etc.

}

void SpExBranchEvent::incrementMapPosition(double ink)
{
    double mapstart = nodeptr->getMapStart();
    double mapend = nodeptr->getMapEnd();

    double temp = getMapTime() + ink;

    // Map end is ROOTWARDS
    // Map start is TIPWARDS

    if (temp > mapend) {
        // Goes to next most "ROOTWARDS" branch
        // BUT if there is no next branch (because ancestor of curr branch is ROOT)
        //      go FORWARD down other descendant.

        ink = temp - mapend; // residual...

        if (getEventNode() == treePtr->getRoot()) {
            std::cout << " Root problem in incrementMapPosition: should never get here" << std::endl;
            throw;
        }

        if (getEventNode()->getAnc() == treePtr->getRoot()) {

            ink = -1 * ink;

            if ((treePtr->getRoot()->getLfDesc() == getEventNode()) &
                    treePtr->getRoot()->getRtDesc()->getCanHoldEvent()) {
                setEventNode(treePtr->getRoot()->getRtDesc());

            } else if ((treePtr->getRoot()->getRtDesc() == getEventNode()) &
                     treePtr->getRoot()->getLfDesc()->getCanHoldEvent()) {
                setEventNode(treePtr->getRoot()->getLfDesc());

            } else {
                // If you get here, either the right or left desc of the root CANNOT hold event
                // Do nothing. Just go backwards and reflect off of root.
            }
            setMapTime(getEventNode()->getMapEnd());
            incrementMapPosition(ink);

        } else {

            setEventNode(getEventNode()->getAnc());
            setMapTime(getEventNode()->getMapStart());
            incrementMapPosition(ink);
        }

    } else if (temp < mapstart) {

        ink = temp - mapstart; // ink is now negative...

        if (getEventNode()->getLfDesc() == NULL &&
                getEventNode()->getRtDesc() == NULL) {
            // At terminal branch:
            // Reflect
            ink  = -1 * ink; // make positive....
            setMapTime(mapstart);
            incrementMapPosition(ink);

        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == true) {
            // both desc branches valid
            double ran = ranPtr->uniformRv();
            if (ran <= 0.5) {
                // left branch
                setEventNode(getEventNode()->getLfDesc());
                setMapTime(getEventNode()->getMapEnd());

            } else {
                setEventNode(getEventNode()->getRtDesc());
                setMapTime(getEventNode()->getMapEnd());
            }
            incrementMapPosition(ink);
        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == false &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == true) {
            setEventNode(getEventNode()->getRtDesc());
            setMapTime(getEventNode()->getMapEnd());
            incrementMapPosition(ink);

        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == false) {

            setEventNode(getEventNode()->getLfDesc());
            setMapTime(getEventNode()->getMapEnd());
            incrementMapPosition(ink);
        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == false &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == false) {
            // neither node can hold event: reflect
            ink  = -1 * ink; // make positive....

            // Reflect:
            setMapTime(mapstart);
            incrementMapPosition(ink);

        } else {
            std::cout << "Problem in incrementMapPosition()" << std::endl;
            throw;
        }

    } else {

        setMapTime(temp);
    }

    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));

}

void SpExBranchEvent::moveEventGlobal(void)
{

    oldNodePtr = nodeptr;
    oldMapTime = mapTime;

    double aa = treePtr->getRoot()->getMapStart();
    double bb = treePtr->getTotalMapLength();
    double position = ranPtr->uniformRv(aa, bb);
    setEventNode(treePtr->mapEventToTree(position));
    setMapTime(position);
    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));

}


// test
void SpExBranchEvent::setEventByMapPosition(double x)
{
    setEventNode(treePtr->mapEventToTree(x));
    setMapTime(x);
    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
}

void SpExBranchEvent::revertOldMapPosition(void)
{

    nodeptr = oldNodePtr;
    mapTime = oldMapTime;
    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));

}

