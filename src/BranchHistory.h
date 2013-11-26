#ifndef BRANCH_HISTORY_H
#define BRANCH_HISTORY_H

#include <set>
#include "Log.h"


template <typename EventType>
class BranchHistory
{

    typedef std::set<EventType*, typename EventType::PtrCompare> EventSet;
    typedef typename EventSet::size_type EventSetSizeType;

private:

    EventType* _nodeEvent;             // event describing focal node
    EventType* _ancestralNodeEvent;    // event describing ancestor

    // Set of all events on branch. Of length 0 if no events occurred on branch.
    // Also, if no events occur on branch, then entire branch is described by
    // the event referenced at nodeEvent
    EventSet _eventsOnBranch;

public:

    BranchHistory();

    EventType* getLastEvent();

    // Get last event from a reference event that occurred on branch:
    EventType* getLastEvent(EventType* x);
    EventType* getLastEvent(double ttime);
    EventType* getNextEvent(double ttime);
    int getNumberOfEventsOnInterval(double t1, double t2);
    EventType* getEventByIndexPosition(int i);

    void       setNodeEvent(EventType* x);
    EventType* getNodeEvent();

    void       setAncestralNodeEvent(EventType* x);
    EventType* getAncestralNodeEvent();

    void printBranchHistory();
    void reversePrintBranchHistory();
    void printEvent(EventType* event);

    void popEventOffBranchHistory(EventType* x);
    void addEventToBranchHistory(EventType* x);
    int  getNumberOfBranchEvents();
};


template <typename EventType>
BranchHistory<EventType>::BranchHistory() :
    _nodeEvent(NULL), _ancestralNodeEvent(NULL)
{
}


template <typename EventType>
EventType* BranchHistory<EventType>::getLastEvent()
{
    if (_eventsOnBranch.size() == 0) {
        return NULL;
    }

    typename EventSet::reverse_iterator it = _eventsOnBranch.rbegin();
    return *it;
}


// Return the last event overall for a given event,
// assuming branch history is set correctly.
template <typename EventType>
EventType* BranchHistory<EventType>::getLastEvent(EventType* x)
{
    EventType* theLastEvent = NULL;

    typename EventSet::iterator it;
    for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end(); ++it) {
        if (*it == x) {
            if (it == _eventsOnBranch.begin()) {
                theLastEvent = getAncestralNodeEvent();
            } else {
                theLastEvent = *(--it);
            }

            break;
        }
    }

    if (theLastEvent == NULL) {
        log(Warning) << "Problem in BranchHistory::getLastEvent()\n";
    }

    return theLastEvent;
}


// Returns most recent upstream event from a given absolute time,
// e.g., the most "rootward" event
template <typename EventType>
EventType* BranchHistory<EventType>::getLastEvent(double ttime)
{
    EventType* event = getAncestralNodeEvent();

    if (getNumberOfBranchEvents() > 0) {
        if (ttime > ((*_eventsOnBranch.begin())->getAbsoluteTime()) ) {
            if (_eventsOnBranch.size() == 1) {
                event = (*_eventsOnBranch.begin());
            } else {
                typename EventSet::iterator it;
                for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end();
                        ++it) {
                    if ((*it)->getAbsoluteTime() < ttime) {
                        event = *it;
                    } else {
                        break;
                    }
                }
            }
        } else {
            // Do nothing. Ancestral node event should be returned,
            // because ttime is BEFORE any events on branch
        }
    }

    return event;
}


// Returns most recent downstream event from a given absolute time,
// e.g.. the most "tipward" event.
// Does not check for range violation, e.g., if ttime is not in the focal branch
// If no events on branch, this will just be the EventNode.
template <typename EventType>
EventType* BranchHistory<EventType>::getNextEvent(double ttime)
{
    EventType* event = getNodeEvent();

    if (getNumberOfBranchEvents() > 0) {
        if (ttime < ((*(--_eventsOnBranch.end()))->getAbsoluteTime()) ) {
            if (_eventsOnBranch.size() == 1) {
                event = (*_eventsOnBranch.begin());
            } else {
                typename EventSet::iterator it;
                for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end();
                        ++it) {
                    if ((*it)->getAbsoluteTime() > ttime) {
                        event = *it;
                        break;
                    }
                }
            }
        } else {
            // Do nothing. EventNode should be returned
            // because ttime is after all events on branch
        }
    }

    return event;
}


// t1, t2 must be absolute time
template <typename EventType>
int BranchHistory<EventType>::getNumberOfEventsOnInterval(double t1, double t2)
{
    int n_events = 0;

    typename EventSet::iterator it;
    for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end(); ++it) {
        double atime = (*it)->getAbsoluteTime();
        if ((atime > t1) && (atime <= t2)) {
            n_events++;
        }
    }

    return n_events;
}


template <typename EventType>
EventType* BranchHistory<EventType>::getEventByIndexPosition(int index)
{
    EventSetSizeType i = static_cast<EventSetSizeType>(index);
    if (i < _eventsOnBranch.size() ) {
        typename EventSet::iterator it = _eventsOnBranch.begin();
        for (EventSetSizeType k = 0; k < i; k++)
            it++;
        return *it;
    } else {
        log(Error) << "BranchHistory::getEventByIndexPosition: "
                   << "accessing invalid event\n";
        std::exit(1);
    }

}


template <typename EventType>
void BranchHistory<EventType>::printBranchHistory()
{
    log() << "Node Event: " << _nodeEvent << "\t"
          << "Ancestor Event: " << _ancestralNodeEvent << "\n";

    EventSetSizeType numEvents = _eventsOnBranch.size();
    log() << "Number of events on branch: " << numEvents << "\n";

    typename EventSet::iterator it;
    for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end(); ++it) {
        printEvent(*it);
    }
}


template <typename EventType>
void BranchHistory<EventType>::reversePrintBranchHistory()
{
    typename EventSet::reverse_iterator it;
    for (it = _eventsOnBranch.rbegin(); it != _eventsOnBranch.rend(); ++it) {
        printEvent(*it);
    }
}


template <typename EventType>
void BranchHistory<EventType>::printEvent(EventType* event)
{
    log() << event << "\t\t"
          << event->getMapTime() << "\t\t"
          << event->getAbsoluteTime() << "\n";
}


template <typename EventType>
inline void BranchHistory<EventType>::setNodeEvent(EventType* x)
{
    _nodeEvent = x;
}


template <typename EventType>
inline EventType* BranchHistory<EventType>::getNodeEvent()
{
    return _nodeEvent;
}


template <typename EventType>
inline void BranchHistory<EventType>::setAncestralNodeEvent(EventType* x)
{
    _ancestralNodeEvent = x;
}


template <typename EventType>
inline EventType* BranchHistory<EventType>::getAncestralNodeEvent()
{
    return _ancestralNodeEvent;
}


template <typename EventType>
inline void BranchHistory<EventType>::popEventOffBranchHistory(EventType* x)
{
    _eventsOnBranch.erase(x);
}


template <typename EventType>
inline void BranchHistory<EventType>::addEventToBranchHistory(EventType* x)
{
    _eventsOnBranch.insert(x);
}


template <typename EventType>
inline int BranchHistory<EventType>::getNumberOfBranchEvents()
{
    return (int)_eventsOnBranch.size();
}


#endif
