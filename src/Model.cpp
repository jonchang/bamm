#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>

#include "Model.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "BranchHistory.h"
#include "Settings.h"
#include "Log.h"
// TODO: Do I need all these includes?


double Model::mhColdness = 1.0;


Model::Model(MbRandom* ranptr, Tree* tp, Settings* sp) :
    _ran(ranptr), _treePtr(tp), _settings(sp)
{
    // Reduce weird autocorrelation of values at start by calling RNG
    // a few times. TODO: Why is there a weird autocorrelation?
    for (int i = 0; i < 100; i++)
        ranptr->uniformRv();

    _treePtr->getTotalMapLength(); // total map length (required to set priors)

    // Set parameter values for model object, including priors etc.

    _gen = 0;

    // Scale for event moves on tree.
    _scale = _settings->getUpdateEventLocationScale();

    _updateEventRateScale = _settings->getUpdateEventRateScale();
    _localGlobalMoveRatio =
        _settings->getLocalGlobalMoveRatio();
    
    _poissonRatePrior = _settings->getPoissonRatePrior();
    _eventRate = 1 / _poissonRatePrior;

    _acceptCount = 0;
    _rejectCount = 0;
    _acceptLast = -1;
}


Model::~Model()
{
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it)
        delete *it;
}


void Model::initializeModelFromEventDataFile(void)
{
    std::ifstream infile(_settings->getEventDataInfile().c_str());

    if (!infile.good()) {
        log(Error) << "Bad Filename.\n";
        exit(1);
    }

    log() << "Initializing model from <"
          << _settings->getEventDataInfile() << ">\n";

    std::string species1;
    std::string species2;
    std::string etimeStr;
    double etime;

    int count = 0;
    while (infile) {
        getline(infile, species1, '\t');
        getline(infile, species2, '\t');
        getline(infile, etimeStr, '\t');
        etime = std::atof(etimeStr.c_str());

        // Read the derived class's parameters
        readModelSpecificParameters(infile);

        Node* x = NULL;

        if ((species1 != "NA") && (species2 != "NA")) {
            x = _treePtr->getNodeMRCA(species1.c_str(), species2.c_str());
        } else if ((species1 == "NA") && (species2 == "NA")) {
            x = _treePtr->getNodeByName(species1.c_str());
        } else {
            log(Error) << "One of the species is NA.\n";
            std::exit(1);
        }

        if (x == _treePtr->getRoot()) {
            setRootEventWithReadParameters();
        } else {
            double deltaT = x->getTime() - etime;
            double newmaptime = x->getMapStart() + deltaT;

            BranchEvent *newEvent =
                newBranchEventWithReadParameters(x, newmaptime);
            newEvent->getEventNode()->getBranchHistory()->
                addEventToBranchHistory(newEvent);

            _eventCollection.insert(newEvent);
            forwardSetBranchHistories(newEvent);
            setMeanBranchParameters();
        }

        count++;
    }

    infile.close();

    log() << "Read a total of " << count << " events.\n";
    log() << "Added " << _eventCollection.size()
          << " pre-defined events to tree, plus root event.\n";
}


// Adds event to tree based on reference map value
// - Adds to branch history set
// - Inserts into _eventCollection

void Model::addEventToTree(double x)
{
    BranchEvent* newEvent = newBranchEventWithRandomParameters(x);

    // Add the event to the branch history.
    // Always done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(newEvent);

    _eventCollection.insert(newEvent);
    forwardSetBranchHistories(newEvent);
    setMeanBranchParameters();

    _lastEventModified = newEvent;
}


// Adds event to tree based on uniform RV
// - Adds to branch history set
// - Inserts into _eventCollection

void Model::addEventToTree()
{
    double aa = _treePtr->getRoot()->getMapStart();
    double bb = _treePtr->getTotalMapLength();
    double x = _ran->uniformRv(aa, bb);

    addEventToTree(x);
}


void Model::printEvents(void)
{
    int numEvents = (int)_eventCollection.size();
    log() << "Number of events: " << numEvents << "\n";

    int counter = 1;
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        log() << "Event: " << counter++ << "\t"
              << "Address: " << *it << "\t" << (*it)->getMapTime() << "\t"
              << "Node: " << (*it)->getEventNode() << "\n";
    }
}


BranchEvent* Model::chooseEventAtRandom()
{
    int numEvents = (int)_eventCollection.size();

    if (numEvents == 0) {
        log(Error) << "Number of events is zero.\n";
        std::exit(1);
    }

    int ctr = 0;
    double xx = _ran->uniformRv();
    int chosen = (int)(xx * (double)numEvents);

    EventSet::iterator sit = _eventCollection.begin();

    for (int i = 0; i < chosen; i++) {
        ++sit;
        ctr++;
    }

    return *sit;
}


void Model::eventLocalMove(void)
{
    eventMove(true);
}


void Model::eventGlobalMove(void)
{
    eventMove(false);
}


// If events are on tree: choose event at random,
// move locally (or globally) and forward set branch histories etc.
// Should also store previous event information to revert to previous

// If parameter local == true, does a local move;
// otherwise, it does a global move

void Model::eventMove(bool local)
{
    if (getNumberOfEvents() > 0) {
        // The event to be moved
        BranchEvent* chosenEvent = chooseEventAtRandom();

        // This is the event preceding the chosen event:
        // histories should be set forward from here.
        BranchEvent* previousEvent = chosenEvent->getEventNode()->
            getBranchHistory()->getLastEvent(chosenEvent);

        // Set this history variable in case move is rejected
        _lastEventModified = chosenEvent;

        chosenEvent->getEventNode()->getBranchHistory()->
            popEventOffBranchHistory(chosenEvent);

        if (local) {
            // Get step size for move
            double step = _ran->uniformRv(0, _scale) - 0.5 * _scale;
            chosenEvent->moveEventLocal(step);
        } else {
            chosenEvent->moveEventGlobal();
        }

        chosenEvent->getEventNode()->getBranchHistory()->
            addEventToBranchHistory(chosenEvent);

        // Get last event from the theEventNode, forward set its history.
        // Then go to the "moved" event and forward set its history.

        forwardSetBranchHistories(previousEvent);
        forwardSetBranchHistories(chosenEvent);
    }

    setMeanBranchParameters();
}


// Used to reset position of event if move is rejected

void Model::revertMovedEventToPrevious()
{
    // Get last event from position of event to be removed
    BranchEvent* newLastEvent = _lastEventModified->getEventNode()->
        getBranchHistory()->getLastEvent(_lastEventModified);

    // Pop event off its new position
    _lastEventModified->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_lastEventModified);

    // Reset nodeptr, reset mapTime
    _lastEventModified->revertOldMapPosition();

    // Now reset forward from _lastEventModified (new position)
    // and from newLastEvent, which holds 'last' event before old position
    _lastEventModified->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(_lastEventModified);

    // Forward set from new position
    forwardSetBranchHistories(newLastEvent);

    // Forward set from event immediately rootwards from previous position
    forwardSetBranchHistories(_lastEventModified);

    // Set _lastEventModified to NULL because it has already been reset.
    // Future implementations should check whether this is NULL
    // before attempting to use it to set event

    _lastEventModified = NULL;

    setMeanBranchParameters();
}


// Recursively count the number of events in the branch histories
int Model::countEventsInBranchHistory(Node* p)
{
    int count = p->getBranchHistory()->getNumberOfBranchEvents();

    if (p->getLfDesc() != NULL) {
        count += countEventsInBranchHistory(p->getLfDesc());
    }

    if (p->getRtDesc() != NULL){
        count += countEventsInBranchHistory(p->getRtDesc());
    }

    return count;
}


void Model::deleteEventFromTree(BranchEvent* be)
{
    if (be == _rootEvent) {
        log(Error) << "Can't delete root event.\n";
        exit(1);
    }

    // Erase from branch history:
    Node* currNode = be->getEventNode();

    // Get event downstream of i
    BranchEvent* newLastEvent = currNode->getBranchHistory()->getLastEvent(be);

    _lastDeletedEventMapTime = be->getMapTime();

    setDeletedEventParameters(be);
    _logQRatioJump = calculateLogQRatioJump();

    currNode->getBranchHistory()->popEventOffBranchHistory(be);

    _eventCollection.erase(be);
    delete be;

    forwardSetBranchHistories(newLastEvent);

    setMeanBranchParameters();
}


void Model::deleteRandomEventFromTree()
{
    int numEvents = (int)_eventCollection.size();
    
    // Can only delete event if more than root node present.
    if (numEvents == 0) {
        return;
    }

    int counter = 0;
    double xx = _ran->uniformRv();
    int chosen = (int)(xx * (double)numEvents);

    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        if (counter++ == chosen) {
            deleteEventFromTree(*it);
            break;
        }
    }
}


void Model::restoreLastDeletedEvent()
{
    // Use constructor for speciation and extinction
    BranchEvent* newEvent = newBranchEventFromLastDeletedEvent();

    // Add the event to the branch history.
    // Always done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(newEvent);

    _eventCollection.insert(newEvent);

    // Event is now inserted into branch history;
    // however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    setMeanBranchParameters();
}


void Model::changeNumberOfEventsMH(void)
{
    int numEvents = (int)_eventCollection.size();
    bool gain = (numEvents > 0) ? (_ran->uniformRv() < 0.5) : true;

    if (gain) {
        addEventMH();
    } else {
        removeEventMH();
    }

    incrementGeneration();
}


void Model::addEventMH()
{
    double oldLogLikelihood = getCurrentLogLikelihood();
    double oldLogPrior = computeLogPrior();

    int currState = (int)_eventCollection.size();

    int K = currState;
    double qRatio = (K > 0) ? 1.0 : 0.5;

    addEventToTree();

    setMeanBranchParameters();

    double logLikelihood = computeLogLikelihood();
    double logPrior = computeLogPrior();

    // Prior ratio is eventRate / (k + 1)
    // but now, _eventCollection.size() == k + 1
    // because event has already been added.
    // Here HR is just the prior ratio

    double logHR = computeEventGainLogHR(K, logLikelihood, oldLogLikelihood,
        logPrior, oldLogPrior, qRatio);

    bool acceptMove = false;
    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    if (!std::isinf(logLikelihoodRatio)) {  // TODO: When is this ever inifinte?
        acceptMove = acceptMetropolisHastings(logHR);
    }

    bool isValidConfig = isEventConfigurationValid(_lastEventModified);
    if (acceptMove && isValidConfig) {
        setCurrentLogLikelihood(logLikelihood);
        _acceptCount++;
        _acceptLast = 1;
    } else {
        deleteEventFromTree(_lastEventModified);
        setMeanBranchParameters();
        _rejectCount++;
        _acceptLast = 0;
    }
}


double Model::computeEventGainLogHR(double K, double logLikelihood,
    double oldLogLikelihood, double logPrior, double oldLogPrior, double qRatio)
{
    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    double logPriorRatio = logPrior - oldLogPrior;

    double logHR = std::log(_eventRate) - std::log(K + 1.0);
    logHR += logLikelihoodRatio;
    logHR += logPriorRatio;
    logHR += std::log(qRatio);

    // Now for jumping density of the bijection between parameter spaces
    logHR -= _logQRatioJump;

    return logHR;
}


void Model::removeEventMH()
{
    double oldLogLikelihood = getCurrentLogLikelihood();
    double oldLogPrior = computeLogPrior();

    int currState = (int)_eventCollection.size();

    int K = currState;
    double qRatio = (K != 1) ? 1.0 : 2.0;

    deleteRandomEventFromTree();

    setMeanBranchParameters();

    double logLikelihood = computeLogLikelihood();
    double logPrior = computeLogPrior();

    // This is probability of going from k to k - 1
    // So, prior ratio is (k / eventRate)

    // First get prior ratio:
    double logHR = computeEventLossLogHR(K, logLikelihood, oldLogLikelihood,
            logPrior, oldLogPrior, qRatio);

    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    bool acceptMove = false;
    if (!std::isinf(logLikelihoodRatio)) {
        acceptMove = acceptMetropolisHastings(logHR);
    }

    if (acceptMove) {
        setCurrentLogLikelihood(logLikelihood);
        _acceptCount++;
        _acceptLast = 1;
    } else {
        restoreLastDeletedEvent();
        _rejectCount++;
        _acceptLast = 0;
    }
}


double Model::computeEventLossLogHR(double K, double logLikelihood,
    double oldLogLikelihood, double logPrior, double oldLogPrior, double qRatio)
{
    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    double logPriorRatio = logPrior - oldLogPrior;

    // This is probability of going from K to K - 1
    // So, prior ratio is (K / eventRate)

    double logHR = std::log(K) - std::log(_eventRate);
    logHR += logLikelihoodRatio;
    logHR += logPriorRatio;
    logHR += std::log(qRatio);

    // Now for jumping density of the bijection between parameter spaces
    logHR += _logQRatioJump;

    return logHR;
}


void Model::moveEventMH(void)
{
    // Consider proposal rejected (can't move nonexistent event)
    if (_eventCollection.size() == 0) {
        _rejectCount++;
        _acceptLast = 0;
        incrementGeneration();
        return;
    }

    double localMoveProb = _localGlobalMoveRatio / (1 + _localGlobalMoveRatio);

    bool isLocalMove = _ran->uniformRv() < localMoveProb;

    if (isLocalMove) {
        // Local move, with event drawn at random
        eventLocalMove();
    } else {
        eventGlobalMove();
    }

    setMeanBranchParameters();

    double logLikelihood = computeLogLikelihood();
    double logLikelihoodRatio = logLikelihood - getCurrentLogLikelihood();
    double logHR = logLikelihood;

    bool isValid = isEventConfigurationValid(_lastEventModified);

    bool acceptMove = false;
    if (!std::isinf(logLikelihoodRatio) && isValid) {
        acceptMove = acceptMetropolisHastings(logHR);
    }

    if (acceptMove) {
        setCurrentLogLikelihood(logLikelihood);
        _acceptCount++;
        _acceptLast = 1;
    } else {
        revertMovedEventToPrevious();
        setMeanBranchParameters();
        _rejectCount++;
        _acceptLast = 0;
    }

    incrementGeneration();
}


// Select an event at random.
// If partition is time-constant, flip state to time-variable.
// If partition is time-variable, flip state to time-constant
// TODO: This function does nothing with be

void Model::updateTimeVariablePartitionsMH()
{
    int toUpdate = _ran->sampleInteger(0, (int)_eventCollection.size());
    BranchEvent* be = _rootEvent;

    if (toUpdate > 0) {
        EventSet::iterator it = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            ++it;
        be = *it;
    } else {
        // Event remains as root event
    }
}


// Metropolis-Hastings step to update Poisson event rate.
// Note that changing this rate does not affect the likelihood,
// so the priors and qratio determine acceptance rate.

void Model::updateEventRateMH()
{
    double oldEventRate = getEventRate();

    double cterm = exp( _updateEventRateScale * (_ran->uniformRv() - 0.5) );
    setEventRate(cterm * oldEventRate);

    double LogPriorRatio =
        _ran->lnExponentialPdf(_poissonRatePrior, getEventRate()) -
        _ran->lnExponentialPdf(_poissonRatePrior, oldEventRate);
    double logProposalRatio = std::log(cterm);
    double logHR = LogPriorRatio + logProposalRatio;
    const bool acceptMove = acceptMetropolisHastings(logHR);

    if (acceptMove) {
        _acceptCount++;
        _acceptLast = 1;
    } else {
        setEventRate(oldEventRate);
        _rejectCount++;
        _acceptLast = 0;
    }

    incrementGeneration();
}


double Model::computeLogLikelihood()
{
    return computeSpecificLogLikelihood();
}


double Model::computeLogPrior()
{
    return computeSpecificLogPrior();
}


bool Model::acceptMetropolisHastings(double lnR)
{
    double r = safeExponentiation(Model::mhColdness * lnR);
    return _ran->uniformRv() < r;
}


void Model::initializeBranchHistories(Node* x)
{
    x->getBranchHistory()->setNodeEvent(_rootEvent);

    if (x->getAnc() != NULL)
        x->getBranchHistory()->setAncestralNodeEvent(_rootEvent);

    if (x->getLfDesc() != NULL)
        initializeBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        initializeBranchHistories(x->getRtDesc());
}



void Model::printStartAndEndEventStatesForBranch(Node* x)
{
    if (x != _treePtr->getRoot())
        log() << "Node: " << x << "\t"
              << "Anc: "
              << x->getBranchHistory()->getAncestralNodeEvent() << "\t"
              << "Event: " << x->getBranchHistory()->getNodeEvent() << "\n";

    if (x->getLfDesc() != NULL)
        printStartAndEndEventStatesForBranch(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        printStartAndEndEventStatesForBranch(x->getRtDesc());
}


/*
    If this works correctly, this will take care of the following:
    1. if a new event is created or added to tree,
    this will forward set all branch histories from the insertion point

    2. If an event is deleted, you find the next event rootwards,
    and call forwardSetBranchHistories from that point. It will replace
    settings due to the deleted node with the next rootwards node.
*/

void Model::forwardSetBranchHistories(BranchEvent* x)
{
    // If there is another event occurring more recent (closer to tips)
    // do nothing. Even just sits in BranchHistory but doesn't affect
    // state of any other nodes.

    // This seems circular, but what else to do?
    // Given an event (which references the node defining the branch on which
    // event occurs) you get the corresponding branch history and the last
    // event since the events will have been inserted in the correct order.

    Node* myNode = x->getEventNode();

    if (x == _rootEvent) {
        forwardSetHistoriesRecursive(myNode->getLfDesc());
        forwardSetHistoriesRecursive(myNode->getRtDesc());
    } else if (x == myNode->getBranchHistory()->getLastEvent()) {
        // If TRUE, x is the most tip-wise event on branch.
        myNode->getBranchHistory()->setNodeEvent(x);

        // if myNode is not a tip:
        if (myNode->getLfDesc() != NULL && myNode->getRtDesc() != NULL) {
            forwardSetHistoriesRecursive(myNode->getLfDesc());
            forwardSetHistoriesRecursive(myNode->getRtDesc());
        }
        // else: node is a tip : do nothing.
    }
    //else: there is another more tipwise event on same branch; do nothing
}


void Model::forwardSetHistoriesRecursive(Node* p)
{
    // Get event that characterizes parent node
    BranchEvent* lastEvent = p->getAnc()->getBranchHistory()->getNodeEvent();

    // Set the ancestor equal to the event state of parent node:
    p->getBranchHistory()->setAncestralNodeEvent(lastEvent);

    // If no events on the branch, go down to descendants and do same thing
    // otherwise, process terminates (because it hits another event on branch
    if (p->getBranchHistory()->getNumberOfBranchEvents() == 0) {
        p->getBranchHistory()->setNodeEvent(lastEvent);
        if (p->getLfDesc() != NULL)
            forwardSetHistoriesRecursive(p->getLfDesc());
        if (p->getRtDesc() != NULL)
            forwardSetHistoriesRecursive(p->getRtDesc());
    }
}


void Model::printBranchHistories(Node* x)
{
    if (x != _treePtr->getRoot()) {
        log() << "Node: " << x << "\t"
              << "Number of events: "
              << x->getBranchHistory()->getNumberOfBranchEvents() << "\t"
              << "Start: "
              << x->getBranchHistory()->getAncestralNodeEvent() << "\t"
              << "End: " << x->getBranchHistory()->getNodeEvent() << "\n";
    }

    if (x->getLfDesc() != NULL)
        printBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        printBranchHistories(x->getRtDesc());
}


double Model::getMHacceptanceRate()
{
    double arate = (double)_acceptCount / (_acceptCount + _rejectCount);
    return arate;
}


void Model::resetMHacceptanceParameters()
{
    _acceptCount = 0;
    _rejectCount = 0;
}


BranchEvent* Model::getEventByIndex(int x)
{
    EventSet::iterator it = _eventCollection.begin();
    for (int i = 0; i <= x; i++)
        ++it;

    return *it;
}


// Counts number of time-varying rate partitions

int Model::countTimeVaryingRatePartitions()
{
    int count = 0;
    count += (int)_rootEvent->getIsEventTimeVariable();
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it)
        count += (int)(*it)->getIsEventTimeVariable();
    return count;
}


//  Write event data to file for all events "on" tree
//  at a given point in the MCMC chain

void Model::getEventDataString(std::stringstream& ss)
{
    ss << getGeneration() << ",";

    BranchEvent* be = _rootEvent;
    Node* xl = _treePtr->getRoot()->getRandomLeftTipNode();
    Node* xr = _treePtr->getRoot()->getRandomRightTipNode();
    ss << xl->getName() << ","
       << xr->getName() << ","
       << be->getAbsoluteTime() << ",";

    // Implemented in derived class
    getSpecificEventDataString(ss, be);

    if (_eventCollection.size() > 0) {
        EventSet::iterator it;
        for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
            ss << "\n" << getGeneration() << ",";
            be = *it;
            if (be->getEventNode()->getLfDesc() == NULL) {
                ss << be->getEventNode()->getName() << "," << "NA" << ",";
            } else {
                Node* xl = be->getEventNode()->getRandomLeftTipNode();
                Node* xr = be->getEventNode()->getRandomRightTipNode();
                ss << xl->getName() << "," << xr->getName() << ",";
            }
            ss << be->getAbsoluteTime() << ",";

            // Implemented in derived class
            getSpecificEventDataString(ss, be);
        }
    }
}


bool Model::isEventConfigurationValid(BranchEvent* be)
{
    bool isValidConfig = false;

    if (be->getEventNode() == _treePtr->getRoot()) {
        Node* rt = _treePtr->getRoot()->getRtDesc();
        Node* lf = _treePtr->getRoot()->getLfDesc();
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0 &&
            lf->getBranchHistory()->getNumberOfBranchEvents() > 0) {
            // Events on both descendants of root. This fails.
            isValidConfig = false;
        } else
            isValidConfig = true;

    } else {
        int badsum = 0;

        Node* anc = be->getEventNode()->getAnc();
        Node* lf = anc->getLfDesc();
        Node* rt = anc->getRtDesc();

        // Test ancestor for events on branch

        if (anc == _treePtr->getRoot())
            badsum++;
        else if (anc->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;
        else {
            // nothing
        }

        // Test lf desc
        if (lf->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        // Test rt desc
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        if (badsum == 3)
            isValidConfig = false;
        else if (badsum < 3)
            isValidConfig = true;
        else {
            log(Error) << "Problem in Model::isEventConfigurationValid\n";
            exit(1);
        }
    }

    return isValidConfig;
}


void Model::debugLHcalculation(void)
{
    log() << "This does not currently support anything\n";
}


double safeExponentiation(double x)
{
    if (x > 0.0)
        return 1.0;
    else if (x < -100.0)
        return 0.0;
    else
        return exp(x);
}
