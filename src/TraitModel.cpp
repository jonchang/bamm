// Undefining this macro constrains analysis to NEGATIVE values
// for the beta shift parameter
#define NEGATIVE_SHIFT_PARAM
//#undef NEGATIVE_SHIFT_PARAM

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm>

#include "TraitModel.h"
#include "TraitBranchEvent.h"
#include "BranchHistory.h"

#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "Settings.h"
#include "Log.h"


TraitModel::TraitModel(MbRandom* ranptr, Tree* tp, Settings* sp) :
    Model(ranptr, tp, sp), _lastLH(0.0)
{
    _updateBetaScale = _settings->getUpdateBetaScale();
    _updateBetaShiftScale = _settings->getUpdateBetaShiftScale();
    _updateNodeStateScale = _settings->getUpdateNodeStateScale();

    _betaInitPrior = _settings->getBetaInitPrior();
    _betaShiftPrior = _settings->getBetaShiftPrior();

    setMinMaxTraitPriors();

#ifdef NEGATIVE_SHIFT_PARAM
    // Constrain beta shift to be zero or less than zero.
    if (_settings->getBetaShiftInit() > 0){
        log(Error) << "Initial value of beta shift (betaShiftInit) cannot\n"
                   << "be positive. This parameter is constrained to negative\n"
                   << "values.\n";
        std::exit(1);
    }
#endif
    
    BranchEvent* x = new TraitBranchEvent
        (_settings->getBetaInit(), _settings->getBetaShiftInit(),
            _treePtr->getRoot(), _treePtr, _ran, 0);
    _rootEvent = x;
    _lastEventModified = x;

    TraitBranchEvent* rootEvent = static_cast<TraitBranchEvent*>(_rootEvent);
    log() << "\nRoot beta: " << rootEvent->getBetaInit() << "\t"
          << _settings->getBetaInit() << "\t"
          << "Shift: " << rootEvent->getBetaShift() << "\n";

    // Set NodeEvent of root node equal to the rootEvent
    tp->getRoot()->getBranchHistory()->setNodeEvent(rootEvent);

    // Initialize all branch histories to equal the root event:
    forwardSetBranchHistories(rootEvent);

    _treePtr->setMeanBranchTraitRates();

    if (_settings->getLoadEventData()) {
        log() << "\nLoading model data from file: "
              << _settings->getEventDataInfile() << "\n";
        initializeModelFromEventDataFile();
    }

    setCurrentLogLikelihood(computeLogLikelihood());

    log() << "\nInitial log-likelihood: " << getCurrentLogLikelihood() << "\n";
    if (_settings->getSampleFromPriorOnly())
        log() << "Note that you have chosen to sample from prior only.\n";

    // this parameter only set during model-jumping.
    _logQRatioJump = 0.0;
}


void TraitModel::readModelSpecificParameters(std::ifstream& infile)
{
    std::string tempString;

    getline(infile, tempString, '\t');
    _readBetaInit = atof(tempString.c_str());

    getline(infile, tempString, '\t');
    _readBetaShift = atof(tempString.c_str());
}


void TraitModel::setRootEventWithReadParameters()
{
    TraitBranchEvent* rootEvent = static_cast<TraitBranchEvent*>(_rootEvent);

    rootEvent->setBetaInit(_readBetaInit);
    rootEvent->setBetaShift(_readBetaShift);
}


BranchEvent* TraitModel::newBranchEventWithReadParameters(Node* x, double time)
{
    return new TraitBranchEvent(_readBetaInit, _readBetaShift, x,
        _treePtr, _ran, time);
}


void TraitModel::setMeanBranchParameters()
{
    _treePtr->setMeanBranchTraitRates();
}


// TODO: Originally this was the last statement in
// initializeModelFromEventDataFile(), is it needed in SpExModel, too?
// printEventData();


BranchEvent* TraitModel::newBranchEventWithRandomParameters(double x)
{
    // Sample beta and beta shift from prior
    double newbeta = _ran->exponentialRv(_settings->getBetaInitPrior());
    double newBetaShift = _ran->normalRv(0.0, _settings->getBetaShiftPrior());
    
#ifdef NEGATIVE_SHIFT_PARAM
    newBetaShift = -fabs(newBetaShift);
    double dens_term = log(2.0);
#else
    double dens_term = 0.0;
#endif
    
    _logQRatioJump = 0.0;
    _logQRatioJump += _ran->lnExponentialPdf(_betaInitPrior, newbeta);
    
    // Add log(2) [see dens_term above] because this is truncated
    // normal distribution constrained to negative values
    _logQRatioJump += dens_term +
        _ran->lnNormalPdf(0.0, _betaShiftPrior, newBetaShift);

    return new TraitBranchEvent(newbeta, newBetaShift,
        _treePtr->mapEventToTree(x), _treePtr, _ran, x);
}


void TraitModel::setDeletedEventParameters(BranchEvent* event)
{
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(event);

    _lastDeletedEventBetaInit = be->getBetaInit();
    _lastDeletedEventBetaShift = be->getBetaShift();
}


double TraitModel::calculateLogQRatioJump()
{
    double logQRatioJump = 0.0;

    logQRatioJump +=
        _ran->lnExponentialPdf(_betaInitPrior, _lastDeletedEventBetaInit);
    logQRatioJump +=
        _ran->lnNormalPdf(0.0, _betaShiftPrior, _lastDeletedEventBetaShift);

    return logQRatioJump;
}


BranchEvent* TraitModel::newBranchEventFromLastDeletedEvent()
{
    TraitBranchEvent* newEvent = new TraitBranchEvent(0.0, 0.0,
        _treePtr->mapEventToTree(_lastDeletedEventMapTime), _treePtr, _ran,
            _lastDeletedEventMapTime);

    newEvent->setBetaInit(_lastDeletedEventBetaInit);
    newEvent->setBetaShift(_lastDeletedEventBetaShift);

    return newEvent;
}


void TraitModel::updateBetaMH(void)
{
    int toUpdate = _ran->sampleInteger(0, (int)_eventCollection.size());

    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        EventSet::iterator it = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            ++it;

        be = static_cast<TraitBranchEvent*>(*it);
    }

    double oldRate = be->getBetaInit();
    double cterm = exp( _updateBetaScale * (_ran->uniformRv() - 0.5) );
    be->setBetaInit(cterm * oldRate);
    _treePtr->setMeanBranchTraitRates();

    double PropLnLik = computeLogLikelihood();

    double LogPriorRatio = _ran->lnExponentialPdf(_settings->getBetaInitPrior(),
                           be->getBetaInit());
    LogPriorRatio -= _ran->lnExponentialPdf(_settings->getBetaInitPrior(), oldRate);

    double LogProposalRatio = log(cterm);

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = likeRatio + LogPriorRatio + LogProposalRatio;

    const bool acceptMove = acceptMetropolisHastings(logHR);

//  std::cout << getGeneration() << "\tL1: " << startLH << "\tL2: " << getCurrentLogLikelihood() << std::endl;

    if (acceptMove == true) {
        //std::cout << "accept: " << oldRate << "\t" << be->getBetaInit() << std::endl;
        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {

        // revert to previous state
        _lastLH = PropLnLik;


        be->setBetaInit(oldRate);
        _treePtr->setMeanBranchTraitRates();
        _acceptLast = 0;
        _rejectCount++;

    }

    /*if (!acceptMove){
        std::cout << std::endl;
        std::cout << startLL << "\tCurr: " << getCurrentLogLikelihood() << "\tcalc: " << computeLogLikelihood() << std::endl;
    }*/

    incrementGeneration();
}


void TraitModel::updateBetaShiftMH(void)
{
    int toUpdate = _ran->sampleInteger(0, (int)_eventCollection.size());

    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        EventSet::iterator it = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            ++it;

        be = static_cast<TraitBranchEvent*>(*it);
    }

    double oldShift = be->getBetaShift();
    double newShift = oldShift + _ran->normalRv((double)0.0, _updateBetaShiftScale);
 
    // Convert to negative via reflection:
    newShift = -fabs(newShift);
 
    be->setBetaShift(newShift);
    _treePtr->setMeanBranchTraitRates();

    double PropLnLik = computeLogLikelihood();

    // Normal prior on shift parameter:
    // We ignore the log(2) term which cancels out.
    
    double LogPriorRatio = _ran->lnNormalPdf((double)0.0,
                                            _settings->getBetaShiftPrior(), newShift);
    LogPriorRatio -= _ran->lnNormalPdf((double)0.0, _settings->getBetaShiftPrior(),
                                      oldShift);


    double LogProposalRatio = 0.0;

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = likeRatio + LogPriorRatio + LogProposalRatio;

    const bool acceptMove = acceptMetropolisHastings(logHR);

    if (acceptMove == true) {

        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {

        // revert to previous state
        be->setBetaShift(oldShift);
        _treePtr->setMeanBranchTraitRates();
        _acceptLast = 0;
        _rejectCount++;

    }

    incrementGeneration();
}


void TraitModel::updateNodeStateMH(void)
{
    Node* xnode = _treePtr->chooseInternalNodeAtRandom();

    double oldTriadLogLik = computeTriadLikelihoodTraits(xnode);

    double oldstate = xnode->getTraitValue();
    double newstate = oldstate + _ran->uniformRv((-1.0 *
                      _settings->getUpdateNodeStateScale()), _settings->getUpdateNodeStateScale());
    xnode->setTraitValue(newstate);

    // Compute triad likelihood
    double newTriadLoglik = computeTriadLikelihoodTraits(xnode);
    double PropLnLik = getCurrentLogLikelihood() - oldTriadLogLik + newTriadLoglik;

    // set to zero for now...flat prior, unifor (see below).
    double LogPriorRatio = 0.0;

    double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = LogPriorRatio + logProposalRatio + likeRatio;
    bool acceptMove = acceptMetropolisHastings(logHR);

    // Here we do prior calculation to avoid computing infinity...
    if (newstate > _settings->getTraitPriorMax() ||
            newstate < _settings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {
        xnode->setTraitValue(oldstate);

        _rejectCount++;
        _acceptLast = 0;
    }

    incrementGeneration();
}


void TraitModel::updateNodeStateMH(Node* xnode)
{
    double oldTriadLogLik = computeTriadLikelihoodTraits(xnode);

    double oldstate = xnode->getTraitValue();
    double newstate = oldstate + _ran->uniformRv((-1.0 *
                      _settings->getUpdateNodeStateScale()), _settings->getUpdateNodeStateScale());
    xnode->setTraitValue(newstate);

    // Compute triad likelihood
    double newTriadLoglik = computeTriadLikelihoodTraits(xnode);
    double PropLnLik = getCurrentLogLikelihood() - oldTriadLogLik + newTriadLoglik;

    // set to zero for now...flat prior, unifor (see below).
    double LogPriorRatio = 0.0;

    double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = LogPriorRatio + logProposalRatio + likeRatio;
    bool acceptMove = acceptMetropolisHastings(logHR);

    // Here we do prior calculation to avoid computing infinity...
    if (newstate > _settings->getTraitPriorMax() ||
            newstate < _settings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {
        xnode->setTraitValue(oldstate);

        _rejectCount++;
        _acceptLast = 0;
    }

    incrementGeneration();
}


void TraitModel::updateDownstreamNodeStatesMH(Node* xnode)
{
    // Get list of internal node descendants from node * x
    // update each (or some fraction thereof).

    _treePtr->setTempInternalNodeArray(xnode);
    for (int i = 0; i < 100; i++)
        updateNodeStateMH(_treePtr->getRandomNodeFromTempArray());

    _treePtr->clearTempNodeArray();
}


double TraitModel::computeSpecificLogLikelihood(void)
{
    double LnL = 0.0;

    //Node * tmpnode = _treePtr->getRoot()->getLfDesc();

    if (_settings->getSampleFromPriorOnly())
        return 0.0;

    int numNodes = _treePtr->getNumberOfNodes();

    // iterate over non-root nodes and compute LnL

    for (int i = 0; i < numNodes; i++) {
        Node* xnode = _treePtr->getNodeFromDownpassSeq(i);
        if ( (xnode != _treePtr->getRoot()) && (xnode->getCanHoldEvent() == true) ) {
            double var = xnode->getBrlen() * xnode->getMeanBeta();

            // change in phenotype:
            double delta = xnode->getTraitValue() - xnode->getAnc()->getTraitValue();

            LnL += _ran->lnNormalPdf(0, var, delta);

        }

    }

    return LnL;
}


double TraitModel::computeTriadLikelihoodTraits(Node* x)
{
    if (_settings->getSampleFromPriorOnly())
        return 0.0;

    double logL = 0.0;

    // Can only use this likelihood if node contributes to
    // likelihood of observed data

    if (x->getCanHoldEvent() == true) {

        // computation for left descendant branch:

        if (x->getLfDesc()->getCanHoldEvent() == true) {
            double delta = x->getLfDesc()->getTraitValue() - x->getTraitValue();
            logL += _ran->lnNormalPdf(0,
                                     (x->getLfDesc()->getBrlen() * x->getLfDesc()->getMeanBeta()), delta);
        }


        if (x->getRtDesc()->getCanHoldEvent() == true) {
            // computation for right descendant branch
            double delta = x->getRtDesc()->getTraitValue() - x->getTraitValue();
            logL += _ran->lnNormalPdf(0,
                                     (x->getRtDesc()->getBrlen() * x->getRtDesc()->getMeanBeta()), delta);


        }


        // computation for ancestral branch (unless == root)

        if (x != _treePtr->getRoot()) {

            double delta = x->getTraitValue() - x->getAnc()->getTraitValue();
            logL += _ran->lnNormalPdf(0, (x->getBrlen() * x->getMeanBeta()), delta);

        }
    }

    return logL;
}


double TraitModel::computeSpecificLogPrior(void)
{
#ifdef NEGATIVE_SHIFT_PARAM
    double dens_term = log(2.0);
#else
    double dens_term = 0.0;
#endif

    TraitBranchEvent* rootEvent = static_cast<TraitBranchEvent*>(_rootEvent);
    
    double logPrior = 0.0;
    logPrior += _ran->lnExponentialPdf(_settings->getBetaInitPrior(),
                                      rootEvent->getBetaInit());
    logPrior += dens_term + _ran->lnNormalPdf((double)0.0, _settings->getBetaShiftPrior(),
                                 rootEvent->getBetaShift());
    EventSet::iterator i;
    for (i = _eventCollection.begin(); i != _eventCollection.end(); ++i) {
        TraitBranchEvent* event = static_cast<TraitBranchEvent*>(*i);
        logPrior += _ran->lnExponentialPdf(_settings->getBetaInitPrior(),
            event->getBetaInit());

        // Add log(2) to make truncated normal density but only if NEGATIVE_SHIFT_PARAM
        logPrior += dens_term + _ran->lnNormalPdf((double)0.0, _settings->getBetaShiftPrior(),
                                     event->getBetaShift());
    }

    // and prior on number of events:

    logPrior += _ran->lnExponentialPdf(_poissonRatePrior , getEventRate());

    return logPrior;
}


void TraitModel::getSpecificEventDataString
    (std::stringstream& ss, BranchEvent* event)
{
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(event);

    ss << be->getBetaInit() << ","
       << be->getBetaShift();
}


void TraitModel::setMinMaxTraitPriors(void)
{
    int nnodes = _treePtr->getNumberOfNodes();
    std::vector<double> tvec;
    for (int i = 0; i < nnodes; i++) {
        Node* xnode = _treePtr->getNodeFromDownpassSeq(i);
        if (xnode->getTraitValue() != 0)
            tvec.push_back(xnode->getTraitValue());
    }

    std::sort(tvec.begin(), tvec.end());

    // Default here will be to use observed range +/- 20%
    double rg = tvec[(tvec.size() - 1)] - tvec[0];
    double minprior = tvec[0] - (0.2 * rg);
    double maxprior = tvec[(tvec.size() - 1)] + (0.2 * rg);

    log() << "\nMin and max phenotype limits set using observed data: \n"
          << "\t\tMin: " << minprior << "\tMax: " << maxprior << "\n";
    _settings->setTraitPriorMin(minprior);
    _settings->setTraitPriorMax(maxprior);
}
