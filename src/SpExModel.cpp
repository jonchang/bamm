#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>

#include "SpExModel.h"
#include "BranchEvent.h"
#include "SpExBranchEvent.h"
#include "BranchHistory.h"

#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"

#include "Settings.h"
#include "Log.h"


#define MAX_E_PROB 0.999
#define JUMP_VARIANCE_NORMAL 0.05


SpExModel::SpExModel(MbRandom* ranptr, Tree* tp, Settings* sp) :
    Model(ranptr, tp, sp)
{
    // Initialize MCMC proposal/tuning parameters
    _updateLambdaInitScale = _settings->getUpdateLambdaInitScale();
    _updateMuInitScale = _settings->getUpdateMuInitScale();
    _updateLambdaShiftScale = _settings->getUpdateLambdaShiftScale();
    _updateMuShiftScale = _settings->getUpdateMuShiftScale();

    // Initial values
    _lambdaInit0 = _settings->getLambdaInit0();
    _lambdaShift0 = _settings->getLambdaShift0();
    _muInit0 = _settings->getMuInit0();
    _muShift0 = _settings->getMuShift0();

    // For Poisson process
    _lambdaInitPrior = _settings->getLambdaInitPrior();
    _lambdaShiftPrior = _settings->getLambdaShiftPrior();
    _muInitPrior = _settings->getMuInitPrior();
    _muShiftPrior = _settings->getMuShiftPrior();

    // Parameter for splitting branch into pieces for numerical computation
    _segLength = _settings->getSegLength();

    BranchEvent* x = new SpExBranchEvent
        (_lambdaInit0, _lambdaShift0, _muInit0, _muShift0,
            _treePtr->getRoot(), _treePtr, _ran, 0);
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the _rootEvent
    tp->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event
    forwardSetBranchHistories(_rootEvent);

    _treePtr->setMeanBranchSpeciation();
    _treePtr->setMeanBranchExtinction();

    // Initialize by previous event histories.
    if (_settings->getLoadEventData()) {
        log() << "\nLoading model data from file.\n";
        initializeModelFromEventDataFile();
    }
	
    setCurrentLogLikelihood(computeLogLikelihood());

    log() << "\nInitial log-likelihood: " << getCurrentLogLikelihood() << "\n";
    if (_settings->getSampleFromPriorOnly())
        log() << "Note that you have chosen to sample from prior only.\n";

    // This parameter only set during model-jumping.
    _logQRatioJump = 0.0;
}


void SpExModel::readModelSpecificParameters(std::ifstream& infile)
{
    std::string tempString;

    getline(infile, tempString, '\t');
    _readLambdaInit = atof(tempString.c_str());

    getline(infile, tempString, '\t');
    _readLambdaShift = atof(tempString.c_str());

    getline(infile, tempString, '\t');
    _readMuInit = atof(tempString.c_str());

    getline(infile, tempString, '\n');
    _readMuShift = atof(tempString.c_str());
}


void SpExModel::setRootEventWithReadParameters()
{
    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    rootEvent->setLamInit(_readLambdaInit);
    rootEvent->setLamShift(_readLambdaShift);
    rootEvent->setMuInit(_readMuInit);
    rootEvent->setMuShift(_readMuShift);
}


BranchEvent* SpExModel::newBranchEventWithReadParameters(Node* x, double time)
{
    return new SpExBranchEvent(_readLambdaInit, _readLambdaShift,
        _readMuInit, _readMuShift, x, _treePtr, _ran, time);
}


void SpExModel::setMeanBranchParameters()
{
    _treePtr->setMeanBranchSpeciation();
    _treePtr->setMeanBranchExtinction();
}


BranchEvent* SpExModel::newBranchEventWithRandomParameters(double x)
{
    double newLam = _ran->exponentialRv(_lambdaInitPrior) ;
    double newLambdaShift = _ran->normalRv(0.0, _lambdaShiftPrior);
    double newMu = _ran->exponentialRv(_muInitPrior);
    double newMuShift = 0.0;

    _logQRatioJump = 0.0;    // Set to zero to clear previous values
    _logQRatioJump = _ran->lnExponentialPdf(_lambdaInitPrior, newLam);
    _logQRatioJump += _ran->lnExponentialPdf(_muInitPrior, newMu);
    _logQRatioJump += _ran->lnNormalPdf
        (0.0, _lambdaShiftPrior, newLambdaShift);
    _logQRatioJump += _ran->lnNormalPdf(0.0, _muShiftPrior, newMuShift);

    return new SpExBranchEvent(newLam, newLambdaShift, newMu,
        newMuShift, _treePtr->mapEventToTree(x), _treePtr, _ran, x);
}


void SpExModel::setDeletedEventParameters(BranchEvent* event)
{
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(event);

    _lastDeletedEventLambdaInit = be->getLamInit();
    _lastDeletedEventLambdaShift = be->getLamShift();
    _lastDeletedEventMuInit = be->getMuInit();
    _lastDeletedEventMuShift = be->getMuShift();
}


double SpExModel::calculateLogQRatioJump()
{
    double logQRatioJump = 0.0;    // Set to zero to clear previous values

    logQRatioJump +=
        _ran->lnExponentialPdf(_lambdaInitPrior, _lastDeletedEventLambdaInit);
    logQRatioJump +=
        _ran->lnExponentialPdf(_muInitPrior, _lastDeletedEventMuInit);
    logQRatioJump +=
        _ran->lnNormalPdf(0.0, _lambdaShiftPrior, _lastDeletedEventLambdaShift);
    logQRatioJump +=
        _ran->lnNormalPdf(0.0, _muShiftPrior, _lastDeletedEventMuShift);

    return logQRatioJump;
}


BranchEvent* SpExModel::newBranchEventFromLastDeletedEvent()
{
    SpExBranchEvent* newEvent = new SpExBranchEvent(0.0, 0.0, 0.0, 0.0,
        _treePtr->mapEventToTree(_lastDeletedEventMapTime), _treePtr, _ran,
            _lastDeletedEventMapTime);

    newEvent->setLamInit(_lastDeletedEventLambdaInit);
    newEvent->setLamShift(_lastDeletedEventLambdaShift);
    newEvent->setMuInit(_lastDeletedEventMuInit);
    newEvent->setMuShift(_lastDeletedEventMuShift);

    return newEvent;
}


void SpExModel::updateLambdaInitMH(void)
{

    //int n_events = _eventCollection.size() + 1;
    int toUpdate = _ran->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        EventSet::iterator it = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            ++it;

        be = static_cast<SpExBranchEvent*>(*it);
    }

    double oldRate = be->getLamInit();
    double cterm = exp( _updateLambdaInitScale * (_ran->uniformRv() - 0.5) );

    be->setLamInit(cterm * oldRate);

    _treePtr->setMeanBranchSpeciation();
    _treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLogLikelihood();

    double logPriorRatio = _ran->lnExponentialPdf(_lambdaInitPrior,
        be->getLamInit()) - _ran->lnExponentialPdf(_lambdaInitPrior, oldRate);

    double LogProposalRatio = log(cterm);

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
    } else {

        // revert to previous state
        be->setLamInit(oldRate);

        _treePtr->setMeanBranchSpeciation();
        _treePtr->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
    }

    incrementGeneration();
}


void SpExModel::updateLambdaShiftMH(void)
{
    //int n_events = _eventCollection.size() + 1;
    int toUpdate = _ran->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        EventSet::iterator it = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            ++it;

        be = static_cast<SpExBranchEvent*>(*it);
    }

    double oldLambdaShift = be->getLamShift();
    double newLambdaShift = oldLambdaShift + _ran->normalRv((double)0.0,
                            _updateLambdaShiftScale);
    be->setLamShift(newLambdaShift);

    _treePtr->setMeanBranchSpeciation();
    _treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLogLikelihood();

    double  logPriorRatio = _ran->lnNormalPdf((double)0.0,
        _settings->getLambdaShiftPrior(), newLambdaShift);
    logPriorRatio -= _ran->lnNormalPdf((double)0.0, _settings->getLambdaShiftPrior(),
                                      oldLambdaShift);

    double LogProposalRatio = 0.0;

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);



    if (acceptMove == true) {

        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
    } else {

        // revert to previous state
        be->setLamShift(oldLambdaShift);

        _treePtr->setMeanBranchSpeciation();
        _treePtr->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
    }

    incrementGeneration();
}


void SpExModel::updateMuInitMH(void)
{
    int toUpdate = _ran->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        EventSet::iterator it = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            ++it;

        be = static_cast<SpExBranchEvent*>(*it);
    }

    double oldRate = be->getMuInit();
    double cterm = exp( _updateMuInitScale * (_ran->uniformRv() - 0.5) );

    be->setMuInit(cterm * oldRate);

    _treePtr->setMeanBranchSpeciation();
    _treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLogLikelihood();

    double logPriorRatio = _ran->lnExponentialPdf(_muInitPrior,
        be->getMuInit()) - _ran->lnExponentialPdf(_muInitPrior, oldRate);



    double LogProposalRatio = log(cterm);

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
    } else {

        // revert to previous state
        be->setMuInit(oldRate);

        _treePtr->setMeanBranchSpeciation();
        _treePtr->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
    }

    incrementGeneration();
}


void SpExModel::updateMuShiftMH(void)
{
    int toUpdate = _ran->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        EventSet::iterator it = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            ++it;

        be = static_cast<SpExBranchEvent*>(*it);
    }

    double oldMuShift = be->getMuShift();
    double newMuShift = oldMuShift + _ran->normalRv((double)0.0,
                        _updateMuShiftScale);

    be->setMuShift(newMuShift);

    _treePtr->setMeanBranchSpeciation();
    _treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLogLikelihood();

    double  logPriorRatio = _ran->lnNormalPdf((double)0.0,
                            _settings->getMuShiftPrior(), newMuShift);
    logPriorRatio -= _ran->lnNormalPdf((double)0.0, _settings->getMuShiftPrior(),
                                      oldMuShift);


    double LogProposalRatio = 0.0;

    double likeRatio = PropLnLik - getCurrentLogLikelihood();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;



    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {

        // revert to previous state
        be->setMuShift(oldMuShift);

        _treePtr->setMeanBranchSpeciation();
        _treePtr->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
    }

    incrementGeneration();
}


double SpExModel::computeSpecificLogLikelihood()
{
    return computeLikelihoodBranchesByInterval();
}


double SpExModel::computeLikelihoodBranchesByInterval()
{
    double LnL = 0.0;

    if (_settings->getSampleFromPriorOnly())
        return 0.0;

    int numNodes = _treePtr->getNumberOfNodes();

    // LEft and right extinction probabilities for root node
    double rootEleft = 0.0;
    double rootEright = 0.0;

    for (int i = 0; i < numNodes; i++) {
        Node* xnode = _treePtr->getNodeFromDownpassSeq(i);
        if (xnode->getLfDesc() != NULL && xnode->getRtDesc() != NULL) {
            // NOT tip, but MUST ultimately be root.

            // Do left descendant:
            Node* ldesc = xnode->getLfDesc();

            double lDinit = ldesc->getDinit();
            double lEinit = ldesc->getEinit();
            double starttime = ldesc->getBrlen();
            double endtime = ldesc->getBrlen();

            double LtotalT = 0.0; // what is this doing?

            double meanRateAtBranchBase = 0.0;
            double curLam = 0.0;

            while (starttime > 0) {
                //std::cout << starttime << "\t" << endtime << std::endl;
                starttime -= _segLength;
                if (starttime < 0)
                    starttime = 0.0;
                double deltaT = endtime - starttime;

                LtotalT += deltaT;

                curLam = ldesc->computeSpeciationRateIntervalRelativeTime(starttime, endtime);

                double curMu = ldesc->computeExtinctionRateIntervalRelativeTime(starttime,
                               endtime);

                double numL = 0.0;
                double denomL = 0.0;


                numL = (exp( deltaT * (curMu - curLam)) * lDinit * ((curLam - curMu) *
                        (curLam - curMu) ) );
                denomL = ( curLam - (lEinit * curLam) + (exp(deltaT * (curMu - curLam)) *
                           (lEinit * curLam - curMu)));

                lDinit = (numL / (denomL * denomL));
                LnL += log(lDinit);
                lDinit = 1.0;


                double Enum = (1 - lEinit) * (curLam - curMu);
                double Edenom =  (1 - lEinit) * curLam - (exp((curMu - curLam) * (deltaT))) *
                                 (curMu - curLam * lEinit);

                lEinit = 1.0 - (Enum / Edenom);


                endtime = starttime; // reset starttime to old endtime
            }

            // this is to get node speciation rate using approximations
            //   to correspond to fact that branch likelihoods themselves are computed
            //      using approximations.

            meanRateAtBranchBase = curLam / 2;
            curLam = 0.0;

            // Setting extinction prob at root node IF xnode is the root
            if (xnode == _treePtr->getRoot())
                rootEleft = lEinit;

            // Compute speciation for right descendant
            // Do right descendant:
            Node* rdesc = xnode->getRtDesc();

            double rDinit = rdesc->getDinit();
            double rEinit = rdesc->getEinit();

            starttime = rdesc->getBrlen();
            endtime = rdesc->getBrlen();

            double RtotalT = 0.0;

            while (starttime > 0) {

                starttime -= _segLength;
                if (starttime < 0)
                    starttime = 0.0;
                double deltaT = endtime - starttime;

                RtotalT += deltaT;

                curLam = rdesc->computeSpeciationRateIntervalRelativeTime(starttime, endtime);

                double curMu = rdesc->computeExtinctionRateIntervalRelativeTime(starttime,
                               endtime);

                double numL = 0.0;
                double denomL = 0.0;

                numL = (exp( deltaT * (curMu - curLam)) * rDinit * ((curLam - curMu) *
                        (curLam - curMu) ) );
                denomL = ( curLam - (rEinit * curLam) + (exp(deltaT * (curMu - curLam)) *
                           (rEinit * curLam - curMu)));

                rDinit = (numL / (denomL * denomL));
                LnL += log(rDinit);
                rDinit = 1.0;

                double Enum = 0.0;
                double Edenom = 0.0;

                Enum = (1 - rEinit) * (curLam - curMu);
                Edenom =  (1 - rEinit) * curLam - (exp((curMu - curLam) * (deltaT))) *
                          (curMu - curLam * rEinit);


                rEinit = 1.0 - (Enum / Edenom);



                endtime = starttime; // reset starttime to old endtime


            }

            meanRateAtBranchBase += curLam / 2;

            // ########################## What to use as Einit for start of NEXT downstream branch?
            // At this point, lEinit is actually the lEinit for the parent node:
            //  Here, we will  (as of 9.1.2012) arbitrarily take this to be the LEFT extinction value:
            xnode->setEinit(lEinit);

            // ######## But alternatively, we could do:
            // Like the above, but now we randomly resolve this. We choose at RANDOM whether to use the right or left Evalue.
            // as we don't know which descendant represents the "parent" state.

            /*          ***************
            if (ran->uniformRv() <= 0.5){
                xnode->setEinit(lEinit);
            }else{
                xnode->setEinit(rEinit);
            }
                        ****************        */


            // Clearly a problem if extinction values approaching/equaling 1.0
            // If so, set to -Inf, leading to automatic rejection of state

            if ((lEinit > MAX_E_PROB) || (rEinit > MAX_E_PROB)) {
                //std::cout << xnode << "\t" << lEinit << "\t" << rEinit << std::endl;
                return -INFINITY;
            }


            if (xnode == _treePtr->getRoot())
                rootEright = rEinit;
            // rDinit at this point should be FINAL value:
            // save as new variable, to keep clear:

            /* SHould be abele to ignore all of these calculations for the root node:
             Must also compute speciation rate for focal node. THis is a critical step.

             Since we are using approximations for the calculations on branches, we should set node speciation
             rate to be equivalent. Currently, I am not doing this - just computing exact rates
             at nodes.
            */

            if (xnode != _treePtr->getRoot()) {

                // Does not include root node, so it is conditioned on basal speciation event occurring:

                LnL  += log(xnode->getNodeLambda());

                //LnL += log(meanRateAtBranchBase);



                xnode->setDinit(1.0);


            }



        } // IF not tip

    } // FOR each node in set



    // 09.15.2012
    // To CONDITION, uncomment the line below:
    // Or, if UNCOMMENTED, comment the line to NOT condition on survival
    LnL -= (log(1 - rootEleft) + log(1 -
                                     rootEright)); // replacement to above for condiioning.


    return LnL;
}


double SpExModel::computeSpecificLogPrior(void)
{
    double logPrior = 0.0;

    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    logPrior += _ran->lnExponentialPdf(_settings->getLambdaInitPrior(),
                                      rootEvent->getLamInit());

    logPrior += _ran->lnNormalPdf((double)0.0, _settings->getLambdaShiftPrior(),
                                 rootEvent->getLamShift());

    logPrior += _ran->lnExponentialPdf(_settings->getMuInitPrior(),
                                      rootEvent->getMuInit());

    logPrior += _ran->lnNormalPdf((double)0.0, _settings->getMuShiftPrior(),
                                 rootEvent->getMuShift());

    int ctr = 0;


    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        SpExBranchEvent* event = static_cast<SpExBranchEvent*>(*it);

        logPrior += _ran->lnExponentialPdf(_settings->getLambdaInitPrior(),
            event->getLamInit());

        logPrior += _ran->lnNormalPdf(0.0, _settings->getLambdaShiftPrior(),
            event->getLamShift());

        logPrior += _ran->lnExponentialPdf(_settings->getMuInitPrior(),
            event->getMuInit());

        logPrior += _ran->lnNormalPdf(0.0, _settings->getMuShiftPrior(),
            event->getMuShift());

        ctr++;

    }

    // Here's prior density on the event rate:
    logPrior += _ran->lnExponentialPdf(_settings->getPoissonRatePrior(),
        getEventRate());

    return logPrior;
}


void SpExModel::getSpecificEventDataString
    (std::stringstream& ss, BranchEvent* event)
{
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(event);

    ss << be->getLamInit() << ","
       << be->getLamShift() << ","
       << be->getMuInit() << ","
       << be->getMuShift();
}
