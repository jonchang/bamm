// defining this macro constrains analysis to NEGATIVE values
// for the beta shift parameter

//#define NEGATIVE_SHIFT_PARAM
#undef NEGATIVE_SHIFT_PARAM

#undef DEBUG  // This is a problem.


#include "TraitModel.h"
#include "Random.h"
#include "Settings.h"
#include "Tree.h"
#include "Node.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "TraitBranchEvent.h"
#include "BetaInitProposal.h"
#include "BetaShiftProposal.h"
#include "BetaTimeModeProposal.h"
#include "NodeStateProposal.h"
#include "JumpProposal.h"

#include "global_macros.h"


//#include "JumpVarianceProposal.h"

#include "Log.h"
#include "Prior.h"
#include "Stat.h"
#include "Tools.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <cstdlib>
#include <sstream>
#include <cmath>


TraitModel::TraitModel(Random& random, Settings& settings) :
    Model(random, settings)
{
#ifdef NEGATIVE_SHIFT_PARAM
    // Constrain beta shift to be zero or less than zero.
    if (_settings.getBetaShiftInit() > 0) {
        log(Error) << "Initial value of beta shift (betaShiftInit) cannot be\n"
            << "positive. This parameter is constrained to negative values\n";
        std::exit(1);
    }
#endif
    
    _sampleFromPriorOnly = _settings.get<bool>("sampleFromPriorOnly");
    
    _jumpVariance = _settings.get<double>("jumpVariancePrior");
    
    
    double betaInit = _settings.get<double>("betaInit");
    double betaShiftInit = _settings.get<double>("betaShiftInit");

    bool isTimeVariable = false;
    double timeVarPrior = _settings.get<double>("betaIsTimeVariablePrior");
    if (timeVarPrior == 0.0) {
        isTimeVariable = false;
        if (betaShiftInit != 0.0) {
            exitWithError("betaShiftInit needs to be 0.0 if "
                "betaIsTimeVariablePrior is also 0.0");
        }
    } else if (timeVarPrior == 1.0) {
        isTimeVariable = true;
    } else {
        isTimeVariable = betaShiftInit != 0.0;
    }

    BranchEvent* x = new TraitBranchEvent(betaInit, betaShiftInit,
        isTimeVariable, false, 0.0, _tree->getRoot(), _tree, _random, 0);
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the _rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event:
    forwardSetBranchHistories(_rootEvent);

    _tree->setMeanBranchTraitRates();

    // Initialize by previous event histories (or from initial event number)
    if (_settings.get<bool>("loadEventData")) {
        initializeModelFromEventDataFile(_settings.get("eventDataInfile"));
    } else {
        int initialNumberOfEvents = _settings.get<int>("initialNumberEvents");
        for (int i = 0; i < initialNumberOfEvents; i++) {
            addRandomEventToTree();
        }
    }
    setCurrentLogLikelihood(computeLogLikelihood());

    // TODO: Code duplication with SpExModel
    if (std::isinf(getCurrentLogLikelihood())) {
        log(Error) << "Initial log-likelihood is infinity.\n"
            << "Please check your initial parameter values.\n";
        std::exit(1);
    }

    log() << "\nInitial log-likelihood: " << getCurrentLogLikelihood() << "\n";
    if (_sampleFromPriorOnly) {
        log() << "Note that you have chosen to sample from prior only.\n";
    }

    // Add proposals
    _proposals.push_back(new BetaInitProposal(random, settings, *this, _prior));
    _proposals.push_back
        (new BetaShiftProposal(random, settings, *this, _prior));
    _proposals.push_back(new NodeStateProposal(random, settings, *this));
    _proposals.push_back(new BetaTimeModeProposal(random, settings, *this));
    _proposals.push_back(new JumpProposal(random, settings, *this, _prior));
    
    // Unimplemented hierarchical model
    //_proposals.push_back(new JumpVarianceProposal(random, settings, *this, _prior));
    
    Model::calculateUpdateWeights();
}


void TraitModel::setRootEventWithReadParameters
    (const std::vector<std::string>& parameters)
{
    TraitBranchEvent* rootEvent = static_cast<TraitBranchEvent*>(_rootEvent);

    rootEvent->setBetaInit(betaInitParameter(parameters));
    rootEvent->setBetaShift(betaShiftParameter(parameters));
}


BranchEvent* TraitModel::newBranchEventWithReadParameters
    (Node* x, double time, const std::vector<std::string>& parameters)
{
    double betaInit = betaInitParameter(parameters);
    double betaShift = betaShiftParameter(parameters);

    // TODO: Return true for now for time-variable
    return new TraitBranchEvent(betaInit, betaShift, true,
                false, 0.0,
                 x, _tree, _random, time);
}




double TraitModel::betaInitParameter(const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[0]);
}


double TraitModel::betaShiftParameter
    (const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[1]);
}


void TraitModel::setMeanBranchParameters()
{
    // DEBUG
    // std::cout << "TraitModel::setMeanBranchParams(): ROOT: \t" << _tree->getRoot()->getTraitValue() << std::endl;
    _tree->setMeanBranchTraitRates();

}

// update for single node only...
void TraitModel::setMeanBranchParameters(Node* x)
{
    
    _tree->computeMeanTraitRatesByNode(x);
    
}



BranchEvent* TraitModel::newBranchEventWithRandomParameters(double x)
{
 
    // Sample whether event is jump from prior:
    bool isNewEventJump = _prior.generateIsEventJumpFromPrior();
    
    // choose new event as jump from uniform
    //bool isNewEventJump = (_random.uniform() < 0.5);
    
    bool newIsTimeVariable = _prior.generateBetaIsTimeVariableFromPrior();
    
    
    
    double newbeta = 0.0;
    double newBetaShift = 0.0;
    double newjump = 0.0;
    
    _logQRatioJump = 0.0;
    
    if (isNewEventJump){

        // for hierarchical model:
        //newjump = _prior.generateJumpFromPrior(_jumpVariance);
        //_logQRatioJump += _prior.jumpPrior(newjump, _jumpVariance);
        
        // For fixed jump variance
        newjump = _prior.generateJumpFromPrior();
        _logQRatioJump += _prior.jumpPrior(newjump);
        
    }else{
 
        
        // Sample beta and beta shift from prior
        newbeta = _prior.generateBetaInitFromPrior();
        newBetaShift = _prior.generateBetaShiftFromPrior();
        newIsTimeVariable = _prior.generateBetaIsTimeVariableFromPrior();
 
        
#ifdef NEGATIVE_SHIFT_PARAM
        newBetaShift = -fabs(newBetaShift);
        double dens_term = std::log(2.0);
#else
        double dens_term = 0.0;
#endif    
        
        /*
        if (isNewEventJump){
            _logQRatioJump += std::log(_prior.isEventJumpPrior());
        }else{
            _logQRatioJump +=   std::log(1 - _prior.isEventJumpPrior());
        }
        */
        
        _logQRatioJump += _prior.betaInitPrior(newbeta);
        
        if (newIsTimeVariable) {
            _logQRatioJump += dens_term + _prior.betaShiftPrior(newBetaShift);
        }
        
    }
 

    BranchEvent* zz =  new TraitBranchEvent(newbeta, newBetaShift, newIsTimeVariable, isNewEventJump,
                                          newjump, _tree->mapEventToTree(x), _tree, _random, x);
   
    
    return zz;
}


// TODO: test this : has not been checked
//     also not implemented in constructor for TraitModel yet

// Should not need this for jump-type events.
BranchEvent* TraitModel::newBranchEventWithParametersFromSettings(double x)
{
    
    // x is map time
    double newbeta = _settings.get<double>("betaInit0");
    double newBetaShift = _settings.get<double>("betaShift0");
    bool newIsTimeVariable = _prior.generateBetaIsTimeVariableFromPrior();

    
    // TODO: This needs to be refactored somewhere else
    // Computes the jump density for the addition of new parameters.
#ifdef NEGATIVE_SHIFT_PARAM
    newBetaShift = -fabs(newBetaShift);
    double dens_term = std::log(2.0);
#else
    double dens_term = 0.0;
#endif
    
    _logQRatioJump = 0.0;
    
    _logQRatioJump += _prior.betaInitPrior(newbeta);
    if (newIsTimeVariable) {
        _logQRatioJump += dens_term + _prior.betaShiftPrior(newBetaShift);
    }
    
    return new TraitBranchEvent(newbeta, newBetaShift, newIsTimeVariable,
                                false, 0.0,
                                _tree->mapEventToTree(x), _tree, _random, x);
    
}



void TraitModel::setDeletedEventParameters(BranchEvent* be)
{
    TraitBranchEvent* event = static_cast<TraitBranchEvent*>(be);

    _lastDeletedEventBetaInit = event->getBetaInit();
    _lastDeletedEventBetaShift = event->getBetaShift();
    _lastDeletedEventTimeVariable = event->isTimeVariable();
    _lastDeletedEventIsJump = event->isJump();
    _lastDeletedEventJump = event->getJump();
}


double TraitModel::calculateLogQRatioJump()
{
    _logQRatioJump = 0.0;

    _logQRatioJump = _prior.betaInitPrior(_lastDeletedEventBetaInit);
    _logQRatioJump += _prior.betaShiftPrior(_lastDeletedEventBetaShift);
    if (_lastDeletedEventIsJump){
        
        // Hierarchical model (neeeds work)
        //_logQRatioJump = _prior.jumpPrior(_lastDeletedEventJump, _jumpVariance);
        
        _logQRatioJump = _prior.jumpPrior(_lastDeletedEventJump);
        //_logQRatioJump += std::log(_prior.isEventJumpPrior());
 
    }else{
        _logQRatioJump = _prior.betaInitPrior(_lastDeletedEventBetaInit);
        
        if(_lastDeletedEventTimeVariable){
            _logQRatioJump += _prior.betaShiftPrior(_lastDeletedEventBetaShift);
        }
        //_logQRatioJump +=   std::log(1 - _prior.isEventJumpPrior());
    }
    return _logQRatioJump;
}


BranchEvent* TraitModel::newBranchEventFromLastDeletedEvent()
{
    return new TraitBranchEvent(_lastDeletedEventBetaInit,
        _lastDeletedEventBetaShift, _lastDeletedEventTimeVariable,
        _lastDeletedEventIsJump, _lastDeletedEventJump,
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _random,
        _lastDeletedEventMapTime);
}


double TraitModel::computeLogLikelihood()
{

    if (_likelihoodPower < 0.000000001)
        return 0.0;
    
    double LnL = 0.0;

    //Node * tmpnode = _tree->getRoot()->getLfDesc();

    if (_sampleFromPriorOnly)
        return 0.0;

#ifdef NO_DATA
    LnL = 0.0;
#else
    int numNodes = _tree->getNumberOfNodes();

    // iterate over non-root nodes and compute LnL

    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        Node* xnode = postOrderNodes[i];
        if ( (xnode != _tree->getRoot()) && (xnode->getCanHoldEvent() == true) ) {


            double var = xnode->getBrlen() * xnode->getMeanBeta();
            
            // std::cout << xnode << "\tMeanBeta: " << xnode->getMeanBeta() << std::endl;
            
            // change in phenotype:
            double delta = xnode->getTraitValue() - xnode->getAnc()->getTraitValue();

            // The expected change in phenotype based on the jump process:
            double expectedDelta = xnode->getNetJump();
            
            LnL += Stat::lnNormalPDF(delta, expectedDelta, std::sqrt(var));
        }

    }

#endif
        
    return LnL * _likelihoodPower ;

}

double TraitModel::computeTriadLikelihoodTraits(Node* x)
{
    if (_likelihoodPower < 0.000000001)
        return 0.0;

    // std::cout << "TraitModel::computeTriadLikelihoodTraits ROOT: \t" << _tree->getRoot()->getTraitValue() <<std::endl;
    if (_sampleFromPriorOnly)
        return 0.0;

#ifdef DEBUG
    std::cout << "Enter computeTriadLikelihood: Node : " << x << std::endl;
#endif

    double logL = 0.0;

    // Can only use this likelihood if node contributes to
    // likelihood of observed data

    if (x->getCanHoldEvent() == true) {

        // computation for left descendant branch:

        if (x->getLfDesc()->getCanHoldEvent() == true) {
            double delta = x->getLfDesc()->getTraitValue() - x->getTraitValue();
            double expectedDelta = x->getLfDesc()->getNetJump();
            
            double var = x->getLfDesc()->getBrlen() *
                x->getLfDesc()->getMeanBeta();
            
            // std::cout << "Lf: delt/E[delt]/var :" << delta << "\t" << expectedDelta << "\t" << var << std::endl;
            
            logL += Stat::lnNormalPDF(delta, expectedDelta, std::sqrt(var));
        }


        if (x->getRtDesc()->getCanHoldEvent() == true) {
            // computation for right descendant branch
            double delta = x->getRtDesc()->getTraitValue() - x->getTraitValue();
            double expectedDelta = x->getRtDesc()->getNetJump();
            
            double var = x->getRtDesc()->getBrlen() *
                x->getRtDesc()->getMeanBeta();
            
            // std::cout << "Rt: delt/E[delt]/var :" << delta << "\t" << expectedDelta << "\t" << var << std::endl;

            
            logL += Stat::lnNormalPDF(delta, expectedDelta, std::sqrt(var));
        }


        // computation for ancestral branch (unless == root)

        if (x != _tree->getRoot()) {

            double delta = x->getTraitValue() - x->getAnc()->getTraitValue();
            double expectedDelta = x->getNetJump();
            double var = x->getBrlen() * x->getMeanBeta();
            logL += Stat::lnNormalPDF(delta, expectedDelta, std::sqrt(var));
        }
    }

#ifdef DEBUG
    std::cout << "Leaving computeTriadLikelihood: Node : " << x << std::endl;
#endif
    // std::cout << "TraitModel::computeTriadLikelihood: " << logL << std::endl;

    return logL * _likelihoodPower;

}


double TraitModel::computeLogPrior()
{
#ifdef NEGATIVE_SHIFT_PARAM
    double dens_term = std::log(2.0);
#else
    double dens_term = 0.0;
#endif

    double logPrior = 0.0;

    // TODO: Remove the lines below, was double-counting root event
    /*
    TraitBranchEvent* re = static_cast<TraitBranchEvent*>(_rootEvent);

    logPrior += _prior.betaInitRootPrior(re->getBetaInit());
    if (re->isTimeVariable()) {
        logPrior += dens_term + _prior.betaShiftRootPrior(re->getBetaShift());
    }
    */
     

    for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
         i != _eventCollection.end(); ++i) {

        TraitBranchEvent* event = static_cast<TraitBranchEvent*>(*i);

        if (event->isJump()){
 
            logPrior += _prior.jumpPrior(event->getJump());
            
            // hierarchical
            //logPrior += _prior.jumpPrior(event->getJump(), _jumpVariance);

            // Uncomment this if new events have equal prob of being jump or not
            //logPrior += std::log(_prior.isEventJumpPrior());
            
        }else{
            logPrior += _prior.betaInitPrior(event->getBetaInit());
            if (event->isTimeVariable()) {
                logPrior += dens_term +
                    _prior.betaShiftPrior(event->getBetaShift());
            }
            // Uncomment this if new events have equal prob of being jump or not
            //logPrior += std::log((1 - _prior.isEventJumpPrior() ) );
        }
    }

    // and prior on number of events:

    logPrior += _prior.poissonRatePrior(getEventRate());

    // For unimplemented hierarchical model:
    //logPrior += _prior.jumpVariancePrior(_jumpVariance);
    //std::cout << "logPrior: " << logPrior << std::endl;
    
    return logPrior;

}


void TraitModel::getSpecificEventDataString
    (std::stringstream& ss, BranchEvent* event)
{
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(event);

    ss << be->getBetaInit() << ","
    << be->getBetaShift();
}




void TraitModel::checkModel()
{
    std::cout << "Begin model check..." << std::endl;
    
    // Print Event Data
    std::cout << "Event\tnode\ttime\tbetainit\tbetashift\tisjump\tjump" << std::endl;
     
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it){
        //bool isValidSingle = isEventConfigurationValid((*it));
        
        TraitBranchEvent* be = static_cast<TraitBranchEvent*>((*it));
        
        std::cout << (*it) << "\t" << (*it)->getEventNode() << "\t";
        std::cout << (*it)->getAbsoluteTime() << "\t" << be->getBetaInit() << "\t";
        std::cout <<  be->getBetaShift() << "\t" << be->isJump() << "\t";
        std::cout <<  be->getJump() << std::endl;
    }
    
    BranchEvent* rr = _rootEvent;
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>((rr));
    std::cout << (rr) << "\t" << (rr)->getEventNode() << "\t";
    std::cout << (rr)->getAbsoluteTime() << "\t" << be->getBetaInit() << "\t";
    std::cout <<  be->getBetaShift() << "\t" << be->isJump() << "\t";
    std::cout <<  be->getJump() << std::endl;
    
    // Print Node State Data
    _tree->debugPrintNodeData();
    

}

int TraitModel::getNumberOfJumpEvents()
{
    int n_jumps = 0;
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it){
        TraitBranchEvent* be = static_cast<TraitBranchEvent*>((*it));
        if (be->isJump() == true){
            n_jumps++;
        }
    }
    return n_jumps;
}


int TraitModel::getNumberOfRateShiftEvents()
{
    // This ignores the root event, which is invariant.
    int n_events = 0;
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it){
        
        TraitBranchEvent* be = static_cast<TraitBranchEvent*>((*it));
        if (be->isJump() == false){
            n_events++;
        }
    }
    return n_events;
    
}

double TraitModel::getRootState()
{
    return _tree->getRoot()->getTraitValue();
}


void TraitModel::revertLikelihoodNodeParams()
{

}


void TraitModel::printEventData()
{


}


