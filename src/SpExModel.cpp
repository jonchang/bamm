#include "global_macros.h"

#include "SpExModel.h"
#include "Model.h"
#include "Random.h"
#include "Settings.h"
#include "Tree.h"
#include "Node.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "SpExBranchEvent.h"
#include "LambdaInitProposal.h"
#include "LambdaShiftProposal.h"
#include "MuInitProposal.h"
#include "MuShiftProposal.h"
#include "LambdaTimeModeProposal.h"
#include "PreservationRateProposal.h"

#include "Log.h"
#include "Prior.h"
#include "Tools.h"
#include "Fossil.h"

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

// if undefined,
//  will run (slightly) faster as will not use
//  fossil class
#define ENABLE_FOSSIL

#define JUMP_VARIANCE_NORMAL 0.05

#define NEVER_RECOMPUTE_E0



SpExModel::SpExModel(Random& random, Settings& settings) :
    Model(random, settings), _fossil(new Fossil(settings, *_tree))
{
    // Initial values
    _lambdaInit0 = _settings.get<double>("lambdaInit0");
    _lambdaShift0 = _settings.get<double>("lambdaShift0");
    _muInit0 = _settings.get<double>("muInit0");
    _muShift0 = _settings.get<double>("muShift0");

    if (_muInit0 < 0.0000001 & _settings.get<double>("updateRateMu0") > 0.0000001){
        std::cout << "Invalid starting value for parameter << muInit0 >>" << std::endl;
        std::cout << "Initial value must be greater than zero if updateRateMu0 > 0" << std::endl;
        exit(0);
    }
	
	
	_alwaysRecomputeE0 = _settings.get<std::string>("combineExtinctionAtNodes") == "recompute";
    
    
    _combineExtinctionAtNodes = _settings.get<std::string>("combineExtinctionAtNodes");
    
    // Move this to a separate function at some point

    if (_combineExtinctionAtNodes == "random"){
        
        // This is currently incompatible with MC3
        // so throw exception if called:
        if (_settings.get<int>("numberOfChains") != 1){
            std::cout << "Can't use option 'random' for combineExtinctionAtNodes\n";
            std::cout << "with Metropolis-coupled MCMC. Only single chain analysis";
            std::cout << "\npermitted at present" << std::endl;
            exit(0);
        }
        
        int numNodes = _tree->getNumberOfNodes();
        const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
        
        for (int i = 0; i < numNodes; i++) {
            Node* node = postOrderNodes[i];
            bool left = _random.uniform() <= 0.5;
            node->setInheritFromLeft(left);
            //std::cout << left << std::endl;
        }
    
    }
    // End block for "random" _combineExtinctionAtNodes
    
    int cs = _settings.get<int>("conditionOnSurvival");
    
    std::string fossil_pres_model = _settings.get<std::string>("preservationModel");
    
    if (fossil_pres_model  != "NONE"){
        
        // This assumes that all trees with paleo data should NOT
        // be conditioned on survival.
        _conditionOnSurvival = false;
        
    }
    
    if (cs == 0){
        _conditionOnSurvival = false;
    }else if (cs == 1){
        _conditionOnSurvival = true;
    }else{
        exitWithError("Invalid initial value for parameter <<conditionOnSurvival>>");
    }

    // Initialize fossil preservation rate:
    //      will not be relevant if this is not paleo data.
    _preservationRate = _settings.get<double>("preservationRateInit");
    
     
    double timeVarPrior = _settings.get<double>("lambdaIsTimeVariablePrior");
    if (timeVarPrior == 0.0) {
        _initialLambdaIsTimeVariable = false;
        if (_lambdaShift0 != 0.0) {
            exitWithError("lambdaShift0 needs to be 0.0 if "
                "lambdaIsTimeVariablePrior is also 0.0");
        }
        
        _lambdaIsTimeVariable = false;
        
    } else if (timeVarPrior == 1.0) {
        _initialLambdaIsTimeVariable = true;
        _lambdaIsTimeVariable = true;
    } else {
        _initialLambdaIsTimeVariable = _lambdaShift0 != 0.0;
        _lambdaIsTimeVariable = true;
    }

    _sampleFromPriorOnly = _settings.get<bool>("sampleFromPriorOnly");

    // Parameter for splitting branch into pieces for numerical computation
    _segLength =
        _settings.get<double>("segLength") * _tree->maxRootToTipLength();

    
    //// Change from BranchEvent to SpExBranchEvent:
    BranchEvent* x =  new SpExBranchEvent(_lambdaInit0, _lambdaShift0,
        _muInit0, _muShift0, _initialLambdaIsTimeVariable,
        _tree->getRoot(), _tree, _random, 0);
    
    
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the_rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event
    forwardSetBranchHistories(_rootEvent);

    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();
    
    _extinctionProbMax = _settings.get<double>("extinctionProbMax");

    
    // Initialize by previous event histories (or from initial event number)
    if (_settings.get<bool>("loadEventData")) {
        initializeModelFromEventDataFile(_settings.get("eventDataInfile"));
    } else {
        int initialNumberOfEvents = _settings.get<int>("initialNumberEvents");
        for (int i = 0; i < initialNumberOfEvents; i++) {
            
            // TODO: this adds event to tree with parameters sampled from the prior
            //       should give option to pull parameters from root event (lambdaInit)
            //       etc such that you can add many events but the initial likelihood is
            //       still the same as with a single event
            //addRandomEventToTree();
        
            // This (below) function adds parameters to random locations
            //  but fixes parameters to the initial values specified in control file
            addFixedParameterEventToRandomLocation();
            
            
            // DEBUG
            //std::cout << "Event: \t" << i << "\t" << computeLogLikelihood() << std::endl;
            
        }
    }

    if (_settings.get<bool>("validateEventConfiguration")){
        bool isValid = testEventConfigurationComprehensive();
        if (!isValid){
            std::cout << "\nInitial event configuration is invalid\n";
            std::cout << "Please check your output carefully" << std::endl;
            printEventValidStatus();
        }
    }
    
    setCurrentLogLikelihood(computeLogLikelihood());
#ifdef USE_FAST
    revertLikelihoodNodeParams();
#endif
    
    if (std::isinf(getCurrentLogLikelihood())) {
        log(Error) << "Initial log-likelihood is infinity.\n"
            << "Please check your initial parameter values.\n";
        std::exit(1);
    }

    log() << "\nInitial log-likelihood: " << getCurrentLogLikelihood() << "\n";
    if (_sampleFromPriorOnly)
        log() << "Note that you have chosen to sample from prior only.\n";

    // Add proposals
    _proposals.push_back
        (new LambdaInitProposal(random, settings, *this, _prior));
    _proposals.push_back
        (new LambdaShiftProposal(random, settings, *this, _prior));
    _proposals.push_back(new MuInitProposal(random, settings, *this, _prior));
    _proposals.push_back(new MuShiftProposal(random, settings, *this, _prior));
    _proposals.push_back(new LambdaTimeModeProposal(random, settings, *this));

    //if (_settings.get<std::string>>("preservationModel") != "NONE"){
    if (fossil_pres_model != "NONE"){
        // Cannot set this parameter unless you have paleo data....
        _proposals.push_back(new PreservationRateProposal(random, settings, *this, _prior));
    }


    Model::calculateUpdateWeights();

}

void SpExModel::setRootEventWithReadParameters
    (const std::vector<std::string>& parameters)
{
    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    rootEvent->setLamInit(lambdaInitParameter(parameters));
    rootEvent->setLamShift(lambdaShiftParameter(parameters));
    rootEvent->setMuInit(muInitParameter(parameters));
    rootEvent->setMuShift(muShiftParameter(parameters));
    // TODO: Add time variable/constant
}


BranchEvent* SpExModel::newBranchEventWithReadParameters
    (Node* x, double time, const std::vector<std::string>& parameters)
{
    double lambdaInit = lambdaInitParameter(parameters);
    double lambdaShift = lambdaShiftParameter(parameters);
    double muInit = muInitParameter(parameters);
    double muShift = muShiftParameter(parameters);

    // TODO: Fix reading of parameters (for now, send true for time-variable)
    return new SpExBranchEvent(lambdaInit, lambdaShift,
        muInit, muShift, true, x, _tree, _random, time);
}


double SpExModel::lambdaInitParameter
    (const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[0]);
}


double SpExModel::lambdaShiftParameter
    (const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[1]);
}


double SpExModel::muInitParameter(const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[2]);
}


double SpExModel::muShiftParameter(const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[3]);
}


void SpExModel::setMeanBranchParameters()
{

#ifdef USE_FAST
    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();    
#else

    // TODO: why was BAMM using this
    // before 2016-05-19? Did not hurt
    // but no need to reset mean branch
    // parameters, only node parameters. Wasteful.
    
    _tree->setMeanBranchSpeciation();
     _tree->setMeanBranchExtinction();
#endif

    //_tree->setNodeSpeciationParameters();
    //_tree->setNodeExtinctionParameters();
    
}

void SpExModel::setMeanBranchParameters(Node* x)
{
    //_tree->setMeanBranchSpeciation();
    //_tree->setMeanBranchExtinction();
}


BranchEvent* SpExModel::newBranchEventWithParametersFromSettings(double x)
{
    
    // x is map time
    
    
    double newLam = _settings.get<double>("lambdaInit0");
    double newMu = _settings.get<double>("muInit0");
    double newMuShift = _settings.get<double>("muShift0");
    bool newIsTimeVariable = _prior.generateLambdaIsTimeVariableFromPrior();
    double newLambdaShift = _settings.get<double>("lambdaShift0");
    
    // TODO: This needs to be refactored somewhere else
    // Computes the jump density for the addition of new parameters.
    _logQRatioJump = 0.0;    // Set to zero to clear previous values
    _logQRatioJump = _prior.lambdaInitPrior(newLam);
    if (newIsTimeVariable) {
        _logQRatioJump += _prior.lambdaShiftPrior(newLambdaShift);
    }
    _logQRatioJump += _prior.muInitPrior(newMu);
    _logQRatioJump += _prior.muShiftPrior(newMuShift);
    
    return new SpExBranchEvent(newLam, newLambdaShift, newMu,
                               newMuShift, newIsTimeVariable, _tree->mapEventToTree(x),
                               _tree, _random, x);

}


BranchEvent* SpExModel::newBranchEventWithRandomParameters(double x)
{
    double newLam = _prior.generateLambdaInitFromPrior();
    double newMu = _prior.generateMuInitFromPrior();
    double newMuShift = _prior.generateMuShiftFromPrior();
    bool newIsTimeVariable = _prior.generateLambdaIsTimeVariableFromPrior();

    double newLambdaShift = 0.0;
    if (newIsTimeVariable) {
        newLambdaShift = _prior.generateLambdaShiftFromPrior();
    }

    // TODO: This needs to be refactored somewhere else
    // Computes the jump density for the addition of new parameters.
    _logQRatioJump = 0.0; // Set to zero to clear previous values
    _logQRatioJump = _prior.lambdaInitPrior(newLam);
    if (std::fabs(newLambdaShift) > 1E-6) {
        _logQRatioJump += _prior.lambdaShiftPrior(newLambdaShift);
    }
    _logQRatioJump += _prior.muInitPrior(newMu);
    //_logQRatioJump += _prior.muShiftPrior(newMuShift);

    return new SpExBranchEvent(newLam, newLambdaShift, newMu,
        newMuShift, newIsTimeVariable, _tree->mapEventToTree(x),
        _tree, _random, x);
}


void SpExModel::setDeletedEventParameters(BranchEvent* be)
{
    SpExBranchEvent* event = static_cast<SpExBranchEvent*>(be);

    _lastDeletedEventLambdaInit = event->getLamInit();
    _lastDeletedEventLambdaShift = event->getLamShift();
    _lastDeletedEventMuInit = event->getMuInit();
    _lastDeletedEventMuShift = event->getMuShift();
    _lastDeletedEventTimeVariable = event->isTimeVariable();
}


double SpExModel::calculateLogQRatioJump()
{
    _logQRatioJump = 0.0;

    _logQRatioJump = _prior.lambdaInitPrior(_lastDeletedEventLambdaInit);
 
    // DLR: add check to ensure that model jumping proposal ratio
    //      only multiplies by prior prob of lambdaShift if the
    //      proposed event is time-constant
    
    if (std::fabs(_lastDeletedEventLambdaShift) > 1E-6){
        _logQRatioJump += _prior.lambdaShiftPrior(_lastDeletedEventLambdaShift);
    }
 
    _logQRatioJump += _prior.muInitPrior(_lastDeletedEventMuInit);
    
    //_logQRatioJump += _prior.muShiftPrior(_lastDeletedEventMuShift);
    //std::cout << "logQ in SpExModel::calculate: " << _logQRatioJump << std::endl;
    
    return _logQRatioJump;
}


BranchEvent* SpExModel::newBranchEventFromLastDeletedEvent()
{
    return new SpExBranchEvent(_lastDeletedEventLambdaInit,
        _lastDeletedEventLambdaShift, _lastDeletedEventMuInit,
        _lastDeletedEventMuShift, _lastDeletedEventTimeVariable,
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _random,
        _lastDeletedEventMapTime);
}


// TODO: Not transparent, but this is where
//  Di for internal nodes is being set to 1.0
 
double SpExModel::computeLogLikelihood()
{
    if (_likelihoodPower < 0.000000001){
        return 0.0;
    }
    
    if (_sampleFromPriorOnly)
        return 0.0;
    
    double logLikelihood = 0.0;

    int numNodes = _tree->getNumberOfNodes();

    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    
    // Set all nodes to initial values for NO downstream rate shift.
    for (int i = 0; i < numNodes; i++){
        Node* node = postOrderNodes[i];
        node->setHasDownstreamRateShift(false);
    }    
 
    
// TODO: cleanup EARLY_REJECT / delete this (why was this here?)
// #ifdef EARLY_REJECT
//  double previousLikelihood = 0.0;
//    double bterm = 0.0;
//    double nterm = 0.0;
// #endif
    
    for (int i = 0; i < numNodes; i++) {
        Node* node = postOrderNodes[i];
        
        if (node->isInternal()) {
 
            double LL = computeSpExProbBranch(node->getLfDesc());
            double LR = computeSpExProbBranch(node->getRtDesc());

// TODO: cleanup EARLY_REJECT / delete this
//#ifdef EARLY_REJECT
//            previousLikelihood += node->getLfDesc()->getLogDiCurrent();
//            previousLikelihood += node->getRtDesc()->getLogDiCurrent();
//#endif
 
   

#ifdef USE_FAST
            double E_left = node->getLfDesc()->getExProbProposed();
            double E_right = node->getRtDesc()->getExProbProposed();
#else
            double E_left = node->getLfDesc()->getExtinctionEnd();
            double E_right = node->getRtDesc()->getExtinctionEnd();
#endif
 
            bool left_shift = node->getLfDesc()->getHasDownstreamRateShift();
            bool right_shift = node->getRtDesc()->getHasDownstreamRateShift();
 
			if ( _alwaysRecomputeE0 ){
				
				// Under recompute, values from right and left descendant branches are identical at this point
				//    so can arbitrarily choose the left
				
				node->setEinit(E_left);
				
			} else if (_combineExtinctionAtNodes == "random"){
                
                if (node->getInheritFromLeft() == true){
                    node->setEinit(E_left);
                 
                }else{
                    node->setEinit(E_right);
                 }
                
                
            }else if(_combineExtinctionAtNodes  == "if_different"){
                
                double delta = std::fabs(E_left - E_right);
                
                if (delta < 0.001){
                    node->setEinit(E_left);
                }else{
                    E_left *= E_right;
                    node->setEinit(E_left);
                }
                
			} else if (_combineExtinctionAtNodes == "favor_shift"){
                if (left_shift == true & right_shift == true){
                    node->setEinit( E_left * E_right );
                }else if (left_shift == true & right_shift == false){
                    node->setEinit(E_left);
                }else if (left_shift == false & right_shift == true){
                    node->setEinit(E_right);
                }else if (left_shift == false & right_shift == false){
                    node->setEinit(E_left);
                }else{
                    std::cout << "problem in computeLogLikelihood()" << std::endl;
                    std::cout << "Error in _combineExtinctionAtNodes option" << std::endl;
                    exit(0);
                }
            }else if (_combineExtinctionAtNodes == "left"){
                node->setEinit(E_left);
            }else if (_combineExtinctionAtNodes == "right"){
                node->setEinit(E_right);
            }else{
                std::cout << "unsupported option for combining extinction probabilities" << std::endl;
                exit(0);
            }
			
            
            logLikelihood += (LL + LR);

            
            // Does not include root node, so it is conditioned
            // on basal speciation event occurring:
            if (node != _tree->getRoot()) {
                logLikelihood  += log(node->getNodeLambda());
                
                node->setDinit(1.0);
				
// TODO: clean up  / remove this
//#ifdef EARLY_REJECT
//                previousLikelihood += log(node->getPreviousNodeLambda());
//                if (logLikelihood < (previousLikelihood - (double)5.0) && std::fabs(LL) > 0.00001){
//                     return -INFINITY;
//                }
//#endif
				
            }
        }
    }

    
#ifdef ENABLE_FOSSIL
    
    logLikelihood += _fossil->computeOccurrenceLogLikelihood(_preservationRate);
    
#endif
    
    
    return logLikelihood;
}


double SpExModel::computeSpExProbBranch(Node* node)
{
 
    
#ifdef USE_FAST
    if (node->getProposedUpdate() == false){
        node->setLogDiProposed(node->getLogDiCurrent());
        node->setExProbProposed(node->getExProbCurrent());
        return node->getLogDiProposed();
    }
 
#endif
 
    int n_events = node->getBranchHistory()->getNumberOfBranchEvents();
    
    if (n_events > 0){
        node->setHasDownstreamRateShift(true);
    }
    
    if (node->getHasDownstreamRateShift() == true){
        node->getAnc()->setHasDownstreamRateShift(true);
    }
  
    
    double logLikelihood = 0.0;

    double D0 = node->getDinit();    // Initial speciation probability
    double E0 = node->getEinit();    // Initial extinction probability
    
    bool local_recompute_E0 = false;

    // The event governing tipwards end of branch:
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(node->getBranchHistory()->getLastEvent(node->getTime()));
    double observation_time = _fossil->getObservationTime();
    
    // 3 scenarios:
    //   i. node is extant tip
    //   ii. node is fossil last occurrence, an unsampled or extinct tip
    //   iii. node is internal. Will now treat separately.
    //  case i and iii can be treated the same
    
    // Test if node is extant
    
    // TODO: This test for "extant" vs "non-extant" can be a problem,
    //      depending on numerical error etc. There must be some sort of check.
    //      If lineages that are EXTANT are looped through here,
    //      you will have massively depressed log-likelihoods as it will compute
    //      and add extinction likelihoods for extant lineages.
    //      Problems were observed with simulated trees when the tolerance parameter
    //      was set to 0.00001, as it was flagging many extant taxa as extinct.
    
    bool isExtant = (std::abs(node->getTime() - observation_time)) < 0.01;
    
    // case 1: node is fossil tip
    if (node->isInternal() == false & isExtant == false){
        
        // Get speciation extinction parameters for governing event
        double lam_init = be->getLamInit();
        double lam_shift = be->getLamShift();
        double mu_init  = be->getMuInit();
        double mu_shift = be->getMuShift();
        
        double event_abs_time = be->getAbsoluteTime();
        
        double interval_start_time = observation_time - event_abs_time;
        double interval_stop_time = interval_start_time;
        double tip_time_from_process = node->getTime() - event_abs_time;
        
        while ( interval_start_time > tip_time_from_process ){
            interval_start_time -= _segLength;
            if (interval_start_time < tip_time_from_process ){
                interval_start_time = tip_time_from_process;
            }
            
            double curLam = computeMeanExponentialRateForInterval
            (lam_init, lam_shift, interval_start_time, interval_stop_time);
            
            double curMu  = computeMeanExponentialRateForInterval
            (mu_init, mu_shift, interval_start_time, interval_stop_time);
            
            double abs_time = interval_start_time + event_abs_time;
            double curPsi = _fossil->getCurrentPreservationRate(node, abs_time, _preservationRate);
            
            double spProb = 0.0;
            double exProb = 0.0;
            double deltaT = interval_stop_time - interval_start_time;
            
            computeSpExProb(spProb, exProb, curLam, curMu, curPsi, D0, E0, deltaT);
            
            if (exProb > _extinctionProbMax ){
                
                return -INFINITY;
            }
            
            E0 = exProb;
            
            interval_stop_time = interval_start_time;
            
        }
        
        // Prob that lineage went extinct before present
        // E0 could be the new D0 for the next calculation
        //  however, we will factor this out and start with 1.0.
        
        logLikelihood += log(E0);
 
        
        D0 = 1.0;
        // current value of E0 can now be passed on for
        // calculation down remainder of branch (eg, the observed segement)
 
    } // End case (i) Node is a fossil tip
 
    
    
    double startTime = node->getBrlen();
    double endTime = node->getBrlen();
 
    //SpExBranchEvent* be = static_cast<SpExBranchEvent*>(node->getBranchHistory()->getLastEvent(node->getTime()));
 
    while (startTime > 0) {
        startTime -= _segLength;
        if (startTime < 0) {
            startTime = 0.0;
        }

        double abs_start_time = node->getAnc()->getTime() + startTime;
        double abs_end_time   = node->getAnc()->getTime() + endTime;
        
        // Get speciation extinction parameters for governing event
        double lam_init = be->getLamInit();
        double lam_shift   = be->getLamShift();
        double mu_init  = be->getMuInit();
        double mu_shift = be->getMuShift();

        double absolute_time_event = be->getAbsoluteTime();


        if (be->getAbsoluteTime() >= abs_start_time & be->getAbsoluteTime() < abs_end_time & be != _rootEvent){
            // event on segment.
 
            startTime = be->getAbsoluteTime() - node->getAnc()->getTime();
            
            // Reset start time to absolute time of event if we pass an event on branch
            abs_start_time = node->getAnc()->getTime() + startTime;
            be = static_cast<SpExBranchEvent*>(node->getBranchHistory()->getLastEvent(be));
            
            // set local_recompute_E0 if you switch to new process.
			// which will only be used if combineExtinctionAtNodes = recompute
            
            local_recompute_E0 = true;
            
        }
  
        // Get time relative to the focal event for computation of mean
        //    branch parameters.
        double event_t_start = abs_start_time - absolute_time_event;
        double event_t_end   = abs_end_time   - absolute_time_event;
        
        // ###
        // Have local function compute speciation/extinction rates based on interval
        //  AND the parameters directly.
        //  No more averaging over events on segments.
        
        double deltaT = endTime - startTime;
        double curLam = computeMeanExponentialRateForInterval(lam_init, lam_shift, event_t_start, event_t_end);
        double curMu  = computeMeanExponentialRateForInterval(mu_init, mu_shift, event_t_start, event_t_end);
        
#ifdef ENABLE_FOSSIL
        
        double curPsi = _fossil->getCurrentPreservationRate(node, abs_start_time, _preservationRate);
        
#else
        
        double curPsi = 0.0;
        
#endif
        
        double spProb = 0.0;
        double exProb = 0.0;
        
        
        // Compute speciation and extinction probabilities and store them
        // in spProb and exProb (through reference passing)
        computeSpExProb(spProb, exProb, curLam, curMu, curPsi, D0, E0, deltaT);
 
        if (exProb > _extinctionProbMax ){
        
            return -INFINITY;
        }
        
        logLikelihood += std::log(spProb);
        
        D0 = 1.0;
		
		if (! _alwaysRecomputeE0 ){
			
			// Here we always use the existing E0 for next calculation
			//    but do not recompute.
			E0 = exProb;
		
		}else{
			
			if ( ! local_recompute_E0 ){
				// There is no event on interval, so E0 is the value at the
				//    end of the previous segment
				E0 = exProb;
			
			}else{
				// otherwise, we have to recompute E0 for the focal interval
				// going all the way forward to the present (or observation time)
			
				local_recompute_E0 = false;
				
				
				// Get parameters for next process.
				// SpExBranchEvent* be is now toggled to parent process
				// Get speciation extinction parameters for governing event
				
				double lam_init = be->getLamInit();
				double lam_shift = be->getLamShift();
				double mu_init  = be->getMuInit();
				double mu_shift = be->getMuShift();
				
				
				// The interval of computation should be defined in terms of
				// start_time relative to start time of process,
				// and end_time relative to age of process.
				// NOTE however that we are dealing now with the next upstream process
				// not the process just used to compute D(t).
				// Because only get here if changing process, the difference in age of the
				// last process and next process is the relevant start time.
				
				E0 = node->getEtip();
				
				// get current time relative to age of process:
				double start_rel_to_process = abs_start_time - be->getAbsoluteTime();
				
				// end time relative to age of current process
				double end_rel_to_process = observation_time - be->getAbsoluteTime();
				
				E0 = recomputeE0(start_rel_to_process, end_rel_to_process, lam_init, lam_shift,
								 mu_init, mu_shift, E0);
				
				
			}
			
			
		}

        endTime = startTime;

    
    
    } // while (startTime > 0)
    
    
    Node * parent = node->getAnc();
	
// Addition 2016-11-05
// *should* work for recompute and non-recompute cases
	node->setExtinctionEnd(E0);
	
	
	
// 2016-11-05
// Previous code:
//#ifdef NEVER_RECOMPUTE_E0
    // set extinction end value for branch for current node.
//    node->setExtinctionEnd(E0);
    
    // but do not set parent -- this will happen in the calling function
    // when right and left descendants are computed.
// #else
//    Node * parent = node->getAnc();
    // Should be exactly equal coming from right or left descendant branch at this point.
//    parent->setEinit(E0);

//#endif
    
    if (parent == _tree->getRoot() & _conditionOnSurvival == true){
 
         logLikelihood -= std::log(1.0 - E0);
    }
	
#ifdef USE_FAST
    
    node->setLogDiProposed(logLikelihood);
    node->setExProbProposed(E0);
 
 
    
#endif
    
 
    return logLikelihood;
}


// t_start : start of interval in time units from start of process e.g., where
//  the value of the function rate(t) = rate_init
//  t_end : end of interval in time units from start of process

double SpExModel::computeMeanExponentialRateForInterval(double rate_init, double rate_shift, double t_start, double t_end)
{
    double delta_T = t_end - t_start;
    double integrated = 0.0;
    
    if (rate_shift < 0) {
        integrated = (rate_init / rate_shift) * (std::exp(rate_shift * t_end) - std::exp(rate_shift * t_start));
    } else if (rate_shift > 0) {
        integrated = rate_init * (2 * delta_T + (1.0 / rate_shift) *
                       (std::exp(-rate_shift * t_end) - std::exp(-rate_shift * t_start)));
    } else {
        integrated = rate_init * delta_T;
    }
    
    integrated /= delta_T;
 
    return integrated;

}

double SpExModel::recomputeE0(double start_time, double end_time, double lam_init, double lam_shift,
                   double mu_init, double mu_shift, double Etip)
{
    double E0 = Etip;
    
    // the decrement variable; will be cut by segLength each time...
    double decrementer = end_time;
    
    while (decrementer > start_time){
        decrementer -= _segLength;
        if (decrementer < start_time){
            decrementer = start_time;
        }
        
        double deltaT = end_time - decrementer;
        
        double curLam = computeMeanExponentialRateForInterval(lam_init, lam_shift, decrementer, end_time);
        double curMu  = computeMeanExponentialRateForInterval(mu_init, mu_shift, decrementer, end_time);
        
        double cpsi = _preservationRate;
        
        double sprob = 0.0;
        double eprob = 0.0;
        
        computeSpExProb(sprob, eprob, curLam, curMu, cpsi, (double)1.0, E0, deltaT);
        
        end_time  = decrementer;
        
        E0 = eprob;
		
        // E0 cannot exceed max probability of extinction
        //   necessary to avoid overflow/underflow issues
        //   as E0 approaches 1.0
        if (E0 > _extinctionProbMax){
            return -INFINITY;
        }
            
    }
    return E0;
    
}


//  Notes for the fossil process:
//
//
//
//

////////////  Notes for the non-fossil process
//
//  If we let:
//
//     M = mu - lam
//     L = lam * (1.0 - E0)
//     E = e^(M * deltaT)
//     m = E * M
//     d = L * (1 - E) - m
//
// then the speciation equation from PLoS paper
//
//                    e^((mu - lam) * deltaT) * D0 * (lam - u)^2
//     D(t) = --------------------------------------------------------------
//            [lam - lam * E0 + e^((mu - lam) * deltaT) * (E0 * lam - mu)]^2
//
// becomes
//
//            D0 * m
//     D(t) = ------
//              d^2
//
// and the extinction equation
//
//                               (1 - E0) * (lam - mu)
//     E(t) = 1 - ----------------------------------------------------------
//                (1 - E0) * lam - e^(-(lam - mu) * deltaT) * (mu - lam * E0)
//
// becomes
//
//                (1 - E0) * M
//     E(t) = 1 + ------------
//

// TODO:: Full document equations.
//          equations above are for reconstructed process only.
//          equations as implemented allow for fossilized process.

void SpExModel::computeSpExProb(double& spProb, double& exProb,
                                double lambda, double mu, double psi, double D0, double E0, double deltaT)
{
    
    double FF = lambda - mu - psi;
    double c1 = std::abs(std::sqrt( FF * FF  + (4.0 * lambda * psi) ));
    double c2 = (-1.0) * (FF - 2.0 * lambda * (1.0 - E0)) / c1;
    
    double A = std::exp((-1.0) * c1 * deltaT) * (1.0 - c2);
    double B = c1 * (A - (1 + c2)) / (A + (1.0 + c2));
    
    exProb = (lambda + mu + psi + B)/ (2.0 * lambda);
    
    // splitting up the speciation calculation denominator:
    
    double X = std::exp(c1 * deltaT) * (1.0 + c2)*(1.0 + c2);
    double Y = std::exp((-1.0) * c1 * deltaT) * (1.0 - c2) * (1.0 - c2);
    
    spProb = (4.0 * D0) / ( (2.0 * ( 1 - (c2 * c2)) ) + X + Y );

    // DEBUG
    //std::cout << "\n\tFF\t" << FF << "\tc1\t" << c1 << "\tc2\t" << c2;
    //std::cout << "\n\tA\t" << A << "\tB\t" << B << "\tX\t" << X << "\tY\t" << Y << "\n" << std::endl;
    
    
    
}

/*
void SpExModel::computeSpExProb(double& spProb, double& exProb,
    double lambda, double mu, double D0, double E0, double deltaT)
{
    double M = mu - lambda;
    double L = lambda * (1.0 - E0);
    double E = std::exp(M * deltaT);
    double m = E * M;
    double d = L * (1.0 - E) - m;

    spProb = (D0 * m * M) / sqr(d);
    exProb = 1.0 + (1.0 - E0) * M / d;
}
*/







double SpExModel::computeLogPrior()
{
    double logPrior = 0.0;

// TODO: delete code block here
// Bug: was double-counting root event, by failing
// to account for fact that root event was already in eventCollection
// TODO: removing this deprecates a specific "root" prior
// Since they are no longer computed under root-specific parameters
/*
    
    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    logPrior += _prior.lambdaInitRootPrior(rootEvent->getLamInit());
    if (rootEvent->isTimeVariable()) {
        logPrior += _prior.lambdaShiftRootPrior(rootEvent->getLamShift());
        //logPrior += std::log(_prior.lambdaIsTimeVariablePrior());
    }else{
        //logPrior += std::log(1.0 - _prior.lambdaIsTimeVariablePrior());
    }
    logPrior += _prior.muInitRootPrior(rootEvent->getMuInit());
    //logPrior += _prior.muShiftRootPrior(rootEvent->getMuShift());

*/
 
 
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
 
        SpExBranchEvent* event = static_cast<SpExBranchEvent*>(*it);

        logPrior += _prior.lambdaInitPrior(event->getLamInit());
        
        //if (event->isTimeVariable()) {
        
        if (std::fabs(_lastDeletedEventLambdaShift) > 1E-6){
             logPrior += _prior.lambdaShiftPrior(event->getLamShift());
        //    logPrior += std::log(_prior.lambdaIsTimeVariablePrior());
        }else{
        //    logPrior += std::log(1.0 - _prior.lambdaIsTimeVariablePrior());
        }
        
        logPrior += _prior.muInitPrior(event->getMuInit());
        //logPrior += _prior.muShiftPrior(event->getMuShift());
    }

    // Here's prior density on the event rate
    logPrior += _prior.poissonRatePrior(_eventRate);
 
    // Prior density on the preservation rate, if paleo data:
    // TODO: ultimately move this to class Fossil to allow more complex
    //       preservation scenarios?
#ifdef ENABLE_FOSSIL
    
    if (_settings.get<std::string>("preservationModel") != "NONE"){
        logPrior += _prior.preservationRatePrior(_preservationRate);
    }
    
#endif
 
    //std::cout << "N_events: " << _eventCollection.size() << "\tlogPrior: " << logPrior << std::endl;
    
    
    return logPrior;
}

// TODO: remove this function entirely;
// redundant with class Fossil.
/*
double SpExModel::computePreservationLogProb()
{
    double logLik = (double)_numberOccurrences * std::log(_preservationRate);

    return logLik;
}
*/

/************************/
// DEBUG FUNCTIONS BELOW

void SpExModel::printNodeProbs()
{
    int numNodes = _tree->getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        Node* node = postOrderNodes[i];
        std::cout << node << "\t" << node->getName() << "\t";
        std::cout << node->getEinit() << std::endl;
    
    }

}



void SpExModel::initializeDebugVectors()
{
    int numNodes = _tree->getNumberOfNodes();
    for (int i = 0; i < numNodes; i++){
        _Dinitvec.push_back(0.0);
        _Dfinalvec.push_back(0.0);
        _Einitvec.push_back(0.0);
        _Efinalvec.push_back(0.0);
    }

}

int SpExModel::nodeIndexLookup(Node* node)
{
    int goodnode = -1;
    int numNodes = _tree->getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        if (node == postOrderNodes[i]){
            goodnode = i;
        }
    }
    if (goodnode == -1){
        std::cout << "failed to match node." << std::endl;
        exit(0);
    }
    return goodnode;
}

void SpExModel::outputDebugVectors()
{
    
    std::string outname = "debug_LH.txt";
    std::ofstream outStream;
    outStream.open(outname.c_str());
    
    // Write header:
    outStream << "Rdesc,Ldesc,blen,D_init,D_final,E_init,E_final\n";
    //outStream.close();
 
    int numNodes = _tree->getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        Node* x = postOrderNodes[i];
        outStream << x->getRandomRightDesc() << ",";
        outStream << x->getRandomLeftDesc() << ",";
        outStream << x->getBrlen() << ",";
        outStream << _Dinitvec[i] << ",";
        outStream << _Dfinalvec[i] << ",";
        outStream << _Einitvec[i] << ",";
        outStream << _Efinalvec[i] << "\n";
        
        // x->getRandomLeftDesc() & x->getRandomRightDesc()
        
    }
    
    outStream.close();
    
    
}

void SpExModel::checkModel()
{
    std::cout << "SpExModel::checkModel() " << std::endl;
 
    int numNodes = _tree->getNumberOfNodes();
    
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    
    double llcurrent = 0.0;
    //double llproposed = 0.0;
    
    for (int i = 0; i < numNodes; i++){
        Node* node = postOrderNodes[i];
        llcurrent += node->getLogDiCurrent();
        //llproposed += node->getLogDiProposed();
        
        if (node->getLfDesc() != NULL && node != _tree->getRoot()){
            llcurrent += std::log(node->getNodeLambda());
        }
    
    }
    
    double el = std::log(1.0 - _tree->getRoot()->getLfDesc()->getExProbCurrent());
    double er = std::log(1.0 - _tree->getRoot()->getRtDesc()->getExProbCurrent());
    
    llcurrent -= el;
    llcurrent -= er;
    std::cout << "Summed logL_current: " << llcurrent << "\tel: " << el << "\ter :" << er << std::endl;
    
}


// These are also in Model class
// but reimplemented here to allow USE_FAST compilation
// Flag node params etc...

void SpExModel::forwardSetBranchHistories(BranchEvent* x)
{
 
    
    // If there is another event occurring more recent (closer to tips),
    // do nothing. Even just sits in BranchHistory but doesn't affect
    // state of any other nodes.
    
    // This seems circular, but what else to do?
    // given an event (which references the node defining the branch on which
    // event occurs) you get the corresponding branch history and the last
    // event since the events will have been inserted in the correct order.
    
    if (x->getIsEventValidForNode() == false){
        std::cout << "ERROR forwardSetBranchHistories(BranchEvent* x) / passed jump node" << std::endl;
        exit(0);
    }
    
    Node* myNode = x->getEventNode();
    
    myNode->setProposedUpdate(true);
    
    if (x == _rootEvent) {
        forwardSetHistoriesRecursive(myNode->getLfDesc());
        forwardSetHistoriesRecursive(myNode->getRtDesc());
    } else if (x == myNode->getBranchHistory()->getLastEvent()) {
        // If true, x is the most tip-wise event on branch.
        myNode->getBranchHistory()->setNodeEvent(x);
        
        // If myNode is not a tip
        if (myNode->getLfDesc() != NULL && myNode->getRtDesc() != NULL) {
            forwardSetHistoriesRecursive(myNode->getLfDesc());
            forwardSetHistoriesRecursive(myNode->getRtDesc());
        }
        // Else: node is a tip; do nothing
    }
    // Else: there is another more tipwise event on the same branch; do nothing

    // Set param flags back to root...
    
    _tree->recursiveSetAreParamsCurrentToRoot(myNode);
    
}


/*
 If this works correctly, this will take care of the following:
 1. if a new event is created or added to tree,
 this will forward set all branch histories from the insertion point
 2. If an event is deleted, you find the next event rootwards,
 and call forwardSetBranchHistories from that point. It will replace
 settings due to the deleted node with the next rootwards node.
 */

void SpExModel::forwardSetHistoriesRecursive(Node* p)
{
    // Get event that characterizes parent node
    BranchEvent* lastEvent = p->getAnc()->getBranchHistory()->getNodeEvent();
    
    p->setProposedUpdate(true);
    
    if (lastEvent->getIsEventValidForNode() == false){
        std::cout << "ERROR forwardSetHistoriesRecursive(Node* p) / passed jump node" << std::endl;
        exit(0);
    }
    
    // Set the ancestor equal to the event state of parent node:
    p->getBranchHistory()->setAncestralNodeEvent(lastEvent);
    
    // Ff no events on the branch, go down to descendants and do same thing;
    // otherwise, process terminates (because it hits another event on branch
    if (p->getBranchHistory()->getNumberOfBranchEvents() == 0) {
        p->getBranchHistory()->setNodeEvent(lastEvent);
        
        if (p->getLfDesc() != NULL) {
            forwardSetHistoriesRecursive(p->getLfDesc());
        }
        
        if (p->getRtDesc() != NULL) {
            forwardSetHistoriesRecursive(p->getRtDesc());
        }
    }
}

// only called if accept...
void SpExModel::revertLikelihoodNodeParams()
{
    //std::cout << "SpExModel::revertLikelihoodNodeParams()" << std::endl;
    int numNodes = _tree->getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    
    for (int i = 0; i < numNodes; i++) {
    
        Node* node = postOrderNodes[i];
    
        node->setLogDiCurrent(node->getLogDiProposed());
        node->setExProbCurrent(node->getExProbProposed());
        node->setProposedUpdate(false);
		
// TODO: remove
// #ifdef EARLY_REJECT
//         node->setPreviousNodeLambda(node->getNodeLambda());
// #endif
 
    }

}


void SpExModel::printEventData()
{

    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        SpExBranchEvent* xx = static_cast<SpExBranchEvent*>(*it);
        std::cout << (*it)->getEventNode() << "\t" << xx->getLamInit();
        std::cout << "\t" << xx->getAbsoluteTime() << std::endl;
        
    }
    SpExBranchEvent* xx = static_cast<SpExBranchEvent*>(_rootEvent);
    std::cout << "rootnode:\t"  << xx->getLamInit() << "\t";
    std::cout << "\t" << xx->getAbsoluteTime() << std::endl;
    
}




