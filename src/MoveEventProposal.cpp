#include "MoveEventProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Node.h"
#include "BranchHistory.h"
#include "Tree.h"

#include "global_macros.h"


#include <algorithm>


MoveEventProposal::MoveEventProposal
    (Random& random, Settings& settings, Model& model) :
        _random(random), _settings(settings), _model(model)
{
    _weight = _settings.get<double>("updateRateEventPosition");

    _localToGlobalMoveRatio = _settings.get<double>("localGlobalMoveRatio");
    _scale = _settings.get<double>("updateEventLocationScale") *
        _model.getTreePtr()->maxRootToTipLength();

    _validateEventConfiguration =
        _settings.get<bool>("validateEventConfiguration");
}


void MoveEventProposal::propose()
{
    _currentEventCount = _model.getNumberOfEvents();
    if (_currentEventCount == 0) {
        return;
    }

    //_model.completeDebugPrint("MoveEventProposa::propose() START");
    
    _event = _model.chooseEventAtRandom();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();

    if (_event->getIsEventValidForNode() == true){
        
        // histories should be set forward from here
        BranchEvent* previousEvent = _event->getEventNode()->getBranchHistory()->
        getLastEvent(_event);
        
        _event->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_event);
        
        double localMoveProb = _localToGlobalMoveRatio /
        (1 + _localToGlobalMoveRatio);
        
        // Choose to move locally or globally
        if (_random.trueWithProbability(localMoveProb)) {
            double step = _random.uniform(0, _scale) - 0.5 * _scale;
            _event->moveEventLocal(step);
        } else {
            _event->moveEventGlobal();
        }
        
        _event->getEventNode()->getBranchHistory()->addEventToBranchHistory(_event);
        
        _model.forwardSetBranchHistories(previousEvent);
        _model.forwardSetBranchHistories(_event);
        _model.setMeanBranchParameters();
    
        //std::cout << "MEP::propose() 1 " << std::endl;
        
    }else{

        _event->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_event);
        
        double localMoveProb = _localToGlobalMoveRatio /
        (1 + _localToGlobalMoveRatio);
        
        // Choose to move locally or globally
        if (_random.trueWithProbability(localMoveProb)) {
            double step = _random.uniform(0, _scale) - 0.5 * _scale;
            _event->moveEventLocal(step);
        } else {
            _event->moveEventGlobal();
        }
        
        _event->getEventNode()->getBranchHistory()->addEventToBranchHistory(_event);
 
        // what is problem w this?
        _model.setMeanBranchParameters(_event->getEventNode());
        
        _model.setMeanBranchParameters();

    }
    
    // This is the event preceding the chosen event;


    _proposedLogLikelihood = _model.computeLogLikelihood();
    //_model.globalSetAllNodesNewUpdate();
 
    //std::cout << "MoveEventProposal::propose() _proposed LogL " << _proposedLogLikelihood;
    //std::cout << "\tAfterREset: " << _model.computeLogLikelihood();
    //std::cout << "\tCurrStored: " << _currentLogLikelihood << std::endl;
    //_model.completeDebugPrint("MoveEventProposa::propose() END");
}


void MoveEventProposal::accept()
{
    if (_currentEventCount == 0) {
        return;
    }

    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
    
    //_model.completeDebugPrint("MoveEventProposal::accept // after setCurrentLogLikelihood");
}


void MoveEventProposal::reject()
{
    if (_currentEventCount == 0) {
          return;
    }

    //std::cout << "start node status in MEP::reject" << std::endl;
    //std::cout << "intitial LogL " << _model.computeLogLikelihood() << std::endl;

    //std::cout << "after reflag all events" << std::endl;
    //_model.globalSetAllNodesNewUpdate();
    //std::cout << "intitial LogL " << _model.computeLogLikelihood() << std::endl;
    
    //_model.printEventData();
    //_model.printNodeUpdateStatus();
    
    if (_event->getIsEventValidForNode() == true){
        // Get last event from position of event to be removed
        BranchEvent* lastEvent = _event->getEventNode()->getBranchHistory()->
        getLastEvent(_event);
        
        // Pop event off its new location
        _event->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_event);
        
        // Reset nodeptr, reset mapTime
        _event->revertOldMapPosition();
        
        // Now reset forward from _lastEventChanged (new position)
        // and from newLastEvent, which holds 'last' event before old position
        _event->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(_event);
        
        _model.forwardSetBranchHistories(lastEvent);
        _model.forwardSetBranchHistories(_event);
        _model.setMeanBranchParameters();
        
        
        
        //std::cout << "MEP::reject / curstored: " << _model.getCurrentLogLikelihood() << "\tComputed: ";
        //std::cout << "after branch updates" << _model.computeLogLikelihood() << std::endl;
        //std::cout << "after reflag all events" << std::endl;
        //_model.globalSetAllNodesNewUpdate();
        //std::cout << "intitial LogL " << _model.computeLogLikelihood() << std::endl;
        //_model.completeDebugPrint("MoveEventProposal::REJECT // after resetting branch histories");
        
    
    }else{
 
        // Pop event off its new location
        _event->getEventNode()->getBranchHistory()->popEventOffBranchHistory(_event);
        
        // Reset nodeptr, reset mapTime
        _event->revertOldMapPosition();
        
        // Now reset forward from _lastEventChanged (new position)
        // and from newLastEvent, which holds 'last' event before old position
        _event->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(_event);
        _model.setMeanBranchParameters();
     }
    //std::cout << "end node data in MEP::reject" << std::endl;
    //_model.printNodeUpdateStatus();
    //_model.printEventData();
}


double MoveEventProposal::acceptanceRatio()
{
    if (_currentEventCount == 0) {
        return 0.0;
    }

    if (_validateEventConfiguration &&
            !_model.isEventConfigurationValid(_event)) {
        return 0.0;
    }

    double logLikelihoodRatio = computeLogLikelihoodRatio();

    double t = _model.getTemperatureMH();
    double logRatio = t * logLikelihoodRatio;

    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }
}


double MoveEventProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}

bool MoveEventProposal::checkIsLikelihoodValid()
{
    double logL_compute = _model.computeLogLikelihood();
    double logL_stored = _model.getCurrentLogLikelihood();
    double delta = std::fabs(logL_compute - logL_stored);
    if (delta > 0.0000001){
        return false;
    }else{
        return true;
    }
}
