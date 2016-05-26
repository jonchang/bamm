#include "EventParameterProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"

#include <algorithm>
#include <cmath>

#include "global_macros.h"

EventParameterProposal::EventParameterProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        _random(random), _settings(settings), _model(model), _prior(prior),
        _tree(_model.getTreePtr())
{
}


void EventParameterProposal::propose()
{
    //_model.completeDebugPrint("EventParameterProposal::proposed() START");

    
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
 
    _event = _model.chooseEventAtRandom(true);
    _currentParameterValue = getCurrentParameterValue();
    _proposedParameterValue = computeNewParameterValue();
   
    //std::cout << "Node: " << _event->getEventNode() << "\tCurrPar: " << _currentParameterValue;
    //std::cout << "\tProposed par: " << _proposedParameterValue << std::endl;
    
    setProposedParameterValue();
    updateParameterOnTree();

    _proposedLogLikelihood = _model.computeLogLikelihood();
    //_model.completeDebugPrint("EventParameterProposal::proposed() END");

}


void EventParameterProposal::accept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
    
    //_model.completeDebugPrint("EventParameterProposal::reject // after setCurrentLogLikelihood");
    
}


void EventParameterProposal::reject()
{
    
    //std::cout << "EventParameterProposal::reject() // Event data before reverting::" << std::endl;
    // _model.printEventData();
    
    _model.computeLogLikelihood();
    
    
    revertToOldParameterValue();
    updateParameterOnTree();
    
    //_model.completeDebugPrint("EventParameterProposal::reject // after reverting");


    // std::cout << "EventParameterProposal::reject() // Event data AFTER reverting::" << std::endl;
    // _model.printEventData();
    // _model.computeLogLikelihood();
    
}


double EventParameterProposal::acceptanceRatio()
{
    double logLikelihoodRatio = computeLogLikelihoodRatio();
    double logPriorRatio = computeLogPriorRatio();
    double logQRatio = computeLogQRatio();

    double t = _model.getTemperatureMH();
    double logRatio = t * (logLikelihoodRatio + logPriorRatio) + logQRatio;

    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }
}


double EventParameterProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double EventParameterProposal::computeLogPriorRatio()
{
    if (_event == _model.getRootEvent()) {
        return computeRootLogPriorRatio();
    } else {
        return computeNonRootLogPriorRatio();
    }
}


double EventParameterProposal::computeRootLogPriorRatio()
{
    return 0.0;
}


double EventParameterProposal::computeNonRootLogPriorRatio()
{
    return 0.0;
}


double EventParameterProposal::computeLogQRatio()
{
    return 0.0;
}
