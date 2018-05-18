#include "EventNumberProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"

#include "global_macros.h"


#include <algorithm>


EventNumberProposal::EventNumberProposal
    (Random& random, Settings& settings, Model& model) :
        _random(random), _model(model)
{
    _weight = settings.get<double>("updateRateEventNumber");

    _validateEventConfiguration =
        settings.get<bool>("validateEventConfiguration");
}


void EventNumberProposal::propose()
{
    _currentEventCount = _model.getNumberOfEvents();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
    _currentLogPrior = _model.computeLogPrior();

    bool shouldAddEvent = (_currentEventCount == 0) ||
        _random.trueWithProbability(0.5);

    if (shouldAddEvent) {
        _lastEventChanged = _model.addRandomEventToTree();
        _lastProposal = AddEvent;
    } else {
        _lastEventChanged = _model.removeRandomEventFromTree();
        _lastProposal = RemoveEvent;
    }
    _model.setMeanBranchParameters();
    _proposedEventCount = _model.getNumberOfEvents();
    _proposedLogLikelihood = _model.computeLogLikelihood();
    _proposedLogPrior = _model.computeLogPrior();
    
    // std::cout << "Current from node current " << _model.sumNodeLikelihoods() << std::endl;
    // std::cout << "EventNumberProposal::propose() proposed LogL : " << _proposedLogLikelihood << std::endl;
    
}


void EventNumberProposal::accept()
{
    if (_lastProposal == RemoveEvent) {
        if (_lastEventChanged != NULL) {
            delete _lastEventChanged;
            _lastEventChanged = NULL;
        }
    }

    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void EventNumberProposal::reject()
{
    if (_lastProposal == AddEvent) {
        
        //std::cout << " EventNumberProposal::reject() // about to remove event from tree" << std::endl;
        
        _model.removeEventFromTree(_lastEventChanged);
        _model.setMeanBranchParameters();
        delete _lastEventChanged;
        _lastEventChanged = NULL;
    } else if (_lastProposal == RemoveEvent) {
        // std::cout << " EventNumberProposal::reject() // about to add previous event" << std::endl;
        
        _model.addEventToTree(_lastEventChanged);
    }
    
    // double ll = _model.computeLogLikelihood();
    
    // std::cout << "done with EventNumberProposal::reject()" << "\tcurr LogL" << _model.getCurrentLogLikelihood();
    // std::cout << "\tcomputed: " << ll << std::endl;
}


double EventNumberProposal::acceptanceRatio()
{
    if (_validateEventConfiguration && _lastProposal == AddEvent &&
            !_model.isEventConfigurationValid(_lastEventChanged)) {
        return 0.0;
    }

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


double EventNumberProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double EventNumberProposal::computeLogPriorRatio()
{
    if (_lastProposal == AddEvent) {
        return _proposedLogPrior - _currentLogPrior +
            std::log(_model.getEventRate()) - std::log(_currentEventCount + 1);
    } else {
        return _proposedLogPrior - _currentLogPrior +
            std::log(_currentEventCount) - std::log(_model.getEventRate());
    }
}


double EventNumberProposal::computeLogQRatio()
{
    if (_lastProposal == AddEvent) {
        // -0.6931... is ln 0.5
        double logQRatio = (_currentEventCount > 0) ? 0.0 : -0.69314718055995;
        return logQRatio - _model.logQRatioJump();
    } else {
        // 0.6931... is ln 2.0
        double logQRatio = (_currentEventCount != 1) ? 0.0 : 0.69314718055995;
        return logQRatio + _model.logQRatioJump();
    }
}

bool EventNumberProposal::checkIsLikelihoodValid()
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
