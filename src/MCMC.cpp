#include "MCMC.h"
#include "Random.h"
#include "Model.h"
#include "ModelFactory.h"

#include "global_macros.h"


#include <climits>

#define DEBUG_MCMC 
#undef DEBUG_MCMC

#define FAIL_CHECK
#undef FAIL_CHECK

// Choose a random number up to INT_MAX - 1, not INT_MAX,
// because MbRandom adds 1 internally, causing an overflow
MCMC::MCMC(Random& seeder, Settings& settings, ModelFactory& modelFactory) :
    _random(seeder.uniformInteger(0, INT_MAX - 1))
{
    _model = modelFactory.createModel(_random, settings);
 
}


MCMC::~MCMC()
{
    delete _model;
}


void MCMC::run(int generations)
{
    for (int g = 0; g < generations; g++) {
         step();
     }
}


void MCMC::step()
{
#ifdef DEBUG_MCMC
  
    std::cout << "\tstored: " << _model->getCurrentLogLikelihood();
    std::cout << "\tActual: " << _model->computeLogLikelihood();
    std::cout << std::endl;

    std::cout << "STEP BEGIN" << std::endl;
    _model->checkModel();


#endif
    //std::cout << "pre-step" << std::endl;
    //_model->checkModel();
    
#ifdef FAIL_CHECK
    double logL_compute = _model->computeLogLikelihood();
    double delta = std::fabs(logL_compute - _model->getCurrentLogLikelihood());
    if (delta > 0.0000001){
        std::cout << "\tNoMatch1: " << logL_compute << "\t" << _model->getCurrentLogLikelihood() << std::endl;
        exit(0);
    }
    
    
#endif
     _model->proposeNewState();

    double acceptanceRatio = _model->acceptanceRatio();
    if (_random.trueWithProbability(acceptanceRatio)) {
    //std::cout << "MCMC::step() ACCEPT " << std::endl;
        _model->acceptProposal();
    } else {
        
        _model->rejectProposal();
    //    std::cout << "MCMC::step() REJECT " << _model->getCurrentLogLikelihood() << std::endl;
    }

    //std::cout << "post-step" << std::endl;
    //_model->checkModel();
    
#ifdef DEBUG_MCMC
 
    
    
#endif
    
    
    
    
#ifdef FAIL_CHECK
    
    logL_compute = _model->computeLogLikelihood();
    delta = std::fabs(logL_compute - _model->getCurrentLogLikelihood());
    
     //_model->printEventData();
    //std::cout << "\nupdateStatus after MCMC update" << std::endl;
    //_model->printNodeUpdateStatus();
    //std::cout << "MCMC logLike after: " << _model->getCurrentLogLikelihood() << std::endl;
    //std::cout << "Last parm updated: " << _model->getLastParameterUpdated() << "\tAccept: ";
    //std::cout << _model->getAcceptLastUpdate() << std::endl;
    if (delta > 0.0000001){
        std::cout << "\tMCMC::step() NoMatch1: " << logL_compute << "\t" << _model->getCurrentLogLikelihood() << std::endl;
        std::cout << "Last parm updated: " << _model->getLastParameterUpdated() << "\tAccept: ";
        std::cout << _model->getAcceptLastUpdate() << std::endl;
        exit(0);
    }
#endif 
}
