#include "JumpProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "TraitBranchEvent.h"
#include "TraitModel.h"


JumpProposal::JumpProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(random, settings, model, prior)
{
    _weight = _settings.get<double>("updateRateJump");
    _updateJumpScale = _settings.get<double>("updateJumpScale");
}


double JumpProposal::acceptanceRatio()
{
    if ( static_cast<TraitBranchEvent*>(_event)->isJump() ) {
        return EventParameterProposal::acceptanceRatio();
    } else {
        return 0.0;
    }
}


double JumpProposal::getCurrentParameterValue()
{
    //std::cout << static_cast<TraitBranchEvent*>(_event)->getJump() << std::endl;
    
    return static_cast<TraitBranchEvent*>(_event)->getJump();

}


double JumpProposal::computeNewParameterValue()
{
    return _currentParameterValue + _random.normal(0.0, _updateJumpScale);
}


void JumpProposal::setProposedParameterValue()
{
    static_cast<TraitBranchEvent*>(_event)->
        setJump(_proposedParameterValue);
}


void JumpProposal::revertToOldParameterValue()
{
    static_cast<TraitBranchEvent*>(_event)->
        setJump(_currentParameterValue);

    //if (std::fabs(static_cast<TraitBranchEvent*>(_event)->getJump() < 0.000001)){
    // std::cout << "New: " << _proposedParameterValue << "\tOld: " << _currentParameterValue << std::endl;
    // }
}


void JumpProposal::updateParameterOnTree()
{
    
    // TODO: This should be handled elsewhere
    //   as jump proposal always just affects a single branch
    _tree->setMeanBranchTraitRates();
}


double JumpProposal::computeNonRootLogPriorRatio()
{
    // HIERARCHICAL
    //double jumpvariance = static_cast<TraitModel*>(&_model)->getJumpVariance();
    
    
    return _prior.jumpPrior(_proposedParameterValue) -
           _prior.jumpPrior(_currentParameterValue);
}


double JumpProposal::computeRootLogPriorRatio()
{
    return 0.0;
}
