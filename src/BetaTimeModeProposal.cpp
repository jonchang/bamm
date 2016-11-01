#include "BetaTimeModeProposal.h"

#include "Settings.h"
#include "Prior.h"
#include "TraitBranchEvent.h"

#include "global_macros.h"


class Random;
class Model;


BetaTimeModeProposal::BetaTimeModeProposal
    (Random& random, Settings& settings, Model& model)
    : TimeModeProposal(random, settings, model)
{
    _weight = settings.get<double>("updateRateBetaTimeMode");
}

double BetaTimeModeProposal::acceptanceRatio()
{
    if ( static_cast<TraitBranchEvent*>(_event)->isJump() == false ) {
        return TimeModeProposal::acceptanceRatio();
    } else {
        return 0.0;
    }
}

double BetaTimeModeProposal::initialParameter(BranchEvent* event)
{
    return static_cast<TraitBranchEvent*>(event)->getBetaInit();
}


double BetaTimeModeProposal::rateParameter(BranchEvent* event)
{
    return static_cast<TraitBranchEvent*>(event)->getBetaShift();
}


bool BetaTimeModeProposal::isTimeVariable(BranchEvent* event)
{
    return static_cast<TraitBranchEvent*>(event)->isTimeVariable();
}


void BetaTimeModeProposal::setEventParameters(BranchEvent* event,
    double initParam, double rateParam, bool isTimeVariable)
{
    TraitBranchEvent* traitEvent = static_cast<TraitBranchEvent*>(event);

    traitEvent->setBetaInit(initParam);
    traitEvent->setBetaShift(rateParam);
    traitEvent->setTimeVariable(isTimeVariable);
}


void BetaTimeModeProposal::setModelParameters()
{
    //std::cout << "DEBUG / in beta time mode proposal" << std::endl;

    _tree->setMeanBranchTraitRates();
}


double BetaTimeModeProposal::rateParameterFromPrior()
{
    return _prior.generateBetaShiftFromPrior();
}
