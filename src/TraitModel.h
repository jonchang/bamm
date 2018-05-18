#ifndef TRAIT_MODEL_H
#define TRAIT_MODEL_H


#include "Model.h"

#include <iosfwd>
#include <vector>
#include <string>

class Node;
class Random;
class Settings;
class BranchEvent;
class Proposal;


class TraitModel : public Model
{

public:

    TraitModel(Random& rng, Settings& settings);

    virtual double computeLogLikelihood();
    virtual double computeTriadLikelihoodTraits(Node* x);

    virtual double computeLogPrior();

    virtual void checkModel();
    
    double getJumpVariance();
    // void setJumpVariance(double x);
    
    int getNumberOfJumpEvents();
    int getNumberOfRateShiftEvents();
    double getRootState();

private:

    virtual void setRootEventWithReadParameters
        (const std::vector<std::string>& parameters);
    virtual BranchEvent* newBranchEventWithReadParameters
        (Node* x, double time, const std::vector<std::string>& parameters);

    double betaInitParameter(const std::vector<std::string>& parameters);
    double betaShiftParameter(const std::vector<std::string>& parameters);

    virtual BranchEvent* newBranchEventWithRandomParameters(double x);
    virtual BranchEvent* newBranchEventWithParametersFromSettings(double x);
    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    virtual void setMeanBranchParameters();
    virtual void setMeanBranchParameters(Node* x);
    virtual void setDeletedEventParameters(BranchEvent* be);

    virtual double calculateLogQRatioJump();

    virtual void getSpecificEventDataString
        (std::stringstream& ss, BranchEvent* event);

    bool _sampleFromPriorOnly;

    double _lastDeletedEventBetaInit;;
    double _lastDeletedEventBetaShift;
    bool _lastDeletedEventTimeVariable;
    
    bool _lastDeletedEventIsJump;
    double _lastDeletedEventJump;

    // Here are several variables that track the previous
    // state. At some point, these should have their own class

    // Ultimately, initializations should be handled in TraitModel
    // and Model classes
    // NOT in class TREE!!
    void   initializeTraitParamsForNodes();

    double _readBetaInit;
    double _readBetaShift;

    // hierarchical model parameter for jumps
    double _jumpVariance;
};


/*
inline void TraitModel::setJumpVariance(double x)
{
    _jumpVariance = x;
}
*/

inline double TraitModel::getJumpVariance()
{
    return _jumpVariance;
}


#endif
