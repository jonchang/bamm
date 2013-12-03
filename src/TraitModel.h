#ifndef TRAIT_MODEL_H
#define TRAIT_MODEL_H

#include "Model.h"
#include <iosfwd>

class MbRandom;
class Tree;
class Settings;


class TraitModel : public Model
{

public:

    TraitModel(MbRandom* ranptr, Tree* tp, Settings* sp);

    // Trait evolution stuff
    void updateBetaMH();
    void updateNodeStateMH();
    void updateNodeStateMH(Node* xnode);
    void updateBetaShiftMH();
    void updateDownstreamNodeStatesMH(Node* xnode);
    void setMinMaxTraitPriors();

    double getLastLH();

private:

    virtual void readModelSpecificParameters(std::ifstream& infile);
    virtual void setRootEventWithReadParameters();
    virtual BranchEvent* newBranchEventWithReadParameters
        (Node* x, double time);
    virtual BranchEvent* newBranchEventWithRandomParameters(double x);
    virtual BranchEvent* newBranchEventFromLastDeletedEvent();
    virtual void setDeletedEventParameters(BranchEvent* event);
    virtual double computeSpecificLogLikelihood();
    virtual double computeSpecificLogPrior();
    virtual double calculateLogQRatioJump();
    virtual void setMeanBranchParameters();
    virtual void getSpecificEventDataString
        (std::stringstream& ss, BranchEvent* event);

    double computeTriadLikelihoodTraits(Node* x);

    double _updateBetaScale;
    double _updateBetaShiftScale;
    double _updateNodeStateScale;

    double _betaInitPrior;
    double _betaShiftPrior;

    double _lastDeletedEventBetaInit;;
    double _lastDeletedEventBetaShift;

    double _lastLH;

    double _readBetaInit;
    double _readBetaShift;
};


inline double TraitModel::getLastLH()
{
    return _lastLH;
}


#endif
