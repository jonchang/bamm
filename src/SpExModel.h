#ifndef SP_EX_MODEL_H
#define SP_EX_MODEL_H

#include "Model.h"
#include <iosfwd>

class Node;
class MbRandom;
class Settings;
class Prior;
class BranchEvent;


class SpExModel : public Model
{

public:

    SpExModel(MbRandom* rng, Settings* settings, Prior* prior);

    virtual double computeLogLikelihood();
    
    virtual double computeLogPrior();

    void updateLambdaInitMH();
    void updateLambdaShiftMH();
    void updateMuInitMH();
    void updateMuShiftMH();

	// Methods for auto-tuning
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);

private:

    double computeLogLikelihoodByInterval();
    
    virtual void initializeTree();

    virtual void readModelSpecificParameters(std::ifstream& inputFile);
    virtual void setRootEventWithReadParameters();

    virtual BranchEvent* newBranchEventWithReadParameters(Node* x, double time);
    virtual BranchEvent* newBranchEventWithRandomParameters(double x);
    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    virtual void setMeanBranchParameters();
    virtual void setDeletedEventParameters(BranchEvent* be);

    virtual double calculateLogQRatioJump();

    virtual void getSpecificEventDataString
        (std::stringstream& ss, BranchEvent* event);

    double _updateLambdaInitScale;
    double _updateLambdaShiftScale;
    double _updateMuInitScale;
    double _updateMuShiftScale;

    // Root event parameters
    double _lambdaInit0;
    double _lambdaShift0;
    double _muInit0;
    double _muShift0;

    // Priors
    double _lambdaInitPrior;
    double _lambdaShiftPrior;
    double _muInitPrior;
    double _muShiftPrior;

    double _lastDeletedEventLambdaInit;
    double _lastDeletedEventLambdaShift;
    double _lastDeletedEventMuInit;
    double _lastDeletedEventMuShift;

    double _segLength;

    double _readLambdaInit;
    double _readLambdaShift;
    double _readMuInit;
    double _readMuShift;
};


inline void SpExModel::setUpdateLambdaInitScale(double x)
{
	_updateLambdaInitScale = x;
}


inline void SpExModel::setUpdateMuInitScale(double x)
{
	_updateMuInitScale = x;
}


inline void SpExModel::setUpdateLambdaShiftScale(double x)
{
	_updateLambdaShiftScale = x;
}


#endif
