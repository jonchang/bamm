#ifndef MODEL_H
#define MODEL_H

#include "BranchEvent.h"
#include <set>
#include <iosfwd>

class MbRandom;
class Tree;
class Settings;


double safeExponentiation(double x);


class Model
{

public:

    typedef std::set<BranchEvent*, BranchEvent::PtrCompare> EventSet;

    static double mhColdness;

    Model(MbRandom* ranptr, Tree* tp, Settings* sp);
    virtual ~Model();

    void initializeModelFromEventDataFile();

    void setCurrentLogLikelihood(double x);
    double getCurrentLogLikelihood();

    double computeLogLikelihood();
    double computeLogPrior();

    void incrementGeneration();
    int  getGeneration();
    void resetGeneration();

    // Event handling
    void addEventToTree(double x);
    void addEventToTree();
    void deleteEventFromTree(BranchEvent* be);
    void deleteRandomEventFromTree();

    int getNumberOfEvents();

    BranchEvent* getRootEvent();

    // These functions take a branch event and recursively update
    // branch histories for all nodes going toward the tips
    void forwardSetBranchHistories(BranchEvent* x);
    void forwardSetHistoriesRecursive(Node* p);

    int countEventsInBranchHistory(Node* p);

    // Initialize all branch histories to the root node
    void initializeBranchHistories(Node* x);

    void printStartAndEndEventStatesForBranch(Node* x);

    void printEvents();

    // Move random event
    void eventLocalMove();
    void eventGlobalMove();

    BranchEvent* chooseEventAtRandom();

    // Propose addition or deletion; accept/reject move
    void changeNumberOfEventsMH();
    void moveEventMH();
    void revertMovedEventToPrevious();

    void updateEventRateMH();
    void updateTimeVariablePartitionsMH();

    void restoreLastDeletedEvent();

    bool acceptMetropolisHastings(double lnR);

    // Troubleshooting
    void printBranchHistories(Node* x);

    Tree* getTreePtr();

    double getMHacceptanceRate();
    void   resetMHacceptanceParameters();

    // 0 = last was rejected; 1 = accepted; -1 = not set
    void setAcceptLastUpdate(int x);
    int  getAcceptLastUpdate();

    void   setPoissonRatePrior(double x);
    double getPoissonRatePrior();

    void   setEventRate(double x);
    double getEventRate();

    BranchEvent* getEventByIndex(int x);

    int countTimeVaryingRatePartitions();

    void getEventDataString(std::stringstream& ss);

    bool isEventConfigurationValid(BranchEvent* be);

    // Methods for auto-tuning
    void setMoveSizeScale(double x);
    void setUpdateEventRateScale(double x);

protected:

    void eventMove(bool local);

    void addEventMH();
    void removeEventMH();

    double computeEventGainLogHR(double K, double logLikelihood,
        double oldLogLikelihood, double logPrior, double oldLogPrior,
            double qRatio);

    double computeEventLossLogHR(double K, double logLikelihood,
        double oldLogLikelihood, double logPrior, double oldLogPrior,
            double qRatio);

    void debugLHcalculation();

    virtual void readModelSpecificParameters(std::ifstream& infile) = 0;
    virtual void setRootEventWithReadParameters() = 0;
    virtual BranchEvent* newBranchEventWithReadParameters
        (Node* x, double time) = 0;
    virtual BranchEvent* newBranchEventWithRandomParameters(double x) = 0;
    virtual BranchEvent* newBranchEventFromLastDeletedEvent() = 0;
    virtual void setDeletedEventParameters(BranchEvent* event) = 0;
    virtual double computeSpecificLogLikelihood() = 0;
    virtual double computeSpecificLogPrior() = 0;
    virtual double calculateLogQRatioJump() = 0;
    virtual void setMeanBranchParameters() = 0;
    virtual void getSpecificEventDataString
        (std::stringstream& ss, BranchEvent* event) = 0;

    double _logLikelihood;
    double _eventRate;

    // This parameter hold the density of the new parameters proposed
    // during jump moves. If the parameters are sampled from the prior,
    // these should exactly cancel.
    double _logQRatioJump;

    double _gen;

    MbRandom* _ran;
    Tree* _treePtr;
    Settings* _settings;

    int _acceptCount;
    int _rejectCount;

    EventSet _eventCollection;
    BranchEvent* _rootEvent;

    // True if last generation was accept; false otherwise
    int _acceptLast;

    double _poissonRatePrior;

    // Parameters for MCMC proposals
    double _scale;
    double _updateEventRateScale;
    double _localGlobalMoveRatio;

    double _lastDeletedEventMapTime;

    // Last event modified, whether it is moved, or has value updated
    BranchEvent* _lastEventModified;
};


inline void Model::setCurrentLogLikelihood(double x)
{
    _logLikelihood = x;
}


inline double Model::getCurrentLogLikelihood()
{
    return _logLikelihood;
}


inline void Model::setEventRate(double x)
{
    _eventRate = x;
}


inline double Model::getEventRate()
{
    return _eventRate;
}


inline void Model::incrementGeneration()
{
    _gen++;
}


inline int Model::getGeneration()
{
    return _gen;
}


// Used ater TraitPreBurnin
inline void Model::resetGeneration()
{
    _gen = 0;
}


inline int Model::getNumberOfEvents()
{
    return (int)_eventCollection.size();
}


inline BranchEvent* Model::getRootEvent()
{
    return _rootEvent;
}


inline Tree* Model::getTreePtr()
{
    return _treePtr;
}


inline void Model::setAcceptLastUpdate(int x)
{
    _acceptLast = x;
}


inline int Model::getAcceptLastUpdate()
{
    return _acceptLast;
}


inline void Model::setPoissonRatePrior(double x)
{
    _poissonRatePrior = x;
}


inline double Model::getPoissonRatePrior()
{
    return _poissonRatePrior;
}


inline void Model::setMoveSizeScale(double x)
{
    _scale = x;
}


inline void Model::setUpdateEventRateScale(double x)
{
    _updateEventRateScale = x;
}


#endif
