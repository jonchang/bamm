#ifndef SP_EX_MODEL_H
#define SP_EX_MODEL_H

#include <set>
#include <vector>

#include "BranchHistory.h"
#include "SpExBranchEvent.h"

//Forward declarations
class Tree;
class Node;
class MbRandom;
class Settings;

double safeExponentiation(double x);
double proportionalShrink(double x, double scale);


class SpExModel
{

    typedef BranchHistory<SpExBranchEvent> SpExBranchHistory;
    typedef std::set<SpExBranchEvent*, SpExBranchEvent::PtrCompare> EventSet;

public:

    static double mhColdness;

    SpExModel(MbRandom* ranptr, Tree* tp, Settings* sp);
    ~SpExModel();

    // Full likelihood will be lnLikTraits + lnLikBranches
    void   setCurrLnLTraits(double x);
    double getCurrLnLTraits();

    void   setCurrLnLBranches(double x);
    double getCurrLnLBranches();

    // Full likelihood function


    // Likelihood functions for branches data
    double computeLikelihoodBranches();
    double computeLikelihoodBranchesByInterval();

    double computeLogPrior();

    double getEventRate();
    void setEventRate(double x);

    double proportionalShrink(double x, double scale);
    bool   acceptMetropolisHastings(const double lnR);

    void incrementGeneration();
    int  getGeneration();
    void resetGeneration();

    // Stuff for event handling:
    void addEventToTree();
    void addEventToTree(double x); // add to specific map position...
    void deleteRandomEventFromTree();
    void deleteEventFromTree(SpExBranchEvent* be);

    void printEvents();

    int getNumberOfEvents();

    SpExBranchEvent* getRootEvent();

    // These functions take a branch event
    // and recursively update branch histories for all nodes
    // going towards the tips
    void forwardSetBranchHistories(SpExBranchEvent* x);
    void forwardSetHistoriesRecursive(Node* p);

    int countEventsInBranchHistory(Node* p);

    // Initialize all branch histories to the root node.
    void initializeBranchHistories(Node* x);

    void printStartAndEndEventStatesForBranch(Node* x);

    void eventLocalMove(SpExBranchEvent* x);   // move specific event
    void eventGlobalMove(SpExBranchEvent* x);  // move specific event
    void eventLocalMove();             // move random event
    void eventGlobalMove();            // move random event
    void revertEventToPreviousPosition();

    // Return random event, or NULL if no events on tree (other than root)
    SpExBranchEvent* chooseEventAtRandom();

    // MCMC:

    // Propose addition or deletion; accept/reject move.
    void changeNumberOfEventsMH();
    void moveEventMH();
    void revertMovedEventToPrevious();

    /*  ***************** */
    // lambda/mu related stuff:
    void updateLambdaInitMH();
    void updateLambdaShiftMH();

    void updateMuInitMH();
    void updateMuShiftMH();

    void updateNodeStateMH();
    void updateNodeStateMH(Node* xnode);
    void updateDownstreamNodeStatesMH(Node* xnode);
    void updateEventRateMH();
    void updateTimeVariablePartitionsMH();

    /*  ********************  */
    void printEventData();

    void restoreLastDeletedEvent();

    // Troubleshooting
    void printBranchHistories(Node* x);

    Tree* getTreePtr();

    // More output: acceptance rates
    double getMHacceptanceRate();
    void   resetMHacceptanceParameters();

    int  getAcceptLastUpdate();
    void setAcceptLastUpdate(int x);

    // 0 = last was rejected; 1 = accepted; -1 = not set.

    void setPoissonRatePrior(double x);
    double getPoissonRatePrior();

    SpExBranchEvent*  getEventByIndex(int x);

    void printExtinctionParams();
    int countTimeVaryingRatePartitions();

    // Generate string with event data:
    void getEventDataString(std::stringstream& ss);

    bool isEventConfigurationValid(SpExBranchEvent* be);

    void initializeModelFromEventDataFile();

    void debugLHcalculation();


	// Functions for auto-tuning
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);
	void setMoveSizeScale(double x);
	void setUpdateEventRateScale(double x);

private:

    // Parameters of the model:

    double lnLikTraits;
    double lnLikBranches;

    // new parameters: March 23 2012
    double _updateLambdaInitScale;
    double _updateLambdaShiftScale;

    double _updateMuInitScale;
    double _updateMuShiftScale;

    // This parameter holds the density of the new
    // parameters proposed during jump moves.
    // If the parameters are sampled from the prior,
    //    these should exactly cancel.

    double _logQratioJump;

    // Root event parameters:

    double _lambdaInit0;
    double _lambdaShift0;

    double _muInit0;
    double _muShift0;

    // Parameters for MCMC proposals
    double _scale; // scale for moving event
    double _updateEventRateScale;
    double _localGlobalMoveRatio;

    double eventLambda; // Poisson rate

    // Priors
    double _poissonRatePrior; // exponential  /* keep this in Model */
    double _lambdaInitPrior;
    double _lambdaShiftPrior;

    double _muInitPrior;
    double _muShiftPrior;

    // Other private variables

    int gen;
    MbRandom* ran;
    Tree* treePtr;
    Settings* sttings;

    int acceptCount;
    int rejectCount;

    // Event collection does not contain the root event
    EventSet eventCollection;
    SpExBranchEvent* rootEvent; // branch event at root node; can't be modified

    double lastDeletedEventMapTime; // map time of last deleted event

    double _lastDeletedEventLambdaInit;
    double _lastDeletedEventLambdaShift;

    double _lastDeletedEventMuInit;
    double _lastDeletedEventMuShift;

    // Here are several variables that track the previous
    // state. At some point, these should have their own class

    // this is a pointer to the last event modified, whether
    // it is moved, or has lambda updated, or whatever.
    SpExBranchEvent*    lastEventModified;

    // General acceptreject flag:
    int acceptLast; // true if last generation was accept; false otherwise

    double _segLength; // for splitting branches
};


inline void SpExModel::setCurrLnLTraits(double x)
{
    lnLikTraits = x;
}


inline double SpExModel::getCurrLnLTraits()
{
    return lnLikTraits;
}


inline void SpExModel::setCurrLnLBranches(double x)
{
    lnLikBranches = x;
}


inline double SpExModel::getCurrLnLBranches()
{
    return lnLikBranches;
}


inline double SpExModel::getEventRate()
{
    return eventLambda;
}


inline void SpExModel::setEventRate(double x)
{
    eventLambda = x;
}


inline void SpExModel::incrementGeneration()
{
    gen++;
}


inline int SpExModel::getGeneration()
{
    return gen;
}


inline void SpExModel::resetGeneration()
{
    gen = 0;   // to be used after TraitPreBurnin
}


inline int SpExModel::getNumberOfEvents()
{
    return (int)eventCollection.size();
}


inline SpExBranchEvent* SpExModel::getRootEvent()
{
    return rootEvent;
}


inline Tree* SpExModel::getTreePtr()
{
    return treePtr;
}


inline int SpExModel::getAcceptLastUpdate()
{
    return acceptLast;
}


inline void SpExModel::setAcceptLastUpdate(int x)
{
    acceptLast = x;
}


inline void SpExModel::setPoissonRatePrior(double x)
{
    _poissonRatePrior  = x;
}


inline double SpExModel::getPoissonRatePrior()
{
    return _poissonRatePrior;
}



inline void SpExModel::setUpdateLambdaInitScale(double x){
	_updateLambdaInitScale = x;
}

inline void SpExModel::setUpdateMuInitScale(double x){
	_updateMuInitScale = x;
}

inline void SpExModel::setUpdateLambdaShiftScale(double x){
	_updateLambdaShiftScale = x;
}

inline void SpExModel::setMoveSizeScale(double x){
	_scale = x;
}

inline void SpExModel::setUpdateEventRateScale(double x){
	_updateEventRateScale = x;
}



#endif
