/*
 *  Node.h
 *  proj7
 *
 *  Created by Dan Rabosky on 9/16/11.
  *
 */

#ifndef NODE_H
#define NODE_H

#include <string>

class branchEvent;
class eventSet;
class Phenotype;

class BranchHistory;


class Node
{

private:

    void init(int x = 0);

    Node*  _lfDesc;
    Node*  _rtDesc;
    Node*  _anc;
    std::string _name;

    std::string _cladeName;

    int    _index;
    double _time;
    double _brlen;
    double _branchTime;

    int  _tipDescCount;
    bool _isExtant;
    bool _isTip;
    bool _isConstant;
    bool _isLivingTip;

    // specific for mapping:
    double _mapStart;
    double _mapEnd;

    BranchHistory* _history;

    // For phenotypes:
    double _trait; // trait value
    double _meanBeta; // mean phenotypic rate
    double _nodeBeta; // exact value at node.
    bool   _isTraitFixed; // is trait value a free parameter?
    double _netJump;

    //specific stuff for compound poisson rate model
    double _meanSpeciationRate;
    double _meanExtinctionRate;

    // Node rates for time-varying models:
    double _nodeLambda;
    double _nodeMu;

    double _di;   // initial value of speciation probability (at node)
    double _ei;   // initial value of extinction probability (at node)
    double _etip; // initial value (sampling frac.) at tip descended from node
    
    // value of log-likelihood at end of branch
    // (after combining w speciation event):
    //  Zero if node is terminal.
    double _nodeLikelihood;

    // Flag for whether node can or cannot define branch that can hold event:
    bool _canHoldEvent;
    
    // HOlds value of Extinction probability at end of branch
    double _eEnd; // remove
    
    // determines whether the descendants of a node have
    // a rates shift somewhere in their history
    
    bool _hasDownstreamRateShift;
    
    bool _inheritFromLeft;
    
    // For experimental speedup: USE_FAST macro
    double _logDi_End;   // remove
    double _logDi_End_last; // remove
    
    bool   _nodeParamsAreCurrent;  // remove
    bool   _lastNodeParamsAreCurrent;  // remove
    
    double _eEndLast; // remove
    
    double _logDi_proposed;
    double _logDi_current;
    
    double _exProb_proposed;
    double _exProb_current;
    
    bool _proposedUpdate;
    
    double _previousLambda;

public:

    Node();
    Node(int x);

    void  setLfDesc(Node* x);
    Node* getLfDesc();

    void  setRtDesc(Node* x);
    Node* getRtDesc();

    void nullifyLfDesc();
    void nullifyRtDesc();
    void nullifyAnc();

    void  setAnc(Node* x);
    Node* getAnc();

    void   setName(std::string x);
    std::string getName();

    void setIndex(int x);
    int  getIndex();

    void   setTime(double x);
    double getTime();

    void   setBrlen(double x);
    double getBrlen();

    void setTipDescCount(int x);
    int  getTipDescCount();

    int getDescCount();

    void setExtantStatus(bool x);
    bool getExtantStatus();

    void setIsTip(bool x);
    bool getIsTip();

    bool isInternal();

    void setIsConstant(bool x);
    bool getIsConstant();

    void setIsLivingTip(bool x);
    bool getIsLivingTip();

    // Get a random TIP node descended from lfdesc of a given node
    Node* getRandomLeftTipNode();

    // Get random RIGHT TIP node
    Node* getRandomRightTipNode();

    double pathLengthToRoot();

    // Specific for treemap:
    void   setMapStart(double x);
    double getMapStart();

    void   setMapEnd(double x);
    double getMapEnd();

    //Need to includet this
    BranchHistory* getBranchHistory();

    void   setMeanSpeciationRate(double x);
    double getMeanSpeciationRate();

    void   setMeanExtinctionRate(double x);
    double getMeanExtinctionRate();

    void computeNodeBranchSpeciationParams();
    void computeNodeBranchExtinctionParams();
    
    void computeAndSetNodeSpeciationParams();
    void computeAndSetNodeExtinctionParams();

    void   setNodeLambda(double x);
    double getNodeLambda();

    void   setNodeMu(double x);
    double getNodeMu();

    void   setNodeBeta(double x);
    double getNodeBeta();

    /*********/

    // Phenotypic stuff:
    void setTraitValue(double x);
    double getTraitValue();

    // Jump model
    void setNetJump(double x);
    double getNetJump();

    void setMeanBeta(double x);
    double getMeanBeta();

    void setIsTraitFixed(bool x);
    bool getIsTraitFixed();

    // Speciation-extinction calculations
    void   setEinit(double x);
    double getEinit();

    void   setDinit(double x);
    double getDinit();

    void   setEtip(double x);
    double getEtip();

    void setNodeLikelihood(double x);
    double getNodeLikelihood();

    bool getCanHoldEvent();
    void setCanHoldEvent(bool x);

    void   setBranchTime(double x);
    double getBranchTime();

    void setCladeName(std::string x);
    std::string getCladeName();

    double computeSpeciationRateIntervalRelativeTime
        (double tstart, double tstop);
    double computeSpeciationRateIntervalAbsoluteTime
        (double tstart, double tstop);
    double computeExtinctionRateIntervalRelativeTime
        (double tstart, double tstop);
    double computeSpeciationRateIntervalRelativeTime(double t_init, double tstart, double tstop);
    double computeSpeciationRateIntervalAbsoluteTime(double t_init, double tstart, double tstop);
    double computeExtinctionRateIntervalRelativeTime(double t_init, double tstart, double tstop);
    double getPointExtinction(double branchtime);
    
    double integrateExponentialRateFunction(double par_init, double shift, double t1, double t2);
    double getExponentialRate(double par_init, double shift, double tm);


    std::string getRandomRightDesc();
    std::string getRandomLeftDesc();
    
    double getExtinctionEnd(void);
    void setExtinctionEnd(double x);
    
    bool getHasDownstreamRateShift(void);
    void setHasDownstreamRateShift(bool x);
    
    bool getInheritFromLeft(void);
    void setInheritFromLeft(bool x);
    
    // USE_FAST macro options
    void setLogDiEnd(double x);
    double getLogDiEnd();
    
    void setLogDiEndLast(double x);
    double getLogDiEndLast();
    
    double getExtinctionEndLast(void);
    void setExtinctionEndLast(double x);
    
    bool getLastNodeParamsAreCurrent();
    void setLastNodeParamsAreCurrent(bool x);
    
    bool getNodeParamsAreCurrent();
    void setNodeParamsAreCurrent(bool x);
    
    bool getProposedUpdate();
    void setProposedUpdate(bool);
    
    double getLogDiProposed();
    void setLogDiProposed(double);
    
    double getLogDiCurrent();
    void setLogDiCurrent(double);
    
    double getExProbProposed();
    void setExProbProposed(double);
    
    double getExProbCurrent();
    void setExProbCurrent(double);
    
    double getPreviousNodeLambda();
    void setPreviousNodeLambda(double x);
    
};

// USE_FAST macro functions

inline double Node::getPreviousNodeLambda()
{
    return _previousLambda;
}

inline void Node::setPreviousNodeLambda(double x)
{
    _previousLambda = x;
}

inline bool Node::getProposedUpdate()
{
    return _proposedUpdate;
}

inline void Node::setProposedUpdate(bool x)
{
    _proposedUpdate = x;
}

inline double Node::getLogDiProposed()
{
    return _logDi_proposed;
}

inline void Node::setLogDiProposed(double x)
{
    _logDi_proposed = x;
}

inline double Node::getLogDiCurrent()
{
    return _logDi_current;
}

inline void Node::setLogDiCurrent(double x)
{
    _logDi_current = x;
}

inline double Node::getExProbProposed()
{
    return _exProb_proposed;
}

inline void Node::setExProbProposed(double x)
{
    _exProb_proposed = x;
}

inline double Node::getExProbCurrent()
{
    return _exProb_current;
}

inline void Node::setExProbCurrent(double x)
{
    _exProb_current = x;
}

// Preceding block of functions OK... following block deprecated


inline void Node::setLogDiEnd(double x)
{
    _logDi_End = x;
}

inline double Node::getLogDiEnd()
{
    return _logDi_End;
}

inline void Node::setLogDiEndLast(double x)
{
    _logDi_End_last = x;
}

inline double Node::getLogDiEndLast()
{
    return _logDi_End_last;
}


inline void Node::setLastNodeParamsAreCurrent(bool x)
{
    _lastNodeParamsAreCurrent = x;
}

inline bool Node::getLastNodeParamsAreCurrent()
{
    return _lastNodeParamsAreCurrent;
}

inline void Node::setNodeParamsAreCurrent(bool x)
{
    _nodeParamsAreCurrent = x;
}

inline bool Node::getNodeParamsAreCurrent()
{
    return _nodeParamsAreCurrent;
}

inline void Node::setExtinctionEndLast(double x)
{
    _eEndLast = x;
}

inline double Node::getExtinctionEndLast(void)
{
    return _eEndLast;
}


// preceding block deprecated



inline void Node::setLfDesc(Node* x)
{
    _lfDesc = x;
}


inline void Node::setRtDesc(Node* x)
{
    _rtDesc = x;
}


inline void Node::nullifyLfDesc()
{
    _lfDesc = NULL;
}


inline void Node::nullifyRtDesc()
{
    _rtDesc = NULL;
}


inline void Node::nullifyAnc()
{
    _anc = NULL;
}


inline void Node::setAnc(Node* x)
{
    _anc = x;
}


inline void Node::setName(std::string x)
{
    _name = x;
}


inline void Node::setIndex(int x)
{
    _index = x;
}


inline void Node::setTime(double x)
{
    _time = x;
}


inline void Node::setBrlen(double x)
{
    _brlen = x;
}


inline void Node::setTipDescCount(int x)
{
    _tipDescCount = x;
}


inline Node* Node::getLfDesc()
{
    return _lfDesc;
}


inline Node* Node::getRtDesc()
{
    return _rtDesc;
}


inline Node* Node::getAnc()
{
    return _anc;
}


inline std::string Node::getName()
{
    return _name;
}


inline int Node::getIndex()
{
    return _index;
}


inline double Node::getTime()
{
    return _time;
}


inline double Node::getBrlen()
{
    return _brlen;
}


inline int Node::getTipDescCount()
{
    return _tipDescCount;
}


inline void Node::setExtantStatus(bool x)
{
    _isExtant = x;
}


inline bool Node::getExtantStatus()
{
    return _isExtant;
}


inline void Node::setIsTip(bool x)
{
    _isTip = x;
}


inline bool Node::getIsTip()
{
    return _isTip;
}


inline void Node::setIsConstant(bool x)
{
    _isConstant = x;
}


inline bool Node::getIsConstant()
{
    return _isConstant;
}


inline void Node::setIsLivingTip(bool x)
{
    _isLivingTip = x;
}


inline bool Node::getIsLivingTip()
{
    return _isLivingTip;
}


inline void Node::setMapStart(double x)
{
    _mapStart = x;
}


inline double Node::getMapStart()
{
    return _mapStart;
}


inline void Node::setMapEnd(double x)
{
    _mapEnd = x;
}


inline double Node::getMapEnd()
{
    return _mapEnd;
}


inline BranchHistory* Node::getBranchHistory()
{
    return _history;
}


inline void Node::setMeanSpeciationRate(double x)
{
    _meanSpeciationRate = x;
}


inline double Node::getMeanSpeciationRate()
{
    return _meanSpeciationRate;
}


inline void Node::setMeanExtinctionRate(double x)
{
    _meanExtinctionRate = x;
}


inline double Node::getMeanExtinctionRate()
{
    return _meanExtinctionRate;
}


inline void Node::setNodeLambda(double x)
{
    _nodeLambda = x;
}


inline double Node::getNodeLambda()
{
    return _nodeLambda;
}


inline void Node::setNodeMu(double x)
{
    _nodeMu = x;
}


inline double Node::getNodeMu()
{
    return _nodeMu;
}


inline void Node::setNodeBeta(double x)
{
    _nodeBeta = x;
}


inline double Node::getNodeBeta()
{
    return _nodeBeta;
}


inline void Node::setTraitValue(double x)
{
    _trait = x;
}


inline double Node::getTraitValue()
{
    return _trait;
}


inline void Node::setMeanBeta(double x)
{
    _meanBeta = x;
}


inline double Node::getMeanBeta()
{
    return _meanBeta;
}


inline void Node::setIsTraitFixed(bool x)
{
    _isTraitFixed = x;
}


inline bool Node::getIsTraitFixed()
{
    return _isTraitFixed;
}


inline void Node::setEinit(double x)
{
    _ei = x;
}


inline double Node::getEinit()
{
    return _ei;
}


inline void Node::setDinit(double x)
{
    _di = x;
}


inline double Node::getDinit()
{
    return _di;
}


inline void Node::setEtip(double x)
{
    _etip = x;
}


inline double Node::getEtip()
{
    return _etip;
}


inline void Node::setNodeLikelihood(double x)
{
    _nodeLikelihood = x;
}


inline double Node::getNodeLikelihood()
{
    return _nodeLikelihood;
}


inline bool Node::getCanHoldEvent()
{
    return _canHoldEvent;
}


inline void Node::setCanHoldEvent(bool x)
{
    _canHoldEvent = x;
}


inline void Node::setBranchTime(double x)
{
    _branchTime = x;
}


inline double Node::getBranchTime()
{
    return _branchTime;
}


inline void Node::setCladeName(std::string x)
{
    _cladeName = x;
}


inline std::string Node::getCladeName()
{
    return _cladeName;
}

inline void Node::setExtinctionEnd(double x)
{
    _eEnd = x;
}

inline double Node::getExtinctionEnd(void)
{
    return _eEnd;
}


inline void Node::setHasDownstreamRateShift(bool x)
{
    _hasDownstreamRateShift = x;
}

inline bool Node::getHasDownstreamRateShift(void)
{
    return _hasDownstreamRateShift;
}

inline bool Node::getInheritFromLeft(void)
{
    return _inheritFromLeft;
}

inline void Node::setInheritFromLeft(bool x)
{
    _inheritFromLeft = x;
}

inline void Node::setNetJump(double x)
{
    _netJump = x;
}

inline double Node::getNetJump()
{
    return _netJump;
}


#endif
