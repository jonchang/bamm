//
//  Fossil.h
//  
//
//  Created by Dan Rabosky on 12/2/15.
//
//

#ifndef ____Fossil__
#define ____Fossil__

#include <stdio.h>
#include <vector>
#include <string>

class Settings;
class Node;
class Tree;

class Fossil
{
public:
    Fossil(Settings& settings, Tree& tree);
    ~Fossil();
    
    double getCurrentPreservationRate(Node* x, double abstime, double psi);
    double computeOccurrenceLogLikelihood(double psi);
    
    /* double computePreservationLogPrior(double psi)  */
    // could be added to Fossil in the future
    // but for now, with a single psi parameter, we will retain in class Prior.
 
    
    
    std::string getPreservationModel();
    
    inline double getObservationTime();
    
    
private:
    
    std::string _preservationModel;
    
    Settings&  _settings;
    Tree&      _tree;
 
    double     _observationTime;
    double     _numberOccurrences;
    

    
    // options: NONE, one_rate, stage_specific, binomial
    // NONE is for no paleo data.
    
    //Vectors are used to hold stage-specific
    // preservation information. They should be moved to their
    // own class at some point.
    std::vector<std::string> _stagenames;
    std::vector<double> _startTime;
    std::vector<double> _endTime;
    std::vector<double> _relPresRate;
    std::vector<double> _fossilCount;
    
    // Functions
 
    void testFossilParametersFromSettings();
    void errorInvalidParameterValue(std::string xx);
    void invalidParameterValueIgnore(std::string xx);
    void getFossilDataFromFile(std::string fileName);
};


inline double Fossil::getObservationTime()
{
    return _observationTime;
}



#endif /* defined(____Fossil__) */









