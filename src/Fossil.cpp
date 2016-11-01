//
//  Fossil.cpp
//  
//
//  Created by Dan Rabosky on 12/2/15.
//
//

#include "Fossil.h"
#include "Settings.h"
#include "Tree.h"
#include "Node.h"

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

Fossil::Fossil(Settings& settings, Tree& tree):
    _settings(settings), _tree(tree)
{
 
    _preservationModel = _settings.get<std::string>("preservationModel");

    _observationTime = _settings.get<double>("observationTime");
    if (_observationTime <= 0){
        _observationTime = _tree.getAge();
    }else if ( _observationTime < _tree.getAge() ){
        std::cout << "WARNING: invalid initial observation time" << std::endl;
        std::cout << "\t... setting observationTime to tree MAX TIME" << std::endl;
        _observationTime = _tree.getAge();
    }
    
    
    testFossilParametersFromSettings();
    
    if (_preservationModel == "NONE"){
        if (! _tree.isUltrametric()){
            std::cout << "Tree must be ultrametric if no fossil data" << std::endl;
            std::cout << "Exiting...." << std::endl;
            exit(0);
        }
    }else if (_preservationModel == "one_rate"){
    
        _numberOccurrences = _settings.get<double>("numberOccurrences");
        
    }else if (_preservationModel == "stage_specific"){
    
        getFossilDataFromFile(_settings.get("fossilDataFile"));
    
    }else if (_preservationModel == "binomial"){
        std::cout << "binomial sampling model not yet implemented" << std::endl;
        exit(0);
    }else{
        std::cout << "Invalid preservationModel << " << _preservationModel;
        std::cout << ">>" << std::endl;
        exit(0);
    }
    
 
    
}


void Fossil::testFossilParametersFromSettings()
{
    
    double updateRatePreservationRate = _settings.get<double>("updateRatePreservationRate");
    double preservationRateInit =
    _settings.get<double>("preservationRateInit");

    if (_preservationModel == "NONE"){
        if (updateRatePreservationRate > 0)
            invalidParameterValueIgnore("updateRatePreservationRate");
        if (preservationRateInit > 0)
            invalidParameterValueIgnore("preservationRateInit");
    }
    
    if (_preservationModel == "one_rate"){
        if (preservationRateInit <= 0)
            errorInvalidParameterValue("preservationRateInit");
    }
    
    if (_preservationModel == "stage_specific"){
        if (preservationRateInit <= 0)
            errorInvalidParameterValue("preservationRateInit");

        if (_settings.get("fossilDataFile") == ""){
            std::cout << "Invalid input file for fossil occurrences\n";
            std::cout << "This must be specified if <<preservationModel>>\n";
            std::cout << "is stage_specific" << std::endl;
            exit(0);
        }
    
    }

}

void Fossil::errorInvalidParameterValue(std::string xx)
{
    std::cout << "ERROR: parameter << " << xx;
    std::cout << " >> has invalid value for the specified";
    std::cout << "\npreservation model. Exiting." << std::endl;
    exit(0);
}

void Fossil::invalidParameterValueIgnore(std::string xx)
{
    std::cout << "WARNING: parameter << " << xx;
    std::cout << " >> has invalid value for the specified";
    std::cout << "\npreservation model" << std::endl;
    std::cout << "This parameter will be ignored" << std::endl;
}


// Evalates settings and tree to determine if
// it is a valid instance of a tree with some paleontological data.
// Sets the _hasPaleoData parameter.
/*
void Fossil::initializeHasPaleoData()
{
    
    double updateRatePreservationRate = _settings.get<double>("updateRatePreservationRate");
    
    if (_preservationModel == "one_rate"){
        
        _numberOccurrences = _settings.get<int>("numberOccurrences");
        
        if (_numberOccurrences > 0 & _settings.get<double>("preservationRateInit") < 0.000000001){
            std::cout << "Invalid initial settings " << std::endl;
            std::cout << " cannot have <<numberOccurrences>> greater than 0 and " << std::endl;
            std::cout << " <<preservationRateInit>> equal to zero.";
            std::cout << "Check control file " << std::endl;
            exit(0);
            
        }
    }
    
    if (_numberOccurrences == 0){
        _hasPaleoData = false;
        _observationTime = _tree->getAge();
        
        if (getTreePtr()->isUltrametric() == false){
            std::cout << "Tree must be ultrametric if no fossil data" << std::endl;
            std::cout << "Exiting...." << std::endl;
            exit(0);
        }
        
        if (updateRatePreservationRate > 0.00000001){
            std::cout << "Attempt to set preservation rate for non-paleo data" << std::endl;
            std::cout << "This parameter will be ignored..." << std::endl;
        }
        
    }else if (_numberOccurrences > 0){
        
        _hasPaleoData = true;
 
        _observationTime = _settings.get<double>("observationTime");
        if (_observationTime <= 0){
            _observationTime = _tree->getAge();
        }else if ( _observationTime < _tree->getAge() ){
            std::cout << "WARNING: invalid initial observation time" << std::endl;
            std::cout << "\t... setting observationTime to tree MAX TIME" << std::endl;
            _observationTime = _tree->getAge();
        }
        
    }else if (_numberOccurrences == -1){
   
        
        _hasPaleoData = true;
        _observationTime = _settings.get<double>("observationTime");
        if (_observationTime <= 0){
            _observationTime = _tree->getAge();
        }else if ( _observationTime < _tree->getAge() ){
            std::cout << "WARNING: invalid initial observation time" << std::endl;
            std::cout << "\t... setting observationTime to tree MAX TIME" << std::endl;
            _observationTime = _tree->getAge();
        }
        
        if (_settings.get("fossilDataFile") == ""){
            std::cout << "Invalid input file for fossil occurrences" << std::endl;
            std::cout << "This must be specified if <<numberOccurrences>> = -1" << std::endl;
            exit(0);
        }else{
            
            getFossilDataFromFile(_settings.get("fossilDataFile"));
            
        }
        
        
    }else{
        std::cout << "Invalid number of occurrences in controlfile" << std::endl;
        exit(0);
    }
    
    _hasMassExtinctionData = _settings.get<bool>("hasMassExtinctions");
    
    if (_hasMassExtinctionData){
     
        std::cout << "Error in SpExModel::initializeHasPaleoData" << std::endl;
        std::cout << "Parameter <<_hasMassExtinctionData>> has been deprecated" << std::endl;
        
        getMassExtinctionDataFromFile();
    }
    
    
}
*/

void Fossil::getFossilDataFromFile(std::string fileName)
{
    
    
    std::ifstream inputFile(fileName.c_str());
    
    if (!inputFile) {
        log(Error) << "Could not read data from file "
        << "<<" << fileName << ">>.\n";
        std::exit(1);
    }
    std::cout << "in class Fossil: " << std::endl;
    log() << "Reading fossil data from file <<" << fileName << ">>.\n";
    
    while (inputFile){
        std::string tempstring;
        getline(inputFile, tempstring, '\t');
        _stagenames.push_back(tempstring);
        getline(inputFile, tempstring, '\t');
        _startTime.push_back(atof(tempstring.c_str()));
        getline(inputFile, tempstring, '\t');
        _endTime.push_back(atof(tempstring.c_str()));
        getline(inputFile, tempstring, '\t');
        _relPresRate.push_back(atof(tempstring.c_str()));
        getline(inputFile, tempstring, '\n');
        _fossilCount.push_back(atof(tempstring.c_str()));
 
        if (inputFile.peek() == EOF) {
            break;
        }
    }
    
    inputFile.close();
    
    if (_endTime[(int)_endTime.size()-1] < _observationTime){
        _endTime[(int)_endTime.size()-1] = _observationTime;
        
        std::cout << "WARNING: preservation file should have final time bin" << std::endl;
        std::cout << "consistent with <<observationTime>>" << std::endl;
        std::cout << "Resetting end time of final bin to <<observationTime>> " << std::endl;
        
    }
    
    
}



double Fossil::computeOccurrenceLogLikelihood(double psi)
{

    double loglik = 0.0;
    
    if (_preservationModel == "NONE"){
        // no change to loglik;
        
    }else if (_preservationModel == "one_rate"){
        
        loglik = _numberOccurrences * std::log(psi);
        
    }else if (_preservationModel == "stage_specific"){
    
        for (int i = 0; i < (int)_relPresRate.size(); i++){
            if (_fossilCount[i] > 0){
                
                double logprate = std::log(psi) + std::log(_relPresRate[i]);
                loglik += logprate * (double)_fossilCount[i];
            
            }
        }
        
        
    }else{
        std::cout << "Problem in Fossil::computeOccurrenceLogLikelihood()" << std::endl;
        exit(0);
    }
    
    return loglik;

}



double Fossil::getCurrentPreservationRate(Node* x, double abstime, double psi)
{
    double presrate = 0.0;
    
    if (_preservationModel == "NONE"){
        presrate = 0.0;
    }else if (_preservationModel == "one_rate"){
        
        presrate = psi;
        
    }else if (_preservationModel == "stage_specific"){
        
        for (int i = 0; i < (int)_startTime.size(); i++){
            if (abstime >= _startTime[i] & abstime <= _endTime[i]){
                presrate = _relPresRate[i] * psi;
                break;
            }
        }
        
        
    }else{
        std::cout << "Error in Fossil::getCurrentPreservationRate" << std::endl;
        exit(0);
    }
 
    return presrate;
    
}
















