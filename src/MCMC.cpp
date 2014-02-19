#include "MCMC.h"
#include "MbRandom.h"
#include "Model.h"
#include "Settings.h"
#include "Log.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <string>


MCMC::MCMC(MbRandom* rng, Model* model, Settings* settings, int chain = 0) :
    _rng(rng), _model(model), _settings(settings), _chain(chain)
{
    // MCMC-related output
    _mcmcOutputFileName = _settings->getMCMCoutfile();
    _mcmcOutputFreq = _settings->getMCMCwriteFreq();

    // Event-related output
    _eventDataOutputFileName = _settings->getEventDataOutfile();
    _eventDataOutputFreq = _settings->getEventDataWriteFreq();

    // Append the chain number (if > 0) to the output file names
    if (_chain > 0) {
        _mcmcOutputFileName += ("-" + std::to_string(chain));
        _eventDataOutputFileName += ("-" + std::to_string(chain));
    }

    // Open streams for output
    _mcmcOutputStream.open(_mcmcOutputFileName.c_str());
    _eventDataOutputStream.open(_eventDataOutputFileName.c_str());

    for (int i = 0; i < _settings->getInitialNumberEvents(); i++) {
        _model->addEventToTree();
    }
}


void MCMC::finishConstruction()
{
    setUpdateWeights();
    outputHeaders();
}


MCMC::~MCMC()
{
    _mcmcOutputStream.close();
    _eventDataOutputStream.close();
}


void MCMC::run(int nGenerations)
{
    for (int gen = 0; gen < nGenerations; gen++) {
        int parameterToUpdate = chooseRandomParameter();
        updateState(parameterToUpdate);
        outputData(gen);
    }
}


void MCMC::setUpdateWeights()
{
    setUpParameterWeights();

    double sumWeights = _parameterWeights[0];
    for (SizeType i = 1; i < _parameterWeights.size(); i++) {
        sumWeights += _parameterWeights[i];
        _parameterWeights[i] += _parameterWeights[i - 1];
    }

    for (SizeType i = 0; i < _parameterWeights.size(); i++) {
        _parameterWeights[i] /= sumWeights;
    }

    // Define vectors to hold accept/reject data:
    for (SizeType i = 0; i < _parameterWeights.size(); i++) {
        _acceptCount.push_back(0);
        _rejectCount.push_back(0);
    }
}


void MCMC::setUpParameterWeights()
{
    _parameterWeights.push_back(_settings->getUpdateRateEventNumber());
    _parameterWeights.push_back(_settings->getUpdateRateEventPosition());
    _parameterWeights.push_back(_settings->getUpdateRateEventRate());

    // Defined by concrete subclass
    setUpSpecificParameterWeights();
}


int MCMC::chooseRandomParameter()
{
    double r = _rng->uniformRv();

    for (SizeType i = 0; i < _parameterWeights.size(); i++) {
        if (r < _parameterWeights[i]) {
            return i;
        }
    }

    return -1;
}


void MCMC::updateState(int parameter)
{
    if (parameter == 0) {
        _model->changeNumberOfEventsMH();
    } else if (parameter == 1) {
        _model->moveEventMH();
    } else if (parameter == 2) {
        _model->updateEventRateMH();
    } else if (parameter > 2) {
        // Defined in concrete subclass
        updateSpecificState(parameter);
    } else {
        // Should never get here
        log(Error) << "Bad parameter to update\n";
        std::exit(1);
    }

    if (_model->getAcceptLastUpdate() == 1) {
        _acceptCount[parameter]++;
    } else if (_model->getAcceptLastUpdate() == 0) {
        _rejectCount[parameter]++;
    } else if (_model->getAcceptLastUpdate() == -1) {
        log(Error) << "Failed somewhere in MH step, parameter "
                   << parameter << "\n";
        std::exit(1);
    } else {
        log(Error) << "Invalid accept/reject flag in model object\n";
        std::exit(1);
    }

    // Reset to unmodified value
    _model->setAcceptLastUpdate(-1);
}


void MCMC::outputHeaders()
{
    outputMCMCHeaders();
    outputEventDataHeaders();
}


void MCMC::outputMCMCHeaders()
{
    _mcmcOutputStream << "generation,N_shifts,logPrior,logLik,"
                      << "eventRate,acceptRate\n";
}


void MCMC::outputEventDataHeaders()
{
    _eventDataOutputStream << "generation,leftchild,rightchild,abstime";
    outputSpecificEventDataHeaders();
}


void MCMC::outputData(int generation)
{
    if (generation % _mcmcOutputFreq == 0) {
        outputMCMCData();
    }

    if (generation % _eventDataOutputFreq == 0) {
        outputEventData();
    }

    // Defined in concrete subclass
    outputSpecificData(generation);
}


void MCMC::outputMCMCData()
{
    _mcmcOutputStream << _model->getGeneration()            << ","
                      << _model->getNumberOfEvents()        << ","
                      << _model->computeLogPrior()          << ","
                      << _model->getCurrentLogLikelihood()  << ","
                      << _model->getEventRate()             << ","
                      << _model->getMHAcceptanceRate()
                      << std::endl;
    _model->resetMHAcceptanceParameters();
}


// TODO: Perhaps the model should be printing this directly
void MCMC::outputEventData()
{
    std::stringstream eventData;
    _model->getEventDataString(eventData);
    _eventDataOutputStream << eventData.str() << std::endl;
}
