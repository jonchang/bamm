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


MCMC::MCMC(MbRandom* rng, Model* model, Settings* settings) :
    _rng(rng), _model(model), _settings(settings)
{
    for (int i = 0; i < _settings->getInitialNumberEvents(); i++) {
        _model->addEventToTree();
    }

    setUpdateWeights();
}


MCMC::~MCMC()
{
    _mcmcOutputStream.close();
    _eventDataOutputStream.close();
}


void MCMC::run()
{
    log() << "\nRunning MCMC chain for "
          << _numGenerations << " generations.\n";

    log() << "\n";
    outputHeaders();

    int parameterToUpdate = 0;
    for (int gen = 0; gen < _numGenerations; gen++) {
        parameterToUpdate = chooseRandomParameter();
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
    outputStdOutHeaders();
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


void MCMC::outputStdOutHeaders()
{
    log() << std::setw(15) << "Generation"
          << std::setw(15) << "LogLikelihood"
          << std::setw(15) << "NumShifts"
          << std::setw(15) << "LogPrior"
          << std::setw(15) << "AcceptRate" << "\n";
}


void MCMC::outputData(int generation)
{
    if (generation % _mcmcOutputFreq == 0) {
        outputMCMCData();
    }

    if (generation % _eventDataOutputFreq == 0) {
        outputEventData();
    }

    if (generation % _stdOutFreq == 0) {
        outputStdOutData();

        // Reset acceptance rates when state data are printed to the screen.
        // This could lead to NANs if outputEventData() is called after this.
        _model->resetMHAcceptanceParameters();
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
}


// TODO: Perhaps the model should be printing this directly
void MCMC::outputEventData()
{
    std::stringstream eventData;
    _model->getEventDataString(eventData);
    _eventDataOutputStream << eventData.str() << std::endl;
}


void MCMC::outputStdOutData()
{
    log() << std::setw(15) << _model->getGeneration()
          << std::setw(15) << _model->getCurrentLogLikelihood()
          << std::setw(15) << _model->getNumberOfEvents()
          << std::setw(15) << _model->computeLogPrior()
          << std::setw(15) << _model->getMHAcceptanceRate()
          << std::endl;
}
