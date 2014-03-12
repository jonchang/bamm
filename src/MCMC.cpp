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
    _numGenerations = _settings->getNGENS();

    // Output file names
    _mcmcOutputFileName      = _settings->getMCMCoutfile();
    _eventDataOutputFileName = _settings->getEventDataOutfile();

    // Acceptance info
    _outputAcceptanceInfo = _settings->getOutputAcceptanceInfo();
    _acceptanceInfoFileName = _settings->getAcceptanceInfoFileName();

    // Output frequencies
    _mcmcOutputFreq      = _settings->getMCMCwriteFreq();
    _eventDataOutputFreq = _settings->getEventDataWriteFreq();
    _stdOutFreq          = _settings->getPrintFreq();

    // Open streams for writing
    _mcmcOutputStream.open(_mcmcOutputFileName.c_str());
    _eventDataOutputStream.open(_eventDataOutputFileName.c_str());

    if (_outputAcceptanceInfo) {
        _acceptanceInfoStream.open(_acceptanceInfoFileName.c_str());
    }
}


MCMC::~MCMC()
{
    _mcmcOutputStream.close();
    _eventDataOutputStream.close();

    if (_outputAcceptanceInfo) {
        _acceptanceInfoStream.close();
    }
}


void MCMC::run()
{
    log() << "\nRunning MCMC chain for "
          << _numGenerations << " generations.\n";

    log() << "\n";
    outputHeaders();

    for (_generation = 0; _generation < _numGenerations; _generation++) {
        _model->proposeNewState();

        outputAcceptanceInfo();
        outputData(_generation);
    }
}


void MCMC::outputHeaders()
{
    outputMCMCHeaders();
    outputEventDataHeaders();
    outputStdOutHeaders();
    outputAcceptanceInfoHeaders();
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


void MCMC::outputAcceptanceInfoHeaders()
{
    _acceptanceInfoStream << "proposedParam,accepted\n";
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
    _mcmcOutputStream << _generation                        << ","
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
    _model->getEventDataString(eventData, _generation);
    _eventDataOutputStream << eventData.str() << std::endl;
}


void MCMC::outputStdOutData()
{
    log() << std::setw(15) << _generation
          << std::setw(15) << _model->getCurrentLogLikelihood()
          << std::setw(15) << _model->getNumberOfEvents()
          << std::setw(15) << _model->computeLogPrior()
          << std::setw(15) << _model->getMHAcceptanceRate()
          << std::endl;
}


void MCMC::outputAcceptanceInfo()
{
    int param = _model->getLastParameterUpdated();
    int accepted = _model->getAcceptLastUpdate();
    _model->setAcceptLastUpdate(-1);

    _acceptanceInfoStream << param << "," << (accepted ? 1 : 0) << std::endl;
}
