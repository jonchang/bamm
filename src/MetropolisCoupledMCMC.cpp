#include "MetropolisCoupledMCMC.h"

#include "MbRandom.h"
#include "Settings.h"
#include "MCMC.h"
#include "Model.h"

#include <vector>
#include <sstream>


MetropolisCoupledMCMC::MetropolisCoupledMCMC(MbRandom* rng, Settings* settings,
    Prior* prior, int nChains, double deltaT, int swapPeriod)
        : _rng(rng), _settings(settings), _prior(prior), _nChains(nChains),
          _deltaT(deltaT), _swapPeriod(swapPeriod), _coldChainIndex(0)
{
    // Total number of steps/generations to run chains
    _nGenerations = _settings->getNGENS();

    // MCMC-related output
    _mcmcOutputFileName = _settings->getMCMCoutfile();
    _mcmcOutputFreq = _settings->getMCMCwriteFreq();

    // Event-related output
    _eventDataOutputFileName = _settings->getEventDataOutfile();
    _eventDataOutputFreq = _settings->getEventDataWriteFreq();

    // Output to the screen (stdandard output)
    _stdOutFreq = _settings->getPrintFreq();

    // Open streams for writing
    _mcmcOutputStream.open(_mcmcOutputFileName.c_str());
    _eventDataOutputStream.open(_eventDataOutputFileName.c_str());
}


MetropolisCoupledMCMC::~MetropolisCoupledMCMC()
{
    _mcmcOutputStream.close();
    _eventDataOutputStream.close();

    for (int i = 0; i < (int)_chains.size(); i++) {
        // This class created both the model and the MCMC
        delete _chains[i]->getModel();
        delete _chains[i];
    }
}


MCMC* MetropolisCoupledMCMC::createMCMC(int chainIndex) const
{
    Model* model = createModel(chainIndex);
    return createSpecificMCMC(chainIndex, model);
}


Model* MetropolisCoupledMCMC::createModel(int chainIndex) const
{
    Model* model = createSpecificModel();
    double beta = calculateTemperature(chainIndex, _deltaT);
    model->setTemperatureMH(beta);
    return model;
}


double MetropolisCoupledMCMC::calculateTemperature
    (int i, double deltaT) const
{
    return 1.0 / (1.0 + deltaT * i);
}


void MetropolisCoupledMCMC::run()
{
    // Create chains (the 0th chain is the cold chain at first)
    for (int i = 0; i < _nChains; i++) {
        _chains.push_back(createMCMC(i));
    }

    log() << "\nRunning " << _chains.size() << " chains for "
          << _nGenerations << " generations.\n\n";

    outputHeaders();

    for (int gen = 0; gen < _nGenerations; gen += _swapPeriod) {
        // Run each chain up to _swapPeriod steps
        for (int i = 0; i < (int)_chains.size(); i++) {
            _chains[i]->run(_swapPeriod);
        }

        // If there is more than one chain, swap temperatures
        if (_chains.size() > 1) {
            int chain_1 = _rng->discreteUniformRv(0, (int)_chains.size() - 1);
            int chain_2 = _rng->discreteUniformRv(0, (int)_chains.size() - 1);

            // TODO: Can the random chains be the same? I'm assuming no.
            while (chain_1 == chain_2) {
                chain_2 = _rng->discreteUniformRv(0, (int)_chains.size() - 1);
            }

            if (acceptChainSwap(chain_1, chain_2)) {
                swapTemperature(chain_1, chain_2);
            }
        }

        // The current generation is gen + _swapPeriod
        outputData(gen + _swapPeriod);
    }
}


void MetropolisCoupledMCMC::outputHeaders()
{
    outputMCMCHeaders();
    outputEventDataHeaders();
    outputStdOutHeaders();
}


void MetropolisCoupledMCMC::outputMCMCHeaders()
{
    _mcmcOutputStream << "cold_chain,chain,temp,generation,N_shifts,"
        << "logPrior,logLik,eventRate,acceptRate\n";
}


void MetropolisCoupledMCMC::outputEventDataHeaders()
{
    _eventDataOutputStream << "generation,leftchild,rightchild,abstime";
    outputSpecificEventDataHeaders();
}


void MetropolisCoupledMCMC::outputStdOutHeaders()
{
    log() << "cold_chain,chain,temp,generation,N_shifts,logPrior,logLik,"
          << "eventRate,acceptRate\n";
}


bool MetropolisCoupledMCMC::acceptChainSwap(int chain_1, int chain_2) const
{
    return trueWithProbability(chainSwapProbability(chain_1, chain_2));
}


bool MetropolisCoupledMCMC::trueWithProbability(double p) const
{
    return _rng->uniformRv() < p;
}


double MetropolisCoupledMCMC::chainSwapProbability
    (int chain_1, int chain_2) const
{
    Model* model_1 = _chains[chain_1]->getModel();
    Model* model_2 = _chains[chain_2]->getModel();

    double beta_1 = model_1->getTemperatureMH();
    double beta_2 = model_2->getTemperatureMH();

    double log_post_1 = calculateLogPosterior(model_1);
    double log_post_2 = calculateLogPosterior(model_2);

    double swapPosteriorRatio = std::exp(logSwapPosteriorRatio
        (beta_1, beta_2, log_post_1, log_post_2));

    return std::min(1.0, swapPosteriorRatio);
}


double MetropolisCoupledMCMC::calculateLogPosterior(Model* model) const
{
    return model->getCurrentLogLikelihood() + model->computeLogPrior();
}


double MetropolisCoupledMCMC::logSwapPosteriorRatio
    (double beta_1, double beta_2, double log_post_1, double log_post_2) const
{
    return (beta_1 * log_post_2 + beta_2 * log_post_1) /
           (beta_1 * log_post_1 + beta_2 * log_post_2);
}


void MetropolisCoupledMCMC::swapTemperature(int chain_1, int chain_2)
{
    double beta_1 = _chains[chain_1]->getModel()->getTemperatureMH();
    double beta_2 = _chains[chain_2]->getModel()->getTemperatureMH();

    _chains[chain_1]->getModel()->setTemperatureMH(beta_2);
    _chains[chain_2]->getModel()->setTemperatureMH(beta_1);

    // Properly keep track of the cold chain
    if (chain_1 == _coldChainIndex) {
        _coldChainIndex = chain_2;
    } else if (chain_2 == _coldChainIndex) {
        _coldChainIndex = chain_1;
    }
}


void MetropolisCoupledMCMC::outputData(int generation)
{
    if (generation % _mcmcOutputFreq == 0) {
        outputMCMCData();
    }

    if (generation % _eventDataOutputFreq == 0) {
        outputEventData();
    }

    if (generation % _stdOutFreq == 0) {
        outputStdOutData();
    }
}


void MetropolisCoupledMCMC::outputMCMCData()
{
    for (int i = 0; i < (int)_chains.size(); i++) {
        Model* model = _chains[i]->getModel();
        _mcmcOutputStream << _coldChainIndex + 1 << ","
            << (i + 1)                           << ","
            << model->getTemperatureMH()         << ","
            << model->getGeneration()            << ","
            << model->getNumberOfEvents()        << ","
            << model->computeLogPrior()          << ","
            << model->getCurrentLogLikelihood()  << ","
            << model->getEventRate()             << ","
            << model->getMHAcceptanceRate()      << std::endl;
    }
}


void MetropolisCoupledMCMC::outputEventData()
{
    std::stringstream eventData;

    // Output only the cold chain's data
    _chains[_coldChainIndex]->getModel()->getEventDataString(eventData);
    _eventDataOutputStream << eventData.str() << std::endl;
}


void MetropolisCoupledMCMC::outputStdOutData()
{
    for (int i = 0; i < (int)_chains.size(); i++) {
        Model* model = _chains[i]->getModel();
        log() << _coldChainIndex + 1 << ","
            << (i + 1)                           << ","
            << model->getTemperatureMH()         << ","
            << model->getGeneration()            << ","
            << model->getNumberOfEvents()        << ","
            << model->computeLogPrior()          << ","
            << model->getCurrentLogLikelihood()  << ","
            << model->getEventRate()             << ","
            << model->getMHAcceptanceRate()      << std::endl;
    }
}
