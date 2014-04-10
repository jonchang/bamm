#ifndef METROPOLIS_COUPLED_MCMC_H
#define METROPOLIS_COUPLED_MCMC_H


#include "StdOutDataWriter.h"
#include "MCMCDataWriter.h"
#include "AcceptanceDataWriter.h"
#include "ChainSwapDataWriter.h"
#include <vector>

class MbRandom;
class Settings;
class ModelFactory;
class MCMC;
class Model;
class EventDataWriter;


class MetropolisCoupledMCMC
{
public:

    MetropolisCoupledMCMC
        (MbRandom& rng, Settings& settings, ModelFactory* modelFactory);
    ~MetropolisCoupledMCMC();

    void run();

private:

    void createChains();
    MCMC* createMCMC(int chainIndex) const;
    double calculateTemperature(int i, double deltaT) const;

    void runChains(int genStart, int genEnd);
    void tryChainSwap(int generation);

    void chooseTwoNumbers(int* x, int* y, int from, int to);
    bool acceptChainSwap(int chain_1, int chain_2) const;
    bool trueWithProbability(double p) const;
    double chainSwapProbability(int chain_1, int chain_2) const;
    double calculateLogPosterior(Model& model) const;
    double logSwapPosteriorRatio(double beta_1, double beta_2,
        double log_post_1, double log_post_2) const;
    void swapTemperature(int chain_1, int chain_2);

    MbRandom& _rng;
    Settings& _settings;
    ModelFactory* _modelFactory;

    int _nGenerations;

    // Holds a variable number of Markov chains
    std::vector<MCMC*> _chains;
    int _nChains;

    // From Altekar, et al. 2004: delta T (> 1) is a temparature
    // increment parameter chosen such that swaps are accepted
    // between 20 and 60% of the time.
    double _deltaT;

    // Number of steps/generations in which chair swapping occurs
    int _swapPeriod;

    // Current index of the cold chain (it changes when a swap occurs)
    int _coldChainIndex;

    StdOutDataWriter _stdOutDataWriter;
    MCMCDataWriter _mcmcDataWriter;
    EventDataWriter* _eventDataWriter;
    AcceptanceDataWriter _acceptanceDataWriter;
    ChainSwapDataWriter _chainSwapDataWriter;

    int _acceptanceResetFreq;
};


#endif