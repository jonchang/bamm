#ifndef METROPOLIS_COUPLED_MCMC
#define METROPOLIS_COUPLED_MCMC


#include <vector>
#include <fstream>

class MbRandom;
class Settings;
class Tree;
class Prior;
class MCMC;
class Model;


class MetropolisCoupledMCMC
{
public:

    MetropolisCoupledMCMC(MbRandom* rng, Settings* settings, Prior* prior,
        int nChains, double deltaT, int swapPeriod);
    virtual ~MetropolisCoupledMCMC();

    void run();

protected:

    MCMC* createMCMC(int chainIndex) const;
    Model* createModel(int chainIndex) const;

    double calculateTemperature(int i, double deltaT) const;
    void swapTemperature(int chain_1, int chain_2);

    void outputHeaders();
    void outputMCMCHeaders();
    void outputEventDataHeaders();
    void outputStdOutHeaders();

    bool acceptChainSwap(int chain_1, int chain_2) const;
    bool trueWithProbability(double p) const;
    double chainSwapProbability(int chain_1, int chain_2) const;

    double calculateLogPosterior(Model* model) const;
    double logSwapPosteriorRatio(double beta_1, double beta_2,
        double log_post_1, double log_post_2) const;

    void outputData(int generation);
    void outputMCMCData();
    void outputEventData();
    void outputStdOutData();

    virtual MCMC* createSpecificMCMC(int chainIndex, Model* model) const = 0;
    virtual Model* createSpecificModel() const = 0;

    virtual void outputSpecificEventDataHeaders() = 0;

    MbRandom* _rng;
    Settings* _settings;
    Tree* _tree;
    Prior* _prior;

    int _nGenerations;

    int _nChains;

    // Holds a variable number of Markov chains
    std::vector<MCMC*> _chains;

    // From Altekar, et al. 2004: delta T (> 1) is a temparature
    // increment parameter chosen such that swaps are accepted
    // between 20 and 60% of the time.
    double _deltaT;

    // Number of steps/generations in which chair swapping occurs
    int _swapPeriod;

    // Current index of the cold chain (it changes when a swap occurs)
    int _coldChainIndex;

    std::string _mcmcOutputFileName;
    int _mcmcOutputFreq;

    std::string _eventDataOutputFileName;
    int _eventDataOutputFreq;

    int _stdOutFreq;

    std::ofstream _mcmcOutputStream;
    std::ofstream _eventDataOutputStream;
};


#endif
