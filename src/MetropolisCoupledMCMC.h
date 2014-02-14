#ifndef METROPOLIS_COUPLED_MCMC
#define METROPOLIS_COUPLED_MCMC


#include <vector>

class MbRandom;
class Settings;
class Tree;
class Prior;
class MCMC;
class Model;


class MetropolisCoupledMCMC
{
public:

    MetropolisCoupledMCMC(MbRandom* rng, Settings* settings, Tree* tree,
        Prior* prior, int nChains, double deltaT, int swapPeriod);
    virtual ~MetropolisCoupledMCMC();

    void run();

protected:

    MCMC* createMCMC(int chainIndex) const;
    Model* createModel(int chainIndex) const;

    double calculateTemperature(int i, double deltaT) const;
    void swapTemperature(int chainIndex);

    virtual MCMC* createSpecificMCMC(Model* model) const = 0;
    virtual Model* createSpecificModel() const = 0;

    MbRandom* _rng;
    Settings* _settings;
    Tree* _tree;
    Prior* _prior;

    int _nGenerations;

    // Holds a variable number of Markov chains
    std::vector<MCMC*> _chains;

    // From Altekar, et al. 2004: delta T (> 1) is a temparature
    // increment parameter chosen such that swaps are accepted
    // between 20 and 60% of the time.
    double _deltaT;

    // Number of steps/generations in which chair swapping occurs
    int _swapPeriod;
};


#endif
