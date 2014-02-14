#ifndef TRAIT_METROPOLIS_COUPLED_MCMC
#define TRAIT_METROPOLIS_COUPLED_MCMC


#include "MetropolisCoupledMCMC.h"
#include "TraitMCMC.h"
#include "TraitModel.h"

class MbRandom;
class Settings;
class Tree;
class MCMC;
class Model;


class TraitMetropolisCoupledMCMC : public MetropolisCoupledMCMC
{
public:

    TraitMetropolisCoupledMCMC(MbRandom* rng, Settings* settings, Tree* tree,
        Prior* prior, int nChains, double deltaT, int swapPeriod);

protected:

    virtual MCMC* createSpecificMCMC(Model* model) const;
    virtual Model* createSpecificModel() const;
};


TraitMetropolisCoupledMCMC::TraitMetropolisCoupledMCMC(MbRandom* rng,
    Settings* settings, Tree* tree, Prior* prior, int nChains, double deltaT,
        int swapPeriod) : MetropolisCoupledMCMC(rng, settings, tree, prior,
            nChains, deltaT, swapPeriod)
{
}


inline MCMC* TraitMetropolisCoupledMCMC::createSpecificMCMC(Model* model) const
{
    return new TraitMCMC(_rng, model, _settings);
}


inline Model* TraitMetropolisCoupledMCMC::createSpecificModel() const
{
    return new TraitModel(_rng, _tree, _settings, _prior);
}


#endif
