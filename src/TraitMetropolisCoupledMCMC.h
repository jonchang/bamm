#ifndef TRAIT_METROPOLIS_COUPLED_MCMC
#define TRAIT_METROPOLIS_COUPLED_MCMC


#include "MetropolisCoupledMCMC.h"
#include "TraitMCMC.h"
#include "TraitModel.h"

class MbRandom;
class Settings;
class MCMC;
class Model;


class TraitMetropolisCoupledMCMC : public MetropolisCoupledMCMC
{
public:

    TraitMetropolisCoupledMCMC(MbRandom* rng, Settings* settings,
        Prior* prior, int nChains, double deltaT, int swapPeriod);

protected:

    virtual MCMC* createSpecificMCMC(int chainIndex, Model* model) const;
    virtual Model* createSpecificModel() const;

    virtual void outputSpecificEventDataHeaders();
};


TraitMetropolisCoupledMCMC::TraitMetropolisCoupledMCMC(MbRandom* rng,
    Settings* settings, Prior* prior, int nChains, double deltaT,
        int swapPeriod) : MetropolisCoupledMCMC(rng, settings, prior,
            nChains, deltaT, swapPeriod)
{
}


inline MCMC* TraitMetropolisCoupledMCMC::createSpecificMCMC
    (int chainIndex, Model* model) const
{
    return new TraitMCMC(_rng, model, _settings, chainIndex + 1);
}


inline Model* TraitMetropolisCoupledMCMC::createSpecificModel() const
{
    return new TraitModel(_rng, _settings, _prior);
}


inline void TraitMetropolisCoupledMCMC::outputSpecificEventDataHeaders()
{
    _eventDataOutputStream << ",betainit,betashift\n";
}


#endif
