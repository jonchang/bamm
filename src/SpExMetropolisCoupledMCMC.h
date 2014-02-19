#ifndef SP_EX_METROPOLIS_COUPLED_MCMC
#define SP_EX_METROPOLIS_COUPLED_MCMC


#include "MetropolisCoupledMCMC.h"
#include "SpExMCMC.h"
#include "SpExModel.h"

class MbRandom;
class Settings;
class MCMC;
class Model;


class SpExMetropolisCoupledMCMC : public MetropolisCoupledMCMC
{
public:

    SpExMetropolisCoupledMCMC(MbRandom* rng, Settings* settings,
        Prior* prior, int nChains, double deltaT, int swapPeriod);

protected:

    virtual MCMC* createSpecificMCMC(int chainIndex, Model* model) const;
    virtual Model* createSpecificModel() const;

    virtual void outputSpecificEventDataHeaders();
};


SpExMetropolisCoupledMCMC::SpExMetropolisCoupledMCMC(MbRandom* rng,
    Settings* settings, Prior* prior, int nChains, double deltaT,
        int swapPeriod) : MetropolisCoupledMCMC(rng, settings, prior,
            nChains, deltaT, swapPeriod)
{
}


inline MCMC* SpExMetropolisCoupledMCMC::createSpecificMCMC
    (int chainIndex, Model* model) const
{
    return new SpExMCMC(_rng, model, _settings, chainIndex + 1);
}


inline Model* SpExMetropolisCoupledMCMC::createSpecificModel() const
{
    return new SpExModel(_rng, _settings, _prior);
}


inline void SpExMetropolisCoupledMCMC::outputSpecificEventDataHeaders()
{
    _eventDataOutputStream << ",lambdainit,lambdashift,muinit,mushift\n";
}


#endif
