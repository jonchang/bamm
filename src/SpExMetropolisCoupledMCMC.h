#ifndef SP_EX_METROPOLIS_COUPLED_MCMC
#define SP_EX_METROPOLIS_COUPLED_MCMC


#include "MetropolisCoupledMCMC.h"
#include "SpExMCMC.h"
#include "SpExModel.h"

class MbRandom;
class Settings;
class Tree;
class MCMC;
class Model;


class SpExMetropolisCoupledMCMC : public MetropolisCoupledMCMC
{
public:

    SpExMetropolisCoupledMCMC(MbRandom* rng, Settings* settings, Tree* tree,
        Prior* prior, int nChains, double deltaT, int swapPeriod);

protected:

    virtual MCMC* createSpecificMCMC(Model* model) const;
    virtual Model* createSpecificModel() const;
};


SpExMetropolisCoupledMCMC::SpExMetropolisCoupledMCMC(MbRandom* rng,
    Settings* settings, Tree* tree, Prior* prior, int nChains, double deltaT,
        int swapPeriod) : MetropolisCoupledMCMC(rng, settings, tree, prior,
            nChains, deltaT, swapPeriod)
{
}


inline MCMC* SpExMetropolisCoupledMCMC::createSpecificMCMC(Model* model) const
{
    return new SpExMCMC(_rng, model, _settings);
}


inline Model* SpExMetropolisCoupledMCMC::createSpecificModel() const
{
    return new SpExModel(_rng, _tree, _settings, _prior);
}


#endif
