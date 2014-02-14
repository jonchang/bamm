#include "MetropolisCoupledMCMC.h"

#include "MbRandom.h"
#include "Settings.h"
#include "MCMC.h"
#include "Model.h"

#include <vector>


MetropolisCoupledMCMC::MetropolisCoupledMCMC(MbRandom* rng, Settings* settings,
    Tree* tree, Prior* prior, int nChains, double deltaT, int swapPeriod)
        : _rng(rng), _settings(settings), _tree(tree), _prior(prior),
            _deltaT(deltaT), _swapPeriod(swapPeriod)
{
    _nGenerations = _settings->getNGENS();

    // Create chains (the 0th chain is assumed to be the cold chain)
    for (int i = 0; i < nChains; i++) {
        _chains.push_back(createMCMC(i));
    }
}


MCMC* MetropolisCoupledMCMC::createMCMC(int chainIndex) const
{
    Model* model = createModel(chainIndex);
    return createSpecificMCMC(model);
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


MetropolisCoupledMCMC::~MetropolisCoupledMCMC()
{
    for (int i = 0; i < (int)_chains.size(); i++) {
        // This class created both the model and the MCMC
        delete _chains[i]->getModel();
        delete _chains[i];
    }
}


void MetropolisCoupledMCMC::run()
{
    for (int gen = 0; gen < _nGenerations; gen += _swapPeriod) {
        // Run each chain up to _swapPeriod steps
        for (int i = 0; i < (int)_chains.size(); i++) {
            // TODO: Implement this run method in MCMC
            // _chains[i].run(_swapPeriod);
        }

        if (_chains.size() > 1) {
            int chainToSwap =
                _rng->discreteUniformRv(1, (int)_chains.size() - 1);
            swapTemperature(chainToSwap);
        }
    }
}


// Swaps the temperatures between the cold chain and the specified chain
void MetropolisCoupledMCMC::swapTemperature(int chainIndex)
{
    double coldBeta = _chains[0]->getModel()->getTemperatureMH();
    double hotBeta = _chains[chainIndex]->getModel()->getTemperatureMH();

    _chains[0]->getModel()->setTemperatureMH(hotBeta);
    _chains[chainIndex]->getModel()->setTemperatureMH(coldBeta);
}
