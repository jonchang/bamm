#ifndef MCMC_H
#define MCMC_H


#include "Random.h"

class Settings;
class Model;
class ModelFactory;


class MCMC
{
public:

    MCMC(Random& seeder, Settings& settings, ModelFactory& modelFactory);
    ~MCMC();

    void run(int generations);
    void step();

    Model& model();

protected:

    int _generation;

    // MCMC has its own random generator, using the seeder
    // (another random generator) to seed it
    Random _random;
    Model* _model;
};


inline Model& MCMC::model()
{
    return *_model;
}


#endif
