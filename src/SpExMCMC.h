#ifndef SP_EX_MCMC_H
#define SP_EX_MCMC_H

#include <stdlib.h>
#include <string>
#include <vector>
#include <iosfwd>

class MbRandom;
class Model;
class Settings;


class SpExMCMC
{

public:

    SpExMCMC(MbRandom* ran, Model* mymodel, Settings* sp);
    ~SpExMCMC();

    void writeStateToFile();
    void printStateData();
    void writeBranchSpeciationRatesToFile();
    void writeBranchExtinctionRatesToFile();
    void writeEventDataToFile();

    int  pickParameterClassToUpdate();
    void updateState(int parm);

    void setUpdateWeights();

private:

    void writeHeaderToStream(std::ostream& outStream);
    void writeStateToStream(std::ostream& outStream);

    bool anyOutputFileExists();
    bool fileExists(const std::string& filename);
    void writeHeadersToOutputFiles();
    void exitWithErrorOutputFileExists();

    MbRandom* ranPtr;
    Model*    ModelPtr;
    Settings* sttings;

    std::vector<double> parWts;

    std::vector<int> acceptCount;
    std::vector<int> rejectCount;

    std::string _mcmcOutFilename;
    std::string _lambdaOutFilename;
    std::string _muOutFilename;
    std::string _eventDataOutFilename;

    std::ofstream _mcmcOutStream;
    std::ofstream _lambdaOutStream;
    std::ofstream _muOutStream;
    std::ofstream _eventDataOutStream;

    bool _writeMeanBranchLengthTrees;

    int _treeWriteFreq;
    int _eventDataWriteFreq;
    int _mcmcWriteFreq;
    int _acceptWriteFreq;
    int _printFreq;
    int _NGENS;
};


#endif
