#ifndef TRAIT_MODEL_PARAMS_DATA_WRITER_H
#define TRAIT_MODEL_PARAMS_DATA_WRITER_H

#include <stdio.h>
#include <string>
#include <fstream>

class Settings;
class TraitModel;

class TraitModelParamsDataWriter
{
public:
    TraitModelParamsDataWriter(Settings& settings);
    ~TraitModelParamsDataWriter();
    
    void writeData(int generation, TraitModel& model);
    
private:

    void writeHeader();
    std::string header();
    bool _headerWritten;
    
    void initializeStream();
    std::string _outputFilename;
    int _outputFreq;
    
    std::ofstream _outputStream;

};




#endif /* defined(TRAIT_MODEL_PARAMS_DATA_WRITER_H) */
