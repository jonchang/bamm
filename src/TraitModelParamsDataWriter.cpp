//
//  TraitModelParamsDataWriter.cpp
//  
//
//  Created by Dan Rabosky on 5/15/16.
//
//

#include "TraitModelParamsDataWriter.h"
#include "Settings.h"
#include "TraitModel.h"


TraitModelParamsDataWriter::TraitModelParamsDataWriter(Settings& settings) :
    _outputFilename(settings.get("traitModelFilename")),
    _outputFreq(settings.get<int>("mcmcWriteFreq")),
    _headerWritten(false)
{
    if (_outputFreq > 0){
        initializeStream();
    }

}

void TraitModelParamsDataWriter::initializeStream()
{
    _outputStream.open(_outputFilename.c_str());
}

TraitModelParamsDataWriter::~TraitModelParamsDataWriter()
{
    if (_outputFreq > 0){
        _outputStream.close();
    }
}


void TraitModelParamsDataWriter::writeData(int generation, TraitModel& model)
{
    if (_outputFreq == 0 || generation % _outputFreq != 0){
        return;
    }
    
    if (_headerWritten == false){
        writeHeader();
        _headerWritten = true;
    }
    
    _outputStream   <<    generation                          << ","
                    <<    model.getNumberOfRateShiftEvents()  << ","
                    <<    model.getNumberOfJumpEvents()       << ","
                    <<    model.getRootState()                << ","
                    <<    model.getJumpVariance()
                    <<    std::endl;
}


void TraitModelParamsDataWriter::writeHeader()
{
    _outputStream << header() << std::endl;

}


std::string TraitModelParamsDataWriter::header()
{
    return "generation,rateshifts,jumps,rootstate,jumpvariance";
}




