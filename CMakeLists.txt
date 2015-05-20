CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(BAMM C CXX)

SET(BAMM_VERSION 2.3.0)
SET(BAMM_VERSION_DATE 2015-03-09)

# Specify executable and source files
SET(CPP_LIB
    src/AcceptanceDataWriter.cpp
    src/BetaInitProposal.cpp
    src/BetaShiftProposal.cpp
    src/BetaTimeModeProposal.cpp
    src/BranchEvent.cpp
    src/BranchHistory.cpp
    src/ChainSwapDataWriter.cpp
    src/CommandLineProcessor.cpp
    src/EventCountLog.cpp
    src/EventDataWriter.cpp
    src/EventNumberForBranchProposal.cpp
    src/EventNumberProposal.cpp
    src/EventParameterProposal.cpp
    src/EventRateProposal.cpp
    src/FastSimulatePrior.cpp
    src/LambdaInitProposal.cpp
    src/LambdaShiftProposal.cpp
    src/LambdaTimeModeProposal.cpp
    src/Log.cpp
    src/MCMC.cpp
    src/MCMCDataWriter.cpp
    src/MbRandom.cpp
    src/MetropolisCoupledMCMC.cpp
    src/Model.cpp
    src/ModelDataWriter.cpp
    src/MoveEventProposal.cpp
    src/MuInitProposal.cpp
    src/MuShiftProposal.cpp
    src/NewickTreeReader.cpp
    src/Node.cpp
    src/NodeStateDataWriter.cpp
    src/NodeStateProposal.cpp
    src/PreservationRateProposal.cpp
    src/Prior.cpp
    src/Proposal.cpp
    src/Random.cpp
    src/Settings.cpp
    src/SettingsParameter.cpp
    src/SpExBranchEvent.cpp
    src/SpExDataWriter.cpp
    src/SpExEventDataWriter.cpp
    src/SpExModel.cpp
    src/Stat.cpp
    src/StdOutDataWriter.cpp
    src/TimeModeProposal.cpp
    src/Tools.cpp
    src/TraitBranchEvent.cpp
    src/TraitDataWriter.cpp
    src/TraitEventDataWriter.cpp
    src/TraitModel.cpp
    src/Tree.cpp
)

ADD_LIBRARY(objlib OBJECT ${CPP_LIB})
set_property(TARGET ${objlib} PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(bamm_shared SHARED $<TARGET_OBJECTS:objlib>)
add_library(bamm_static STATIC $<TARGET_OBJECTS:objlib>)
ADD_EXECUTABLE(bamm src/main.cpp)
target_link_libraries(bamm bamm_static)

# Specify flags according to compiler
IF(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    SET(CMAKE_CXX_FLAGS "-g -Wall -Wextra -O3 -std=c++11 -stdlib=libc++")
    SET(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++")
ELSEIF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
    IF(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.6.3)
        MESSAGE(FATAL_ERROR "GCC version must be >= 4.6.3")
    ELSEIF(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.7)
        SET(CMAKE_CXX_FLAGS "-g -Wall -Wextra -O3 -std=c++0x")
    ELSE()
        SET(CMAKE_CXX_FLAGS "-g -Wall -Wextra -O3 -std=c++11")
    ENDIF()
    FIND_PACKAGE(Threads REQUIRED)
    IF(Threads_FOUND)
        TARGET_LINK_LIBRARIES (bamm ${CMAKE_THREAD_LIBS_INIT})
    ENDIF()
ELSEIF(${CMAKE_CXX_COMPILER_ID} MATCHES MSVC)
    SET(CMAKE_CXX_FLAGS "/W4")
ENDIF()

# Provide BAMM version to the compiler
ADD_DEFINITIONS(-DBAMM_VERSION=\"${BAMM_VERSION}\")
ADD_DEFINITIONS(-DBAMM_VERSION_DATE=\"${BAMM_VERSION_DATE}\")

# Set Git commit id
EXECUTE_PROCESS(COMMAND "git" "rev-parse" "HEAD"
    OUTPUT_VARIABLE GIT_COMMIT_ID
    OUTPUT_STRIP_TRAILING_WHITESPACE)
ADD_DEFINITIONS(-DGIT_COMMIT_ID=\"${GIT_COMMIT_ID}\")

INSTALL(TARGETS bamm RUNTIME DESTINATION bin)
