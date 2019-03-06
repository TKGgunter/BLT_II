// =============================================================================
// A simple analysis on Bacon ntuples
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => number of entries
//   argv[3] => ...
//
// Users should inherit from BLTSelector and implement the three functions:
//   Begin()
//   Process()
//   Terminate()
// =============================================================================


#ifndef DEMOANALYZER_HH
#define DEMOANALYZER_HH

// Analysis tools 
// TODO: Moved to BLT_II directory over
#include "BLT_II/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT_II/BLTAnalysis/interface/Parameters.hh"
#include "BLT_II/BLTAnalysis/interface/TriggerSelector.hh"



//CMSSW libraries
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TRandom.h>

// C++ headers
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>


class DemoAnalyzer: public BLTSelector {
public:
    DemoAnalyzer();
    ~DemoAnalyzer();

    void    Begin(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();

    void    ReportPostBegin();
    void    ReportPostTerminate();

    TFile       *outFile;
    TTree       *outTree;
    TTree       *debugTree;

    std::string  outFileName;
    std::string  outTreeName;


		//rochcor2016 *muonCorr;
		//RoccoR *muonCorr;
		RunLumiRangeMap lumiMask;
    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<baconhep::TTrigger>    triggerSelector;
    //std::unique_ptr<EnergyScaleCorrection_class> eGammaCorrections; 
};


#endif  // DEMOANALYZER_HH
