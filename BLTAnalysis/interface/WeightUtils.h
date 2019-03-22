/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _WeightUtils_H
#define _WeightUtils_H

// c++ libraries
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// ROOT libraries
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

// Calibration libraries
#include "BTagCalibrationStandalone.h"
#include "CSVReader.hh"



struct WeightAndSigma{
		float nominal;
		float up;
		float down;
};



class WeightUtils {
    public:
        WeightUtils();
        //virtual ~WeightUtils() {};

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(std::string);
        void    SetSelection(std::string);

				//Pileup weights
        WeightAndSigma   GetPUWeight(float);
        WeightAndSigma   GetPUWeightMF(float);

				//Muon weights
        WeightAndSigma  GetMuonTriggerEff(std::string, TLorentzVector) const;
        WeightAndSigma  GetMuonRecoEff(TLorentzVector) const; 

				//Electron weights
        WeightAndSigma  GetElectronTriggerEff(TLorentzVector);
        WeightAndSigma  GetElectronRecoEff(TLorentzVector) const;
        WeightAndSigma  GetElectronTightRecoEff(TLorentzVector);


				WeightAndSigma GetBJetWeight(float pt, float eta, int flavor, float csv_score, float csv_cut);

				//Dilepnton weights
        WeightAndSigma  GetDiElectronTrig(TLorentzVector);
        WeightAndSigma  GetDiMuonTrig(TLorentzVector);
        WeightAndSigma  GetMuonElectronTrig(TLorentzVector, TLorentzVector);
				WeightAndSigma  GetDiLeptonTrig(TLorentzVector, TLorentzVector, std::string);

				WeightAndSigma  TEMP_GetDiLeptonTrig(TLorentzVector, TLorentzVector, std::string);


        //input parameters
        std::string _dataPeriod;
        std::string _sampleName;
        std::string _selection;
        bool				_isRealData;

        TGraph  _puReweight, _puReweight_up, _puReweight_down;
        TGraph  _puReweight_mf, _puReweight_up_mf, _puReweight_down_mf;
        
        //Muon Stuff
        //https://calderon.web.cern.ch/calderon/MuonPOG/2016dataReRecoEfficiencies/
        //muon reconstruction scalefactors: Measurement of inclusive W and Z boson cross section in pp collisions at \sqrt{s}=13 TeV
        //TGraphAsymmErrors *_eff_IsoMu22_DATA[4]; 
        TH2D _muSF_BCDEF_TRIG;
        TH2D _muSF_GH_TRIG;
        TGraphAsymmErrors _muSF_BCDEF_ID_DATA[4], _muSF_BCDEF_ID_MC[4]; 
        TGraphAsymmErrors _muSF_GH_ID_DATA[4],		_muSF_GH_ID_MC[4];
				TGraphAsymmErrors _muSF_TRACK; 

        TH2D _muSF_BCDEF_ISO;
        TH2D _muSF_GH_ISO;

	
				BTagCalibration       _btag_calibrator;
				BTagCalibrationReader _btagreader; 

				//https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
				TH2D _elSF_id_iso;
				TH2D _elSF_reco;
				TH2D _elSF_trigger; //https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2017Feb-A-Popov_TriggerSF_Run2016All_v1.root

				//ADD up and down
				
				CSV _trigger_dimuon_highleg;
				CSV _trigger_dimuon_lowleg;

				CSV _trigger_dielectron_highleg;
				CSV _trigger_dielectron_lowleg;

				CSV _trigger_muonelectron_highleg;
				CSV _trigger_muonelectron_lowleg;

				CSV _trigger_dilepton;

				CSV _tight_electron;
				
};

#endif

