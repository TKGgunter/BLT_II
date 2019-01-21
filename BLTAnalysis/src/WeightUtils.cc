#include "BLT_II/BLTAnalysis/interface/WeightUtils.h"
#include "BLT_II/BLTAnalysis/interface/tgUtility.hh"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>


WeightUtils::WeightUtils() //std::string dataPeriod, std::string selection, bool isRealData)
{

		srand ( 124 );

    const std::string cmssw_base = getenv("CMSSW_BASE");

    // PU METFILTER Weights  
		{
				std::string puFileName_mf = cmssw_base + "/src/BLT_II/BLTAnalysis/data/pileup_sf_2016_full_69216_bins75_metFilter_mar2018.root"; 
				std::string puFileName_up_mf = cmssw_base + "/src/BLT_II/BLTAnalysis/data/pileup_sf_2016_full_72386_bins75_metFilter_mar2018.root"; 
				std::string puFileName_down_mf = cmssw_base + "/src/BLT_II/BLTAnalysis/data/pileup_sf_2016_full_66013_bins75_metFilter_mar2018.root"; 
		

				OPEN_TFILE(puFileName_mf, puFile_mf);
				_puReweight_mf = *(TGraph*)puFile_mf->Get("pileup_sf");

				OPEN_TFILE(puFileName_up_mf, puFile_up_mf);
				_puReweight_up_mf = *(TGraph*)puFile_up_mf->Get("pileup_sf");

				OPEN_TFILE(puFileName_down_mf, puFile_down_mf);
				_puReweight_down_mf = *(TGraph*)puFile_down_mf->Get("pileup_sf");
		}


		// PU Weights
		{
				std::string puFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/pileup_sf_2016_full_69216_bins75_mar2018.root"; 
				std::string puFileName_up = cmssw_base + "/src/BLT_II/BLTAnalysis/data/pileup_sf_2016_full_72386_bins75_mar2018.root"; 
				std::string puFileName_down = cmssw_base + "/src/BLT_II/BLTAnalysis/data/pileup_sf_2016_full_66013_bins75_mar2018.root"; 

				OPEN_TFILE(puFileName, puFile);
				_puReweight = *(TGraph*)puFile->Get("pileup_sf");

				OPEN_TFILE(puFileName_up, puFile_up);
				_puReweight_up = *(TGraph*)puFile_up->Get("pileup_sf");

				OPEN_TFILE(puFileName_down, puFile_down);
				_puReweight_down = *(TGraph*)puFile_down->Get("pileup_sf");
		}
		///////////// 
		//muon Trigger
		{
				std::string recoFileName_ = cmssw_base + "/src/BLT_II/BLTAnalysis/data/MuonRecoEfficienciesAndSF_RunBtoF.root";
				OPEN_TFILE(recoFileName_, f_muSF_reco);
				_muSF_BCDEF_TRIG    = *(TH2D*) f_muSF_reco->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
				if (_muSF_BCDEF_TRIG.IsZombie()){
						std::cout << "muon trigger shit is nulled" << std::endl;
				}
		}
		{
				std::string recoFileName_ = cmssw_base + "/src/BLT_II/BLTAnalysis/data/MuonRecoEfficienciesAndSF_Period4.root";
				OPEN_TFILE(recoFileName_, f_muSF_reco);
				_muSF_GH_TRIG    = *(TH2D*) f_muSF_reco->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
		}
		///////////// 

    // tight muon ID sf
		{
				std::string idFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/MuonID_EfficienciesAndSF_BCDEF.root";
				OPEN_TFILE(idFileName, f_muRecoSF_ID );
				
				std::string graphPath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
				_muSF_BCDEF_ID_DATA[0] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin0_DATA").c_str());
				_muSF_BCDEF_ID_DATA[1] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin1_DATA").c_str());
				_muSF_BCDEF_ID_DATA[2] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin2_DATA").c_str());
				_muSF_BCDEF_ID_DATA[3] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin3_DATA").c_str());
				
				graphPath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
				_muSF_BCDEF_ID_MC[0] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin0_MC").c_str());
				_muSF_BCDEF_ID_MC[1] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin1_MC").c_str());
				_muSF_BCDEF_ID_MC[2] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin2_MC").c_str());
				_muSF_BCDEF_ID_MC[3] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin3_MC").c_str());
		}


		// Runs GH
		//Muon ID
		{
				std::string idFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/MuonID_EfficienciesAndSF_GH.root";
				OPEN_TFILE(idFileName, f_muRecoSF_ID );
				
				std::string graphPath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
				_muSF_GH_ID_DATA[0] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin0_DATA").c_str());
				_muSF_GH_ID_DATA[1] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin1_DATA").c_str());
				_muSF_GH_ID_DATA[2] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin2_DATA").c_str());
				_muSF_GH_ID_DATA[3] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin3_DATA").c_str());
				
				graphPath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
				_muSF_GH_ID_MC[0] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin0_MC").c_str());
				_muSF_GH_ID_MC[1] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin1_MC").c_str());
				_muSF_GH_ID_MC[2] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin2_MC").c_str());
				_muSF_GH_ID_MC[3] = *(TGraphAsymmErrors*)f_muRecoSF_ID->Get((graphPath + "pt_PLOT_abseta_bin3_MC").c_str());
		}


		///////////// 
		//muon ISO 
		{
				std::string isoFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/IsoEfficienciesAndSF_BCDEF.root";

				OPEN_TFILE(isoFileName, f_muSF_ISO );
				_muSF_BCDEF_ISO = *(TH2D*) f_muSF_ISO->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
		}
		{
				std::string isoFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/IsoEfficienciesAndSF_GH.root";

				OPEN_TFILE(isoFileName, f_muSF_ISO );
				_muSF_GH_ISO    = *(TH2D*) f_muSF_ISO->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
		}

		///////////// 
		//muon Tracking
		{
				std::string trackingFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/MuonTracking_EfficienciesAndSF_BCDEFGH.root";

				OPEN_TFILE(trackingFileName, f_muSF_TRACK );
				_muSF_TRACK = *(TGraphAsymmErrors*)f_muSF_TRACK->Get("ratio_eff_aeta_dr030e030_corr");
		}
	

		//Electron ID
		{
				std::string el_idFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/egammaEffi_id_iso.root";

				OPEN_TFILE(el_idFileName, f_elSF_ID );
				_elSF_id_iso = *(TH2D*)f_elSF_ID->Get("EGamma_SF2D");
		}
		//Electron Reco
		{
				std::string el_idFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/egammaEffi_reco.root";

				OPEN_TFILE(el_idFileName, f_elSF_reco );
				_elSF_reco = *(TH2D*)f_elSF_reco->Get("EGamma_SF2D");
		}
		//Electron Trigger
		//https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2 link to root file at bottom of page.
		{
				std::string el_idFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/2017Feb-A-Popov_TriggerSF_Run2016All_v1.root";

				OPEN_TFILE(el_idFileName, f_elSF_trigger );
				_elSF_trigger = *(TH2D*)f_elSF_trigger->Get("Ele27_WPTight_Gsf");
		}



		//BTagging stuff
		std::ifstream file( cmssw_base + "/src/BLT_II/BLTAnalysis/data/CSVv2_Moriond17_B_H.csv");

		if ( file )
		{
				std::stringstream buffer;
				buffer << file.rdbuf();

				file.close();	

				_btag_calibrator = BTagCalibration("csv");
				_btag_calibrator.readCSV(buffer);

				//_btagreader = BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
				_btagreader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
				_btagreader.load(_btag_calibrator, BTagEntry::FLAV_B,		"comb");
				_btagreader.load(_btag_calibrator, BTagEntry::FLAV_C,		"comb");
				_btagreader.load(_btag_calibrator, BTagEntry::FLAV_UDSG,"incl");
		}

		if (false) {
				
				read_csv(cmssw_base + "/src/BLT_II/BLTAnalysis/data/", &_trigger_dimuon_lowleg);
				read_csv(cmssw_base + "/src/BLT_II/BLTAnalysis/data/", &_trigger_dimuon_highleg);

				read_csv(cmssw_base + "/src/BLT_II/BLTAnalysis/data/", &_trigger_dielectron_lowleg);
				read_csv(cmssw_base + "/src/BLT_II/BLTAnalysis/data/", &_trigger_dielectron_highleg);

				read_csv(cmssw_base + "/src/BLT_II/BLTAnalysis/data/", &_trigger_muonelectron_lowleg);
				read_csv(cmssw_base + "/src/BLT_II/BLTAnalysis/data/", &_trigger_muonelectron_highleg);

				read_csv(cmssw_base + "/src/BLT_II/BLTAnalysis/data/ditriggers", &_trigger_dilepton);
			
		}
}


void WeightUtils::SetDataBit(bool isRealData)
{
    _isRealData = isRealData;
}

void WeightUtils::SetDataPeriod(std::string dataPeriod)
{
    _dataPeriod = dataPeriod;
}

void WeightUtils::SetSelection(std::string selection)
{
    _selection = selection;
}

WeightAndSigma WeightUtils::GetPUWeight(float nPU)
{
    return {(float)_puReweight.Eval(nPU), (float)_puReweight_up.Eval(nPU), (float)_puReweight_down.Eval(nPU)}; 
}

WeightAndSigma WeightUtils::GetPUWeightMF(float nPU)
{
    return {(float)_puReweight_mf.Eval(nPU), (float)_puReweight_up_mf.Eval(nPU), (float)_puReweight_down_mf.Eval(nPU)}; 
}

WeightAndSigma WeightUtils::GetMuonTriggerEff(std::string triggerName, TLorentzVector muon) const
{
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    float binningPt[] = {25., 30., 40, 50, 60, 120, 200};
    int ptBin = 0;
    for (int i = 0; i < 6; ++i) {
        if (fabs(muon.Pt()) > binningPt[i] && fabs(muon.Pt()) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }
		if (muon.Pt() > 200) ptBin = 6;
    
    float weight   = 1;
    float error = 0;
    if (triggerName == "HLT_IsoMu24_v*") {
		if (float(rand() % 100) / 100. > (1.26E08 + 5.69E07 + 8.31E07 + 7.67E07 + 5.59E07) /
																		 (1.26E08 + 5.69E07 + 8.31E07 + 7.67E07 + 5.59E07 + 1.32E08 + 1.53E08)) 
		{ //Replace with actual ratio:  BF/ (BF + GH) 
				weight *= _muSF_BCDEF_TRIG.GetBinContent(etaBin + 1, ptBin + 1);
				error		= _muSF_BCDEF_TRIG.GetBinError(etaBin + 1, ptBin + 1) + weight * 0.005;
			}
			else {
				weight *= _muSF_GH_TRIG.GetBinContent(etaBin + 1, ptBin + 1);
				error		= _muSF_GH_TRIG.GetBinError(etaBin + 1, ptBin + 1) + weight * 0.005;
			}
    }

    return {weight, error, error};
}

WeightAndSigma WeightUtils::GetMuonRecoEff(TLorentzVector muon) const
{

		//Reco scale factors
    float weight_reco = 1;
		float error_reco = 0;
		{
			float binningEta[] =  {-2.4,		-2.1,		-1.2,		-0.90,	-0.3,		-0.2,		0.,			0.2,		0.3,		0.9,		1.2,		2.1,		2.4};
			float pt_25_40[]=			{0.9845, 0.9929, 0.9883, 0.9914, 0.9171, 0.9947, 0.9958, 0.9156, 0.9905, 0.9857, 0.9960, 0.9899};
			float pt_40[]		=			{0.9867, 0.9965, 0.9912, 0.9904, 0.9172, 0.9951, 0.9981, 0.9151, 0.9899, 0.9863, 0.9972, 0.9883};
			float err_pt_25_40[]= {0.0018, 0.0011, 0.0011, 0.0008, 0.0022, 0.0004, 0.0012, 0.0026, 0.0008, 0.0014, 0.0009, 0.0015};
			float err_pt_40[]		= {0.0015, 0.0004, 0.0000, 0.0000, 0.0023, 0.0007, 0.0000, 0.0023, 0.0005, 0.0008, 0.0004, 0.0010};

			//eta
			for (int i = 0; i < 12; i++){
        if (muon.Eta() > binningEta[i] && muon.Eta() <= binningEta[i+1]) {
					
					if (muon.Pt() < 40){
						weight_reco *= pt_25_40[i];
						error_reco = err_pt_25_40[i];
					}
					else{
						weight_reco *= pt_40[i];
						error_reco = err_pt_40[i];
					}
				} 
			}
		}


    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    float binningPt[] = {20., 25, 30, 40, 50, 60, 200};
    int ptBin = 0;
    for (int i = 0; i < 6; ++i) {
        if (fabs(muon.Pt()) > binningPt[i] && fabs(muon.Pt()) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }
		if (muon.Pt() > 200) ptBin = 5;


    float weight = 1;
		float high = 0;
		float low = 0;

		auto error = [](float a, float b, float a_, float b_){
			return pow( pow(1./b * a_, 2) + pow( a/b * b_, 2) , .5);
		};
		auto quad_add = [](float a, float b){
			return float(pow( pow(a, 2.0) + pow(b, 2.0), 0.5 ));
		};
		auto square = [](float a){
			return pow(a, 2.0);
		};

		if (float(rand() % 100) / 100. > (1.26E08 + 5.69E07 + 8.31E07 + 7.67E07 + 5.59E07) /
																		 (1.26E08 + 5.69E07 + 8.31E07 + 7.67E07 + 5.59E07 + 1.32E08 + 1.53E08)) 
		{ //Replace with actual ratio:  BF/ (BF + GH) 
			float w_data = _muSF_BCDEF_ID_DATA[etaBin].Eval(muon.Pt());
			float w_mc   = _muSF_BCDEF_ID_MC[etaBin].Eval(muon.Pt()) ; 
			weight   *= w_data/w_mc ;

			high = quad_add(error(w_data, w_mc, _muSF_BCDEF_ID_DATA[etaBin].GetErrorYhigh(ptBin), _muSF_BCDEF_ID_MC[etaBin].GetErrorYhigh(ptBin)), 0.01 * w_data/w_mc);
			low  = quad_add(error(w_data, w_mc, _muSF_BCDEF_ID_DATA[etaBin].GetErrorYlow(ptBin),  _muSF_BCDEF_ID_MC[etaBin].GetErrorYlow(ptBin)), 0.01 * w_data/w_mc);
			//ISO
			float _weight = _muSF_BCDEF_ISO.GetBinContent(etaBin + 1, ptBin + 1);
			weight   *= _weight;
			high      = pow(square(high) + square(quad_add(_muSF_BCDEF_ISO.GetBinError(etaBin + 1, ptBin + 1), 0.005 * _weight)), .5);
			low       = pow(square(low)  + square(quad_add(_muSF_BCDEF_ISO.GetBinError(etaBin + 1, ptBin + 1), 0.005 * _weight)), .5);
		}
		else{
			float w_data = _muSF_GH_ID_DATA[etaBin].Eval(muon.Pt());
			float w_mc   = _muSF_GH_ID_MC[etaBin].Eval(muon.Pt());
			weight	*=	w_data/w_mc;

			high = quad_add(error(w_data, w_mc, _muSF_GH_ID_DATA[etaBin].GetErrorYhigh(ptBin), _muSF_GH_ID_MC[etaBin].GetErrorYhigh(ptBin)), 0.01 * w_data/w_mc);
			low  = quad_add(error(w_data, w_mc, _muSF_GH_ID_DATA[etaBin].GetErrorYlow(ptBin),  _muSF_GH_ID_MC[etaBin].GetErrorYlow(ptBin)), 0.01 * w_data/w_mc);
			//ISO
			float _weight = _muSF_GH_ISO.GetBinContent(etaBin + 1, ptBin + 1);
			weight   *= _weight;
			high      = pow(square(high) + square(quad_add(_muSF_GH_ISO.GetBinError(etaBin + 1, ptBin + 1), 0.005 * _weight)), .5);
			low       = pow(square(low)  + square(quad_add(_muSF_GH_ISO.GetBinError(etaBin + 1, ptBin + 1), 0.005 * _weight)), .5);
    }

		float weight_track = _muSF_TRACK.Eval(abs(muon.Eta()));
		//std::cout << " muon weight: " << weight << " " << high << " " << low << "ISO " << _muSF_BCDEF_ISO->GetBinContent(etaBin + 1, ptBin + 1) << std::endl;
    return {weight * weight_reco * weight_track, quad_add(high, error_reco) , quad_add(low, error_reco)};
}

float th2_interpolate(const TH2* hist, Double_t x, Double_t y);
WeightAndSigma WeightUtils::GetElectronTriggerEff(TLorentzVector electron1)
{
		//Electron Trigger SF
		//`/tthome/ksung/cms/Analysis/12/CMSSW_8_0_20/src/DMSAna/ttDM/limits/LeptonEffUtils.hh`
		float eta = (float)electron1.Eta();
		float pt =  (float)electron1.Pt();
		

		float trigger_sf = 1;

		
		/*{ //Originally from Kevin
				float numerator = 1;
				float _denominator = 1;
				//Data Trigger
				{
					if(eta > -2.1 && eta < -1.6) {
						if(pt > 30 && pt < 32) { numerator *= 0.695; }
						else if(pt < 35)       { numerator *= 0.753; }
						else if(pt < 40)       { numerator *= 0.801; }
						else if(pt < 50)       { numerator *= 0.839; }
						else if(pt < 60)       { numerator *= 0.857; }
						else if(pt < 120)      { numerator *= 0.869; }
						else if(pt > 120)      { numerator *= 0.884; }

					} else if(eta < -1.4) {
						if(pt > 30 && pt < 32) { numerator *= 0.721; }
						else if(pt < 35)       { numerator *= 0.774; }
						else if(pt < 40)       { numerator *= 0.824; }
						else if(pt < 50)       { numerator *= 0.865; }
						else if(pt < 60)       { numerator *= 0.884; }
						else if(pt < 120)      { numerator *= 0.892; }
						else if(pt > 120)      { numerator *= 0.912; }

					} else if(eta < -0.8) {
						if(pt > 30 && pt < 32) { numerator *= 0.814; }
						else if(pt < 35)       { numerator *= 0.861; }
						else if(pt < 40)       { numerator *= 0.895; }
						else if(pt < 50)       { numerator *= 0.924; }
						else if(pt < 60)       { numerator *= 0.940; }
						else if(pt < 120)      { numerator *= 0.954; }
						else if(pt > 120)      { numerator *= 0.966; }

					} else if(eta < -0.4) {
						if(pt > 30 && pt < 32) { numerator *= 0.791; }
						else if(pt < 35)       { numerator *= 0.853; }
						else if(pt < 40)       { numerator *= 0.886; }
						else if(pt < 50)       { numerator *= 0.922; }
						else if(pt < 60)       { numerator *= 0.940; }
						else if(pt < 120)      { numerator *= 0.956; }
						else if(pt > 120)      { numerator *= 0.967; }

					} else if(eta < 0) {
						if(pt > 30 && pt < 32) { numerator *= 0.710; }
						else if(pt < 35)       { numerator *= 0.784; }
						else if(pt < 40)       { numerator *= 0.830; }
						else if(pt < 50)       { numerator *= 0.883; }
						else if(pt < 60)       { numerator *= 0.911; }
						else if(pt < 120)      { numerator *= 0.938; }
						else if(pt > 120)      { numerator *= 0.964; }

					} else if(eta < 0.4) {
						if(pt > 30 && pt < 32) { numerator *= 0.713; }
						else if(pt < 35)       { numerator *= 0.782; }
						else if(pt < 40)       { numerator *= 0.828; }
						else if(pt < 50)       { numerator *= 0.877; }
						else if(pt < 60)       { numerator *= 0.907; }
						else if(pt < 120)      { numerator *= 0.929; }
						else if(pt > 120)      { numerator *= 0.958; }

					} else if(eta < 0.8) {
						if(pt > 30 && pt < 32) { numerator *= 0.801; }
						else if(pt < 35)       { numerator *= 0.853; }
						else if(pt < 40)       { numerator *= 0.885; }
						else if(pt < 50)       { numerator *= 0.918; }
						else if(pt < 60)       { numerator *= 0.936; }
						else if(pt < 120)      { numerator *= 0.950; }
						else if(pt > 120)      { numerator *= 0.965; }

					} else if(eta < 1.4) {
						if(pt > 30 && pt < 32) { numerator *= 0.810; }
						else if(pt < 35)       { numerator *= 0.863; }
						else if(pt < 40)       { numerator *= 0.898; }
						else if(pt < 50)       { numerator *= 0.927; }
						else if(pt < 60)       { numerator *= 0.944; }
						else if(pt < 120)      { numerator *= 0.956; }
						else if(pt > 120)      { numerator *= 0.971; }

					} else if(eta < 1.6) {
						if(pt > 30 && pt < 32) { numerator *= 0.698; }
						else if(pt < 35)       { numerator *= 0.766; }
						else if(pt < 40)       { numerator *= 0.815; }
						else if(pt < 50)       { numerator *= 0.864; }
						else if(pt < 60)       { numerator *= 0.887; }
						else if(pt < 120)      { numerator *= 0.888; }
						else if(pt > 120)      { numerator *= 0.900; }

					} else if(eta < 2.1) {
						if(pt > 30 && pt < 32) { numerator *= 0.687; }
						else if(pt < 35)       { numerator *= 0.750; }
						else if(pt < 40)       { numerator *= 0.801; }
						else if(pt < 50)       { numerator *= 0.842; }
						else if(pt < 60)       { numerator *= 0.860; }
						else if(pt < 120)      { numerator *= 0.878; }
						else if(pt > 120)      { numerator *= 0.891; }
					}
				}
				trigger_sf *= numerator;
				//MC Trigger
				{
					if(eta > -2.1 && eta < -1.6) {
						if(pt > 30 && pt < 32) { _denominator *=0.801; }
						else if(pt < 35)       { _denominator *=0.838; }
						else if(pt < 40)       { _denominator *=0.879; }
						else if(pt < 50)       { _denominator *=0.908; }
						else if(pt < 60)       { _denominator *=0.923; }
						else if(pt < 120)      { _denominator *=0.937; }
						else if(pt > 120)      { _denominator *=0.966; }

					} else if(eta < -1.4) {
						if(pt > 30 && pt < 32) { _denominator *=0.795; }
						else if(pt < 35)       { _denominator *=0.836; }
						else if(pt < 40)       { _denominator *=0.869; }
						else if(pt < 50)       { _denominator *=0.904; }
						else if(pt < 60)       { _denominator *=0.922; }
						else if(pt < 120)      { _denominator *=0.932; }
						else if(pt > 120)      { _denominator *=0.957; }

					} else if(eta < -0.8) {
						if(pt > 30 && pt < 32) { _denominator *=0.848; }
						else if(pt < 35)       { _denominator *=0.888; }
						else if(pt < 40)       { _denominator *=0.904; }
						else if(pt < 50)       { _denominator *=0.929; }
						else if(pt < 60)       { _denominator *=0.949; }
						else if(pt < 120)      { _denominator *=0.970; }
						else if(pt > 120)      { _denominator *=0.983; }

					} else if(eta < -0.4) {
						if(pt > 30 && pt < 32) { _denominator *=0.829; }
						else if(pt < 35)       { _denominator *=0.866; }
						else if(pt < 40)       { _denominator *=0.883; }
						else if(pt < 50)       { _denominator *=0.912; }
						else if(pt < 60)       { _denominator *=0.940; }
						else if(pt < 120)      { _denominator *=0.963; }
						else if(pt > 120)      { _denominator *=0.979; }

					} else if(eta < 0) {
						if(pt > 30 && pt < 32) { _denominator *=0.791; }
						else if(pt < 35)       { _denominator *=0.824; }
						else if(pt < 40)       { _denominator *=0.848; }
						else if(pt < 50)       { _denominator *=0.883; }
						else if(pt < 60)       { _denominator *=0.917; }
						else if(pt < 120)      { _denominator *=0.949; }
						else if(pt > 120)      { _denominator *=0.983; }

					} else if(eta < 0.4) {
						if(pt > 30 && pt < 32) { _denominator *=0.796; }
						else if(pt < 35)       { _denominator *=0.819; }
						else if(pt < 40)       { _denominator *=0.847; }
						else if(pt < 50)       { _denominator *=0.884; }
						else if(pt < 60)       { _denominator *=0.920; }
						else if(pt < 120)      { _denominator *=0.947; }
						else if(pt > 120)      { _denominator *=0.979; }

					} else if(eta < 0.8) {
						if(pt > 30 && pt < 32) { _denominator *=0.832; }
						else if(pt < 35)       { _denominator *=0.863; }
						else if(pt < 40)       { _denominator *=0.882; }
						else if(pt < 50)       { _denominator *=0.913; }
						else if(pt < 60)       { _denominator *=0.943; }
						else if(pt < 120)      { _denominator *=0.966; }
						else if(pt > 120)      { _denominator *=0.986; }

					} else if(eta < 1.4) {
						if(pt > 30 && pt < 32) { _denominator *=0.856; }
						else if(pt < 35)       { _denominator *=0.887; }
						else if(pt < 40)       { _denominator *=0.907; }
						else if(pt < 50)       { _denominator *=0.930; }
						else if(pt < 60)       { _denominator *=0.954; }
						else if(pt < 120)      { _denominator *=0.969; }
						else if(pt > 120)      { _denominator *=0.988; }

					} else if(eta < 1.6) {
						if(pt > 30 && pt < 32) { _denominator *=0.798; }
						else if(pt < 35)       { _denominator *=0.836; }
						else if(pt < 40)       { _denominator *=0.865; }
						else if(pt < 50)       { _denominator *=0.903; }
						else if(pt < 60)       { _denominator *=0.923; }
						else if(pt < 120)      { _denominator *=0.934; }
						else if(pt > 120)      { _denominator *=0.945; }

					} else if(eta < 2.1) {
						if(pt > 30 && pt < 32) { _denominator *=0.808; }
						else if(pt < 35)       { _denominator *=0.852; }
						else if(pt < 40)       { _denominator *=0.889; }
						else if(pt < 50)       { _denominator *=0.913; }
						else if(pt < 60)       { _denominator *=0.930; }
						else if(pt < 120)      { _denominator *=0.939; }
						else if(pt > 120)      { _denominator *=0.971; }
					}
				}
				denominator *= _denominator;
		}
		trigger_sf = trigger_sf  / denominator;
		*/
		trigger_sf = th2_interpolate(&_elSF_trigger, pt, eta); //.Interpolate(pt, eta);
		if( trigger_sf == 0 ) trigger_sf = 1;
		return {trigger_sf, 0.0, 0.0};
}

WeightAndSigma WeightUtils::GetElectronRecoEff(TLorentzVector electron) const
{

		//Electron RECO SF
		//http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/Moriond17/approval/RECO/passingRECO/
		float reco_sf = 1;
		float reco_error = 1;
		{
				float binningEta[] =  {-2.5, -2.45, -2.4, -2.3, -2.2, -2.0, -1.8, -1.630, -1.566, -1.444, -1.2, -1.0, -0.6, -0.4, -0.2,
														0., 0.2, 0.4, 0.6, 1.0, 1.2, 1.444, 1.566, 1.630, 1.8, 2., 2.2, 2.3, 2.4, 2.45, 2.5};

				float etaBin = 0;
				for (int i = 0; i < 31; i++)
				{
						if (fabs(electron.Eta()) > binningEta[i] && fabs(electron.Eta()) <= binningEta[i+1]) {
								etaBin = i+1;
								break;
						}
							
				}
				reco_sf   *= _elSF_reco.GetBinContent(etaBin, 1);
				reco_error = _elSF_reco.GetBinError(etaBin, 1);
		}
/*
		{
		//`/tthome/ksung/cms/Analysis/12/CMSSW_8_0_20/src/DMSAna/ttDM/limits/LeptonEffUtils.hh`
		//
				float eta = electron.Eta();
				if     (eta <-2.5   ) { reco_sf *= 1; }
				else if(eta <-2.45  ) { reco_sf *= 1.3176; }
				else if(eta <-2.4   ) { reco_sf *= 1.11378; }
				else if(eta <-2.3   ) { reco_sf *= 1.02463; }
				else if(eta <-2.2   ) { reco_sf *= 1.01364; }
				else if(eta <-2     ) { reco_sf *= 1.00728; }
				else if(eta <-1.8   ) { reco_sf *= 0.994819; }
				else if(eta <-1.63  ) { reco_sf *= 0.994786; }
				else if(eta <-1.566 ) { reco_sf *= 0.991632; }
				else if(eta <-1.4442) { reco_sf *= 0.963128; }
				else if(eta <-1.2   ) { reco_sf *= 0.989701; }
				else if(eta <-1     ) { reco_sf *= 0.985656; }
				else if(eta <-0.6   ) { reco_sf *= 0.981595; }
				else if(eta <-0.4   ) { reco_sf *= 0.984678; }
				else if(eta <-0.2   ) { reco_sf *= 0.981614; }
				else if(eta < 0     ) { reco_sf *= 0.980433; }
				else if(eta < 0.2   ) { reco_sf *= 0.984552; }
				else if(eta < 0.4   ) { reco_sf *= 0.988764; }
				else if(eta < 0.6   ) { reco_sf *= 0.987743; }
				else if(eta < 1     ) { reco_sf *= 0.987743; }
				else if(eta < 1.2   ) { reco_sf *= 0.987743; }
				else if(eta < 1.4442) { reco_sf *= 0.98768; }
				else if(eta < 1.566 ) { reco_sf *= 0.967598; }
				else if(eta < 1.63  ) { reco_sf *= 0.989627; }
				else if(eta < 1.8   ) { reco_sf *= 0.992761; }
				else if(eta < 2     ) { reco_sf *= 0.991761; }
				else if(eta < 2.2   ) { reco_sf *= 0.99794; }
				else if(eta < 2.3   ) { reco_sf *= 1.00104; }
				else if(eta < 2.4   ) { reco_sf *= 0.989507; }
				else if(eta < 2.45  ) { reco_sf *= 0.970519; }
				else if(eta < 2.5   ) { reco_sf *= 0.906667; }
				else                  { reco_sf *= 1; }
		}
*/


    std::vector<float> binningEta =  {-2.5, -2.0, -1.556, -1.442, -0.8, 0., 0.8, 1.442, 1.556, 2., 2.5};
    int etaBin = 0;
    for (int i = 0; i < int(binningEta.size())-1 ; ++i) {
        if (electron.Eta() > binningEta[i] && electron.Eta() <= binningEta[i+1]) {
            etaBin = i+1;
            break;
        }
    }

    std::vector<float> binningPt = { 10, 20, 35, 50, 90, 2000 };
    int ptBin = 0;
    for (int i = 0; i < int(binningPt.size()) - 1; ++i) {
        if (fabs(electron.Pt()) > binningPt[i] && fabs(electron.Pt()) <= binningPt[i+1]) {
            ptBin = i+1;
            break;
        }
    }

		if(electron.Pt() > 2000)
		{
				ptBin = binningPt.size() - 1;
		}

    float weight = 1;
		float high = 0;
		float low = 0;

		weight   *= _elSF_id_iso.GetBinContent(etaBin, ptBin);

		high			= pow( pow( reco_error, 2.0 ) + pow(_elSF_id_iso.GetBinError(etaBin, ptBin), 2.0), 0.5 );
		low				= high;


    return {weight * reco_sf, high, low};
}

WeightAndSigma WeightUtils::GetBJetWeight(float pt, float eta, int flavor, float csv_score, float csv_cut)
{
	int binningPt_low[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
	int binningPt_high[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1000};

	float btag_eff[] = {
  0.5542576310070667,
  0.5606102430319312,
  0.5914605031773699,
  0.6180488071650029,
  0.6344665209139178,
  0.6479497556938635,
  0.6608028574597213,
  0.6644499981441845,
  0.6721201396126378,
  0.6394896941827578,
  0.616956077630235,
  0.5740658362989324,
  0.5450061652281134,
  0.4839797639123103,
  0.437125748502994,
  0.4032258064516129,
  0.25};

	float cmistag_eff[] = {
  0.098510927868065,
  0.098510927868065,
  0.09679576483700195,
  0.09573161369568556,
  0.09555001039717197,
  0.09905397885364496,
  0.11177688201329121,
  0.10397376543209877,
  0.12468683754095201,
  0.11891679748822606,
  0.11082024432809773,
  0.11702127659574468,
  0.11559139784946236,
  0.1115702479338843,
  0.1,
  0.19047619047619047,
  0.11764705882352941};

	float lmistag_eff[] = {
  0.010224535206605491,
  0.010224535206605491,
  0.0090752110514198,
  0.00806412198817262,
  0.008831463864830328,
  0.008478437431905073,
  0.008668121599606553,
  0.009386300787066582,
  0.011345999464811346,
  0.010502471169686986,
  0.013101476893758932,
  0.018541930046354824,
  0.01507537688442211,
  0.030280649926144758,
  0.015306122448979591,
  0.030054644808743168,
  0.012345679012345678};






	float mc_eff = 0;
	BTagEntry::JetFlavor flavor_tag;

	if (flavor == 5)			flavor_tag = BTagEntry::JetFlavor::FLAV_B;
	else if (flavor == 4) flavor_tag = BTagEntry::JetFlavor::FLAV_C;
	else									flavor_tag = BTagEntry::JetFlavor::FLAV_UDSG;

	for(int i = 0; i < 16; i++){

		if(pt > binningPt_low[i] && pt <= binningPt_high[i]){
			if(flavor_tag == BTagEntry::JetFlavor::FLAV_B) mc_eff = btag_eff[i];
			else if(flavor_tag == BTagEntry::JetFlavor::FLAV_C) mc_eff = cmistag_eff[i];
			else mc_eff = lmistag_eff[i];
		}

	 if( pt < binningPt_high[0] ){
			if(flavor_tag == BTagEntry::JetFlavor::FLAV_B )			mc_eff = btag_eff[0];
			if(flavor_tag == BTagEntry::JetFlavor::FLAV_C )			mc_eff = cmistag_eff[0];
			if(flavor_tag == BTagEntry::JetFlavor::FLAV_UDSG )	mc_eff = lmistag_eff[0];
		}
		if( pt > binningPt_high[15] ){
			if(flavor_tag == BTagEntry::JetFlavor::FLAV_B )			mc_eff = btag_eff[15];
			if(flavor_tag == BTagEntry::JetFlavor::FLAV_C )			mc_eff = cmistag_eff[15];
			if(flavor_tag == BTagEntry::JetFlavor::FLAV_UDSG )	mc_eff = lmistag_eff[15];
		}
	}

    float weight       = 1;
    float weight_up    = 1;
    float weight_down  = 1;


		bool tagged = false;
		if(csv_score > csv_cut){ 
				tagged = true;
		}
		float sf      = _btagreader.eval_auto_bounds("central",	flavor_tag, eta, pt, 0);
		float sf_up   = _btagreader.eval_auto_bounds("up",			flavor_tag, eta, pt, 0);
		float sf_down = _btagreader.eval_auto_bounds("down",		flavor_tag, eta, pt, 0);


		if(tagged == false){
				if(sf != 0){ 
            weight      = (1 - sf * mc_eff) / (1 - mc_eff);
            weight_up   = (1 - sf_up * mc_eff) / (1 - mc_eff);
            weight_down = (1 - sf_down * mc_eff) / (1 - mc_eff);
        }
				else{ 
            weight      = 1; 
            weight_up   = 1; 
            weight_down = 1; 
		    }
		}
		if (tagged == true) {
				if(sf != 0){
            weight      = sf;
            weight_up   = sf_up;
            weight_down = sf_down;
        }
				else{ 
            weight      = 1;
            weight_up   = 1;
            weight_down = 1;
		    }		
		}		
		return {weight, weight_up, weight_down};
}


//TODO
//check to see if this stuff works
//NEED TO AS LEGS 
/*
WeightAndSigma  GetDiElectronTrig(TLorentzVector p4){
		_trigger_dielectron;
		float pt = p4.Pt()
		float eta = p4.Eta()

		auto ptmin  =  get_column(&_trigger_dimuon, "ptmin");
		auto ptmax  =  get_column(&_trigger_dimuon, "ptmax");
		auto etamin =  get_column(&_trigger_dimuon, "etamin");
		auto etamax =  get_column(&_trigger_dimuon, "etamax");

		auto eff    =  get_column(&_trigger_dielectron, "eff");
		auto sigma  =  get_column(&_trigger_dielectron, "sigma");
		
		WeightAndSigma w;
		for (int i = 0; i < _trigger_dielectron.nentries; i++) {
				if( ptmin[i]  < pt  && ptmax[i]  >= pt && 
						etamin[i] < eta && etamax[i] >= eta  ){
						
						w.nominal = eff[i];
						w.up 			= sigma[i];
						w.low 		= sigma[i];
				}
		}

		return w;
}
WeightAndSigma  GetDiMuonTrig(TLorentzVector p4){
		_trigger_dimuon;
		float pt = p4.Pt()
		float eta = p4.Eta()

		auto ptmin  =  get_column(&_trigger_dielectron, "ptmin");
		auto ptmax  =  get_column(&_trigger_dielectron, "ptmax");
		auto etamin =  get_column(&_trigger_dielectron, "etamin");
		auto etamax =  get_column(&_trigger_dielectron, "etamax");

		auto eff       =  get_column(&_trigger_dielectron, "eff");
		auto deff_high =  get_column(&_trigger_dielectron, "deff_high");
		auto deff_low  =  get_column(&_trigger_dielectron, "deff_low");

		
		WeightAndSigma w;
		for (int i = 0; i < _trigger_dielectron.nentries; i++) {
				if( ptmin[i]  < pt  && ptmax[i]  >= pt && 
						etamin[i] < eta && etamax[i] >= eta  ){
						
						w.nominal = eff[i];
						w.up 			= deff_high[i];
						w.low 		= deff_low[i];
				}
		}


		WeightAndSigma w;
		return w;
}
*/
//WeightAndSigma  GetMuonElectronTrig(TLorentzVector p4, TLorentzVector q4){
/*
WeightAndSigma  WeightUtils::GetDiLeptonTrig(TLorentzVector p4, TLorentzVector q4, std::string flavor){
		WeightAndSigma w;

		float pt1 = p4.Pt();
		float eta1 = p4.Eta();

		float pt2  = q4.Pt();
		float eta2 = q4.Eta();
		
		WeightAndSigma w_high;
		CSV* csv_high;
		CSV* csv_low;
		if (flavor =="dimuon"){
				csv_high = &_trigger_dimuon_highleg;
				csv_low  = &_trigger_dimuon_lowleg;
		}
		else if (flavor == "dielectron"){
				csv_high = &_trigger_dielectron_highleg;
				csv_low  = &_trigger_dielectron_lowleg;
		} else if (flavor == "emu"){
				csv_high = &_trigger_muonelectron_highleg;
				csv_low  = &_trigger_muonelectron_lowleg;
		} else{
				printf("Getdilepton trigger got a flavor it did not expect: %s", flavor.c_str());
		}
		
		//TODO
		//overflow index protection needed
		{//high pt leg
				auto ptmin  =  *get_column(csv_high, "ptmin");
				auto ptmax  =  *get_column(csv_high, "ptmax");
				auto etamin =  *get_column(csv_high, "etamin");
				auto etamax =  *get_column(csv_high, "etamax");

				auto eff    =  *get_column(csv_high, "eff");
				auto sigma  =  *get_column(csv_high, "sigma");

				for (int i = 0; i < csv_high->nentries; i++) {
						if( ptmin[i]  < pt1  && ptmax[i]  >= pt1 && 
								etamin[i] < eta1 && etamax[i] >= eta1  ){
								
								w_high.nominal= eff[i];
								w_high.up 		= sigma[i];
								w_high.down 		= sigma[i];
						}
				}
		}
		WeightAndSigma w_low;
		{//low pt leg
				auto ptmin  =  *get_column(csv_low, "ptmin");
				auto ptmax  =  *get_column(csv_low, "ptmax");
				auto etamin =  *get_column(csv_low, "etamin");
				auto etamax =  *get_column(csv_low, "etamax");

				auto eff    =  *get_column(csv_low, "eff");
				auto sigma  =  *get_column(csv_low, "sigma");

				for (int i = 0; i < csv_high->nentries; i++) {
						if( ptmin[i]  < pt2  && ptmax[i]  >= pt2 && 
								etamin[i] < eta2 && etamax[i] >= eta2  ){
								
								w_low.nominal = eff[i];
								w_low.up 			= sigma[i];
								w_low.down 		= sigma[i];
						}
				}
		}

		return w;
}
*/
WeightAndSigma  WeightUtils::TEMP_GetDiLeptonTrig(TLorentzVector p4, TLorentzVector q4, std::string flavor){
		WeightAndSigma w;

		float pt1 = p4.Pt();
		float eta1 = p4.Eta();

		float pt2  = q4.Pt();
		float eta2 = q4.Eta();

    auto pt1min  =  *get_column(&_trigger_dilepton, "pt1min");
    auto pt1max  =  *get_column(&_trigger_dilepton, "pt1max");
    auto eta1min =  *get_column(&_trigger_dilepton, "eta1min");
    auto eta1max =  *get_column(&_trigger_dilepton, "eta1max");

    auto pt2min  =  *get_column(&_trigger_dilepton, "pt2min");
    auto pt2max  =  *get_column(&_trigger_dilepton, "pt2max");
    auto eta2min =  *get_column(&_trigger_dilepton, "eta2min");
    auto eta2max =  *get_column(&_trigger_dilepton, "eta2max");
	
		std::vector<float> eff;

		if (flavor == "dimuon"){
				auto eff    =  *get_column(&_trigger_dilepton, "mumu_ratio");
		}
		else if (flavor == "dielectron"){
				auto eff    =  *get_column(&_trigger_dilepton, "ee_ratio");
		}
		else if (flavor == "emu"){
				auto eff    =  *get_column(&_trigger_dilepton, "emu_ratio");
		}

		else{
		    printf("Getdilepton trigger got a flavor it did not expect: %s", flavor.c_str());
		}
		for (int i = 0; i < _trigger_dilepton.nentries; i++) {
				if( pt1min[i]  < pt1  && pt1max[i]  >= pt1  && 
						eta1min[i] < eta1 && eta1max[i] >= eta1 &&

				    pt2min[i]  < pt2  && pt2max[i]  >= pt2  && 
						eta2min[i] < eta2 && eta2max[i] >= eta2  ){
						
						w.nominal = eff[i];
						//TODO
						//we need to get real uncertainties but these can come later
						//Jan 20, 2019
						w.up 			= 0.05;
						w.down 		= 0.05;
				}
		}
		return w;
}
 
//NOTE
//Stolen from https://root.cern.ch/doc/master/TH2_8cxx_source.html#l01328
//with minor alterations
float th2_interpolate(const TH2* hist, Double_t x, Double_t y)
{
   float f=0;
   float x1=0,x2=0,y1=0,y2=0;
   float dx,dy;
   int bin_x = hist->GetXaxis()->FindBin(x);
   int bin_y = hist->GetYaxis()->FindBin(y);


   if( bin_x < 1 )bin_x = 1;
	 if( bin_x > hist->GetNbinsX()) bin_x = hist->GetNbinsX();
	 
	 if( bin_y < 1 ) bin_y = 1;
	 if( bin_y > hist->GetNbinsY()) bin_y = hist->GetNbinsY();



   Int_t quadrant = 0; // CCW from UR 1,2,3,4
   // which quadrant of the bin (bin_P) are we in?
   dx = hist->GetXaxis()->GetBinUpEdge(bin_x)-x;
   dy = hist->GetYaxis()->GetBinUpEdge(bin_y)-y;
   if (dx <= hist->GetXaxis()->GetBinWidth(bin_x)/2 && dy <= hist->GetYaxis()->GetBinWidth(bin_y)/2)
   quadrant = 1; // upper right
   if (dx > hist->GetXaxis()->GetBinWidth(bin_x)/2 && dy <= hist->GetYaxis()->GetBinWidth(bin_y)/2)
   quadrant = 2; // upper left
   if (dx > hist->GetXaxis()->GetBinWidth(bin_x)/2 && dy > hist->GetYaxis()->GetBinWidth(bin_y)/2)
   quadrant = 3; // lower left
   if (dx <= hist->GetXaxis()->GetBinWidth(bin_x)/2 && dy > hist->GetYaxis()->GetBinWidth(bin_y)/2)
   quadrant = 4; // lower right
   switch(quadrant) {
   case 1:
      x1 = hist->GetXaxis()->GetBinCenter(bin_x);
      y1 = hist->GetYaxis()->GetBinCenter(bin_y);
      x2 = hist->GetXaxis()->GetBinCenter(bin_x+1);
      y2 = hist->GetYaxis()->GetBinCenter(bin_y+1);
      break;
   case 2:
      x1 = hist->GetXaxis()->GetBinCenter(bin_x-1);
      y1 = hist->GetYaxis()->GetBinCenter(bin_y);
      x2 = hist->GetXaxis()->GetBinCenter(bin_x);
      y2 = hist->GetYaxis()->GetBinCenter(bin_y+1);
      break;
   case 3:
      x1 = hist->GetXaxis()->GetBinCenter(bin_x-1);
      y1 = hist->GetYaxis()->GetBinCenter(bin_y-1);
      x2 = hist->GetXaxis()->GetBinCenter(bin_x);
      y2 = hist->GetYaxis()->GetBinCenter(bin_y);
      break;
   case 4:
      x1 = hist->GetXaxis()->GetBinCenter(bin_x);
      y1 = hist->GetYaxis()->GetBinCenter(bin_y-1);
      x2 = hist->GetXaxis()->GetBinCenter(bin_x+1);
      y2 = hist->GetYaxis()->GetBinCenter(bin_y);
      break;
   }

   int bin_x1 = hist->GetXaxis()->FindBin(x1);
   if(bin_x1 < 1) bin_x1 = 1;

   int bin_x2 = hist->GetXaxis()->FindBin(x2);
   if(bin_x2 > hist->GetNbinsX()) bin_x2 = hist->GetNbinsX();

   int bin_y1 = hist->GetYaxis()->FindBin(y1);
   if(bin_y1 < 1) bin_y1=1;

   int bin_y2 = hist->GetYaxis()->FindBin(y2);
   if(bin_y2 > hist->GetNbinsY()) bin_y2 = hist->GetNbinsY();

   int bin_q22 = hist->GetBin(bin_x2, bin_y2);
   int bin_q12 = hist->GetBin(bin_x1, bin_y2);
   int bin_q11 = hist->GetBin(bin_x1, bin_y1);
   int bin_q21 = hist->GetBin(bin_x2, bin_y1);
   float q11 = hist->GetBinContent(bin_q11);
   float q12 = hist->GetBinContent(bin_q12);
   float q21 = hist->GetBinContent(bin_q21);
   float q22 = hist->GetBinContent(bin_q22);
   float d = 1.0*(x2-x1)*(y2-y1);
   f = 1.0*q11/d*(x2-x)*(y2-y) + 1.0*q21/d*(x-x1)*(y2-y) + 1.0*q12/d*(x2-x)*(y-y1) + 1.0*q22/d*(x-x1)*(y-y1);

   return f;
}
