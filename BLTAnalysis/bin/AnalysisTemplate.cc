#include "AnalysisTemplate.hh"

//BLT_II dependancies
#include "BLT_II/BLTAnalysis/interface/SelectionCuts.hh"
#include "BLT_II/BLTAnalysis/interface/TGPhysicsObject.hh"
#include "BLT_II/BLTAnalysis/interface/tgUtility.hh"
#include "BLT_II/BLTAnalysis/interface/WeightUtils.h"
#include "BLT_II/BLTAnalysis/interface/RoccoR.h"
#include "BLT_II/BLTAnalysis/src/SDF.cc"

//C++ libraries
#include <algorithm>    // std::find
#include <cmath>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */





std::string variables =	"lep1_pt float;"
												"lep2_pt float;"
												"avg_muon_correction float;"
												"numb_jets int;"	
												"numb_bjets int;"
												"mll float;"
												"met float;"
												"lhe_nominal float;"
												"pu_weight float;"
												"pu_weight_up float;"
												"pu_weight_down float;"
												"process string;"
												;	





//Hist the effect of hard events cuts
//Ends event record ends if cut is true 
static int count_cuts = 0;
#define END_RECORDING(hist, cut) {\
hist.Fill( count_cuts );			\
count_cuts++;									\
	if (cut){										\
		count_cuts=0;							\
		return kTRUE;							\
	}														\
}

std::vector<std::string> else_print_messages;
#define ELSE_PRINT(message)																																		\
	else{																																												\
			if(0 == std::count(else_print_messages.begin(), else_print_messages.end(), message) ){	\
					printf("\e[31m%s\e[0m", message);																										\
					else_print_messages.push_back(message);																							\
					printf("\n");																																				\
			}																																												\
	}																																														\






//////////////////////////////////////
//CLEANUP: 
//WTF is this stuff
std::string process; 
std::string process_decay; 
int n_muons = 0;
int n_elec  = 0;
int gMuons = 0;
int gElectrons = 0;
//////////////////////////////////////



TH1F tot_nEvents;
TH1F nEvents ;
TH1F nPV_plot;

std::vector<std::string> trigger_vector = {
                                      //MuEG
                                      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
                                      "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
                                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",

                                      //DiMuon
                                      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
                                      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",

                                      //DiElectron
                                      "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*"
                                      };






std::map<std::string, float>					vars_float;
std::map<std::string, int32_t>				vars_int;
std::map<std::string, std::string>		vars_string;

SerialDataFormat TKGfile;

WeightUtils weights;
RoccoR			muonCorr;

//TODO
//Move some where more appropriate
ProgressBar progressbar;

using namespace baconhep;

DemoAnalyzer::DemoAnalyzer() : BLTSelector()
{

}

DemoAnalyzer::~DemoAnalyzer()
{

}



void DemoAnalyzer::Begin(TTree *tree)
{
		const std::string cmssw_base = getenv("CMSSW_BASE");


		//***********************************//
		// Parse command line option
		std::string tmp_option = GetOption();
		std::vector<std::string> options;

		bool more_to_parse = true;
		int i = 0;
		int iter = 0;
		std::cout << tmp_option << std::endl;
		while(more_to_parse == true){
			std::string token = tmp_option.substr(i, tmp_option.find(" ", i, 1) - i );
			options.push_back(token);
			iter += 1;
			if( tmp_option.find(" ", i, 1) == std::string::npos ) more_to_parse = false;
			i = tmp_option.find(" ", i, 1) + 1 ;
		}
		//**********************************//
		//**********************************//
		

		enum PARSER_FLAG {
			VAR,
			TYPE,
			ERROR
		};

		
		//Parse variable string and creates a map for the variables you want.
		//variables designated 'float' will be put in float map
		//variables designated 'int' will be put in int map
		std::string buffer_var  = "";
		std::string buffer_type = "";
		PARSER_FLAG pflag = VAR;
		for(auto var_iter = variables.begin(); var_iter!=variables.end(); ++var_iter){
				if(*var_iter == ';' || (var_iter + 1) == variables.end()){
						//complete assumed information if string is just be for end.
						if( (var_iter + 1) == variables.end()  &&  *var_iter != ';'){
								buffer_type += *var_iter;
						}


						if ( buffer_type.compare("float") == 0){
								vars_float[buffer_var] = 0.0;
						}
						else if (buffer_type.compare("int") == 0){ 
								vars_int[buffer_var] = 0;
						}
						else if (buffer_type.compare("string") == 0){ 
								vars_string[buffer_var] = "";
						}
						else{
								std::cout << "Something fucked up in variable parsing!!! TYPE " << buffer_type << " Name " << buffer_var << std::endl;
						}

						pflag= VAR;
						buffer_var = "";
						buffer_type = "";
				}
				else{
						if(*var_iter == ' ' || *var_iter == '\n' ){
								if(buffer_var.length() == 0 && buffer_type.length() == 0) pflag= VAR;
								else if ( buffer_var.length() != 0 && buffer_type.length() == 0) pflag= TYPE;
								else pflag= ERROR;
						}
						else{
								switch(pflag){
										case VAR :
											buffer_var += *var_iter;	
											break;
										case TYPE :
											buffer_type += *var_iter;	
											break;
										case ERROR :
											printf("What ever TYPE you put in has not been implemented.");
											break;
								}
						}
				}	
		}
		
		//*********************************//
		// Prepare the output tree
		outFileName =  options.at(0) + "_" + options.at(2)  +"_" + options.at( options.size() - 1) + ".root";

		outFile = new TFile(outFileName.c_str(),"RECREATE");
		outFile->cd();
		tot_nEvents = TH1F("tot_nEvents", "tot_nEvents", 10, 0, 10);;
		nEvents			= TH1F("nEvents", "nEvents", 10, 0, 10); 
		nPV_plot		= TH1F("nPV", "nPV", 50, -0.05, 50.5); 


		outTree = new TTree("data", "data_vec");

		//Setting up root file with pointers from our map
		for(auto iter = vars_float.begin(); iter!=vars_float.end(); ++iter){
			outTree->Branch(iter->first.c_str(), &iter->second);
		}
		for(auto iter = vars_int.begin(); iter!=vars_int.end(); ++iter){
			outTree->Branch(iter->first.c_str(), &iter->second);
		}
		for(auto iter = vars_string.begin(); iter!=vars_string.end(); ++iter){
			outTree->Branch(iter->first.c_str(), &iter->second);
		}

		//Should contain full physics object
		//But only a sub set of total events
		//Used to debug cuts...maybe 10% of ALL events
					//Every reco object added should we try to 
					//find complimentary gen object, 
					//if non matched return closest
		//debugTree = new TTree("debug", "debug_vec");
		printf("outfile set up\n");


		TKGfile.vars_float = &vars_float;	
		TKGfile.vars_int = &vars_int;	
		TKGfile.vars_str = &vars_string;	
		stage_variables( &TKGfile );
	


		gRandom = new TRandom();
		weights = WeightUtils();

		triggerSelector.reset(new baconhep::TTrigger( cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns" ));

		lumiMask = RunLumiRangeMap();
		std::string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
		lumiMask.AddJSONFile( jsonFileName );

		muonCorr.init(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");

		/////////////////////////
		//PROGRESS BAR Start Up
		{
				progressbar.timer.Start();
				long int max_events = tree->GetEntries();
				progressbar.max_epochs = 100;
				long int epoch = (int) ((float) max_events / (float) progressbar.max_epochs);
				if( epoch  > 10 ) {
						progressbar.epoch = epoch;
				}
		}

		ReportPostBegin();
		printf("branches set up\n");
}




Bool_t DemoAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  
    this->totalEvents++;
		
		tot_nEvents.Fill(1);
		nEvents.Fill( count_cuts );
    const bool isRealData = (fInfo->runNum != 1);
		
		// Set default values for all variables
		for(auto iter = vars_float.begin(); iter!=vars_float.end(); ++iter){
				iter->second = -999.99;
		}
		for(auto iter = vars_int.begin(); iter!=vars_int.end(); ++iter){
				iter->second = -999;
		}
		for(auto iter = vars_string.begin(); iter!=vars_string.end(); ++iter){
				iter->second = "EMPTY";
		}
		//

		

		
		//Do lepton stuff
		std::vector<TGPhysObject> muonList;
		std::vector<TGPhysObject> electronList;
		{
				if(fMuonArr){
					
						muonSelection(muonList, fMuonArr);
						
						//EXCLAIMATION of JOY
						//FOR_IN_fARR(muon, fMuonArr, TMuon){
						//		std::cout << "For Looped MACRO " << muon->pt << std::endl;
						//}
						//FOR NON BELIEVERS
						//for(int i= 0; i < fMuonArr->GetEntries(); i++){
						//		TMuon* muon = (TMuon*) fMuonArr->At(i);
						//		std::cout << "standard " << muon->pt << std::endl;
						//}
				} ELSE_PRINT("fMuonArr is not available.");

				if(fElectronArr){
						electronSelection(electronList, fElectronArr, fInfo->rhoJet);

				} ELSE_PRINT("fElectronArr is not available.");
		}

		//CONCATE LEPTONS
		std::vector<TGPhysObject> leptonList;
		float avg_muon_scalefactor = 0.0;
		FOR_IN(muon, muonList){
				double muon_scalefactor = 1.0;
				if (isRealData){
						muon_scalefactor *= muonCorr.kScaleDT( muon->particle.muon.q, 
																									 muon->particle.muon.pt,
																									 muon->particle.muon.eta,
																									 muon->particle.muon.phi
																									);
				}
				else{
						muon_scalefactor *= muonCorr.kScaleAndSmearMC( muon->particle.muon.q, 
																													 muon->particle.muon.pt,
																													 muon->particle.muon.eta,
																													 muon->particle.muon.phi,
																													 muon->particle.muon.nTkLayers,
																													 gRandom->Rndm(),
																													 gRandom->Rndm()
																										      );
				}
				avg_muon_scalefactor += (float)muon_scalefactor;
				muon->particle.muon.pt *= muon_scalefactor;
				leptonList.push_back(*muon);
		}
		avg_muon_scalefactor /= (float)muonList.size();

		FOR_IN(electron, electronList){
				leptonList.push_back(*electron);
		}
		tgSort(leptonList);


		//Do Jet stuff
		std::vector<TGPhysObject> jetList;
		{
				if(fAK4CHSArr){
						jetSelection( jetList, fAK4CHSArr, leptonList);

				} ELSE_PRINT("fAK4CHArr is not available.");

				if(fAK4PuppiArr){
				} ELSE_PRINT("fAK4PUPPIArr is not available.");
		}

		//Do gen particle && gen jet stuff
		{
				if(fGenParticleArr){
				} ELSE_PRINT("fGenParticleArr is not available.");

				if(fGenJetArr){
				} ELSE_PRINT("fGenJetArr is not available.");
		}


		//LHE && PV stuffs
		float lhe_nominal_weight = 0.0;
		std::vector<float> lheList;
		std::vector<float> qcdList;
		{
				if(fLHEWeightArr){
						ENUMERATE_IN_fARR(i, lhe, fLHEWeightArr, TLHEWeight){

								if( i == 0 ) lhe_nominal_weight = lhe->weight;

								else if( i < 10) qcdList.push_back(lhe->weight);

								else lheList.push_back(lhe->weight);
						}
				} ELSE_PRINT("fLHEWeight is not available.");
				
				if(fPVArr){
				} ELSE_PRINT("fPVArr is not available.")
		}


		//////////////////////////////
		//////////////////////////////
		//Saving Features
    *get_value(&vars_float, "met") = fInfo->pfMET;
    *get_value(&vars_float, "mll") =  leptonList.size() >= 2 ? ((leptonList[0].p4()) + (leptonList[1].p4())).Mag() : -999.99;
    *get_value(&vars_float, "lhe_nominal") = lhe_nominal_weight;

    *get_value(&vars_float, "pu_weight")      = weights.GetPUWeight(fInfo->nPUmean).nominal;
    *get_value(&vars_float, "pu_weight_up")   = weights.GetPUWeight(fInfo->nPUmean).up;
    *get_value(&vars_float, "pu_weight_down") = weights.GetPUWeight(fInfo->nPUmean).down;

		*get_value(&vars_float, "avg_muon_correction") = avg_muon_scalefactor;

			/*
			//TODO save
  		fGenEvtInfo->weight
			fGenEvtInfo->id_1 
			fGenEvtInfo->id_2 
			fGenEvtInfo->x_1
			fGenEvtInfo->x_2
			fGenEvtInfo->scalePDF
			fGenEvtInfo->xs 
			*/


		//////////////////////////////
		//NOTE
		//setting lepton pt, eta, phi, stuffs
		ENUMERATE_IN(i, lepton, leptonList){
				if( i+1 < 3){
						char buffer[10];
						sprintf(buffer, "lep%i_pt", i+1);	
						*get_value(&vars_float, buffer) = lepton->pt();
				}
		}


		//////////////////////////////
		//NOTE
		//jet related stuffs
		*get_value(&vars_int, "numb_jets")  = jetList.size();

		int number_bjets = 0;
		FOR_IN(jet, jetList){
				//NOTE
				//0.8484 should be defined in the Selection headerfile or something
				if( jet->particle.jet.csv == 0.8484 ){
						number_bjets += 1;
				} 
		}
		*get_value(&vars_int, "numb_bjets")	 = number_bjets;


		
	



		////////////////////////////////////////////
		//////////   CUTS        ///////////////////
		////////////////////////////////////////////
	
		progress_update(&progressbar);
		//////////////////////////////
		//Number of Leptons cut
		END_RECORDING(nEvents, leptonList.size() < 2)


		bool pass_trigger = false;
		if (fInfo->triggerBits != 0){
				for(unsigned i = 0; i < trigger_vector.size(); i++){
						pass_trigger |= triggerSelector->pass( trigger_vector[i], fInfo->triggerBits);
				}
		}
		END_RECORDING(nEvents, !pass_trigger)

		if (isRealData){ 
				RunLumiRangeMap::RunLumiPairType rl( fInfo->runNum, fInfo->lumiSec);
				//NOTE:
				// We skip the event if the mask
				// has the specified lumi id.
				END_RECORDING(nEvents, !lumiMask.HasRunLumi(rl))
		}


	
		//NOTE:
		//end and update save structs	
		update_buffers(&TKGfile);
		outTree->Fill();
    this->passedEvents++;
		END_RECORDING(nEvents, true)
}







void DemoAnalyzer::Terminate()
{

    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void DemoAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
//    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void DemoAnalyzer::ReportPostTerminate()
{

    std::cout << "\n  ==== Terminate Job =========================================" << std::endl;
//    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
//    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
		std::cout << " RECO      : Muons "  << n_muons << "\t" << "Electrons " << n_elec << std::endl;
//		std::cout << " GEN       : Muons "  << gMuons << "\t" << "Electrons " << gElectrons << std::endl;
    std::cout << "  ============================================================" << std::endl;


		//write_sdf_to_disk( outFileName + ".tdf", &TKGfile );
		//read_sdf_from_disk( outFileName + ".tdf");

}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
