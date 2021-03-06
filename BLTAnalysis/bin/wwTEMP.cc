#include "AnalysisTemplate.hh"

//BLT_II dependancies
#include "BLT_II/BLTAnalysis/interface/SelectionCuts.hh"
#include "BLT_II/BLTAnalysis/interface/TGPhysicsObject.hh"
#include "BLT_II/BLTAnalysis/interface/tgUtility.hh"
#include "BLT_II/BLTAnalysis/interface/WeightUtils.h"
#include "BLT_II/BLTAnalysis/interface/RoccoR.h"
#include "BLT_II/BLTAnalysis/src/SDF.cc"
#include "BLT_II/BLTAnalysis/interface/CSVReader.hh"
#include "BLT_II/BLTAnalysis/src/RandomForest.cc"
#include "BLT_II/BLTAnalysis/interface/ElectronCorrector.h"

//C++ libraries
#include <algorithm>    // std::find std::max
#include <cmath>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */




                  //**  Lepton Variables  **//
std::string variables = "lep1_pt float;"
                        "lep2_pt float;"
                        "lep3_pt float;"
                        "lep1_phi float;"
                        "lep2_phi float;"
                        "lep3_phi float;"
                        "lep1_eta float;"
                        "lep2_eta float;"
                        "lep3_eta float;"
                        "lep1_q int;"
                        "lep2_q int;"
                        "lep3_q int;"
                        "lep1_type int;"
                        "lep2_type int;"
                        "lep3_type int;"
                        "numb_leptons int;"
                        "numb_loose_leptons int;"
                        "dilepton_type int;"
                        "dilepton_type_debug int;"
                        "dPhill float;"
                        "dPhillmet float;"
                        "mll float;"
                        "mllmet float;"
                        "mllmet_jetres float;"
                        "mllmet_jetscale_up float;"
                        "mllmet_jetscale_down float;"
                        "avg_muon_pt_correction float;"
                        "qt float;"

                        "genlep1_type int;"
                        "genlep2_type int;"

                //** Jet Variables begin here **//
                        "genjet1_pt float;"
                        "genjet1_phi float;"
                        "genjet1_eta float;"

                        "pujet1_pt float;"
                        "pujet1_phi float;"
                        "pujet1_eta float;"

                        "jet1_pt float;"
                        "jet2_pt float;"
                        "jet3_pt float;"
                        "jet4_pt float;"
                        "jet5_pt float;"
                        "jet6_pt float;"
                        "jet1_phi float;"
                        "jet2_phi float;"
                        "jet3_phi float;"
                        "jet4_phi float;"
                        "jet5_phi float;"
                        "jet6_phi float;"
                        "jet1_eta float;"
                        "jet2_eta float;"
                        "jet3_eta float;"
                        "jet4_eta float;"
                        "jet5_eta float;"
                        "jet6_eta float;"
                        "jet1_csv float;"
                        "jet2_csv float;"
                        "jet3_csv float;"
                        "jet4_csv float;"
                        "jet5_csv float;"
                        "jet6_csv float;"
                        "jet1_flv int;"
                        "jet2_flv int;"
                        "jet3_flv int;"
                        "jet4_flv int;"
                        "jet5_flv int;"
                        "jet6_flv int;"
                        "jet1_res float;"
                        "jet2_res float;"
                        "jet3_res float;"
                        "jet4_res float;"
                        "jet5_res float;"
                        "jet6_res float;"
                        "jet1_scale_up float;"
                        "jet2_scale_up float;"
                        "jet3_scale_up float;"
                        "jet4_scale_up float;"
                        "jet5_scale_up float;"
                        "jet6_scale_up float;"
                        "jet1_scale_down float;"
                        "jet2_scale_down float;"
                        "jet3_scale_down float;"
                        "jet4_scale_down float;"
                        "jet5_scale_down float;"
                        "jet6_scale_down float;"
                        "HT float;"
                        "HT_res float;"
                        "HT_scale_up float;"
                        "HT_scale_down float;"

                        "recoil_old float;"
                        "recoil_old_res float;"
                        "recoil_old_scale_up float;"
                        "recoil_old_scale_down float;"

                        "recoil float;"
                        "recoil_res float;"
                        "recoil_scale_up float;"
                        "recoil_scale_down float;"

                        "numb_jets int;"  
                        "numb_jets_res int;"  
                        "numb_jets_scale_up int;" 
                        "numb_jets_scale_down int;" 

                        "numb_jets_24 int;"  
                        "numb_jets_24_res int;"  
                        "numb_jets_24_scale_up int;" 
                        "numb_jets_24_scale_down int;" 

                        "numb_genjets_24 int;"  


                        "numb_bjets int;"
                        "numb_bjets_gen int;"


                        "dPhilljet float;"
                        "dPhimetjet float;"
                        "dPhilljet_res float;"
                        "dPhimetjet_res float;"
                        "dPhilljet_scale_up float;"
                        "dPhimetjet_scale_up float;"
                        "dPhilljet_scale_down float;"
                        "dPhimetjet_scale_down float;"

                //** MET Variables begin here **//

                        "met float;"
                        "met_jetres float;"
                        "met_jetscale_up float;"
                        "met_jetscale_down float;"
                        "met_phi float;"
                        "met_filterflag int;"
                        "met_filterflag_recommended int;"
                        "met_filterflag_recommended_old int;"
                        "met_proj float;"
                        "met_proj_jetres float;"
                        "met_proj_jetscale_up float;"
                        "met_proj_jetscale_down float;"
                        "met_over_sET float;"

                        "lepmet_pt float;"

                //** Weight Variables begin here **//

                        "lhe_nominal float;"
                        "qcd_weight float;"
            
            // muR, muF scale variation mappings:
            // LHEWeight[0] ==>     muR,     muF
            // LHEWeight[1] ==>     muR,   2*miuF
            // LHEWeight[2] ==>     muR, 0.5*muF
            // LHEWeight[3] ==>   2*muR,     muF
            // LHEWeight[4] ==>   2*muR,   2*muF
            // LHEWeight[5] ==>   2*muR, 0.5*muF
            // LHEWeight[6] ==> 0.5*muR,     muF
            // LHEWeight[7] ==> 0.5*muR,   2*muF
            // LHEWeight[8] ==> 0.5*muR, 0.5*muF

                        "qcd1_weight float;"
                        "qcd2_weight float;"
                        "qcd3_weight float;"
                        "qcd4_weight float;"
                        "qcd6_weight float;"
                        "qcd8_weight float;"
                        "qcd9_weight float;"  //<= WTF is this?  There is no 9th weight and it never gets filled.
                        "pdf_weight float;"

                        "pu_weight float;"
                        "pu_weight_up float;"
                        "pu_weight_down float;"

                        "alpha_s float;"

                        "lep1_reco_weight float;"
                        "lep1_reco_weight_up float;"
                        "lep1_reco_weight_down float;"
                        "lep1_trigger_weight float;"
                        "lep1_trigger_weight_up float;"
                        "lep1_trigger_weight_down float;"
                        "lep1_weight float;"

                        "lep2_reco_weight float;"
                        "lep2_reco_weight_up float;"
                        "lep2_reco_weight_down float;"
                        "lep2_trigger_weight float;"
                        "lep2_trigger_weight_up float;"
                        "lep2_trigger_weight_down float;"
                        "lep2_weight float;"

                        "trigger_weight float;"
                        "trigger_weight_up float;"
                        "trigger_weight_down float;"

                        "bjet_weight float;"
                        "bjet_weight_up float;"
                        "bjet_weight_down float;"
                        "bjet_l_weight float;"
                        "bjet_l_weight_up float;"
                        "bjet_l_weight_down float;"
                        "bjet_bc_weight float;"
                        "bjet_bc_weight_up float;"
                        "bjet_bc_weight_down float;"

                        "gen_weight float;"
                        "weight float;"
                        "process string;"
                        "process_decay string;"

                //** WW and TTbar specific Variables begin here **//

                        "wwpt float;"
                        "ttbar1_pt float;"
                        "ttbar2_pt float;"
                        "ttbar_sf float;"
                
                //** Event ID Variables begin here **//

                        "lumiSec int;"
                        "runNumb int;"
                        "eventNumb int;"
                        "nPUmean float;"
                        "npv int;"

                        "_pre_fTT float;"
                        "_pre_fDY float;"

                        "ewkcorr float;"
                        "ewkcorre float;"
                        ; 


struct TGHist{
		int n_bins;
		std::vector<float> bin_edges;
		std::vector<float> bin_contents;
};

TGHist init_pdfhist;
TGHist rf_pdfhist;

void init_tghist( float begin, float end, int n_bins, TGHist* hist){
		hist->n_bins = n_bins;
		float width = (end - begin) / (float) n_bins;
		for(int i = 0 ; i <= n_bins; i++){
				hist->bin_edges.push_back( begin + i*width );
				if( i < n_bins) hist->bin_contents.push_back( 0.0 );
		}
};

void update_tghist(float f, float w, TGHist* hist){
		for(int i = 0; i < hist->n_bins; i++){
				if ( (hist->bin_edges[i] <= f) && (f < hist->bin_edges[i+1] )){
						hist->bin_contents[i] += w;
				}
		}
};

std::string print_tghist( TGHist* hist){
		std::string s;
		s += "Number of bins: " + std::to_string(hist->n_bins);
		s += "Bin edges: \n";
		FOR_IN(it, hist->bin_edges){
				s += std::to_string(*it) + " ";
		}
		s += "\n";

		s += "Bin contents: \n";
		FOR_IN(it, hist->bin_contents){
				s += std::to_string(*it) + " ";
		}
		return s;	
}



struct notes{
    std::vector<std::string> names;
    std::vector<int> counts;
};

notes sdf_notes;

stringSDF construct_note(notes n){
    std::string str_notes = "";
    {
        ENUMERATE_IN(i, it, n.names){
           str_notes += *it + ":" + std::to_string(n.counts[i]) + ";"; 
        }
    }
    stringSDF note;
		for(int i = 0 ; i < 256; i++){
				note.characters[i] = '\0';
		}
    strcpy(note.characters, str_notes.c_str());
    return note;
};

void update_note(std::string str){
    if ( sdf_notes.names.size() > 0 ){
        bool found_name = false;
        ENUMERATE_IN(i, it, sdf_notes.names){
            if( *it == str){
                sdf_notes.counts[i]++;
                found_name = true;
                break;
            }
        }
        if( found_name == false ){
            sdf_notes.names.push_back(str);
            sdf_notes.counts.push_back(0);
        }
    }
    else{
        sdf_notes.names.push_back(str);
        sdf_notes.counts.push_back(0);
    }
};


//NOTE
//Hist the effect of hard events cuts
//Ends event record ends if cut is true 
static int count_cuts = 0;
#define END_RECORDING(hist, cut, cut_label) {\
hist.Fill( count_cuts );            \
update_note(cut_label);  \
count_cuts++;                       \
  if (cut){                         \
    count_cuts=0;                   \
    return kTRUE;                   \
  }                                 \
}

std::vector<std::string> else_print_messages;
#define ELSE_PRINT(message)                                                                   \
  else{                                                                                       \
      if(0 == std::count(else_print_messages.begin(), else_print_messages.end(), message) ){  \
          printf("\e[31m%s\e[0m", message);                                                   \
          else_print_messages.push_back(message);                                             \
          printf("\n");                                                                       \
      }                                                                                       \
  }                                                                                           \





std::string ascii = "#Run number, lumi number, event number, lepton category(mm = 0 ee = 1 and em == 2), event weight, MC normalization, PU, lep1, lep2, trigger, EWK , number of jets\n";

//////////////////////////////////////
//CLEANUP: 
//WTF is this stuff
std::string process; 
std::string process_decay; 
int n_dimuons = 0;
int n_dielec  = 0;
int n_emu  = 0;
int gMuons = 0;
int gElectrons = 0;
//////////////////////////////////////



TH1F tot_nEvents;
TH1F nEvents ;
TH1F nPV_plot;
TH1F init_pdf_plot;
TH1F rf_pdf_plot; 
TH1F dy_pdf_plot; 
TH1F tt_pdf_plot; 
TH1F gen_lep_types;
TH1F rf_gen_lep_types;


TFile* fWWEWKCorrFile;
TH1D *fhDWWEWKCorr;

std::vector<std::string> double_trigger_arr = { 
     //"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
     "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
     "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",

                  //DiMuon
     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
     "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",

                 //DiElectron
     "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", //NEW
     "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
     "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*"
};
std::vector<std::string> single_trigger_arr = { "HLT_IsoMu24_v*",
                                            "HLT_IsoTkMu24_v*",
                                            "HLT_Ele27_WPTight_Gsf_v*"};




std::map<std::string, float>          vars_float;
std::map<std::string, int32_t>        vars_int;
std::map<std::string, std::string>    vars_string;
//NOTE:
//Maybe we want to do things like this instead?
//As of Sept 10, 2018 
//This does not work it can not properly deduce the type of T
template <class T>
T* _get_value(std::string name){
    //Note look through all maps and determine if 
    if ( vars_string.count(name) > 0 ) return &(vars_string.find(name)->second);
    if ( vars_int.count(name) > 0 ) return &(vars_int.find(name)->second);
    if ( vars_float.count(name) > 0 ) return &(vars_float.find(name)->second);
    printf("\e[31mKey %s not found\e[0m\n", name.c_str());
    exit(EXIT_FAILURE);
}


SerialDataFormat TKGfile;
bool initialized_progressbar;
ProgressBar progressbar;


WeightUtils weights;
RoccoR      muonCorr;
EnergyScaleCorrection electronCorr;


Forest TTforest;
Forest DYforest;
		
struct RFScore{
		float dy;
		float tt;
};

RFScore get_rfscore(){
		RFScore s = {0.0,0.0};
		{//TTforest
				std::vector<float> arr;
				FOR_IN(it, TTforest.features){
						if ( vars_int.count(*it) ){
								arr.push_back((float) *get_value(&vars_int, *it));
						}
						else if ( vars_float.count(*it) ){
								arr.push_back(*get_value(&vars_float, *it));
						}
						else {
								printf("\e[31mKey %s not found\e[0m\n", (*it).c_str());
								exit(EXIT_FAILURE);
						}
				}
				s.tt = score_forest( arr, &TTforest);
		}

		{//DYforest
				std::vector<float> arr;
				FOR_IN(it, DYforest.features){
						if ( vars_int.count(*it) ){
								arr.push_back((float) *get_value(&vars_int, *it));
						}
						else if ( vars_float.count(*it) ){
								arr.push_back(*get_value(&vars_float, *it));
						}
						else {
								printf("\e[31mKey %s not found\e[0m\n", (*it).c_str());
								exit(EXIT_FAILURE);
						}
				}
				s.dy = score_forest( arr, &DYforest);
		}
		return s;
};



//////////////////
//    Jet resolution 
//////////////////
float jet_resolution_pt_variation(TLorentzVector jet){
    std::vector<float> etalist  = {0.0,   0.5,    0.8,    1.1,  1.3,    1.7,  1.9,    2.1,  2.3,    2.5,  2.8,    3.0,  3.2,    5.0};
    std::vector<float> sflist   = {1.109, 1.138, 1.114, 1.123, 1.084, 1.082, 1.140, 1.067, 1.177, 1.364, 1.857, 1.328, 1.16};
    std::vector<float> unclist  = {0.008, 0.013, 0.013, 0.024, 0.011, 0.035, 0.047, 0.053, 0.041, 0.039, 0.071, 0.022, 0.029};

    float variation  = 1.0 ;
    for(int eta_bin= 1; eta_bin < int(etalist.size()); eta_bin++){
        if(fabs(jet.Eta()) > etalist[eta_bin-1] && fabs(jet.Eta()) <= etalist[eta_bin]){
            float sf  = sflist[eta_bin - 1];
            float unc = unclist[eta_bin - 1];

            variation = 1.0 + gRandom->Gaus(0, unc) * pow( std::max(pow(sf, 2.0) - 1.0, 0.0), 0.5);
            break;
        }
    }

    return variation;
}


float jet_scale_pt_variation(TLorentzVector jet, float up_down){
    std::vector<float> bins = { 0.03461166,  0.03479799,  0.04091238,  0.02792759,  0.02824131, 0.02832819,  0.02825601,  0.03989046,  0.03549729,  0.03493907,
                                0.03136509,  0.03090756,  0.03668286,  0.02394859,  0.02239104, 0.02249762,  0.02307086,  0.03628306,  0.03031631,  0.03207012,
                                0.02975871,  0.02853599,  0.0352883 ,  0.02129324,  0.01918865, 0.01892064,  0.02077935,  0.03521946,  0.0287719 ,  0.02968091,
                                0.02629722,  0.02502647,  0.03269131,  0.01982166,  0.0166834 , 0.0166285 ,  0.01929461,  0.02990752,  0.02622172,  0.02388362,
                                0.02129543,  0.02493503,  0.02903984,  0.01775597,  0.01487546, 0.01472155,  0.01831201,  0.02790234,  0.02409434,  0.02122127,
                                0.0203598 ,  0.01978004,  0.0313457 ,  0.01490755,  0.01351365, 0.01340035,  0.0164281 ,  0.02502149,  0.0210477 ,  0.04091238,
                                0.04091238,  0.01973345,  0.02326544,  0.01584168,  0.01243062, 0.01259028,  0.01568127,  0.02699364,  0.01805894,  0.04091238,
                                0.04091238,  0.01674347,  0.02343294,  0.0132059 ,  0.01071732, 0.01056681,  0.01349371,  0.02381255,  0.01539776,  0.04091238,
                                0.04091238,  0.0147737 ,  0.01997136,  0.01144881,  0.0084123 , 0.00831597,  0.01163707,  0.02156332,  0.0166752 ,  0.04091238,
                                0.04091238,  0.04091238,  0.01868997,  0.00934076,  0.00733844, 0.00707609,  0.00989171,  0.01639932,  0.01418547,  0.04091238,
                                0.04091238,  0.04091238,  0.01235716,  0.0074562 ,  0.00614868, 0.00607702,  0.00878479,  0.01790062,  0.04091238,  0.04091238,
                                0.04091238,  0.04091238,  0.04091238,  0.00817963,  0.00543559, 0.00530111,  0.00808943,  0.01196179,  0.04091238,  0.04091238,
                                0.04091238,  0.04091238,  0.01233914,  0.00650495,  0.0050472 , 0.00505865,  0.00766776,  0.00985612,  0.04091238,  0.04091238,
                                0.04091238,  0.04091238,  0.04091238,  0.00687309,  0.00457696, 0.00424243,  0.00658787,  0.00960332,  0.04091238,  0.04091238};
    int ncolumns = 10;
    std::vector<float> etalist = {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<float> ptlist  = {30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 10400.0};
    float variation = 0.0;

    for( int eta_bin = 1; eta_bin < int(etalist.size()); eta_bin++){
        if ( jet.Eta() > etalist[eta_bin - 1] && jet.Eta() <= etalist[eta_bin] ){

            for( int pt_bin = 1; pt_bin < int(ptlist.size()); pt_bin++){
                if ( (jet.Pt() > ptlist[pt_bin - 1] && jet.Pt() <= ptlist[pt_bin]) || (jet.Pt() < ptlist[pt_bin] && pt_bin == 1) ){
                    int bin = (eta_bin - 1) + ((pt_bin-1) * ncolumns);
                    float sf = bins[bin];
                    variation = (1.0 + up_down * sf);
                    break;
                }
            }

        }
    }
    return variation;
}




using namespace baconhep;

DemoAnalyzer::DemoAnalyzer() : BLTSelector()
{

}

DemoAnalyzer::~DemoAnalyzer()
{

}



void DemoAnalyzer::Begin(TTree *tree)
{
		
		///////////////////
    const std::string cmssw_base = getenv("CMSSW_BASE");
		///////////////////


    printf("Parsing Commandline Options\n");
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
    printf("Commandline Options Complete\n");
    

    printf("Parse Feature Set\n");
    enum PARSER_FLAG {
      VAR,
      TYPE,
      ERROR
    };

    
    //NOTE:
    //Parse variable string and creates a map for the variables you want.
    //variables designated 'float' will be put in float map
    //variables designated 'int' will be put in int map
    //variables designated 'string' will be put in string map
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
    
    printf("Parse Feature Set Complete\n");
    //*********************************//
    // Prepare the output tree
    printf("Preparing ROOT File\n");
    outFileName =  options.at(0) + "_" + options.at(2)  +"_" + options.at( options.size() - 1);
    std::string outFileName_ROOT =  options.at(0) + "_" + options.at(2)  +"_" + options.at( options.size() - 1) + ".root";

    outFile = new TFile(outFileName_ROOT.c_str(),"RECREATE");
    outFile->cd();
    tot_nEvents = TH1F("tot_nEvents", "tot_nEvents", 10, 0, 10);;
    nEvents     = TH1F("nEvents", "nEvents", 10, 0, 10); 
    nPV_plot    = TH1F("nPV", "nPV", 50, -0.05, 50.5); 

    init_pdf_plot  = TH1F("init_pdf", "init_pdf", 100, -0.1, 99.9); 
    rf_pdf_plot    = TH1F("rf_pdf", "rf_pdf", 100, -0.1, 99.9); 
    dy_pdf_plot    = TH1F("dy_pdf", "dy_pdf", 100, -0.1, 99.9); 
    tt_pdf_plot    = TH1F("tt_pdf", "tt_pdf", 100, -0.1, 99.9); 

    gen_lep_types  = TH1F("gen_lep_types", "gen_lep_types", 40, 0, 40);
    rf_gen_lep_types  = TH1F("rf_gen_lep_types", "rf_gen_lep_types", 40, 0, 40);

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
    printf("ROOT Outfile set up\n");


    TKGfile.vars_float = &vars_float; 
    TKGfile.vars_int = &vars_int; 
    TKGfile.vars_str = &vars_string;  
    stage_variables( &TKGfile );
  
    
    printf("Setup WW ewk");
    fWWEWKCorrFile = TFile::Open((cmssw_base + "/src/BLT_II/BLTAnalysis/data/WW_EWK_Corr.root").c_str());
    fhDWWEWKCorr = (TH1D*)(fWWEWKCorrFile->Get("ratio_Ptlm")); 
    assert(fhDWWEWKCorr); 

    printf("Setup weights/corrections\n");
    gRandom = new TRandom();
    weights = WeightUtils();

    printf("\t Trigger\n");
    triggerSelector.reset(new baconhep::TTrigger( cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns" ));

    printf("\t LumiMask\n");
    lumiMask = RunLumiRangeMap();
    std::string jsonFileName = cmssw_base + "/src/BLT_II/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile( jsonFileName );

    printf("\t Muon Corrections\n");
    muonCorr.init(cmssw_base + "/src/BLT_II/BLTAnalysis/data/rcdata.2016.v3");

    printf("\t Electron Corrections\n");
    electronCorr = EnergyScaleCorrection(cmssw_base + "/src/BLT_II/BLTAnalysis/data/");
    printf("Complete weights/corrections\n");

    printf("Setting process/process decay names");
    process_decay = options.at(0);  
    process = options.at(1) ; 

		//Load Random Forest
		printf("Loading random forest\n");
		load_forest( cmssw_base + "/src/BLT_II/BLTAnalysis/data/rfs/rfTT", &TTforest);
		//printf("Top forest\n");
		//print_forest(&TTforest);	

		load_forest( cmssw_base + "/src/BLT_II/BLTAnalysis/data/rfs/rfDY", &DYforest);
		//printf("Drell Yan forest\n");
		//print_forest(&DYforest);	
    
    //NOTE
    //THIS IS A TEST TO MAKE SURE NEW ELECTRON WEIGHTS LOAD
    //PLZ REMOVE AFTER TEST IS COMPLETE
    //CSV temp_csv;
    //read_csv( cmssw_base + "/src/BLT_II/BLTAnalysis/data/electrons_cut_WP_Tight80X.txt", &temp_csv);
    //print_csv(&temp_csv);
    //exit(0);



		init_tghist(-0.1, 99.9, 100, &init_pdfhist);
		init_tghist(-0.1, 99.9, 100, &rf_pdfhist);

    /////////////////////////
    //PROGRESS BAR Start Up
    //NOTE: THIS IS A !!!BAD!!! HACK
    //The BLTSelector class, what we use to create our ntuples,
    //has member variable "tree" of type TTree* which can be either
    //a TTree or a TChain. A TChain is a child of the TTree class.
    //When loading many files at once, which happens often when
    //created ntuples using condor but not so often with developing code,
    //the TChain is created and cast to a TTree.
    //The problem which is in regard to which virtual functions.
    //In my analysis code I was calling the GetEntries() virtual member
    //function  defined in TTree and redefined in TChain.
    //When trying to use GetEntries on the "tree" that was
    //originally a TChain the segfault occurs. 

    //NOTE:
    //TO BE CLEAR THE FOLLOWING CODE WILL NOT 
    //PREVENT SEGFAULTS.  IF YOU TRY TO LOAD MORE
    //THAN ONE ROOT FILE IT WILL CRASH!!!

    std::string progress_hostname = "ttgrid01";

    char char_hostname[256];
    gethostname(char_hostname, 256);
    printf("Hostname: %s", char_hostname);
    std::string hostname(char_hostname);

    initialized_progressbar = false;
    if (tree && (hostname == progress_hostname)){
        initialized_progressbar = true;
        progressbar.timer.Start();
        long int max_events = tree->GetEntries();
        progressbar.max_epochs = 100;
        long int epoch = (int) ((float) max_events / (float) progressbar.max_epochs);
        if( epoch  > 10 ) {
            progressbar.epoch = epoch;
        }
    }


    std::cout << " FLOAT 16 to 32 I expect 5.0: " << float16to32(float32to16(-23.00)) << std::endl;
    //std::cout << " FLOAT 32 to 16 I expect 2.0: " << float32to16(2.0) << std::endl;

    ReportPostBegin();
    printf("branches set up\n");
}




Bool_t DemoAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  
    this->totalEvents++;
    
    tot_nEvents.Fill(1);
    nEvents.Fill( count_cuts );
    count_cuts++;


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

    const bool isRealData = (fInfo->runNum != 1);
    if (isRealData) {
        process = "Da";
        process_decay = "Da";
    }
      

    
    //Do lepton stuff
    std::vector<TGPhysObject> muonList;
    std::vector<TGPhysObject> electronList;
    std::vector<TGPhysObject> loosemuonList;
    std::vector<TGPhysObject> looseelectronList;
    {
        if(fMuonArr){
            muonSelection(muonList, fMuonArr);
            loosemuonSelection(loosemuonList, fMuonArr);

        } ELSE_PRINT("fMuonArr is not available.");

        if(fElectronArr){
            tightelectronSelection(electronList, fElectronArr, fInfo->rhoJet);
            looseelectronSelection(looseelectronList, fElectronArr, fInfo->rhoJet);

        } ELSE_PRINT("fElectronArr is not available.");
    }
							

    //CONCATE LEPTONS
    std::vector<TGPhysObject> leptonList;
    float avg_muon_pt_scalefactor = 0.0;
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
        avg_muon_pt_scalefactor += (float)muon_scalefactor;
        muon->particle.muon.pt *= muon_scalefactor;
        leptonList.push_back(*muon);
    }
    avg_muon_pt_scalefactor /= (float)muonList.size();

    FOR_IN(electron, electronList){
        
        if(isRealData){ 
            scaleData sdata = electronCorr.GetScaleData( &*electron, fInfo->runNum);
            electron->particle.electron.pt *= sdata.scale;
        } else {
            float sFactor = electronCorr.GetSmearingFactor( &*electron, 0, 0);
            electron->particle.electron.pt *= gRandom->Gaus(1, sFactor);
        }
        
        leptonList.push_back(*electron);
    }
    tgSort(leptonList);
    tgCleanVector(leptonList);


    //Do Jet stuff
    std::vector<TGPhysObject> jetList;
    {
        if(fAK4CHSArr){
            jetSelection( jetList, fAK4CHSArr, leptonList);
            tgSort(jetList);
            tgCleanVector(jetList);

        } ELSE_PRINT("fAK4CHArr is not available.");

        if(fAK4PuppiArr){
        } ELSE_PRINT("fAK4PUPPIArr is not available.");
    }

    //Do gen particle && gen jet stuff
    float gen_event_info_weight = 1.0;
    std::vector<TGPhysObject> wwList;
    std::vector<TGPhysObject> ttbarList;
    std::vector<TGPhysObject> genLepList;
    std::vector<TGPhysObject> genNeutrinoList;
    std::vector<TGPhysObject> genJetList;
		*get_value(&vars_int,   "numb_genjets_24") = 0;
    {
        if(fGenParticleArr){
            wwSelection( wwList, fGenParticleArr);
            wwGenLepSelection( genLepList, fGenParticleArr);
            wwGenNeutrinoSelection( genNeutrinoList, fGenParticleArr);
            FOR_IN_fARR(it, fGenParticleArr, TGenParticle){
                if (abs(it->pdgId) == 6 && it->status == 62) {
                   ttbarList.push_back(it); 
                }
            }

            tgSort(genLepList);
        } ELSE_PRINT("fGenParticleArr is not available.");

        if(fGenJetArr){

            if(fGenParticleArr)
			          *get_value(&vars_int,   "numb_genjets_24") = genJetSelection(fGenParticleArr, fGenJetArr, 30.0, 2.4, genJetList);

        } ELSE_PRINT("fGenJetArr is not available.");

        if(fGenEvtInfo){
            gen_event_info_weight = fGenEvtInfo->weight;
        } ELSE_PRINT("fGenEvtInfo is not available.");
    }


    //LHE && PV stuffs
    float lhe_nominal_weight = 1.0;
    float alphas_scale = 1.0;
    std::vector<float> lheList;
    std::vector<float> qcdList;
    int tot_npv = -999;
    {
        if(fLHEWeightArr){
            ENUMERATE_IN_fARR(i, lhe, fLHEWeightArr, TLHEWeight){

                if( lhe->id == 1001 ) lhe_nominal_weight = lhe->weight;

                else if( lhe->id > 1001 && lhe->id <= 1009) qcdList.push_back(lhe->weight);

                else if (lhe->id <= 2100) lheList.push_back(lhe->weight);

								else  alphas_scale = lhe->weight;
            }
        } ELSE_PRINT("fLHEWeight is not available.");
        
        if(fPVArr){
            tot_npv = fPVArr->GetEntries();
        } ELSE_PRINT("fPVArr is not available.")
    }


    //////////////////////////////
    //////////////////////////////
    //Saving Features

    //////////////////////////////
    //NOTE
    //Environment Vars
    *get_value(&vars_int,   "lumiSec")   = fInfo->lumiSec;
    *get_value(&vars_int,   "runNumb")   = fInfo->runNum;
    *get_value(&vars_int,   "eventNumb") = fInfo->evtNum;
    *get_value(&vars_float, "nPUmean")   = fInfo->nPUmean;

    *get_value(&vars_string, "process")           = process;
    *get_value(&vars_string, "process_decay")     = process_decay;
      
    //NOTE
    *get_value(&vars_int,   "npv")   = tot_npv;



    //NOTE
    //THIS IS UNSAFE AND IS DONE FOR SPEED UP PURPOSES. FEEL FREE TO REMOVE IF YOU ARE NOT COMFORTABLE
    //
    if (leptonList.size() < 2){
        goto _END_NTUPLE_LABEL;
    }

  
    //////////////////////////////
    //NOTE
    //setting MET Vars
    *get_value(&vars_float, "met")            = fInfo->pfMET;
    *get_value(&vars_float, "met_phi")        = fInfo->pfMETphi;
    *get_value(&vars_int,   "met_filterflag") = fInfo->metFilterFailBits;
    {//Determine if event passes met filter flag selection cuts
        uint32_t m_flag = 0;
        //NOTE
        //Recommeneded by 
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
        std::vector<int> flags = { kGoodVerticesFilter, kGlobalSuperTightHalo2016Filter, kHBHENoiseFilter, 
                                   kHBHENoiseIsoFilter, kECALDeadCellTriggerPrimitiveFilter, kMuonBadTrackFilter,
                                   kChargedHadronTrackResolutionFilter, kEEBadScFilter};
        FOR_IN(flag, flags){
            //printf("%x \n", *flag);
            for( int n = 0; n < 32; n++){
                if (*flag == pow(2, n)){
                    m_flag |= (*flag >> n) & ((fInfo->metFilterFailBits >> n) & 1); 
                    break;
                }
            }
        }
        *get_value(&vars_int,   "met_filterflag_recommended") = m_flag;
        
        uint32_t old_flag = 0;
        std::vector<int> flags_old = { kHBHENoiseFilter, kCSCTightHaloFilter, kHCALLaserEventFilter, kECALDeadCellTriggerPrimitiveFilter, kTrackingFailureFilter,
                                   kEEBadScFilter, kECALLaserCorrFilter, kTrkPOGFilter_manystripclus53X,
                                   kTrkPOGFilter_toomanystripclus53X, kTrkPOGFilter_logErrorTooManyClusters};
        FOR_IN(flag, flags){
            //printf("%x \n", *flag);
            for( int n = 0; n < 32; n++){
                if (*flag == pow(2, n)){
                    old_flag |= (*flag >> n) & ((fInfo->metFilterFailBits >> n) & 1); 
                    break;
                }
            }
        }
        *get_value(&vars_int,   "met_filterflag_recommended_old") = old_flag;
    }
    if (leptonList.size() >= 2){ //MET Projected !!!
        
        float dphi_lep1_met = fabs( deltaPhi(leptonList[0].phi() - fInfo->pfMETphi));
        float dphi_lep2_met = fabs( deltaPhi(leptonList[1].phi() - fInfo->pfMETphi));
        float dphi_l_met    = ( dphi_lep1_met < dphi_lep2_met ) ? dphi_lep1_met : dphi_lep2_met;
        float met_proj = 0.0;

        if( fabs(dphi_l_met) < M_PI/2. ){

            met_proj =  fabs(sin(dphi_l_met)*fInfo->pfMET);
        } else{

            met_proj = fInfo->pfMET;
        }
        *get_value(&vars_float, "met_proj") = met_proj;
    } else{
        *get_value(&vars_float, "met_proj") = fInfo->pfMET;
    }
    *get_value(&vars_float, "met_over_sET") = leptonList.size() >= 2 ? fInfo->pfMET / sqrt( pow(leptonList[0].p4().Et(),2.) + pow(leptonList[1].p4().Et(),2.) ): -999.99 ; 

    
    //////////////////////////////
    //NOTE
    //Setting Weights
    {// Setting default weights
        *get_value(&vars_float, "lep1_reco_weight")    = 1.0;
        *get_value(&vars_float, "lep1_trigger_weight") = 1.0;
        *get_value(&vars_float, "lep2_reco_weight")    = 1.0;
        *get_value(&vars_float, "lep2_trigger_weight") = 1.0;
        *get_value(&vars_float, "bjet_weight")         = 1.0;
    }

    {
        auto pu_weights = weights.GetPUWeight(fInfo->nPUmean);
        *get_value(&vars_float, "pu_weight")      = pu_weights.nominal;
        *get_value(&vars_float, "pu_weight_up")   = pu_weights.up;
        *get_value(&vars_float, "pu_weight_down") = pu_weights.down;

    }
    *get_value(&vars_float, "gen_weight")     = gen_event_info_weight;


    //NOTE
    //QCD weights
    *get_value(&vars_float, "lhe_nominal")    = lhe_nominal_weight;
    
    {
        // muR, muF scale variation mappings:
        // qcd_weight[0] ==>     muR,   2*miuF
        // qcd_weight[1] ==>     muR, 0.5*muF
        // qcd_weight[2] ==>   2*muR,     muF
        // qcd_weight[3] ==>   2*muR,   2*muF
        // qcd_weight[4] ==>   2*muR, 0.5*muF
        // qcd_weight[5] ==> 0.5*muR,     muF
        // qcd_weight[6] ==> 0.5*muR,   2*muF
        // qcd_weight[7] ==> 0.5*muR, 0.5*muF
        float qcd_weight = 1.0;
        ENUMERATE_IN( i, qcd, qcdList){
            if ( i == 4 || i == 6) continue;

            float temp_weights = *qcd / lhe_nominal_weight;

            if ( fabs(1 - qcd_weight) < fabs(1 - temp_weights) ){
                qcd_weight = temp_weights;
            }
            char buffer[15];
            sprintf(buffer, "qcd%i_weight", i+1);
            *get_value(&vars_float, buffer) = temp_weights;
        }
        *get_value(&vars_float, "qcd_weight") = qcd_weight;
    }

    //NOTE
    //PDF weights
    {
        float pdf_weight = 0.0;

        ENUMERATE_IN(i, pdf, lheList){
            pdf_weight += pow(lhe_nominal_weight - *pdf, 2.);
        }
        pdf_weight = pow(pdf_weight/(99. * lhe_nominal_weight * lhe_nominal_weight), .5);
        *get_value(&vars_float, "pdf_weight") = pdf_weight;

        *get_value(&vars_float, "alpha_s") = alphas_scale / lhe_nominal_weight;
    }


    if (leptonList.size() >= 2){
        {// Leading lepton weights
            
            //NOTE
            //It is not clear that the defaults should be set to 1.0.
            //This prevents our ability to discover issues.
            //Things are set to 1.0 now because taus and jets that pass our 
            //selection criteria should be counted. 


            WeightAndSigma lep1_reco_weights    = {1.0, 0.0, 0.0};
            WeightAndSigma lep1_trigger_weights = {1.0, 0.0, 0.0};
            if (leptonList[0].type_flag == MUON){
                lep1_reco_weights     = weights.GetMuonRecoEff(leptonList[0].p4());
                lep1_trigger_weights  = weights.GetMuonTriggerEff("HLT_IsoMu24_v*", leptonList[0].p4());
                
            }
            else if (leptonList[0].type_flag == ELECTRON){
                lep1_reco_weights    = weights.GetElectronTightRecoEff(leptonList[0].p4());
                lep1_trigger_weights = weights.GetElectronTriggerEff( leptonList[0].p4());
            } else{}
            *get_value(&vars_float, "lep1_reco_weight")          = lep1_reco_weights.nominal;
            *get_value(&vars_float, "lep1_trigger_weight")       = lep1_trigger_weights.nominal;

            *get_value(&vars_float, "lep1_reco_weight_up")       = 1 + lep1_reco_weights.up;
            *get_value(&vars_float, "lep1_trigger_weight_up")    = 1 + lep1_trigger_weights.up;

            *get_value(&vars_float, "lep1_reco_weight_down")     = 1 - lep1_reco_weights.down;
            *get_value(&vars_float, "lep1_trigger_weight_down")  = 1 - lep1_trigger_weights.down;

            *get_value(&vars_float, "lep1_weight") = lep1_reco_weights.nominal;//  * lep1_trigger_weights.nominal;
        }
        {// Sub lepton weights
          
            //NOTE
            //It is not clear that the defaults should be set to 1.0.
            //This prevents our ability to discover issues.
            //Things are set to 1.0 now because taus and jets that pass our 
            //selection criteria should be counted. 


            WeightAndSigma lep2_reco_weights    = {1.0, 0.0, 0.0};
            WeightAndSigma lep2_trigger_weights = {1.0, 0.0, 0.0};
            if (leptonList[1].type_flag == MUON){
                lep2_reco_weights     = weights.GetMuonRecoEff(leptonList[1].p4());
                lep2_trigger_weights  = weights.GetMuonTriggerEff("HLT_IsoMu24_v*", leptonList[1].p4());
                
            }
            else if (leptonList[1].type_flag == ELECTRON){
                lep2_reco_weights    = weights.GetElectronTightRecoEff(leptonList[1].p4());
                lep2_trigger_weights = weights.GetElectronTriggerEff(leptonList[1].p4() );
            } else{}
            *get_value(&vars_float, "lep2_reco_weight")          = lep2_reco_weights.nominal;
            *get_value(&vars_float, "lep2_trigger_weight")       = lep2_trigger_weights.nominal;

            *get_value(&vars_float, "lep2_reco_weight_up")       = 1 + lep2_reco_weights.up;
            *get_value(&vars_float, "lep2_trigger_weight_up")    = 1 + lep2_trigger_weights.up;

            *get_value(&vars_float, "lep2_reco_weight_down")     = 1 - lep2_reco_weights.down;
            *get_value(&vars_float, "lep2_trigger_weight_down")  = 1 - lep2_trigger_weights.down;

            *get_value(&vars_float, "lep2_weight") = lep2_reco_weights.nominal ;//* lep2_trigger_weights.nominal;
        }

        //TODO 
        //Complete update on trigger weights
        //Currently Oct 24 2018 we are trying to 
        //check the pt against which trigger is fired 
        enum FUTURE_TRIGGER_COMB{
            SINGLE_MUON = 0,
            SINGLE_ELECTRON,
            DOUBLE_MUON,
            DOUBLE_ELECTRON,
            MUON_ELECTRON,
            DEFAULT_ME
        };

				{//Trigger stuffs Jan 8, 2018
						//Get the total set of triggers passed
						//to do this we need to check if we have the trigger bits 
						//and if the proper leptons have fired
						struct trig_info{
								FUTURE_TRIGGER_COMB trigger_flag;
								std::string trigger_name;
								int leg1; //NOTE index of lepton associated with that leg
								int leg2;
						};
						std::vector<trig_info> good_triggers;

        		if (fInfo->triggerBits != 0){

								FOR_IN(it, single_trigger_arr){
										if( triggerSelector->pass(*it, fInfo->triggerBits) == false ) continue;
										
									  ENUMERATE_IN(i, lepton_it, leptonList){

												if( triggerSelector->passObj(*it, 1, lepton_it->hltMatchBits)){
														if(lepton_it->type_flag == MUON){
																good_triggers.push_back({SINGLE_MUON, it->c_str(), i, -999});
														}
														if(lepton_it->type_flag == ELECTRON){
																good_triggers.push_back({SINGLE_ELECTRON, it->c_str(), i, -999});
														}	
												}	
										}	
								}
								FOR_IN(it, double_trigger_arr){
										if( triggerSelector->pass(*it, fInfo->triggerBits) == false ) continue;
										
									  ENUMERATE_IN(i, lepton_it, leptonList){
									  ENUMERATE_IN(j, lepton_jt, leptonList){
												if (i >= j) continue;

												if( triggerSelector->passObj(*it, 1, lepton_it->hltMatchBits) && triggerSelector->passObj(*it, 2, lepton_jt->hltMatchBits) ){
														if(lepton_it->type_flag == lepton_jt->type_flag){

																if (lepton_it->type_flag == MUON) good_triggers.push_back({DOUBLE_MUON, it->c_str(), i, j});

																else if (lepton_it->type_flag == ELECTRON) good_triggers.push_back({DOUBLE_ELECTRON, it->c_str(), i, j});

														}
														else{
																good_triggers.push_back({MUON_ELECTRON, it->c_str(), i, j});
														}	
												}	

										}	
										}	
								}
								//DEBUG
								//Print out the results of getting good triggers to peek at what we've got/
								//printf(" triggertype  leg1   leg2  trigger name\n");	
								//FOR_IN(it, good_triggers){
								//		printf(" %d           %d        %d  %s \n", (int)it->trigger_flag, it->leg1, it->leg2, it->trigger_name.c_str());	
								//}
						}

					// There was some special thing to take into account when it comes to applying the trigger weights. I need to check nates code for the answers. 
						struct weights_and_unc{
								float w;
								float unc_up;
								float unc_down;
						};
						std::vector<weights_and_unc> triggerweights_arr;
						FOR_IN(it, good_triggers){
								if (it->trigger_flag == DOUBLE_ELECTRON){
										auto triggerweight = weights.TEMP_GetDiLeptonTrig(leptonList[it->leg1].p4(), leptonList[it->leg2].p4(), "dielectron");
										triggerweights_arr.push_back({triggerweight.nominal, 0.0, 0.0});
								}
								else if (it->trigger_flag == DOUBLE_MUON){
										auto triggerweight = weights.TEMP_GetDiLeptonTrig(leptonList[it->leg1].p4(), leptonList[it->leg2].p4(), "dimuon");
										triggerweights_arr.push_back({triggerweight.nominal, 0.0, 0.0});
								}
								else if (it->trigger_flag == MUON_ELECTRON){
										auto triggerweight = weights.TEMP_GetDiLeptonTrig(leptonList[it->leg1].p4(), leptonList[it->leg2].p4(), "emu");
										triggerweights_arr.push_back({triggerweight.nominal, 0.0, 0.0});
								}
								else if      (it->trigger_flag == SINGLE_MUON){
                		auto triggerweight = weights.GetMuonTriggerEff( "HLT_IsoMu24_v*", leptonList[it->leg1].p4());
										triggerweights_arr.push_back({triggerweight.nominal, 0.0, 0.0});
								}
								else if (it->trigger_flag == SINGLE_ELECTRON){
                		auto triggerweight = weights.GetElectronTriggerEff( leptonList[it->leg1].p4());
										triggerweights_arr.push_back({triggerweight.nominal, 0.0, 0.0});
								}
								ELSE_PRINT("Housten we have a problem in the new trigger section has a good trigger with no actual trigger info set.");
						}

						//TODO we should look at unc in the future
						//float temp_unc = 999.9;
						float _final_trigger_weight = 0.0;
						FOR_IN( it, triggerweights_arr){
								//TODO
								//this should be used in the future
								//jan 20, 2019
								//float delta = fabs(it->unc_up - it->unc_down) / 2.0;
								//if ( delta < temp_unc){
								//		temp_unc = delta;
								//		_final_trigger_weight = it->w;
								//}
								if (_final_trigger_weight < it->w)  _final_trigger_weight = it->w;
						}
            if ( _final_trigger_weight == 0){
                //printf("\n\n\t\t\ttrigger fail %f %f %f %f %f %d\n\n", _final_trigger_weight, *get_value(&vars_float, "lep1_trigger_weight"), *get_value(&vars_float, "lep2_trigger_weight"), leptonList[0].pt(), leptonList[0].eta(), (int)leptonList.size());
                _final_trigger_weight = *get_value(&vars_float, "lep1_trigger_weight");
            }
            *get_value(&vars_float, "trigger_weight") = _final_trigger_weight;
				}
				
  
    }
   	 
    {// BJet weights
        WeightAndSigma bjet_weight = {1.0, 1.0, 1.0};
        WeightAndSigma bjet_bc_weight = {1.0, 1.0, 1.0};
        WeightAndSigma bjet_l_weight = {1.0, 1.0, 1.0};
        FOR_IN(jet, jetList){
            auto temp_weight     = weights.GetBJetWeight(jet->pt(), jet->eta(), jet->particle.jet.hadronFlavor, jet->particle.jet.csv, 0.8484);
            bjet_weight.nominal *= temp_weight.nominal;
            bjet_weight.up      *= temp_weight.up;
            bjet_weight.down    *= temp_weight.down;

            //Light Flavor jets
            if ( jet->particle.jet.hadronFlavor != 5 && jet->particle.jet.hadronFlavor != 4 ) {
                bjet_l_weight.nominal *= temp_weight.nominal;
                bjet_l_weight.up      *= temp_weight.up;
                bjet_l_weight.down    *= temp_weight.down;
            }
            //B/C Flavor jets
            if ( jet->particle.jet.hadronFlavor == 5 || jet->particle.jet.hadronFlavor == 4 ) {
                bjet_bc_weight.nominal *= temp_weight.nominal;
                bjet_bc_weight.up      *= temp_weight.up;
                bjet_bc_weight.down    *= temp_weight.down;
            }
        }

        *get_value(&vars_float, "bjet_weight")      = bjet_weight.nominal;
        *get_value(&vars_float, "bjet_weight_up")   = bjet_weight.up;
        *get_value(&vars_float, "bjet_weight_down") = bjet_weight.down;

        *get_value(&vars_float, "bjet_l_weight")      = bjet_l_weight.nominal;
        *get_value(&vars_float, "bjet_l_weight_up")   = bjet_l_weight.up;
        *get_value(&vars_float, "bjet_l_weight_down") = bjet_l_weight.down;

        *get_value(&vars_float, "bjet_bc_weight")      = bjet_bc_weight.nominal;
        *get_value(&vars_float, "bjet_bc_weight_up")   = bjet_bc_weight.up;
        *get_value(&vars_float, "bjet_bc_weight_down") = bjet_bc_weight.down;
    }


    {
        float pu_weight      =   *get_value(&vars_float, "pu_weight"); 
        float lep1_weight    =   *get_value(&vars_float, "lep1_weight");
        float lep2_weight    =   *get_value(&vars_float, "lep2_weight");
        float trigger_weight =   *get_value(&vars_float, "trigger_weight");
        *get_value(&vars_float, "weight")  = pu_weight * lep1_weight * lep2_weight * trigger_weight;
        if(isRealData){ 
            *get_value(&vars_float, "weight")  = 1.0;
        }
    }


    //////////////////////////////
    //NOTE
    //setting lepton pt, eta, phi, stuffs
    {
        ENUMERATE_IN(i, lepton, leptonList){
            //TODO
            //the range cut should not be hard coded
            if( i+1 < 4){
                {
                    char buffer[10];
                    sprintf(buffer, "lep%i_pt", i+1); 
                    *get_value(&vars_float, buffer) = lepton->pt();
                }
                {
                    char buffer[10];
                    sprintf(buffer, "lep%i_phi", i+1);  
                    *get_value(&vars_float, buffer) = lepton->phi();
                }
                {
                    char buffer[10];
                    sprintf(buffer, "lep%i_eta", i+1);  
                    *get_value(&vars_float, buffer) = lepton->eta();
                }
                {
                    char buffer[10];
                    sprintf(buffer, "lep%i_q", i+1);
                    *get_value(&vars_int, buffer) = lepton->q();
                }
                {
                    char buffer[10];
                    sprintf(buffer, "lep%i_type", i+1); 
                    int lep_type = lepton->type_flag == MUON     ? 13 : -999; 
                    lep_type     = lepton->type_flag == ELECTRON ? 11 : lep_type; 
                    *get_value(&vars_int, buffer) = lep_type;

                    if (i < 2){
                        if (i == 0 ) *get_value(&vars_int, "dilepton_type") =  lep_type; 
                        else{
                            if (*get_value(&vars_int, "dilepton_type") !=  lep_type){
                                if (lep_type == 11) *get_value(&vars_int, "dilepton_type") = 1;
                                else if (lep_type == 13) *get_value(&vars_int, "dilepton_type") = 2;
                                else *get_value(&vars_int, "dilepton_type") = -999;
                            }
                            else if ( lep_type == 11 && *get_value(&vars_int, "dilepton_type") == 11 )  *get_value(&vars_int, "dilepton_type") = -2;
                            
                            else if ( lep_type == 13 && *get_value(&vars_int, "dilepton_type") == 13 )  *get_value(&vars_int, "dilepton_type") = -1;

                            else *get_value(&vars_int, "dilepton_type") = -999;
                        }
                    }
                }
            }
        }
    }
    *get_value(&vars_int, "numb_leptons")      = leptonList.size();
    *get_value(&vars_int, "numb_loose_leptons")= loosemuonList.size() + looseelectronList.size();
    *get_value(&vars_float, "mll")             = leptonList.size() >= 2 ? ((leptonList[0].p4()) + (leptonList[1].p4())).Mag() : -999.99;
    *get_value(&vars_float, "dPhill")          = leptonList.size() >= 2 ? fabs( deltaPhi(leptonList[0].phi() - leptonList[1].phi())): -999.9;
    *get_value(&vars_float, "dPhillmet")       = leptonList.size() >= 2 ? fabs( deltaPhi( (leptonList[0].p4() + leptonList[1].p4()).Phi() - fInfo->pfMETphi) ): -999.9;
    *get_value(&vars_float, "avg_muon_pt_correction") = avg_muon_pt_scalefactor;
    *get_value(&vars_float, "qt")                     = leptonList.size() >= 2 ? (leptonList[0].p4() + leptonList[1].p4() ).Pt(): -999.99;
    //OLD
    //*get_value(&vars_float, "qt")                     = leptonList.size() >= 2 ? qtCalc(&leptonList[0], &leptonList[1]) : -999.99;
    *get_value(&vars_int, "dilepton_type_debug") =  *get_value(&vars_int, "dilepton_type"); 
    if (*get_value(&vars_int, "dilepton_type") > 1) *get_value(&vars_int, "dilepton_type") =  1;

    {
        TLorentzVector p4; 
        p4.SetPtEtaPhiM(fInfo->pfMET, 0.0, fInfo->pfMETphi, 0.0);
        *get_value(&vars_float, "lepmet_pt") = leptonList.size() >= 2 ? (leptonList[0].p4() + leptonList[1].p4() + p4).Pt(): -999.9;
    }

    //////////////////////////////
    //NOTE
    //jet related stuffs
    {
        int number_bjets = 0;
        int number_bjets_gen = 0;

        int number_jets = 0;
        int number_jets_res = 0;
        int number_jets_scale_up = 0;
        int number_jets_scale_down = 0;

        int number_jets_24 = 0;
        int number_jets_24_res = 0;
        int number_jets_24_scale_up = 0;
        int number_jets_24_scale_down = 0;

        float HT                = 0.0;
        float HT_res            = 0.0;
        float HT_scale_up       = 0.0;
        float HT_scale_down     = 0.0;
        float tot_jetres        = 0.0;
        float tot_jetscale_up   = 0.0;
        float tot_jetscale_down = 0.0;

        *get_value(&vars_float, "recoil")                 = 0;
        *get_value(&vars_float, "recoil_res")             = 0;
        *get_value(&vars_float, "recoil_scale_up")        = 0;
        *get_value(&vars_float, "recoil_scale_down")      = 0;

        TLorentzVector recoil_vec;
        if (leptonList.size() >= 2){
            TLorentzVector metp4;
            metp4.SetPtEtaPhiM(fInfo->pfMET, 0.0, fInfo->pfMETphi, 0.0);
            recoil_vec                              = leptonList[0].p4() + leptonList[1].p4() + metp4;
            *get_value(&vars_float, "recoil_old")   = recoil_vec.Pt();
        } 
        *get_value(&vars_float, "recoil_old_res")         = recoil_vec.Pt();  
        *get_value(&vars_float, "recoil_old_scale_up")    = recoil_vec.Pt();
        *get_value(&vars_float, "recoil_old_scale_down")  = recoil_vec.Pt();

        if (jetList.size() > 0){
            auto _recoil = jetList[0].p4();
            auto _recoil_res = jetList[0].p4();
            auto _recoil_scale_up = jetList[0].p4();
            auto _recoil_scale_down = jetList[0].p4();
            ENUMERATE_IN(i, jet, jetList){
                float _resolution_jet_pt_variation = jet_resolution_pt_variation(jet->p4());
                float _scale_jet_pt_up = jet_scale_pt_variation(jet->p4(), 1.0);
                float _scale_jet_pt_down = jet_scale_pt_variation(jet->p4(), -1.0);
                
                tot_jetres        += jet->pt() * (_resolution_jet_pt_variation - 1.0);
                tot_jetscale_up   += jet->pt() * (_scale_jet_pt_up - 1.0);
                tot_jetscale_down += jet->pt() * (_scale_jet_pt_down - 1.0);
                if (i+1 <= 6){
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_pt", i+1); 
                        *get_value(&vars_float, buffer) = jet->pt();
                    }
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_phi", i+1);  
                        *get_value(&vars_float, buffer) = jet->phi();
                    }
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_eta", i+1);  
                        *get_value(&vars_float, buffer) = jet->eta();
                    }
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_csv", i+1);  
                        *get_value(&vars_float, buffer) = jet->particle.jet.csv;
                    }
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_flv", i+1);  
                        *get_value(&vars_int, buffer) = jet->particle.jet.hadronFlavor;
                    }
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_res", i+1);  
                        *get_value(&vars_float, buffer) = _resolution_jet_pt_variation;
                    }
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_scale_up", i+1); 
                        *get_value(&vars_float, buffer) = _scale_jet_pt_up;
                    }
                    {
                        char buffer[10];
                        sprintf(buffer, "jet%i_scale_down", i+1); 
                        *get_value(&vars_float, buffer) = _scale_jet_pt_down;
                    }

                }
                if( jet->pt() >= 30){
                    number_jets += 1;
                    HT += jet->pt();
                    if (i == 0) {} 
                    else  _recoil += jet->p4();
                    
                    if (fabs(jet->eta()) < 2.4) number_jets_24 += 1;
                }
                if( jet->pt()*_resolution_jet_pt_variation >= 30 ){
                    number_jets_res += 1;
                    HT_res += jet->pt()*_resolution_jet_pt_variation;

                    if (i == 0) _recoil_res *= _resolution_jet_pt_variation;
                    else  _recoil_res += jet->p4()*_resolution_jet_pt_variation;

                    if (fabs(jet->eta()) < 2.4) number_jets_24_res += 1;
                }
                if( jet->pt()*_scale_jet_pt_up >= 30 ){
                    number_jets_scale_up += 1;
                    HT_scale_up += jet->pt()*_scale_jet_pt_up;

                    if (i == 0) _recoil_scale_up *= _scale_jet_pt_up;
                    else  _recoil_scale_up += jet->p4()*_scale_jet_pt_up;

                    if (fabs(jet->eta()) < 2.4) number_jets_24_scale_up += 1;
                }
                if( jet->pt()*_scale_jet_pt_down >= 30 ){
                    number_jets_scale_down += 1;
                    HT_scale_down += jet->pt()*_scale_jet_pt_down;

                    if (i == 0) _recoil_scale_down *= _scale_jet_pt_down;
                    else  _recoil_scale_down += jet->p4()*_scale_jet_pt_down;

                    if (fabs(jet->eta()) < 2.4) number_jets_24_scale_down += 1;
                }
                //NOTE
                //0.8484 is the Mid working point for bjet identification using the cssv alg.
                //Maybe this should be defined in the Selection headerfile or something
                
                //TODO test the following
                //if( jet->particle.jet.csv >= 0.8484 && fabs(jet->particle.jet.eta) <= 2.4){
                if( jet->particle.jet.csv >= 0.8484 ){
                    number_bjets += 1;
                } 
                //if( abs(jet->particle.jet.hadronFlavor) == 5 && fabs(jet->particle.jet.eta) <= 2.4){
                if( abs(jet->particle.jet.hadronFlavor) == 5 ){ //&& fabs(jet->particle.jet.eta) <= 2.4){
                    number_bjets_gen += 1;
                } 
            }
            //NOTE
            //Why do we check for recoil pt being greater than 30GeV we already do a check when we add recoil pts. 
            //Look above.
            *get_value(&vars_float, "recoil")                 = _recoil.Pt();
            *get_value(&vars_float, "recoil_res")             = _recoil_res.Pt();
            *get_value(&vars_float, "recoil_scale_up")        = _recoil_scale_up.Pt();
            *get_value(&vars_float, "recoil_scale_down")      = _recoil_scale_down.Pt();

            *get_value(&vars_float, "recoil_old_res")         = (recoil_vec - _recoil + _recoil_res).Pt();
            *get_value(&vars_float, "recoil_old_scale_up")    = (recoil_vec - _recoil + _recoil_scale_up).Pt();
            *get_value(&vars_float, "recoil_old_scale_down")  = (recoil_vec - _recoil + _recoil_scale_down).Pt();
        }
        *get_value(&vars_int,   "numb_jets")           = number_jets;
        *get_value(&vars_int,   "numb_jets_res")       = number_jets_res;
        *get_value(&vars_int,   "numb_jets_scale_up")  = number_jets_scale_up;
        *get_value(&vars_int,   "numb_jets_scale_down")  = number_jets_scale_down;

        *get_value(&vars_int,   "numb_jets_24")             = number_jets_24;
        *get_value(&vars_int,   "numb_jets_24_res")         = number_jets_24_res;
        *get_value(&vars_int,   "numb_jets_24_scale_up")    = number_jets_24_scale_up;
        *get_value(&vars_int,   "numb_jets_24_scale_down")  = number_jets_24_scale_down;

        *get_value(&vars_int,   "numb_bjets")     = number_bjets;
        *get_value(&vars_int,   "numb_bjets_gen") = number_bjets_gen;
        *get_value(&vars_float, "HT")             = HT;
        *get_value(&vars_float, "HT_res")         = HT_res;
        *get_value(&vars_float, "HT_scale_up")    = HT_scale_up;
        *get_value(&vars_float, "HT_scale_down")  = HT_scale_down;



        if ( leptonList.size() >= 2 && jetList.size() > 0){
            *get_value(&vars_float, "dPhilljet")  = jetList[0].pt() > 30 ? fabs(deltaPhi((leptonList[0].p4() + leptonList[1].p4()).Phi() - jetList[0].phi())) : -999.99;
            *get_value(&vars_float, "dPhimetjet") = jetList[0].pt() > 30 ? fabs(deltaPhi( fInfo->pfMETphi - jetList[0].phi())) : -999.99;

            *get_value(&vars_float, "dPhilljet_res")  = jetList[0].pt() * *get_value(&vars_float, "jet1_res") > 30 ? fabs(deltaPhi((leptonList[0].p4() + leptonList[1].p4()).Phi() - jetList[0].phi())) : -999.99;
            *get_value(&vars_float, "dPhimetjet_res") = jetList[0].pt() * *get_value(&vars_float, "jet1_res") > 30 ? fabs(deltaPhi( fInfo->pfMETphi - jetList[0].phi())) : -999.99;

            *get_value(&vars_float, "dPhilljet_scale_up")  = jetList[0].pt() * *get_value(&vars_float, "jet1_scale_up") > 30 ? fabs(deltaPhi((leptonList[0].p4() + leptonList[1].p4()).Phi() - jetList[0].phi())) : -999.99;
            *get_value(&vars_float, "dPhimetjet_scale_up") = jetList[0].pt() * *get_value(&vars_float, "jet1_scale_up") > 30 ? fabs(deltaPhi( fInfo->pfMETphi - jetList[0].phi())) : -999.99;
            *get_value(&vars_float, "dPhilljet_scale_down")  = jetList[0].pt() * *get_value(&vars_float, "jet1_scale_down") > 30 ? fabs(deltaPhi((leptonList[0].p4() + leptonList[1].p4()).Phi() - jetList[0].phi())) : -999.99;
            *get_value(&vars_float, "dPhimetjet_scale_down") = jetList[0].pt() * *get_value(&vars_float, "jet1_scale_down") > 30 ? fabs(deltaPhi( fInfo->pfMETphi - jetList[0].phi())) : -999.99;
        }
        *get_value(&vars_float, "met_jetres")           = fInfo->pfMET - tot_jetres;
        *get_value(&vars_float, "met_jetscale_up")      = fInfo->pfMET - tot_jetscale_up;
        *get_value(&vars_float, "met_jetscale_down")    = fInfo->pfMET - tot_jetscale_down;
    }






    if (leptonList.size() >= 2){
        TLorentzVector metp4;
        metp4.SetPtEtaPhiM(fInfo->pfMET, 0.0, fInfo->pfMETphi, 0.0);
        *get_value(&vars_float, "mllmet")   = ((leptonList[0].p4() + leptonList[1].p4()) + metp4).Mag();

        metp4.SetPtEtaPhiM(*get_value(&vars_float, "met_jetres"), 0.0, fInfo->pfMETphi, 0.0);
        *get_value(&vars_float, "mllmet_jetres")   = ((leptonList[0].p4() + leptonList[1].p4()) + metp4).Mag();

        metp4.SetPtEtaPhiM(*get_value(&vars_float, "met_jetscale_down"), 0.0, fInfo->pfMETphi, 0.0);
        *get_value(&vars_float, "mllmet_jetscale_down")   = ((leptonList[0].p4() + leptonList[1].p4()) + metp4).Mag();

        metp4.SetPtEtaPhiM(*get_value(&vars_float, "met_jetscale_up"), 0.0, fInfo->pfMETphi, 0.0);
        *get_value(&vars_float, "mllmet_jetscale_up")   = ((leptonList[0].p4() + leptonList[1].p4()) + metp4).Mag();



        float dphi_lep1_met = fabs( deltaPhi(leptonList[0].phi() - fInfo->pfMETphi));
        float dphi_lep2_met = fabs( deltaPhi(leptonList[1].phi() - fInfo->pfMETphi));
        float dphi_l_met    = ( dphi_lep1_met < dphi_lep2_met ) ? dphi_lep1_met : dphi_lep2_met;
        float met_proj_jetscale_up   = 0.0;
        float met_proj_jetscale_down = 0.0;
        float met_proj_jetres = 0.0;

        if( fabs(dphi_l_met) < M_PI/2. ){

            met_proj_jetscale_up   =  fabs(sin(dphi_l_met)* (*get_value(&vars_float, "met_jetscale_up")));
            met_proj_jetscale_down =  fabs(sin(dphi_l_met)* (*get_value(&vars_float, "met_jetscale_down")));
            met_proj_jetres =  fabs(sin(dphi_l_met)* (*get_value(&vars_float, "met_jetres")));
        } else{

            met_proj_jetscale_up   = *get_value(&vars_float, "met_jetscale_up");
            met_proj_jetscale_down = *get_value(&vars_float, "met_jetscale_down");
            met_proj_jetres        = *get_value(&vars_float, "met_jetres");
        }
        *get_value(&vars_float, "met_proj_jetscale_up")   = met_proj_jetscale_up;
        *get_value(&vars_float, "met_proj_jetscale_down") = met_proj_jetscale_down;
        *get_value(&vars_float, "met_proj_jetres") = met_proj_jetres;

    } 

    //NOTE
    //WWpt calculation && genlepton type records
    {
        float wwpt = 0.0;
        if (wwList.size() >=2){
            wwpt = (wwList[0].p4() + wwList[1].p4()).Pt();
            *get_value(&vars_float, "wwpt")  = wwpt;

        }
        
        if (genLepList.size() >=2){
            *get_value(&vars_int, "genlep1_type")  = abs(genLepList[0].particle.genparticle.pdgId);
            *get_value(&vars_int, "genlep2_type")  = abs(genLepList[1].particle.genparticle.pdgId);
        }
        if(genLepList.size()      >= 2 && 
           genNeutrinoList.size() >= 2 &&
           wwList.size()          >= 2){// Electro Weak corrections
            
            //NOTE
            //This could prob be move outside of the selection loop
            fhDWWEWKCorr->SetDirectory(0);


            float theEWKCorr = 1.0; 
            float theEWKCorrE = 1.0;
            float the_rhoWW = 0.0; 
            float theLeptonHT = genLepList[0].pt() + genLepList[1].pt() + genNeutrinoList[0].pt() + genNeutrinoList[1].pt();
            auto the_rhoP4 =  wwList[0].p4() + wwList[1].p4();

            if(theLeptonHT > 0) the_rhoWW = the_rhoP4.Pt()/theLeptonHT;

            auto lepNegGen = genLepList[0].q() == -1 ? genLepList[0]  : genLepList[1];
            if(lepNegGen.pt() > 0 && the_rhoWW <= 0.3){
                Int_t EWKValbin =  fhDWWEWKCorr->GetXaxis()->FindBin(TMath::Min((float)lepNegGen.pt(), (float)499.999));

                if(EWKValbin >= 0) theEWKCorr  =  fhDWWEWKCorr->GetBinContent(EWKValbin);
                if(EWKValbin >= 0) theEWKCorrE =  1.0+ fhDWWEWKCorr->GetBinError(EWKValbin) / fhDWWEWKCorr->GetBinContent(EWKValbin);
            }
            *get_value(&vars_float, "ewkcorr") = theEWKCorr;
            *get_value(&vars_float, "ewkcorre") = theEWKCorrE;

            *get_value(&vars_float, "weight") *= theEWKCorr;
        }



    }

    //NOTE
    //TTbar pt 
    {
        if (ttbarList.size() >= 2 ){
            *get_value(&vars_float, "ttbar1_pt")  = ttbarList[0].pt();
            *get_value(&vars_float, "ttbar2_pt")  = ttbarList[1].pt();
            //TODO
            // This should be removed to the upper level

            double _sf  = exp(0.0615 - 0.0005*ttbarList[0].pt());
            _sf *= exp(0.0615 - 0.0005*ttbarList[1].pt());

            *get_value(&vars_float, "ttbar_sf") =  powf(_sf, 0.5);
            *get_value(&vars_float, "weight")  *=  powf(_sf, 0.5);

        }
    }

   
		//NOTE
		//TGhist pdf weights
		{
				auto scores = get_rfscore();
				ENUMERATE_IN(i, it, lheList){
						float w = (*it) / lhe_nominal_weight;

						update_tghist( i, w, &init_pdfhist);
						init_pdf_plot.Fill(i, w);

						if( (scores.dy > 0.96) && (scores.tt > 0.6)){
								update_tghist( i, w, &rf_pdfhist);
								rf_pdf_plot.Fill(i, w);
						}
						if( (scores.dy < 0.6) && (scores.tt > 0.6)){
								dy_pdf_plot.Fill(i, w);
						}
						if( (scores.dy > 0.6) && (scores.tt < 0.6)){
								tt_pdf_plot.Fill(i, w);
						}

				}

        *get_value(&vars_float, "_pre_fTT") = scores.tt;
        *get_value(&vars_float, "_pre_fDY") = scores.dy;
		}
 


    {//Recording gen lepton types
        int genlep1 = *get_value(&vars_int, "genlep1_type");
        int genlep2 = *get_value(&vars_int, "genlep2_type");
        
        int _result = genlep1  == 15 || genlep2 == 15 ?  20 : 0;
        
        gen_lep_types.Fill( genlep1 + genlep2  - 20 + _result);

				auto scores = get_rfscore();
        if (scores.dy > 0.96 && scores.tt > 0.6 &&
            leptonList.size() >= 2 && 
            *get_value(&vars_float, "lep1_pt") > 25 &&
            *get_value(&vars_float, "lep2_pt") > 20 &&
            *get_value(&vars_float, "mll") > 12 &&
            *get_value(&vars_int, "met_filterflag_recommended") == 0 
            )
             rf_gen_lep_types.Fill( genlep1 + genlep2 - 20 + _result );
    }




    _END_NTUPLE_LABEL: ;
    //////////   PROGRESS BAR        ///////////////////
    if (initialized_progressbar) progress_update(&progressbar);


    ////////////////////////////////////////////
    //////////   CUTS        ///////////////////
    ////////////////////////////////////////////
  
    //////////////////////////////
    //Number of Leptons cut
    END_RECORDING(nEvents, (leptonList.size() < 2), "Nleptons" )

    //////////////////////////////
    //Low mass cut
    {
        float mll = *get_value(&vars_float, "mll");
        END_RECORDING(nEvents, (mll < 12), "low mass")
    }
    //////////////////////////////
    //Z mass cut
    {   //NOTE
        //This was changed recently might not be working as expected stay vigilant
        float mll = *get_value(&vars_float, "mll");
        bool cut = ( fabs(mll - 91.18) < 15 );
        END_RECORDING(nEvents, cut, "Zmass")
    }
    { 
        bool lep1_pt = *get_value(&vars_float, "lep1_pt") < 25.0;
        bool lep2_pt = *get_value(&vars_float, "lep2_pt") < 20.0;
        bool lep = lep1_pt || lep2_pt;
        END_RECORDING(nEvents, lep, "lep pt")
    }

    //////////////////////////////
    //Trigger

    //NOTE
    //It is unclear if this is the best way to handle the case were no 
    //trigger information is recorded. For that situation there will be 
    //no trigger information in the entire file. If you take the information on
    //an event by event bases you are including addition events that would have
    //had no triggers passed.  Though in some sense this situation
    //is inprobable because the reason an event was save is because some trigger was fired.
    bool pass_trigger = true;
    if (fInfo->triggerBits != 0){
        pass_trigger = false;
        for(unsigned i = 0; i < single_trigger_arr.size(); i++){
            pass_trigger |= triggerSelector->pass( single_trigger_arr[i], fInfo->triggerBits);
        }
        for(unsigned i = 0; i < double_trigger_arr.size(); i++){
            pass_trigger |= triggerSelector->pass( double_trigger_arr[i], fInfo->triggerBits);
        }
    }
    END_RECORDING(nEvents, !pass_trigger, "Trigger")

    //////////////////////////////
    //Good data
    if (isRealData){ 
        RunLumiRangeMap::RunLumiPairType rl( fInfo->runNum, fInfo->lumiSec);
        //NOTE:
        // We skip the event if the mask
        // has the specified lumi id.
        END_RECORDING(nEvents, !lumiMask.HasRunLumi(rl), "GoodLumiRuns")
    }

    //NOTE:
    //end and update save structs 
    
    update_buffers(&TKGfile);
    outTree->Fill();
		/*{// Random forest scores
				auto _temp = get_rfscore();
				printf("%f %f\n", _temp.dy, _temp.tt);
		}*/
    if (*get_value(&vars_int, "dilepton_type") == -1) {
        n_dimuons++;
    } else if (*get_value(&vars_int, "dilepton_type") == -2) {
        n_dielec++;
    } else {
        n_emu++;
    }

   //#Run number, lumi number, event number, lepton category(mm = 0 ee = 1 and em == 2), event weight, MC normalization, PU, lep1, lep2, trigger, EWK , number of jets
    char buffer[256];
    for (int i = 0; i < 255; i++){
        buffer[i] = '\0';
    } 
   
    int runnumber = fInfo->runNum;
    int lumisection = fInfo->lumiSec;
    int eventnumber = fInfo->evtNum;
    int lep_category = *get_value(&vars_int, "dilepton_type");
    if (lep_category == -2) {
        lep_category = 1;
    }
    else if ( lep_category == -1) {
        lep_category = 0;
    }
    else{
        lep_category = 2;
    }
    double event_weight = *get_value(&vars_float, "weight");
    int numb_jets = *get_value(&vars_int, "numb_jets");
    float norm_w = (35.9e3 * (118.7 - 3.974) * pow( (3*.1086), 2) / 1998956.);
    float pu_weight = *get_value(&vars_float, "pu_weight");
    float trigger_weight = *get_value(&vars_float, "trigger_weight");
    float lep1_weight = *get_value(&vars_float, "lep1_weight");
    float lep2_weight = *get_value(&vars_float, "lep2_weight");
    float ewk = *get_value(&vars_float, "ewkcorr");
    float bjet_weight = *get_value(&vars_float, "bjet_weight");

    sprintf(buffer, "%d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %d\n", runnumber, lumisection, eventnumber, lep_category, event_weight * norm_w, norm_w, pu_weight, lep1_weight, lep2_weight, trigger_weight, ewk, bjet_weight, numb_jets); 
    ascii.append(buffer);

    this->passedEvents++;
    END_RECORDING(nEvents, true, "COMPLETE")
}







void DemoAnalyzer::Terminate()
{

    //NOTE
    //Can we remove this forever????
    //No not quite
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
    std::cout << " RECO      : diMuons "  << n_dimuons << "\t" << "diElectrons " << n_dielec << "\tEMu " << n_emu << std::endl;
//    std::cout << " GEN       : Muons "  << gMuons << "\t" << "Electrons " << gElectrons << std::endl;
    std::cout << "  ============================================================" << std::endl;

    auto ascii_file = fopen("ascii_thoth.txt", "w");
    fprintf(ascii_file, ascii.c_str());
    fclose(ascii_file);
    
    //TKGfile.header.notes = construct_note(sdf_notes);
    //strcat(TKGfile.header.notes.characters, "TOTAL:");
    //strcat(TKGfile.header.notes.characters, std::to_string(this->totalEvents).c_str());
    //write_sdf_to_disk( outFileName + ".tdf", &TKGfile );
    
    //read_sdf_from_disk( outFileName + ".tdf");
		//printf("Initial pdf histogram\n");
		//std::cout << print_tghist( &init_pdfhist) << std::endl;
		//printf("Random forest selection pdf histogram\n");
		//std::cout << print_tghist( &rf_pdfhist) << std::endl;

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
