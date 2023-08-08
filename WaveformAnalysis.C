#include <iostream>
#include <fstream>
#include "like_ivashkin_wants_it.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TF1.h>
#include <TChain.h>
#include <TString.h>
#include "ChannelEntry.h"
#include "Coeffs.h"
#include "MultiCoreAnalysis.h"
#include "ctime"
#include "Map.h"
#define UsedScatterer 1
using std::cout, std::endl, std::array;
// Амплитуда говно не работает. Запомни, Саня!!!!
const int GATE_BEG = 50;
const int GATE_END = 150;
const Int_t N_zl_intervals = 2;
// const Int_t Gate_end_for_main_scatterers = 75;
// const Int_t Intermediate_gate = 200;
int main(int argc, char** argv)
{
    if(argc == 4) {
        argv[1];
        argv[2];
        argv[3];
        // argv[4];
    }

	TString source_path = (TString)argv[1];
    TString run_name = (TString)argv[2];
    TString temp_DIR = (TString)argv[3];

    const Int_t total_channels = 128;

    ChannelEntry channel_info[total_channels];
    short_energy_ChannelEntry short_channel_info[total_channels];
    // MultiCoreAnalysis mca(channel_info,short_channel_info,6,128);
    // diff_short_energy_ChannelEntry diff_short_channel_info[total_channels];
    Coeffs_struct cal_coeff[total_channels];
    
    TFile *combined_root = new TFile (source_path+temp_DIR+run_name+"_calibrated.root", "RECREATE");
    TTree *combined_tree = new TTree ("adc64_data","adc64_data");
    TTree *coefficients_tree = new TTree ("coefficients","coefficients");

    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
	    combined_tree->Branch((TString::Format("channel_%i", channel_number+1)).Data(), &short_channel_info[channel_number]);
    // for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
    //     //(cal_coeff[channel_number]).CreateBranch(coefficients_tree, channel_number);
	//     coefficients_tree->Branch((TString::Format("channel_%i", channel_number)).Data(), &cal_coeff[channel_number]);
    Float_t interval_width = 1.;
    Int_t gr_point_counter = 0;
    TChain *PMT_tree = new TChain;
    PMT_tree->AddFile( (source_path+run_name+".root" + "/adc64_data").Data() );
    PMT_tree->AddFile( (source_path+"2023_5_5_10_20_2"+".root" + "/adc64_data").Data() );

    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
        (channel_info[channel_number]).SetBranch(PMT_tree,channel_number+1);
    Int_t total_entries = PMT_tree->GetEntries()/1;
    TLeaf *tf = PMT_tree->GetLeaf("time_in_seconds");
///////////////////////        
    Int_t intervals = 1;
    Int_t int_width = total_entries;
    for (intervals = 1; intervals < 1000; intervals++) 
    {
        int_width=total_entries/intervals;
        if (int_width < 2000000) break;
    }
    time_t start_time = time(NULL);
    for (Int_t num_int = 0; num_int < intervals; num_int++)
    {

        Int_t Events_for_cuts_left = num_int*int_width;
        Int_t Events_for_cuts_right = (num_int+1)*int_width;
        if (num_int==intervals-1) Events_for_cuts_right = total_entries;

//////////////////////Fill_clibrated_tree
        for (Long_t entry_num = Events_for_cuts_left;  entry_num< Events_for_cuts_right; entry_num++)
        {
            PMT_tree->GetEntry(entry_num);
            MultiCoreAnalysis mca(channel_info,short_channel_info,8,128);
            // mca.SetChannelsToCalculate(map);
            mca.AnalyzeMultiThread();	
            combined_tree->Fill();
            for (int ch = 0; ch < total_channels; ch++) {channel_info[ch].Initialize();short_channel_info[ch].Initialize();}

            if (entry_num%1000 == 0) 
            {
                std::cout<< u8"\033[2J\033[1;1H"; 
                std::cout << (Float_t)entry_num/(Float_t)total_entries*100 << "%" << std::endl;
                time_t time_left = (time(NULL)-start_time)*(float)(total_entries-entry_num)/(float)(entry_num);
                // std::cout << " time left: " << (time(NULL)-start_time)*(float)(total_entries-entry_num)/(float)(entry_num) << "s" <<std::endl;
                std::cout << " time left: ";
                if (time_left/3600 > 0) cout << time_left/3600 <<"h ";
                if ((time_left%3600)/60 > 0 || time_left/3600 == 0) cout << (time_left%3600)/60 << "m ";
                cout << (time_left%3600)%60<< "s " <<std::endl;

            }
        }
    }
    delete PMT_tree;
////////////////////////
    combined_tree->Write();
    combined_root->Close();
    return 0;
}