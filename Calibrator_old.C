#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TString.h>
#include <TGraph2D.h>
#include "ChannelEntry.h"
#include "CHSH_class.h"
#include <string.h>
// #include "Double_scattering_cuts.h"
#include "like_ivashkin_wants_it.h"
#include "muon_struct.h"
#include "Map.h"
#include "MultiCoreAnalysis.h"

using namespace std;
using namespace MAP_H;

int diff(std::vector<int> in)
{
    int max = -100;
    int min = 100;
    for (int i = 0; i< in.size(); i++) 
    {
        if (max < in[i]) max = in[i];
        if (min > in[i]) min = in[i];

    }
    return max-min;
}
int devergeance(std::vector<int> in)
{
    int max = -100;
    int min = 100;
    for (int i = 0; i< in.size(); i++) 
    {
        if (max < in[i]) max = in[i];
        if (min > in[i]) min = in[i];

    }
    return (max-min)/(max+min);
}

int main(int argc, char** argv)
{


    if(argc == 3) {
        argv[1];
        argv[2];
        // argv[4];
    }
    CreateMap();
    cout << Map[93][2] << endl;

	TString source_path = (TString)argv[1];
    TString run_name = (TString)argv[2];

    const Int_t total_channels = 128;

    TFile *source_file = TFile::Open("/media/doc/DATA/muon_calibration/new_data/tr_tracks.root");
    TTree *source_tree = (TTree*)source_file->Get("Tracks");
    //short_tree->SetBranch(MiniDST_tree);

    TrackInfo *trackinfo = new TrackInfo;

    source_tree->SetBranchAddress("TrackInfo", &trackinfo);

    TCanvas *canvas[8];
    for (int i = 0; i < 8; i++)
    {
        canvas[i] = new TCanvas(Form("c_%i",i),Form("c_%i",i));
        canvas[i]->Divide(4,4);
    }

    TH1F *histall[128];
    for (int i = 1; i <= 128; i++)
    {
        histall[i-1] = new TH1F(Form("ch_%i_x=%i_y=%i_z=%i",i,Map[i][0],Map[i][1],Map[i][2]),Form("ch_%i_x=%i_y=%i_z=%i",i,Map[i][0],Map[i][1],Map[i][2]),70,000,40000);
    }
    Int_t total_entries = source_tree->GetEntries();

    for (int EvNum = 0; EvNum < source_tree->GetEntries(); EvNum++)
    {
        trackinfo->Reset();
        source_tree->GetEntry(EvNum);
        int track_size = trackinfo->GetCurrentTrackSize();
        std::array<std::array<std::vector<int>,8>,4> hitscounter; std::array<std::array<std::vector<int>,8>,4> hitsnumber;
        // std::array<std::vector<int>,4> hitscounter; std::array<std::vector<int>,4> hitsnumber;
        // std::array<std::array<std::vector<int>,8>,4> hitscounter;
        TGraph2D *gr = new TGraph2D();
        vector<Float_t> chgr;
        trackinfo->SearchTrueTrack(3000);
        for (int hitnum = 0; hitnum < track_size; hitnum++)
        {
            if (trackinfo->hits[hitnum].charge > 3000)
            {
                hitscounter[trackinfo->hits[hitnum].X][trackinfo->hits[hitnum].Z].push_back(trackinfo->hits[hitnum].Y); hitsnumber[trackinfo->hits[hitnum].X][trackinfo->hits[hitnum].Z].push_back(hitnum);
                // hitscounter[trackinfo->hits[hitnum].X].push_back(trackinfo->hits[hitnum].Y); hitsnumber[trackinfo->hits[hitnum].X].push_back(hitnum);
                // hitscounter[trackinfo->hits[hitnum].X][trackinfo->hits[hitnum].Y].push_back(trackinfo->hits[hitnum].Z);
            }
        }

        for (int x = 1; x < 4;x++) 
        {
            for (int z = 1; z < 8; z++) 
            {
                auto k = diff(hitscounter[x][z]);//-*min_element(hitscounter[x][z].begin(),hitscounter[x][z].end());
                if (k>=3 && hitscounter[x][z].size()>=4 && devergeance(hitscounter[x][z]) < 0.1) 
                {

                    // cout << k << endl;
                    // cout << "X  = " << x << "; Z = " << z << "; hits = " << hitscounter[x][z].size() <<endl;
                    for (int y = 0; y< hitscounter[x][z].size(); y++)
                    {
                        // cout << hitscounter[x][z][y];
                        Int_t hitnumber = hitsnumber[x][z][y];
                        histall[trackinfo->hits[hitnumber].channel-1]->Fill(trackinfo->hits[hitnumber].charge);

                    } 
                    // cout << endl;
                }                   
            }


            // for (int y = 1; y < 5; y++) 
            // {
                
            //     // cout << "shit" << " x= " << x << " y= " << y << endl;
            //     auto k = diff(hitscounter[x][y]);//-*min_element(hitscounter[x][z].begin(),hitscounter[x][z].end());
            //     // cout << "shit" << " x= " << x << " y= " << y << endl;
            //     // int k = 2;

            //     if (k>=5 && hitscounter[x][y].size() >=5) 
            //     {
            //         cout << k << endl;
            //         cout << "X  = " << x << "; Y = " << y << "; hits = " << hitscounter[x][y].size() <<endl;
            //         for (int z = 0; z< hitscounter[x][y].size(); z++) cout << hitscounter[x][y][z];
            //         cout << endl;
            //     }   
            // }                

                // int sk = 5;
                // auto k = diff(hitscounter[x]);//-*min_element(hitscounter[x][z].begin(),hitscounter[x][z].end());
                // if (k>=sk && hitscounter[x].size() >= sk) 
                // {

                //     for (int y = 0; y< hitscounter[x].size(); y++)
                //     {
                //         // cout << hitscounter[x][z][y];
                //         Int_t hitnumber = hitsnumber[x][y];
                //         histall[trackinfo->hits[hitnumber].channel-1]->Fill(trackinfo->hits[hitnumber].charge);

                //     } 


                // }                   
        }

        
        // for (int ch = 0; ch < total_channels; ch++)
        // {
        //     if (short_channel_info[ch]->amp > 150 && short_channel_info[ch]->amp < 60000
        //         //  && short_channel_info[ch]->charge > 1000 
        //         // && short_channel_info[ch]->time > 1500 && short_channel_info[ch]->time < 1600
        //         && ch !=63 && ch!=127
        //         )
        //     {
        //         HitInfo hitinfo;
        //         hitinfo.charge = short_channel_info[ch]->charge;
        //         hitinfo.time = short_channel_info[ch]->time;
        //         hitinfo.channel = ch+1;
        //         hitinfo.X = Map[ch+1][0];
        //         hitinfo.Y = Map[ch+1][1];
        //         hitinfo.Z = Map[ch+1][2];

        //         trackinfo.AddHit(hitinfo);                
        //     }
        // }

            if (EvNum%10000 == 0) {std::cout<< u8"\033[2J\033[1;1H"; std::cout << (Float_t)EvNum/(Float_t)total_entries*100 << "%" << std::endl;}

    }
    for (int i = 0; i < 128; i++)
    {
        canvas[i/16]->cd(i%16+1);
        histall[i]->Draw();
    }

    TFile *result_root = new TFile ("/media/doc/DATA/muon_calibration/new_data/tracks.root", "RECREATE");

    for (int i = 0; i < 8; i++)
    {
        canvas[i]->Write();
    }
    result_root->Close();
}