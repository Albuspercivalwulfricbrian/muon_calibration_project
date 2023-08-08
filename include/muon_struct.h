#ifndef HIFO
#define HIFO 1
#include "ChannelEntry.h"
class HitInfo
{
    public:
    Int_t X=-1;
    Int_t Y=-1;
    Int_t Z=-1;
    Float_t charge = 0;
    Int_t amp = 0;
    Float_t time = 0;
    Int_t channel = -100;
    Float_t zl_rms = 999;
    IntegralInfo II;

};

class TrackInfo
{
    public:
    std::vector<HitInfo> hits = {};

    void Reset()
    {
        hits = {};
        hits_reduced = {};
    }
    int GetCurrentTrackSize()
    {
        return hits.size();
    }
    void AddHit(const HitInfo& hit)
    {
        hits.push_back(hit);
    }

    std::vector<HitInfo> GetReduced()
    {
        return hits_reduced;
    }

    void SetReducedCharge(int numb, Float_t energy)
    {
        hits_reduced[numb].charge = energy;
    }
    void SearchTrueTrack(int threshold)
    {
        // hits.erase(
        //     std::remove_if(
        //         hits.begin(), 
        //         hits.end(),
        //         [&](const HitInfo& p) { 
        //             return (
        //                 // false
        //                 !(p.time > 1550 && p.time < 1580)
        //             || (p.zl_rms > 30) || (p.II.signal_length > 40) 
        //             || (p.II.end_amplitude < -15) || (p.charge < threshold)
        //             ); 
        //         }
        //     ), hits.end()); 

        for (auto p : hits)
        {
            if (    
                    (p.time > 1550 && p.time < 1580) && 
                    (p.zl_rms < 40) && 
                    (p.II.signal_length < 60) && 
                    (p.II.end_amplitude > -15)  && 
                    (p.charge > threshold)

                    ) 
                    hits_reduced.push_back(p);
        }
        if (hits_reduced.size()>0)    
        {
            for (int i = 0; i < hits_reduced.size()-1; i++)
            {
                int flag = 0;
                for (int j = i+1; j < hits_reduced.size(); j++)
                {
                    if ((abs(hits_reduced[j].X-hits_reduced[i].X)<=1) && (abs(hits_reduced[j].Y-hits_reduced[i].Y)<=1) && (abs(hits_reduced[j].Z-hits_reduced[i].Z)<=1)
                    //  && (2*abs(hits_reduced[j].charge-hits_reduced[i].charge)/abs(hits_reduced[j].charge+hits_reduced[i].charge)<0.3)
                     )
                    {
                        flag = 1;
                        break;
                    }    
                }
                if (flag == 0) {hits_reduced[i].charge = -1000;}
            }

            hits_reduced.erase(
                std::remove_if(
                    hits_reduced.begin(), 
                    hits_reduced.end(),
                    [](HitInfo const & p) { return p.charge == -1000; }
                ), hits_reduced.end()); 
            Float_t avgTime = 0;
            for (int i = 0; i < hits_reduced.size(); i++) avgTime+=(Float_t)hits_reduced[i].time/(Float_t)hits_reduced.size();

            // Float_t charge_slice[6];
            // for (int i = 0; i < 6; i++) charge_slice[i] = 0;
            // for (auto p : hits_reduced)
            // {
            //     charge_slice[p.Y-1] += p.charge;
            // }
            // bool track_selector = 1;
            // for (int i = 1; i < 6; i++)
            // {
            //     if (charge_slice[i-1]!= 0 && charge_slice[i]!=0 && 
            //     (2*abs(charge_slice[i-1]-charge_slice[i])/abs(charge_slice[i-1]+charge_slice[i])>0.3)) track_selector = 0;
            // }
            // if (track_selector==0) hits_reduced.clear(); 
        }  
    }

    bool isVertical()
    {
        int maxX = -100;
        int minX = 100;
        int maxZ = -100;
        int minZ = 100;
        int maxY = -100;
        int minY = 100;
        Float_t avX = 0;
        Float_t avY = 0;
        Float_t avZ = 0;
        Float_t totalCharge = 0;
        for (int i = 0; i < hits_reduced.size(); i++)
        {
            avX +=(Float_t)(hits_reduced[i].X)*(hits_reduced[i].charge);
            avY +=(Float_t)(hits_reduced[i].Y)*(hits_reduced[i].charge);
            avZ +=(Float_t)(hits_reduced[i].Z)*(hits_reduced[i].charge);
            totalCharge+=hits_reduced[i].charge;
            if (hits_reduced[i].X<minX) minX = hits_reduced[i].X;
            if (hits_reduced[i].X>maxX) maxX = hits_reduced[i].X;
            if (hits_reduced[i].Z<minZ) minZ = hits_reduced[i].Z;
            if (hits_reduced[i].Z>maxZ) maxZ = hits_reduced[i].Z;            
            if (hits_reduced[i].Y<minY) minY = hits_reduced[i].Y;
            if (hits_reduced[i].Y>maxY) maxY = hits_reduced[i].Y;
            
        }
        avX/=totalCharge;
        avY/=totalCharge;
        avZ/=totalCharge;

        // cout << avX << " " << avY << " " << avZ << endl;
        Float_t avAngle = 0;
        for (int i = 0; i < hits_reduced.size(); i++)
        {
            avAngle += pow((((Float_t)hits_reduced[i].X-avX)*((Float_t)hits_reduced[i].X-avX)+((Float_t)hits_reduced[i].Z-avZ)*((Float_t)hits_reduced[i].Z-avZ))/
            (((Float_t)hits_reduced[i].X-avX)*((Float_t)hits_reduced[i].X-avX)+((Float_t)hits_reduced[i].Z-avZ)*((Float_t)hits_reduced[i].Z-avZ)+((Float_t)hits_reduced[i].Y-avY)*((Float_t)hits_reduced[i].Y-avY))
            ,0.5)*(hits_reduced[i].charge);
        }
        avAngle/=totalCharge;
        // cout << avAngle << endl;

        ////////check vertically whole
        int flag = 0;
        // for (int pos = minY; pos <= maxY; pos++)
        // {
        //     int flag = 0;
        //     for (int i = 0; i < hits_reduced.size(); i++)
        //     { 
        //         if (hits_reduced[i].Y == pos) {flag = 1; break;}
        //     }       
        //     if (flag==0) return 0;
        // }
        ///////////////////////////

        // if ((maxX-minX <= 1) && (maxZ-minZ <=2) && (maxY-minY >=4 )) return 1;
        // else return 0;
        if (avAngle < 0.55 && avAngle >= 0.25) return 1;
        else return 0;
    }
    private:
    std::vector<HitInfo> hits_reduced = {};

};

#endif