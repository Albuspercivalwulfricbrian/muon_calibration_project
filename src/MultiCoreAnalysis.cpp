#include "MultiCoreAnalysis.h"
// #include "Map.h"
MultiCoreAnalysis::MultiCoreAnalysis(ChannelEntry *Echannel_info, short_energy_ChannelEntry* Eshort_channel_info, int Cores, int Channels) : fCores(Cores),fChannels(Channels)
{

    fCores = Cores;
    fChannels = Channels;
    channel_info=Echannel_info;
    short_channel_info=Eshort_channel_info;

    for (int i = 0; i < 128; i++) {
        if (i!=63 || i!=127) ChannelsToCalculate.push_back(i);
    }
}

MultiCoreAnalysis::~MultiCoreAnalysis()
{

}

void MultiCoreAnalysis::ExecuteThread() {
    int channel = GetTaskSafely();
    while(channel >= 0) {
        Analyze1Thread(channel);
        channel = GetTaskSafely();
    }
}

int MultiCoreAnalysis::GetTaskSafely() {
    const std::lock_guard<std::mutex> lock(ChannelsToCalculate_mutex);
    if (ChannelsToCalculate.size() > 0) {
        // cout << ChannelsToCalculate.size() << endl;
        int to_return = ChannelsToCalculate.back();
        ChannelsToCalculate.pop_back();
        return to_return;
    }
    return -1;
}

void MultiCoreAnalysis::Analyze1Thread(Int_t channel_number)
{

    if (channel_info[channel_number].wf_size > 0)
    {

        channel_info[channel_number].SplineWf();
        short_channel_info[channel_number].Initialize();
        channel_info[channel_number].Set_Zero_Level_Area(20);
        Int_t zero_level = channel_info[channel_number].Get_Zero_Level();
        short_channel_info[channel_number].zl_rms = channel_info[channel_number].Get_Zero_Level_RMS();
        //channel_info[channel_number].SplineWf();
        Int_t pp = channel_info[channel_number].Get_time();
        channel_info[channel_number].AssumeSmartScope();
        short_channel_info[channel_number].amp = channel_info[channel_number].Get_Amplitude();
        short_channel_info[channel_number].charge = channel_info[channel_number].Get_Charge();
        short_channel_info[channel_number].time = channel_info[channel_number].Get_time_gauss();
        short_channel_info[channel_number].zl = zero_level;
        short_channel_info[channel_number].II = channel_info[channel_number].GetIntegralInfo();
    }
}

void MultiCoreAnalysis::AnalyzeMultiThread()
{
    // int ch4thread = fChannels/fCores;
    // if (fChannels%fCores!=0) fch4thread++;
    // for (int i = 0; i < 128;i++) std::cout << channel_info[i].wf_size << endl;
    // for (int i = 0; i < fCores; i++) MCthreads[i].join();
    for (int i = 0; i < fCores; i++) MCthreads.push_back(std::thread(&MultiCoreAnalysis::ExecuteThread,this));
    for (int i = 0; i < fCores; i++) MCthreads[i].join();
    // MCthreads.clear();

}

