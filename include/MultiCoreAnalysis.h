#include "ChannelEntry.h"
#include <mutex>
#include <thread>
#include <vector>
#include <array>
#include <iostream>
using namespace std;



class MultiCoreAnalysis
{
    private:
    int fCores = 8;
    int fChannels = 128;
    int ch4thread = 0;
    ChannelEntry *channel_info = nullptr;
    short_energy_ChannelEntry *short_channel_info = nullptr;
    // std::array<std::array<short_energy_ChannelEntry,fChannels>,fCores> short_channel_info = {};
    std::vector<int> ChannelsToCalculate;
    std::mutex ChannelsToCalculate_mutex;
    std::vector<std::thread> MCthreads;

    void Analyze1Thread(Int_t Channel_number);
    void ExecuteThread();
    int GetTaskSafely();

    public:
    MultiCoreAnalysis(ChannelEntry *Echannel_info, short_energy_ChannelEntry *Eshort_channel_info,int Cores, int Channels);
    ~MultiCoreAnalysis();
    void AnalyzeMultiThread();
    // ChannelEntry channel_info[fChannels];
    // short_energy_ChannelEntry short_channel_info[fChannels];



};