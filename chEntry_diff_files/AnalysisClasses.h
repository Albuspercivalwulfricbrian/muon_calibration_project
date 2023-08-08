//##################
//################
struct mini_tree_nrg
{
    Float_t EdepIntermediate0;
    Float_t EdepIntermediate1;
    Float_t EdepScat0;
    Float_t EdepScat1;
    Float_t EdepDet0;
    Float_t EdepDet1;  
    Short_t DetNum0;
    Short_t DetNum1;
    Short_t EventType; 
    // Int_t SetBranch(TTree *tree);
    Int_t Initialize();
};

struct mini_tree_time
{
    Float_t TimeIntermediate0;
    Float_t TimeIntermediate1;
    Float_t TimeScat0;
    Float_t TimeScat1;
    Float_t TimeDet0;
    Float_t TimeDet1;
    // Int_t SetBranch(TTree *tree);
    Int_t Initialize();
};

Int_t mini_tree_nrg::Initialize()
{
    EdepIntermediate0 = 0;
    EdepIntermediate1 = 0;

    EdepScat0 = 0;
    EdepScat1 = 0;
    EdepDet0 = 0;
    EdepDet1 = 0;  
    DetNum0 = 0;
    DetNum1 = 0; 
    EventType = -10;
    return 0;
}

Int_t mini_tree_time::Initialize()
{
    TimeScat0 = 0;
    TimeScat1 = 0;
    TimeDet0 = 0;
    TimeDet1 = 0;
    TimeIntermediate0 = 0;
    TimeIntermediate1 = 0;

    return 0;
}
