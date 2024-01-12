#include <iostream>
#include <time.h>

#include <bitset>
#include <stdint.h>
#include <vector>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TAxis.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TRandom.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <vector>
#include <fstream>
#include <map>

using namespace std;


typedef struct MRTOFHit{
    Double_t daqtss;
    Double_t ts;
    Double_t tof1;
    Double_t tof2;
    Double_t e;
    Double_t vetot;
    std::vector<float> pulsestart1;
    std::vector<float> pulsestart2;
    std::vector<float> pulsestop;
    std::vector<float> pulsecfdstop;
    void Clear(){
        daqtss = -1;
        tof1 = 0;
        tof2 = 0;
        e = 0;
        vetot=-1;
        ts = 0;
        pulsestart1.clear();
        pulsestart2.clear();
        pulsestop.clear();
        pulsecfdstop.clear();
    }
    void Copy(MRTOFHit& obj){
        obj.daqtss = daqtss;
        obj.tof1 = tof1;
        obj.tof2 = tof2;
        obj.e = e;
        obj.vetot = vetot;
        obj.ts = ts;
        obj.pulsestart1 = pulsestart1;
        obj.pulsestart2 = pulsestart2;
        obj.pulsestop = pulsestop;
        obj.pulsecfdstop = pulsecfdstop;
    }
}MRTOFHit;


MRTOFHit * treedata;


void mcs_make_rootfilev1(char* inputfile, char* outputfile){
    //////////////////////////////////////////////////////////
    //   This file has been automatically generated
    //     (Thu Apr 27 18:26:49 2023 by ROOT version6.22/08)
    //   from TTree tree/tree
    //   found on file: mcs_run249_Tc105.root
    //////////////////////////////////////////////////////////


    //Reset ROOT and connect tree file
    gROOT->Reset();

    //! read start time of ch 0
    TString logfilename = TString(inputfile)+TString(".log");
    ifstream ifs(logfilename.Data());
    Double_t StartTimeSecCh[10];
    memset(StartTimeSecCh,0,sizeof(StartTimeSecCh));
    while (!ifs.eof()){
        std::string tmpstr[10];
        for (Int_t i=0;i<10;i++){
            ifs>>tmpstr[i];
        }
        Int_t ch = atoi(tmpstr[1].c_str());
        Double_t startTime = atof(tmpstr[9].c_str());
        StartTimeSecCh[ch] = startTime;
    }
    for (Int_t i=0;i<10;i++){
        if (StartTimeSecCh[i]>0)
            cout<<"Ch="<<i<<"\tTimeSinceEpoch="<<std::setprecision(15)<<StartTimeSecCh[i]<<endl;
    }


    TTree* tree = NULL;
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputfile);
    if (!f) {
        f = new TFile(inputfile);
    }
    f->GetObject("tree",tree);

    //Declaration of leaves types
    Int_t           ch;
    Int_t           edge;
    Int_t           t;
    Int_t           sweeps;
    Int_t           tag;
    Int_t           data_lost;

    // Set branch addresses.
    tree->SetBranchAddress("ch",&ch);
    tree->SetBranchAddress("edge",&edge);
    tree->SetBranchAddress("t",&t);
    tree->SetBranchAddress("sweeps",&sweeps);
    tree->SetBranchAddress("tag",&tag);
    tree->SetBranchAddress("data_lost",&data_lost);

    Long64_t nentries = tree->GetEntries();

    TH1F* htof = new TH1F("h1","h1",6000,14677730.,14685730.);
    TH1F* htof_ref = new TH1F("h2","h2",6000,12640136.,12665136.);

    Int_t maxsweep = 0;
    Int_t minsweep = 10000000;

//    nentries = 1000000;
    TH2F* hsweep = new TH2F("hsweep","hsweep",200,0,nentries,200,0,65535);
    TH2F* hsweepL = new TH2F("hsweepL","hsweepL",200,0,nentries,200,0,10145817);
//    TH2F* htofvssweep = new TH2F("htofvssweep","htofvssweep",200,0,50000,5000,12640136.,12665136.);
    TH2F* htofvssweep = new TH2F("htofvssweep","htofvssweep",2000,0,150000,300,12.6509e6,12.6513e6);
    TH2F* htofvssweep_corr = new TH2F("htofvssweep_corr","htofvssweep_corr",2000,0,150000,300,12.6509e6,12.6513e6);

    Long64_t nbytes = 0;

    ULong64_t sweepL = 0;
    Int_t prevSweeps = -1;
    ULong64_t prevSweepsL = 0;
    ULong64_t cycle = 0;
    ULong64_t cycleStartNo = 0;

    TFile* file0 = new TFile(outputfile,"recreate");

    treedata = new MRTOFHit;
//    treedata = new MRTOFHit;
    TTree* treeout = new TTree("tree","tree");
    treeout->Branch("b",treedata);

    long long driffcorr_tstart = 0;
    Double_t tofmeanVal_start;
    long long mrtof_TOFcorrTimeLength = 217;
    Double_t mrtof_TOFcorrwindow_window[] = {12.6508e6,12.6513e6};


    int driffcorr_n = 0;

    Long64_t tofentry = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nbytes += tree->GetEntry(jentry);

        if (prevSweeps!=-1&&prevSweeps>sweeps){
            cycleStartNo+=65536;
            cycle++;
        }
        sweepL = cycleStartNo+sweeps;
        hsweepL->Fill(jentry,sweepL);
        if (prevSweepsL>sweepL)
            cout<<"something wrong "<<jentry<<"\t"<<prevSweepsL<<"\t"<<sweepL<<"\t"<<cycleStartNo<<endl;
        prevSweepsL  = sweepL;
        prevSweeps  = sweeps;
        if (sweeps<minsweep)
            minsweep = sweeps;
        if (sweeps>maxsweep)
            maxsweep = sweeps;
        hsweep->Fill(jentry,sweeps);

        if (ch==1&&tag==1)
        htofvssweep->Fill(sweepL,(Double_t)t/10.);
        if (ch==1&&tag==0)
        htof->Fill((Double_t)t/10);
        if (ch!=1&&(tag!=0||tag!=1))
            continue;
        tofentry ++;
//        cout<<jentry<<endl;
        treedata->Clear();
        if (tag==1)//ref ions
            treedata->vetot = 1;
        else
            treedata->vetot = 0;
        treedata->daqtss = StartTimeSecCh[1];
        treedata->ts=sweepL;
        treedata->tof2 = t/10.;

        if (treeout) treeout->Fill();
    }
    cout<<treeout->GetEntries()<<"\t"<<tofentry<<endl;
    cout<<sweepL<<endl;
    cout<<cycle<<endl;
    cout<<minsweep<<"\t"<<maxsweep<<endl;
    TCanvas * c1 =new TCanvas("c1","c1",900,700);
    c1->Divide(1,2);
    c1->cd(1);
    htofvssweep->Draw("colz");
    c1->cd(2);
    htofvssweep_corr->Draw("colz");
    treeout->Write();
    file0->Close();
}
