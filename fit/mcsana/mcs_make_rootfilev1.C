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
#include <vector>
#include <fstream>
#include <map>

using namespace std;

typedef struct MRTOFHit{
    Long64_t sw;
    Double_t tof;
    Double_t tofcorr[3];
    Double_t tag;
    void Clear(){
        tof = 0;
        tag=0;
        sw = 0;
        tofcorr[0] = 0;
        tofcorr[1] = 0;
        tofcorr[2] = 0;
        
    }
    void Copy(MRTOFHit& obj){
        obj.tof = tof;
        obj.tofcorr[0] = tofcorr[0];
        obj.tofcorr[1] = tofcorr[1];
        obj.tofcorr[2] = tofcorr[2];
        obj.tag = tag;
        obj.sw = sw;
    }
}MRTOFHit;

MRTOFHit * treedata;
MRTOFHit* treedata_drift_corr;

std::multimap<ULong64_t,MRTOFHit*> datamap_driffcorr;

void mcs_make_rootfilev1(char* inputfile, char* outputfile){
    //Reset ROOT and connect tree file
    gROOT->Reset();
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


    Int_t maxsweep = 0;
    Int_t minsweep = 10000000;

    Long64_t nbytes = 0;

    ULong64_t sweepL = 0;
    Int_t prevSweeps = -1;
    ULong64_t prevSweepsL = 0;
    ULong64_t cycle = 0;
    ULong64_t cycleStartNo = 0;

    TFile* file0 = new TFile(outputfile,"recreate");

    treedata = new MRTOFHit;
    treedata_drift_corr = new MRTOFHit;
    TTree* tree_drift_corr = new TTree("tree","tree");
    tree_drift_corr->Branch("b",treedata_drift_corr);

    long long driffcorr_tstart = 0;
    Double_t tofmeanVal_start;
    long long mrtof_TOFcorrTimeLength = 217;
    Double_t mrtof_TOFcorrwindow_window[] = {12.6508e6,12.6513e6};
    ifstream parmsfile("parms.txt");
    parmsfile>>mrtof_TOFcorrTimeLength;
    parmsfile>>mrtof_TOFcorrwindow_window[0];
    parmsfile>>mrtof_TOFcorrwindow_window[1];
    cout<<"reading parms.txt: number_of_sweep="<<mrtof_TOFcorrTimeLength<<"\tTOF range: "
    <<mrtof_TOFcorrwindow_window[0]<<" - "<<mrtof_TOFcorrwindow_window[1]<<" ns"<<endl;
    

    int driffcorr_n = 0;

    Long64_t tofentry = 0;

    double medianVal = 0;
    double rms = 0;
    // cout<<"AAAAA"<<nentries<<endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nbytes += tree->GetEntry(jentry);

        if (prevSweeps!=-1&&prevSweeps>sweeps){
            cycleStartNo+=65536;
            cycle++;
        }
        sweepL = cycleStartNo+sweeps;
        if (prevSweepsL>sweepL)
            cout<<"something wrong "<<jentry<<"\t"<<prevSweepsL<<"\t"<<sweepL<<"\t"<<cycleStartNo<<endl;
        prevSweepsL  = sweepL;
        prevSweeps  = sweeps;
        if (sweeps<minsweep)
            minsweep = sweeps;
        if (sweeps>maxsweep)
            maxsweep = sweeps;
        if (ch!=1&&(tag!=0||tag!=1))
            continue;
        tofentry ++;
//        cout<<jentry<<endl;
        treedata->Clear();
        if (tag==1)//ref ions
            treedata->tag= 1;
        else
            treedata->tag= 0;
        treedata->sw=sweepL;
        treedata->tof = t/10.;

        //! DRIFT Correction here!
        MRTOFHit* hitTOF = new MRTOFHit;
        treedata->Copy(*hitTOF);

        if (driffcorr_tstart==0)
            driffcorr_tstart = treedata->sw;
        if (treedata->sw - driffcorr_tstart<mrtof_TOFcorrTimeLength){
            datamap_driffcorr.emplace(treedata->sw,hitTOF);
        }else if (treedata->sw - driffcorr_tstart>=mrtof_TOFcorrTimeLength) {
            driffcorr_tstart = treedata->sw;
            Int_t nevtforcorr = 0;
            std::vector<double> arrayOfTOF;
            arrayOfTOF.clear();
            for (auto itdriffcorr =datamap_driffcorr.begin() ; itdriffcorr != datamap_driffcorr.end(); ++itdriffcorr){
                MRTOFHit* myhitTOF = (MRTOFHit*) itdriffcorr->second;
                if (myhitTOF->tag!=0){
                    if (myhitTOF->tof>mrtof_TOFcorrwindow_window[0]&&myhitTOF->tof<mrtof_TOFcorrwindow_window[1]){
                        arrayOfTOF.push_back(myhitTOF->tof);
                        nevtforcorr++;
                    }
                }
            }
            medianVal = TMath::Median(arrayOfTOF.size(),&arrayOfTOF[0]);
            rms = TMath::RMS(arrayOfTOF.size(),&arrayOfTOF[0]);
            if (driffcorr_n==0){
                tofmeanVal_start = medianVal;
                cout<<" data size = "<<datamap_driffcorr.size()<<"\tNevents for correction = "<<nevtforcorr<<"\tMedian Val = "<<medianVal<<"\tRMS = "<<rms<<endl;
            }
            for (auto itdriffcorr =datamap_driffcorr.begin() ; itdriffcorr != datamap_driffcorr.end(); ++itdriffcorr){
                MRTOFHit* myhitTOF = (MRTOFHit*) itdriffcorr->second;
                myhitTOF->Copy(*treedata_drift_corr);
                treedata_drift_corr->tofcorr[0] = medianVal;
                treedata_drift_corr->tofcorr[1] = tofmeanVal_start;
                treedata_drift_corr->tofcorr[2] = rms;
                if (tree_drift_corr) tree_drift_corr->Fill();
            }

            driffcorr_n++;
            for (auto itdriffcorr =datamap_driffcorr.begin() ; itdriffcorr != datamap_driffcorr.end(); ++itdriffcorr){
                delete itdriffcorr->second;
            }
            datamap_driffcorr.clear();
            datamap_driffcorr.emplace(treedata->sw,hitTOF);
        }

        
    }
    for (auto itdriffcorr =datamap_driffcorr.begin() ; itdriffcorr != datamap_driffcorr.end(); ++itdriffcorr){
            MRTOFHit* myhitTOF = (MRTOFHit*) itdriffcorr->second;
            myhitTOF->Copy(*treedata_drift_corr);
            treedata_drift_corr->tofcorr[0] = medianVal;
            treedata_drift_corr->tofcorr[1] = tofmeanVal_start;
            treedata_drift_corr->tofcorr[2] = rms;
            if (tree_drift_corr) tree_drift_corr->Fill();
    }
    cout<<tree_drift_corr->GetEntries()<<"\t"<<tofentry<<endl;
    tree_drift_corr->Write();
    file0->Close();
}