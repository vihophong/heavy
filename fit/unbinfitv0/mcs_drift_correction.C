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
    Long64_t ts;
    Double_t tof1;
    Double_t tof2;
    Double_t e;
    Double_t vetot;
    std::vector<float> pulsestart1;
    std::vector<float> pulsestart2;
    std::vector<float> pulsestop;
    std::vector<float> pulsecfdstop;
    void Clear(){
        tof1 = 0;
        tof2 = 0;
        e = 0;
        vetot=0;
        ts = 0;
        pulsestart1.clear();
        pulsestart2.clear();
        pulsestop.clear();
        pulsecfdstop.clear();
    }
    void Copy(MRTOFHit& obj){
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
MRTOFHit* treedata_drift_corr;

std::multimap<ULong64_t,MRTOFHit*> datamap_driffcorr;

void mcs_drift_correction(){
    //////////////////////////////////////////////////////////
    //   This file has been automatically generated
    //     (Thu Apr 27 18:26:49 2023 by ROOT version6.22/08)
    //   from TTree tree/tree
    //   found on file: mcs_run249_Tc105.root
    //////////////////////////////////////////////////////////


    //Reset ROOT and connect tree file
    gROOT->Reset();
    TTree* tree = NULL;
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcs_run249_Tc105.root");
    if (!f) {
        f = new TFile("mcs_run249_Tc105.root");
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

    TFile* file0 = new TFile("mcs_drift_corr_tree_refreal.root","recreate");

    treedata = new MRTOFHit;
    treedata_drift_corr = new MRTOFHit;
    TTree* tree_drift_corr = new TTree("tree_drift_corr","tree_drift_corr");
    tree_drift_corr->Branch("b",treedata_drift_corr);

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
        treedata->ts=sweepL;
        treedata->tof2 = t/10.;


        //! DRIFT Correction here!
        MRTOFHit* hitTOF = new MRTOFHit;
        treedata->Copy(*hitTOF);

        if (driffcorr_tstart==0){
            driffcorr_tstart = treedata->ts;
        }else{
            if (treedata->ts - driffcorr_tstart<mrtof_TOFcorrTimeLength){
                datamap_driffcorr.emplace(treedata->ts,hitTOF);
            }else{
                driffcorr_tstart = treedata->ts;
                Int_t nevtforcorr = 0;
                std::vector<double> arrayOfTOF;
                arrayOfTOF.clear();
                for (auto itdriffcorr =datamap_driffcorr.begin() ; itdriffcorr != datamap_driffcorr.end(); ++itdriffcorr){
                    MRTOFHit* myhitTOF = (MRTOFHit*) itdriffcorr->second;
                    if (myhitTOF->vetot!=0){
                        if (myhitTOF->tof2>mrtof_TOFcorrwindow_window[0]&&myhitTOF->tof2<mrtof_TOFcorrwindow_window[1]){
                            arrayOfTOF.push_back(myhitTOF->tof2);
                            nevtforcorr++;
                        }
                    }
                }
                double medianVal = TMath::GeomMean(arrayOfTOF.size(),&arrayOfTOF[0]);
                double rms = TMath::RMS(arrayOfTOF.size(),&arrayOfTOF[0]);
//                            double medianVal = TMath::Median(arrayOfTOF.size(),&arrayOfTOF[0]);
                if (driffcorr_n==0){
                    tofmeanVal_start = medianVal;
                    cout<<" data size = "<<datamap_driffcorr.size()<<"\tNevents for correction = "<<nevtforcorr<<"\tMedian Val = "<<medianVal<<"\tRMS = "<<rms<<endl;
                }
                double correctionFactor = tofmeanVal_start/medianVal;

//                cout<<std::setprecision(10)<<(Double_t)jentry/(Double_t)nentries*100<<"\t"<<treedata->ts<<"\t"<<medianVal<<"\t"<<tofmeanVal_start<<"\t"<<correctionFactor<<endl;

                for (auto itdriffcorr =datamap_driffcorr.begin() ; itdriffcorr != datamap_driffcorr.end(); ++itdriffcorr){
                    MRTOFHit* myhitTOF = (MRTOFHit*) itdriffcorr->second;
                    double tof_corrected = myhitTOF->tof2*correctionFactor;
//                    cout<<std::setprecision(10)<<(Double_t)jentry/(Double_t)nentries*100<<"\t"<<tof_corrected<<"\t"<<myhitTOF->tof2<<"\t"<<tofmeanVal_start<<"\t"<<correctionFactor<<endl;

//                    htof_corr[0]->Fill(tof_corrected);
//                    htof2d_corr[0]->Fill(((Double_t)myhitTOF->ts-(Double_t)tsbeg_tof)/1e9,tof_corrected);
//                    if (myhitTOF->vetot==0){
//                        htof_corr[2]->Fill(tof_corrected);
//                        htof2d_corr[2]->Fill(((Double_t)myhitTOF->ts-(Double_t)tsbeg_tof)/1e9,tof_corrected);
//                    }else{
//                        htof_corr[1]->Fill(tof_corrected);
//                        htof2d_corr[1]->Fill(((Double_t)myhitTOF->ts-(Double_t)tsbeg_tof)/1e9,tof_corrected);
//                    }


                    myhitTOF->Copy(*treedata_drift_corr);
                    treedata_drift_corr->tof1 = tof_corrected;

                    htofvssweep_corr->Fill(treedata_drift_corr->ts,tof_corrected);
                    if ((treedata_drift_corr->tof2>12630000&&treedata_drift_corr->tof2<12680000)||
                            (treedata_drift_corr->tof2>14674000&&treedata_drift_corr->tof2<14690000))
                                if (tree_drift_corr) tree_drift_corr->Fill();
                }


                driffcorr_n++;
                for (auto itdriffcorr =datamap_driffcorr.begin() ; itdriffcorr != datamap_driffcorr.end(); ++itdriffcorr){
                    delete itdriffcorr->second;
                }
                datamap_driffcorr.clear();
                datamap_driffcorr.emplace(treedata->ts,hitTOF);
            }
        }
    }
    cout<<tree_drift_corr->GetEntries()<<"\t"<<tofentry<<endl;
    cout<<sweepL<<endl;
    cout<<cycle<<endl;
    cout<<minsweep<<"\t"<<maxsweep<<endl;
//    htof->Draw();
//    hsweep->Draw();
    TCanvas * c1 =new TCanvas("c1","c1",900,700);
    c1->Divide(1,2);
    c1->cd(1);
    htofvssweep->Draw("colz");
    c1->cd(2);
    htofvssweep_corr->Draw("colz");

    tree_drift_corr->Write();
    file0->Close();
//    htofvssweep->SaveAs("htofcounts_mcs.root");
}
