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

#include <RooPolynomial.h>
#include <RooDataSet.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>

#include <vector>
#include <fstream>
#include "fitF.hh"
using namespace std;


using namespace RooFit;

class fitF;

typedef struct parms_t
{
    Double_t val[MAX_N_PARMS+MAX_NPEAKS];
    Double_t min[MAX_N_PARMS+MAX_NPEAKS];
    Double_t max[MAX_N_PARMS+MAX_NPEAKS];
    Bool_t isFix[MAX_N_PARMS+MAX_NPEAKS];
    Int_t index[MAX_N_PARMS+MAX_NPEAKS];
    Int_t category[MAX_N_PARMS+MAX_NPEAKS];
    Int_t group[MAX_N_PARMS+MAX_NPEAKS];
    Int_t peakNo[MAX_N_PARMS+MAX_NPEAKS];
    int print(Int_t nparms){
        for (Int_t i=0;i<nparms;i++){
            std::cout<<"****"<<i<<std::endl;
             std::cout<<"PARMS = "<<i<<std::endl;
             std::cout<<"category = "<<category[i]<<std::endl;
             std::cout<<"index = "<<index[i]<<std::endl;
             std::cout<<"isFix = "<<isFix[i]<<std::endl;
             std::cout<<"max = "<<max[i]<<std::endl;
             std::cout<<"min = "<<min[i]<<std::endl;
             std::cout<<"val = "<<val[i]<<std::endl;
             std::cout<<"peakNo = "<<peakNo[i]<<std::endl;
             std::cout<<"group = "<<group[i]<<std::endl;
        }
        return 0;
    }
    int reset(){
        for (Int_t i=0;i<MAX_N_PARMS+MAX_NPEAKS;i++){
            val[i] = 10;
            min[i] = 0;
            max[i] = 20;
            isFix[i] = false;
            index[i] = -1;
            category[i] = -2;
            group[i] = -2;
            peakNo[i] = -1;
        }
        return 0;
    }
}parms_t;


typedef struct{
    std::string name;
    int id;
    int tag;
    double peakpos;
    bool is_common_shape_parm;
    double shape1;
    double shape2;
    double shape3;
    double lowerbound;
    double upperbound;
    bool isbkgfix;
    double bkgslope;
    double snratio;
    int group;
    int print(){
        std::cout<<"roi = "<<std::endl;
        std::cout<<"id = "<<id<<std::endl;
        std::cout<<"name = "<<name<<std::endl;
        std::cout<<"tag = "<<tag<<std::endl;
        std::cout<<"peakpos = "<<peakpos<<std::endl;
        std::cout<<"is_common_shape_parm = "<<is_common_shape_parm<<std::endl;
        std::cout<<"shape1 = "<<shape1<<std::endl;
        std::cout<<"shape2 = "<<shape2<<std::endl;
        std::cout<<"shape3 = "<<shape3<<std::endl;
        std::cout<<"lowerbound = "<<lowerbound<<std::endl;
        std::cout<<"upperbound = "<<upperbound<<std::endl;
        std::cout<<"isbkgfix = "<<isbkgfix<<std::endl;
        std::cout<<"bkgslope = "<<bkgslope<<std::endl;
        std::cout<<"snratio = "<<snratio<<std::endl;
        std::cout<<"group = "<<group<<std::endl;
        std::cout<<"--------------"<<std::endl;
        return 0;
    }
}roi_t;

typedef struct treedata_t
{
    UInt_t          fUniqueID;
    UInt_t          fBits;
    Long64_t        ts;
    Double_t        tof1;
    Double_t        tof2;
    Double_t        e;
    Double_t        vetot;
    vector<float>   pulsestart1;
    vector<float>   pulsestart2;
    vector<float>   pulsestop;
    vector<float>   pulsecfdstop;
}treedata_t;

vector<roi_t*> ROIs;

void calculateRange(Double_t &rangelow, Double_t &rangehi, roi_t* roi, Double_t ratio_to_peak=5e-3){
    TF1* f1 = new TF1("fcn1","1./sqrt(1+(x-[0])*(x-[0])/[1]/[1]) * exp(-0.5*([2]+[3]*TMath::ASinH((x-[0])/[1])) *([2]+[3]*TMath::ASinH((x-[0])/[1])))",-roi->shape1*100,roi->shape1*100);
    f1->SetParameter(0,0);
    f1->SetParameter(1,roi->shape1);
    f1->SetParameter(2,roi->shape2);
    f1->SetParameter(3,roi->shape3);
    Double_t maxPos = f1->GetMaximumX(-roi->shape1*100,roi->shape1*100);
    Double_t maxVal = f1->GetMaximum(-roi->shape1*100,roi->shape1*100);

    Double_t pp = maxPos;
    Int_t nsteps = 10000;
    Double_t step = roi->shape1*100/nsteps;
    for (Int_t i=0;i<nsteps;i++){
        pp+=step;
        if(f1->Eval(pp)<maxVal*ratio_to_peak){
            rangehi = pp;
            break;
        }
    }
    pp = maxPos;
    for (Int_t i=0;i<nsteps;i++){
        pp-=step;
        if(f1->Eval(pp)<maxVal*ratio_to_peak){
            rangelow = pp;
            break;
        }
    }
}

void getROI(char* infile,Bool_t iscalculateRange=false)
{

    vector<vector<string>> content;
    vector<string> row;
    string line, word;

    fstream file (infile, ios::in);
    if(file.is_open())
    {
        while(getline(file, line))
        {
            row.clear();

            stringstream str(line);

            while(getline(str, word, ','))
                row.push_back(word);
            content.push_back(row);
        }
    }
    else
        cout<<"Could not open the file\n";

    for(int i=0;i<content.size();i++)//row
    {
        if (i>0){
            roi_t* roi = new roi_t;
            roi->id = atoi(content[i][0].c_str());
            roi->name = content[i][1];
            roi->tag = atoi(content[i][2].c_str());
            roi->peakpos = atof(content[i][3].c_str());
            int tmp = atoi(content[i][4].c_str());
            roi->group = tmp;
            if (tmp==0)
                roi->is_common_shape_parm = false;
            else
                roi->is_common_shape_parm = true;

            roi->shape1 = atof(content[i][5].c_str());
            roi->shape2 = atof(content[i][6].c_str());
            roi->shape3 = atof(content[i][7].c_str());
            roi->lowerbound = atof(content[i][8].c_str());
            roi->upperbound = atof(content[i][9].c_str());
            tmp = atoi(content[i][10].c_str());
            if (tmp==0)
                roi->isbkgfix = false;
            else
                roi->isbkgfix = true;
            roi->bkgslope = atof(content[i][11].c_str());
            roi->snratio = atof(content[i][11].c_str());
            roi->upperbound = atof(content[i][12].c_str());

            if (iscalculateRange){
                calculateRange(roi->lowerbound,roi->upperbound,roi);
                roi->lowerbound = roi->lowerbound+roi->peakpos;
                roi->upperbound = roi->upperbound+roi->peakpos;
            }
            ROIs.push_back(roi);
        }
    }
    for (auto ii=ROIs.begin();ii!=ROIs.end();ii++){
        roi_t* s = *ii;
        s->print();
    }
    //! readall
}


RooAbsPdf* GetModel(RooAbsReal* p[],char* fitname, RooRealVar *x, parms_t* parms,Int_t category,Int_t group,Int_t sizeparms)
{
    RooRealVar* pvar[100];
    RooAbsReal* pp[MAX_N_PARMS];
    RooAbsReal* pbkg;
    RooAbsReal* pbkgratio;
    for (Int_t i=0;i<sizeparms;i++){
        if ((parms->category[i]==category&&parms->group[i]==group)||parms->category[i]==-1){
            pvar[i]=(RooRealVar*) p[i];
            pvar[i]->setMax(parms->max[i]);
            pvar[i]->setMin(parms->min[i]);
            pvar[i]->setVal(parms->val[i]);
            pvar[i]->setConstant(parms->isFix[i]);
            if (parms->index[i]<MAX_N_PARMS)
                pp[parms->index[i]] = p[i];
            if (parms->index[i]==MAX_N_PARMS)
                pbkg = p[i];
            else if (parms->index[i]==MAX_N_PARMS+1)
                pbkgratio = p[i];
        }
    }
    fitF* fitmodel=new fitF(Form("peak_%s",fitname),"peak",*x,pp);
    RooPolynomial* pol1 = new RooPolynomial(Form("bkg_%s",fitname), "Linear background", *x, RooArgSet(*pbkg));
    RooAbsPdf* finalfitmodel = new RooAddPdf(Form("finalmodel_%s",fitname),"finalmodel",RooArgList(*pol1,*fitmodel),*pbkgratio);
    return finalfitmodel;
}

TTree* tree =NULL;
void getTree(treedata_t* treedata,char* infile){
    //////////////////////////////////////////////////////////
    //   This file has been automatically generated
    //     (Tue Mar  7 15:15:24 2023 by ROOT version6.22/08)
    //   from TTree tree/tree
    //   found on file: rootfiles/run5_0.root
    //////////////////////////////////////////////////////////


    //Reset ROOT and connect tree file
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
    if (!f) {
        f = new TFile(infile);
    }
//    f->GetObject("tree",tree);
    f->GetObject("tree_drift_corr",tree);

    // Set branch addresses.
    tree->SetBranchAddress("fUniqueID",&treedata->fUniqueID);
    tree->SetBranchAddress("fBits",&treedata->fBits);
    tree->SetBranchAddress("ts",&treedata->ts);
    tree->SetBranchAddress("tof1",&treedata->tof1);
    tree->SetBranchAddress("tof2",&treedata->tof2);
    tree->SetBranchAddress("e",&treedata->e);
    tree->SetBranchAddress("vetot",&treedata->vetot);
    tree->SetBranchAddress("pulsestart1",&treedata->pulsestart1);
    tree->SetBranchAddress("pulsestart2",&treedata->pulsestart2);
    tree->SetBranchAddress("pulsestop",&treedata->pulsestop);
    tree->SetBranchAddress("pulsecfdstop",&treedata->pulsecfdstop);
//    tree->Print();
}
void anaTOF(char* infile,char* parmsfile, char* outfile,Int_t nbins = 600){
    treedata_t* treedata = new treedata_t;
    getTree(treedata,infile);
    tree->Print();
    Long64_t nentries = tree->GetEntries();

    TFile* f1 = new TFile("test.root","recreate");


    getROI(parmsfile,true);
    Int_t nwindows = 10;

    Int_t Npeaks = ROIs.size();
    parms_t parms;
    parms.reset();

    Int_t pidx = 0;
    Int_t roiidx = 0;

    Double_t peakFitRange[2];
    Int_t sizeparms = 0;
    Double_t htof_all_edge = 200000000000;
    Double_t htof_all_edge_hi = -10000000000000;

    Int_t current_group_index = -1;
    for (auto ii=ROIs.begin();ii!=ROIs.end();ii++){
        roi_t* roi = *ii;
        if (roi->lowerbound<htof_all_edge)
            htof_all_edge = roi->lowerbound;
        if (roi->upperbound>htof_all_edge_hi)
            htof_all_edge_hi = roi->upperbound;
        if (roiidx==0){
            peakFitRange[0] = roi->lowerbound-roi->peakpos;
            peakFitRange[1] = roi->upperbound-roi->peakpos;
            //shape parameter 2
            parms.category[pidx] = -1;
            parms.index[pidx] = 2;
            parms.isFix[pidx] = kFALSE;
            parms.min[pidx] = 0;
            parms.max[pidx] = roi->shape2*10;
            parms.val[pidx] = roi->shape2;
            parms.peakNo[pidx] = roiidx;
            parms.group[pidx] = roi->group;
            pidx++;
            //shape parameter 3
            parms.category[pidx] = -1;
            parms.index[pidx] = 3;
            parms.isFix[pidx] = kFALSE;
            parms.min[pidx] = 0;
            parms.max[pidx] = roi->shape3*10;
            parms.val[pidx] = roi->shape3;
            parms.peakNo[pidx] = roiidx;
            parms.group[pidx] = roi->group;
            pidx++;
        }

        //shape parameter 1
        if (roi->group!=current_group_index){//first group occurence
            parms.category[pidx] = roiidx;
            parms.index[pidx] = 1;
            parms.isFix[pidx] = kFALSE;
            parms.min[pidx] = 0;
            parms.max[pidx] = roi->shape1*10;
            parms.val[pidx] = roi->shape1;
            parms.peakNo[pidx] = roiidx;
            parms.group[pidx] = roi->group;
            pidx++;
        }
        current_group_index = roi->group;

        // peak pos
        parms.category[pidx] = roiidx;
        parms.index[pidx] = 0;
        parms.isFix[pidx] = kFALSE;
        parms.min[pidx] = roi->lowerbound-roi->peakpos;
        parms.max[pidx] = roi->upperbound-roi->peakpos;
        parms.val[pidx] = 0.;
        parms.peakNo[pidx] = roiidx;
        parms.group[pidx] = roi->group;
        pidx++;
        // bkg slope
        parms.category[pidx] = roiidx;
        parms.index[pidx] = 4;
        parms.isFix[pidx] = roi->isbkgfix;
        parms.min[pidx] = -10;
        parms.max[pidx] = 10;
        parms.val[pidx] = roi->bkgslope;
        parms.peakNo[pidx] = roiidx;
        parms.group[pidx] = roi->group;
        pidx++;
        // signal to bkg ratio
        parms.category[pidx] = roiidx;
        parms.index[pidx] = 5;
        parms.isFix[pidx] = roi->isbkgfix;
        parms.min[pidx] = -10;
        parms.max[pidx] = 10;
        parms.val[pidx] = roi->snratio;
        parms.peakNo[pidx] = roiidx;
        parms.group[pidx] = roi->group;
        pidx++;

        roiidx++;
    }
    sizeparms = pidx;
    parms.print(sizeparms);
//    return;
//    cout<<htof_all_edge<<endl;
    Double_t nwindowD = (htof_all_edge_hi - htof_all_edge+2*(peakFitRange[1]-peakFitRange[0]))/(peakFitRange[1]-peakFitRange[0]);
    nwindows = (int) round(nwindowD);

    TH1F* htof_all = new TH1F("htof_all","htof_all",nbins*nwindows,htof_all_edge-(peakFitRange[1]-peakFitRange[0]),htof_all_edge-(peakFitRange[1]-peakFitRange[0])+(peakFitRange[1]-peakFitRange[0])*nwindows);


    Double_t xx[MAX_NPEAKS];
    TTree* treeFit[MAX_NPEAKS];
    for (Int_t i=0;i<Npeaks;i++){
        treeFit[i] = new TTree(Form("treeFit%d",i),Form("treeFit%d",i));
        treeFit[i]->Branch("x",&xx[i],"x/D");
    }

    for (Long64_t i=0; i<nentries;i++) {
        tree->GetEntry(i);
        if (treedata->vetot==0){
            for (Int_t i=0;i<Npeaks;i++){
                xx[i] = treedata->tof1-ROIs[i]->peakpos;
                treeFit[i]->Fill();
            }
            htof_all->Fill(treedata->tof1);
        }
    }

    RooRealVar *x=new RooRealVar("x","x",peakFitRange[0],peakFitRange[1]) ;

    RooCategory sample("sample","sample") ;
    for (Int_t i=0;i<Npeaks;i++){
        sample.defineType(Form("data%d",i));
    }

    RooDataSet* data[MAX_NPEAKS];

    std::map<std::string,RooDataSet*> DataCollection;
    for (Int_t i=0;i<Npeaks;i++){
        data[i] = new RooDataSet(Form("data%d",i),Form("data%d",i),*x,Import(*treeFit[i]));
        DataCollection.emplace(Form("data%d",i),data[i]);
    }
    RooDataSet* combData = new RooDataSet("combData","combined data",*x,Index(sample),Import(DataCollection)) ;

    RooAbsReal* p[100];
    for (Int_t i=0;i<sizeparms;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),10,0,20);
    }

    RooAbsPdf* fitModel[MAX_NPEAKS];
    for (Int_t i=0;i<Npeaks;i++){
        fitModel[i] = GetModel(p,Form("model%d",i),x,&parms,i,ROIs[i]->group,sizeparms);
    }

    RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simultaneous pdf",sample) ;
    for (Int_t i=0;i<Npeaks;i++){
        simPdf->addPdf(*fitModel[i],Form("data%d",i)) ;
    }

    simPdf->fitTo(*combData,NumCPU(16),Save(kTRUE),PrintLevel(3));


    TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
    Int_t npadrow = Npeaks/2+Npeaks%2;
    c1->Divide(2,npadrow);

    RooCurve* fitcurve[MAX_NPEAKS];
    for (Int_t i=0;i<Npeaks;i++){
        c1->cd(i+1);
        RooPlot* xframe0 = x->frame(Title(Form("Peak %d",i))) ;
        combData->plotOn(xframe0,Binning(nbins),Cut(Form("sample==sample::data%d",i))) ;
        simPdf->plotOn(xframe0,Slice(sample,Form("data%d",i)),ProjWData(sample,*combData), RooFit::Name(Form("peak%d",i))) ;
        xframe0->Draw();
        fitcurve[i]=(RooCurve*) xframe0->getCurve(Form("peak%d",i));
    }
    c1->Draw();


    TCanvas* c2 = new TCanvas("c2","c2",900,700);
    c2->cd()->SetLogy();
    htof_all->SetMarkerSize(1);
    htof_all->SetMarkerStyle(20);
    htof_all->Draw("E0");
    TGraph* grPeak[MAX_NPEAKS];
    for (Int_t i=0;i<Npeaks;i++){
        vector<double> xpt;
        vector<double> ypt;
        for (Int_t j=0;j<fitcurve[i]->GetN();j++){
            xpt.push_back(fitcurve[i]->GetX()[j]+ROIs[i]->peakpos);
            ypt.push_back(fitcurve[i]->GetY()[j]);
        }
        grPeak[i] = new TGraph(xpt.size(),&xpt[0],&ypt[0]);
        grPeak[i]->SetName(Form("fit_peak%d",i));
        grPeak[i]->SetLineColor(2);
        grPeak[i]->SetLineWidth(2);
        grPeak[i]->Draw("L same");
    }
    c2->Draw();

    ofstream ofs(outfile);
    ofs<<"Peak positions"<<endl;
    ofs<<"id,value,error"<<std::setprecision(15)<<endl;
    for (Int_t i=0;i<sizeparms;i++){
        if (parms.peakNo[i]>=0 && parms.index[i]==0){
            RooRealVar* pvar=(RooRealVar*) p[i];
            Double_t peakPos = pvar->getVal()+ROIs[parms.peakNo[i]]->peakpos;
            ofs<<parms.peakNo[i]<<","<<peakPos<<","<<pvar->getError()<<endl;
        }
    }
    ofs<<endl;
    ofs<<"Shape parameters"<<endl;
    ofs<<"id,value,error"<<endl;
    for (Int_t i=0;i<sizeparms;i++){
        if (parms.peakNo[i]>=0 && parms.index[i]>0){
            RooRealVar* pvar=(RooRealVar*) p[i];
            ofs<<parms.peakNo[i]<<","<<pvar->getVal()<<","<<pvar->getError()<<endl;
        }
    }

}
