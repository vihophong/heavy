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
 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
using namespace std;
 
int mcs_make_rootfilev0(char* inputfile, char* outputfile)
{
    TFile output_file(outputfile, "recreate");
    Int_t ch, edge, t, sweeps, tag, data_lost;
    Int_t data[6];
    TTree tree("tree", "tree");
    tree.Branch("ch", &data[0], "ch/I");
    tree.Branch("edge", &data[1], "edge/I");
    tree.Branch("t", &data[2], "t/I");
    tree.Branch("sweeps", &data[3], "sweeps/I");
    tree.Branch("tag", &data[4], "tag/I");
    tree.Branch("data_lost", &data[5], "data_lost/I");

    string line, word;
 
    fstream file (inputfile, ios::in);
    Int_t nrows = 0;
    if(file.is_open())
    {
        while(getline(file, line))
        {
            if (nrows>0){
                stringstream str(line);
                Int_t ncols = 0;
                while(getline(str, word, ',')){
                    if (ncols>0)
                        data[ncols-1]=atoi(word.c_str());
                    ncols++;
                }
                tree.Fill();
            }
            nrows++;
        }
    }
    else{
        cout<<"Could not open the file\n";
        }

    tree.Write();
    output_file.Close();
 
    return 0;
}
 
