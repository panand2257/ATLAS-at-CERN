
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TKey.h"
#include "TCanvas.h"
#include "Riostream.h"
#include <vector>
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"




namespace Part{
    enum Part{
        partEBA = 0,
        partLBA = 1,
        partLBC = 2,
        partEBC = 3
    };
}


TH1F* A = new TH1F("A", "Channel A",50,-50.,50.);
TH1F* B = new TH1F("B", "Channel B",50,-50.,50.);
TH1F* D = new TH1F("D", "Channel D",50,-50.,50.);

Double_t meanA[95];
Double_t meanB[95];
Double_t meanD[95];

Double_t rmsA[95];
Double_t rmsB[95];
Double_t rmsD[95];

TDatime da[95];
Double_t x[95];

//Reading files

double Las_constant[4][64][48];
double Las_Vol[4][64][48];






//---------------------------------------
//MODIFY
//---------------------------------------
string extension = ".eps";

void separatePartMod(string &tmp_part, int &tmp_module)
{
    stringstream tmp_part_stream(tmp_part.substr(3,2));
    tmp_part_stream >> tmp_module;
    tmp_part = tmp_part.substr(0,3);
}

void evo()
{    //-------- read the data
    //EBC64 38 0    1.000000  -1.000000


  











    string tmp_part;
    int    tmp_module;
    int    tmp_channel;
    int    tmp_zero;
    float  tmp_constant;
    float  tmp_HV;
    string mask;
    string line;
   
   ifstream f_date ;


   string run;
   string sdate;
   string stime;
   string edate;
   string etime;

   f_date.open("laser_dates.txt");

int n=0;
    while(!f_date.eof())
    {
     
     f_date >> run;
     f_date >> sdate;
     f_date >> stime;
     f_date >> edate;
     f_date >> etime;
    
   string date = sdate + " " + stime ;
      

     
    
     string f = "LaserCombined_";
     f.append(run);
     f.append(".txt");
     string g = "masked_channels_";
     g.append(run);
     g.append(".txt");
     ifstream f_Laser;
     ifstream f_Masked;
        
         f_Laser.open(f);
         f_Masked.open(g);
         int j=0;
         if(f_Laser){
          while(!f_Laser.eof())
           {
            f_Laser >> tmp_part;
	    f_Laser >> tmp_channel;
            mask = tmp_part + " " + std::to_string(tmp_channel);
            separatePartMod(tmp_part, tmp_module);

            f_Laser >> tmp_zero;
            f_Laser >> tmp_constant;
            f_Laser >> tmp_HV;
            
            int myPart_ctr = -1;
            if (tmp_part == "EBA") myPart_ctr = Part::partEBA;
            if (tmp_part == "LBA") myPart_ctr = Part::partLBA;
            if (tmp_part == "LBC") myPart_ctr = Part::partLBC;
            if (tmp_part == "EBC") myPart_ctr = Part::partEBC;
            
          while(std::getline(f_Masked,line))
            {
             if (line.compare(mask) ==0 || tmp_HV < 0.) 
               {  
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
                 j--;
                 break;
               }
              else{ 
             double drift= (tmp_constant-1)*100;
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
f_Masked.close();    
f_Laser.close();
    
    
    // Fill histograms
    
    
    for(int m = 1; m < 65; m++)
    {
        
          A -> Fill(Las_constant[1][m-1][1]);        //A1
           A -> Fill(Las_constant[1][m-1][4]);
            
           A-> Fill(Las_constant[1][m-1][5]);        //A2
           A-> Fill(Las_constant[1][m-1][8]);
            
           A -> Fill(Las_constant[1][m-1][9]);        //A3
           A -> Fill(Las_constant[1][m-1][10]);
            
           A -> Fill(Las_constant[1][m-1][15]);       //A4
           A -> Fill(Las_constant[1][m-1][18]);
            
           A -> Fill(Las_constant[1][m-1][19]);       //A5
           A -> Fill(Las_constant[1][m-1][20]);
             
           A -> Fill(Las_constant[1][m-1][23]);       //A6
           A -> Fill(Las_constant[1][m-1][26]);
           
           A -> Fill(Las_constant[1][m-1][29]);       //A7
           A -> Fill(Las_constant[1][m-1][32]);
            
           A -> Fill(Las_constant[1][m-1][35]);       //A8
           A -> Fill(Las_constant[1][m-1][38]);
            
           A -> Fill(Las_constant[1][m-1][36]);       //A9
           A -> Fill(Las_constant[1][m-1][37]);
            
           A -> Fill(Las_constant[1][m-1][45]);       //A10
           A -> Fill(Las_constant[1][m-1][46]);
            
           A -> Fill(Las_constant[0][m-1][6]);         //A12
           A -> Fill(Las_constant[0][m-1][7]);

           A -> Fill(Las_constant[0][m-1][10]);         //A13
           A -> Fill(Las_constant[0][m-1][11]);
    
           A -> Fill(Las_constant[0][m-1][20]);        //A14
           A -> Fill(Las_constant[0][m-1][21]);

           A -> Fill(Las_constant[0][m-1][31]);        //A15
           A -> Fill(Las_constant[0][m-1][32]);
 
           A -> Fill(Las_constant[0][m-1][40]);        //A16
           A -> Fill(Las_constant[0][m-1][41]);
 
           // Right Modules

           A -> Fill(Las_constant[2][m-1][1]);         //A-1
           A -> Fill(Las_constant[2][m-1][4]);
            
           A -> Fill(Las_constant[2][m-1][5]);         //A-2
           A -> Fill(Las_constant[2][m-1][8]);
            
           A -> Fill(Las_constant[2][m-1][9]);        //A-3
           A -> Fill(Las_constant[2][m-1][10]);
            
           A -> Fill(Las_constant[2][m-1][15]);        //A-4
           A -> Fill(Las_constant[2][m-1][18]);
            
           A -> Fill(Las_constant[2][m-1][19]);       //A-5
           A -> Fill(Las_constant[2][m-1][20]);
            
           A -> Fill(Las_constant[2][m-1][23]);       //A-6
           A -> Fill(Las_constant[2][m-1][26]);
            
           A -> Fill(Las_constant[2][m-1][29]);       //A-7
           A -> Fill(Las_constant[2][m-1][32]);
            
           A -> Fill(Las_constant[2][m-1][35]);      //A-8
           A -> Fill(Las_constant[2][m-1][38]);
            
           A -> Fill(Las_constant[2][m-1][36]);       //A-9
           A -> Fill(Las_constant[2][m-1][37]);
            
           A -> Fill(Las_constant[2][m-1][45]);       //A-10
           A -> Fill(Las_constant[2][m-1][46]);
            
           A -> Fill(Las_constant[3][m-1][6]);         //A-12
           A -> Fill(Las_constant[3][m-1][7]);

           A -> Fill(Las_constant[3][m-1][10]);         //A-13
           A -> Fill(Las_constant[3][m-1][11]);
    
           A -> Fill(Las_constant[3][m-1][20]);        //A-14
           A -> Fill(Las_constant[3][m-1][21]);

           A -> Fill(Las_constant[3][m-1][31]);        //A-15
           A -> Fill(Las_constant[3][m-1][32]);
 
           A -> Fill(Las_constant[3][m-1][40]);        //A-16
           A -> Fill(Las_constant[3][m-1][41]);






          // Left Modules
           B -> Fill(Las_constant[1][m-1][2]);         //B1
           B -> Fill(Las_constant[1][m-1][3]);
            
           B -> Fill(Las_constant[1][m-1][6]);         //B2
           B -> Fill(Las_constant[1][m-1][7]);
            
           B -> Fill(Las_constant[1][m-1][11]);        //B3
           B -> Fill(Las_constant[1][m-1][12]);
            
           B -> Fill(Las_constant[1][m-1][16]);        //B4
           B -> Fill(Las_constant[1][m-1][17]);
            
           B -> Fill(Las_constant[1][m-1][21]);        //B5
           B -> Fill(Las_constant[1][m-1][22]);
           
           B -> Fill(Las_constant[1][m-1][27]);       //B6
           B -> Fill(Las_constant[1][m-1][28]);
            
           B -> Fill(Las_constant[1][m-1][33]);        //B7
           B -> Fill(Las_constant[1][m-1][34]);
            
           B -> Fill(Las_constant[1][m-1][39]);       //B8
           B -> Fill(Las_constant[1][m-1][40]);
            
           B -> Fill(Las_constant[1][m-1][42]);        //B9
           B -> Fill(Las_constant[1][m-1][47]);
            
           B -> Fill(Las_constant[0][m-1][8]);           //B11
           B -> Fill(Las_constant[0][m-1][9]);
            
           B -> Fill(Las_constant[0][m-1][14]);         //B12
           B -> Fill(Las_constant[0][m-1][15]);

           B -> Fill(Las_constant[0][m-1][22]);        //B13
           B -> Fill(Las_constant[0][m-1][23]);
    
           B -> Fill(Las_constant[0][m-1][30]);        //B14
           B -> Fill(Las_constant[0][m-1][35]);

           B -> Fill(Las_constant[0][m-1][36]);        //B15
           B -> Fill(Las_constant[0][m-1][39]);
 


           // Right Modules


           B -> Fill(Las_constant[2][m-1][2]);        //B-1
           B -> Fill(Las_constant[2][m-1][3]);
 
           B -> Fill(Las_constant[2][m-1][6]);        //B-2
           B -> Fill(Las_constant[2][m-1][7]);
            
           B -> Fill(Las_constant[2][m-1][11]);        //B-3
           B -> Fill(Las_constant[2][m-1][12]);
            
           B -> Fill(Las_constant[2][m-1][16]);      //B-4
           B -> Fill(Las_constant[2][m-1][17]);
            
           B -> Fill(Las_constant[2][m-1][21]);      //B-5
           B -> Fill(Las_constant[2][m-1][22]);
            
           B -> Fill(Las_constant[2][m-1][27]);        //B-6
           B -> Fill(Las_constant[2][m-1][28]);
            
           B -> Fill(Las_constant[2][m-1][33]);      //B-7
           B -> Fill(Las_constant[2][m-1][34]);
            
           B -> Fill(Las_constant[2][m-1][39]);      //B-8
           B -> Fill(Las_constant[2][m-1][40]);
            
           B -> Fill(Las_constant[2][m-1][42]);      //B-9
           B -> Fill(Las_constant[2][m-1][47]);
            
           B -> Fill(Las_constant[3][m-1][8]);     //B-11
           B -> Fill(Las_constant[3][m-1][9]);
            
           B -> Fill(Las_constant[3][m-1][14]);      //B-12
           B -> Fill(Las_constant[3][m-1][15]);
            
           B -> Fill(Las_constant[3][m-1][22]);         //B-13
           B -> Fill(Las_constant[3][m-1][23]);

           B -> Fill(Las_constant[3][m-1][30]);         //B-14
           B -> Fill(Las_constant[3][m-1][35]);
    
           B -> Fill(Las_constant[3][m-1][36]);        //B-15
           B -> Fill(Las_constant[3][m-1][39]);

 
           D -> Fill(Las_constant[1][m-1][0]);         //D0
           D -> Fill(Las_constant[1][m-1][0]);
           D -> Fill(Las_constant[2][m-1][0]);         //D0
           D -> Fill(Las_constant[2][m-1][0]);        
            
           D -> Fill(Las_constant[1][m-1][13]);         //D1
           D -> Fill(Las_constant[1][m-1][14]);
            
           D -> Fill(Las_constant[1][m-1][24]);        //D2
           D -> Fill(Las_constant[1][m-1][25]);
            
           D -> Fill(Las_constant[1][m-1][41]);        //D3
           D -> Fill(Las_constant[1][m-1][44]);
            
           D -> Fill(Las_constant[0][m-1][2]);        //D4
           D -> Fill(Las_constant[0][m-1][3]);
           
           D -> Fill(Las_constant[0][m-1][16]);       //D5
           D -> Fill(Las_constant[0][m-1][17]);
            
           D -> Fill(Las_constant[0][m-1][37]);        //D6
           D -> Fill(Las_constant[0][m-1][38]);
            
           D -> Fill(Las_constant[2][m-1][13]);       //D-1
           D -> Fill(Las_constant[2][m-1][14]);
            
           D -> Fill(Las_constant[2][m-1][24]);        //D-2
           D -> Fill(Las_constant[2][m-1][25]);
            
           D -> Fill(Las_constant[2][m-1][41]);           //D-3
           D -> Fill(Las_constant[2][m-1][44]);
            
           D -> Fill(Las_constant[3][m-1][2]);         //D-4
           D -> Fill(Las_constant[3][m-1][3]);

           D -> Fill(Las_constant[3][m-1][16]);        //D-5
           D -> Fill(Las_constant[3][m-1][17]);
    
           D -> Fill(Las_constant[3][m-1][37]);        //D-6
           D -> Fill(Las_constant[3][m-1][38]);

            



       
 //loop over modules
    

}
//cout<<A->GetMean()<<endl;
n++;
meanA[n-1] = (A->GetMean());
meanB[n-1] =(B->GetMean());
meanD[n-1] = (D->GetMean());

rmsA[n-1] = A->GetMeanError();
rmsB[n-1] = B->GetMeanError();
rmsD[n-1] = D->GetMeanError();

const char *cstr = date.c_str();

da[n-1].Set(cstr);
x[n-1] = da[n-1].Convert();

A ->Reset("ICES");
B ->Reset("ICES");
D ->Reset("ICES");
}
}



auto myCanvas = new TCanvas("myCanvas", "",484,86,800,750);
 gROOT->Reset();
   TPad *pad = new TPad("pad","",0,0,1,1);
   pad->SetFillColor(kWhite);
   
   pad->Draw();
   pad->cd();

      // draw a frame to define the range
   TH1F *hr = myCanvas->DrawFrame(-0.4,0,1.2,12);
   hr->SetXTitle("X title");
   hr->SetYTitle("Y title");
   pad->GetFrame()->SetFillColor(kWhite);
   pad->GetFrame()->SetBorderSize(12);



TGraphErrors* grA = new TGraphErrors(95,x,meanA,nullptr,rmsA); 
TGraphErrors* grB = new TGraphErrors(95,x,meanB,nullptr,rmsB); 
TGraphErrors* grD = new TGraphErrors(95,x,meanD,nullptr,rmsD);




    




Double_t _fx1[343] = {
   1.461376e+09,
   1.461376e+09,
   1.461389e+09,
   1.461389e+09,
   1.461476e+09,
   1.461476e+09,
   1.461561e+09,
   1.461561e+09,
   1.461733e+09,
   1.461733e+09,
   1.461895e+09,
   1.461895e+09,
   1.462568e+09,
   1.462568e+09,
   1.462585e+09,
   1.462585e+09,
   1.462611e+09,
   1.462611e+09,
   1.462675e+09,
   1.462675e+09,
   1.46273e+09,
   1.46273e+09,
   1.462746e+09,
   1.462746e+09,
   1.462836e+09,
   1.462836e+09,
   1.462854e+09,
   1.462854e+09,
   1.462928e+09,
   1.462928e+09,
   1.463021e+09,
   1.463021e+09,
   1.46311e+09,
   1.46311e+09,
   1.463198e+09,
   1.463198e+09,
   1.463235e+09,
   1.463235e+09,
   1.463285e+09,
   1.463285e+09,
   1.463396e+09,
   1.463396e+09,
   1.463458e+09,
   1.463458e+09,
   1.463555e+09,
   1.463555e+09,
   1.463601e+09,
   1.463601e+09,
   1.463827e+09,
   1.463827e+09,
   1.464311e+09,
   1.464311e+09,
   1.464344e+09,
   1.464344e+09,
   1.46438e+09,
   1.46438e+09,
   1.464452e+09,
   1.464452e+09,
   1.464475e+09,
   1.464475e+09,
   1.464573e+09,
   1.464573e+09,
   1.464653e+09,
   1.464653e+09,
   1.464692e+09,
   1.464692e+09,
   1.464749e+09,
   1.464749e+09,
   1.464823e+09,
   1.464823e+09,
   1.464898e+09,
   1.464898e+09,
   1.464952e+09,
   1.464952e+09,
   1.465021e+09,
   1.465021e+09,
   1.465127e+09,
   1.465127e+09,
   1.465184e+09,
   1.465184e+09,
   1.465651e+09,
   1.465651e+09,
   1.465754e+09,
   1.465754e+09,
   1.465872e+09,
   1.465872e+09,
   1.465978e+09,
   1.465978e+09,
   1.466021e+09,
   1.466021e+09,
   1.466138e+09,
   1.466138e+09,
   1.46625e+09,
   1.46625e+09,
   1.466292e+09,
   1.466292e+09,
   1.466324e+09,
   1.466324e+09,
   1.466361e+09,
   1.466361e+09,
   1.466479e+09,
   1.466479e+09,
   1.466731e+09,
   1.466731e+09,
   1.466915e+09,
   1.466915e+09,
   1.467071e+09,
   1.467071e+09,
   1.467134e+09,
   1.467134e+09,
   1.467215e+09,
   1.467215e+09,
   1.467243e+09,
   1.467243e+09,
   1.467306e+09,
   1.467306e+09,
   1.467333e+09,
   1.467333e+09,
   1.467436e+09,
   1.467436e+09,
   1.467538e+09,
   1.467538e+09,
   1.467608e+09,
   1.467608e+09,
   1.467645e+09,
   1.467645e+09,
   1.467748e+09,
   1.467748e+09,
   1.467856e+09,
   1.467856e+09,
   1.467968e+09,
   1.467968e+09,
   1.46806e+09,
   1.46806e+09,
   1.468125e+09,
   1.468125e+09,
   1.468225e+09,
   1.468225e+09,
   1.468319e+09,
   1.468319e+09,
   1.468439e+09,
   1.468439e+09,
   1.46855e+09,
   1.46855e+09,
   1.468643e+09,
   1.468643e+09,
   1.468732e+09,
   1.468732e+09,
   1.468804e+09,
   1.468804e+09,
   1.468896e+09,
   1.468896e+09,
   1.469043e+09,
   1.469043e+09,
   1.469127e+09,
   1.469127e+09,
   1.469164e+09,
   1.469164e+09,
   1.469193e+09,
   1.469193e+09,
   1.469264e+09,
   1.469264e+09,
   1.469276e+09,
   1.469276e+09,
   1.46935e+09,
   1.46935e+09,
   1.469371e+09,
   1.469371e+09,
   1.469414e+09,
   1.469414e+09,
   1.469473e+09,
   1.469473e+09,
   1.470019e+09,
   1.470019e+09,
   1.47007e+09,
   1.47007e+09,
   1.470127e+09,
   1.470127e+09,
   1.470292e+09,
   1.470292e+09,
   1.470307e+09,
   1.470307e+09,
   1.470365e+09,
   1.470365e+09,
   1.470415e+09,
   1.470415e+09,
   1.470425e+09,
   1.470425e+09,
   1.470531e+09,
   1.470531e+09,
   1.470578e+09,
   1.470578e+09,
   1.470654e+09,
   1.470654e+09,
   1.470707e+09,
   1.470707e+09,
   1.470786e+09,
   1.470786e+09,
   1.471083e+09,
   1.471083e+09,
   1.471149e+09,
   1.471149e+09,
   1.471236e+09,
   1.471236e+09,
   1.471326e+09,
   1.471326e+09,
   1.471347e+09,
   1.471347e+09,
   1.471378e+09,
   1.471378e+09,
   1.471399e+09,
   1.471399e+09,
   1.471462e+09,
   1.471462e+09,
   1.471568e+09,
   1.471568e+09,
   1.471617e+09,
   1.471617e+09,
   1.471639e+09,
   1.471639e+09,
   1.471744e+09,
   1.471744e+09,
   1.472063e+09,
   1.472063e+09,
   1.472101e+09,
   1.472101e+09,
   1.472172e+09,
   1.472172e+09,
   1.472272e+09,
   1.472272e+09,
   1.472302e+09,
   1.472302e+09,
   1.472338e+09,
   1.472338e+09,
   1.472398e+09,
   1.472398e+09,
   1.47247e+09,
   1.47247e+09,
   1.472569e+09,
   1.472569e+09,
   1.472618e+09,
   1.472618e+09,
   1.472633e+09,
   1.472633e+09,
   1.472695e+09,
   1.472696e+09,
   1.472744e+09,
   1.472744e+09,
   1.472773e+09,
   1.472773e+09,
   1.472872e+09,
   1.472872e+09,
   1.472895e+09,
   1.472895e+09,
   1.472929e+09,
   1.472929e+09,
   1.473027e+09,
   1.473027e+09,
   1.473136e+09,
   1.473136e+09,
   1.473225e+09,
   1.473225e+09,
   1.473391e+09,
   1.473391e+09,
   1.473451e+09,
   1.473451e+09,
   1.474794e+09,
   1.474794e+09,
   1.474812e+09,
   1.474812e+09,
   1.474884e+09,
   1.474884e+09,
   1.474954e+09,
   1.474954e+09,
   1.475033e+09,
   1.475033e+09,
   1.475103e+09,
   1.475103e+09,
   1.475253e+09,
   1.475253e+09,
   1.475309e+09,
   1.475309e+09,
   1.475458e+09,
   1.475458e+09,
   1.475746e+09,
   1.475746e+09,
   1.475918e+09,
   1.475918e+09,
   1.475961e+09,
   1.475961e+09,
   1.475995e+09,
   1.475995e+09,
   1.476071e+09,
   1.476071e+09,
   1.476163e+09,
   1.476163e+09,
   1.476173e+09,
   1.476173e+09,
   1.476248e+09,
   1.476248e+09,
   1.476318e+09,
   1.476318e+09,
   1.476355e+09,
   1.476355e+09,
   1.476436e+09,
   1.476436e+09,
   1.476548e+09,
   1.476548e+09,
   1.476602e+09,
   1.476602e+09,
   1.47669e+09,
   1.47669e+09,
   1.476707e+09,
   1.476707e+09,
   1.476801e+09,
   1.476801e+09,
   1.47682e+09,
   1.47682e+09,
   1.476861e+09,
   1.476861e+09,
   1.476938e+09,
   1.476938e+09,
   1.477033e+09,
   1.477033e+09,
   1.477126e+09,
   1.477126e+09,
   1.477179e+09,
   1.477179e+09,
   1.477237e+09,
   1.477237e+09,
   1.477306e+09,
   1.477306e+09,
   1.477361e+09,
   1.477361e+09,
   1.477417e+09,
   1.477417e+09,
   1.477455e+09,
   1.477455e+09,
   1.477504e+09,
   1.477504e+09,
   1.479e+09,
   1.48025e+09,
   1.48033e+09,
   };
   Double_t _fy1[343] = {
   0,
   0.063,
   0.063,
   0.112,
   0.112,
   0.697,
   0.697,
   1.505,
   1.505,
   2.142,
   2.142,
   6.211,
   6.211,
   6.245,
   6.245,
   7.441,
   7.441,
   9.821,
   9.821,
   16.565,
   16.565,
   38.859,
   38.859,
   40.701,
   40.701,
   49.489,
   49.489,
   58.61,
   58.61,
   85.748,
   85.748,
   142.501,
   142.501,
   194.646,
   194.646,
   227.16,
   227.16,
   257.004,
   257.004,
   355.664,
   355.664,
   387.336,
   387.336,
   460.461,
   460.461,
   499.509,
   499.509,
   499.553,
   499.553,
   766.615,
   766.615,
   776.551,
   776.551,
   776.602,
   776.602,
   783.93,
   783.93,
   890.575,
   890.575,
   898.221,
   898.221,
   1043.39,
   1043.39,
   1252.324,
   1252.324,
   1413.655,
   1413.655,
   1474.106,
   1474.106,
   1722.909,
   1722.909,
   2077.174,
   2077.174,
   2085.999,
   2085.999,
   2399.022,
   2399.022,
   2759.846,
   2759.846,
   2999.194,
   2999.194,
   3045.206,
   3045.206,
   3480.304,
   3480.304,
   3972.778,
   3972.778,
   4389.383,
   4389.383,
   4598.184,
   4598.184,
   4805.598,
   4805.598,
   5215.412,
   5215.412,
   5427.456,
   5427.456,
   5540.237,
   5540.237,
   5692.577,
   5692.577,
   6227.892,
   6227.892,
   6318.362,
   6318.362,
   6869.96,
   6869.96,
   7585.989,
   7585.989,
   7844.155,
   7844.155,
   8211.458,
   8211.458,
   8231.989,
   8231.989,
   8336.1,
   8336.1,
   8463.075,
   8463.075,
   8584.492,
   8584.492,
   9152.897,
   9152.897,
   9524.799,
   9524.799,
   9706.352,
   9706.352,
   10219.01,
   10219.01,
   10760.76,
   10760.76,
   11284.18,
   11284.18,
   11736.94,
   11736.94,
   12084.35,
   12084.35,
   12563.03,
   12563.03,
   12843.46,
   12843.46,
   13360.89,
   13360.89,
   13885.08,
   13885.08,
   14367.28,
   14367.28,
   14868.37,
   14868.37,
   15036.4,
   15036.4,
   15536.03,
   15536.03,
   16130.85,
   16130.85,
   16587.31,
   16587.31,
   16769.97,
   16769.97,
   16908.12,
   16908.12,
   17283.15,
   17283.15,
   17315.61,
   17315.61,
   17761.13,
   17761.13,
   17869.12,
   17869.12,
   18020.97,
   18020.97,
   18081.76,
   18081.76,
   18081.8,
   18081.8,
   18094.8,
   18094.8,
   18427.49,
   18427.49,
   18734.61,
   18734.61,
   18768.7,
   18768.7,
   19087.7,
   19087.7,
   19116.36,
   19116.36,
   19136.89,
   19136.89,
   19532.09,
   19532.09,
   19669.06,
   19669.06,
   20003.43,
   20003.43,
   20287.76,
   20287.76,
   20555.54,
   20555.54,
   20721.75,
   20721.75,
   21123.42,
   21123.42,
   21594.84,
   21594.84,
   21800.51,
   21800.51,
   21896.98,
   21896.98,
   21933.63,
   21933.63,
   22047.1,
   22047.1,
   22434.52,
   22434.52,
   22536.35,
   22536.35,
   22562.66,
   22562.66,
   22577.28,
   22577.28,
   22735.53,
   22735.53,
   22735.56,
   22735.56,
   22764.07,
   22764.07,
   22989.46,
   22989.46,
   23396.97,
   23396.97,
   23593.47,
   23593.47,
   23682.22,
   23682.22,
   24076.35,
   24076.35,
   24455.81,
   24455.81,
   24956.26,
   24956.26,
   25171.54,
   25171.54,
   25251.7,
   25251.7,
   25608.07,
   25608.07,
   25929.6,
   25929.6,
   26001.53,
   26001.53,
   26460.18,
   26460.18,
   26549.99,
   26549.99,
   26783.53,
   26783.53,
   27345.83,
   27345.83,
   27659.39,
   27659.39,
   28094.01,
   28094.01,
   28590.93,
   28590.93,
   28957.14,
   28957.14,
   28963.24,
   28963.24,
   28971.07,
   28971.07,
   29050.46,
   29050.46,
   29517.05,
   29517.05,
   29991.08,
   29991.08,
   30382.01,
   30382.01,
   30782.26,
   30782.26,
   31148.07,
   31148.07,
   31538.72,
   31538.72,
   31562.91,
   31562.91,
   31563.07,
   31563.07,
   31571.82,
   31571.82,
   31609.44,
   31609.44,
   32085.17,
   32085.17,
   32538.67,
   32538.67,
   32583.19,
   32583.19,
   32787.63,
   32787.63,
   32956.33,
   32956.33,
   33209.18,
   33209.18,
   33215.76,
   33215.76,
   33796.66,
   33796.66,
   34151.12,
   34151.12,
   34581.33,
   34581.33,
   34581.38,
   34581.38,
   35080.74,
   35080.74,
   35204.43,
   35204.43,
   35420.99,
   35420.99,
   35614.86,
   35614.86,
   35861.43,
   35861.43,
   36261.96,
   36261.96,
   36535.99,
   36535.99,
   36947.22,
   36947.22,
   37333.87,
   37333.87,
   37649.63,
   37649.63,
   38006.91,
   38006.91,
   38115.67,
   38115.67,
   38479.26,
   38479.26,
   38479.26,
   0};
   TGraph *graph = new TGraph(343,_fx1,_fy1);
  // graph->SetName("");
   graph->SetName("");
   graph->SetTitle("");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cccccc");
   graph->SetFillColor(ci);

   ci = TColor::GetColor("#cccccc");
   graph->SetLineColor(ci);
   graph->Draw("AFB1");
   graph->GetXaxis()->SetLimits(x[0],x[94]);
   graph->GetXaxis()->SetLabelOffset(999);
graph->GetXaxis()->SetLabelSize(0);
graph->GetXaxis()->SetTickSize(0);
graph->GetYaxis()->SetLabelSize(0.03);
graph->GetYaxis()->SetTitleSize(0.035);
graph->GetYaxis()->SetTitleOffset(1.36);
graph->GetYaxis()->SetLabelColor(ci-1);
graph->GetYaxis()->SetTitleColor(ci-1);



 /* graph->GetXaxis()->SetTimeDisplay(1);
graph->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%Y}%F1970-01-01 00:00:00s0");
graph->GetXaxis()->SetNdivisions(509);
graph->GetXaxis()->SetTimeDisplay(1);
graph->GetXaxis()->SetTimeOffset(0,"gmt");
graph->GetXaxis()->SetLabelOffset(0.013);
graph->GetXaxis()->SetLabelSize(0.02);
graph->GetXaxis()->SetTitleSize(0.035);
graph->GetXaxis()->SetTitleOffset(1.35);*/
  
graph->SetTitle(";Time [dd/mm and year];Total Integrated Delivered Luminosity [pb^{-1}] ");
myCanvas->cd();
   TPad *overlay = new TPad("overlay","",0,0,1,1);
   overlay->SetFillStyle(4000);
   overlay->SetFillColor(kWhite);
   overlay->SetFrameFillStyle(4000);
   overlay->Draw();
   overlay->cd();


Double_t xmin = pad->GetUxmin();
   Double_t ymin = pad->GetUymin();
   Double_t xmax = pad->GetUxmax();
   Double_t ymax = pad->GetUymax();
   TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
   hframe->GetXaxis()->SetLabelOffset(99);
   hframe->GetYaxis()->SetLabelOffset(99);
 
grA->SetMarkerStyle(21);
grA->SetMarkerColor(kRed);
grA->SetMarkerSize(1.5);
grA->GetYaxis()->SetTitle("Average Drift [%]");
/*grA->GetXaxis()->SetLabelOffset(999);
grA->GetXaxis()->SetLabelSize(0);
grA->GetXaxis()->SetTickSize(0);*/
grA->GetXaxis()->SetLimits(x[0],x[94]);
grA->GetXaxis()->SetTimeDisplay(1);
grA->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%Y}%F1970-01-01 00:00:00s0");
grA->GetXaxis()->SetNdivisions(509);
grA->GetXaxis()->SetTimeDisplay(1);
grA->GetXaxis()->SetTimeOffset(0,"gmt");
grA->GetXaxis()->SetLabelOffset(0.013);
grA->GetXaxis()->SetLabelSize(0.02);
grA->GetXaxis()->SetTitleSize(0.035);
grA->GetXaxis()->SetTitleOffset(1.36);
grA->SetLineColor(kRed);
grB->SetLineColor(kBlue);
grD->SetLineColor(kGreen);

grA->Draw("AY+PE");
grA->GetYaxis()->SetRangeUser(-7.,14.);
grB->SetMarkerStyle(22);
grB->SetMarkerColor(kBlue);
grB->SetMarkerSize(1.5);
grB->Draw("PE");

grD->SetMarkerStyle(20);
grD->SetMarkerColor(kGreen);
grD->SetMarkerSize(1.5);
grD->Draw("PE");

//grA->GetXaxis()->redraw();
TLegend *leg = new TLegend(0.7401669,0.7092683,0.9398093,0.9092683,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.03);
   //leg->SetHeader("Channels","C"); // option "C" allows to center the header
   leg->AddEntry(grA,"Channel A","p");
   leg->AddEntry(grB,"Channel BC","p");
   leg->AddEntry(grD,"Channel D","p");
         leg->AddEntry(graph,"Luminosity","f");
  leg->Draw();
  
TLatex printInfo;
    printInfo.SetNDC();
    printInfo.SetTextFont(72);
    printInfo.DrawLatex(0.15, 0.85, "#scale[0.9]{ATLAS}");
    printInfo.SetTextFont(42);
    printInfo.DrawLatex(0.25, 0.85, "#scale[0.9]{Internal}");
    printInfo.DrawLatex(0.15, 0.80, "#scale[0.9]{Tile Calorimeter}");
    printInfo.DrawLatex(0.15, 0.75, "#scale[0.9]{Laser Combined Method 2016}");



  

  
   
   

   


   
    
    
 } 

