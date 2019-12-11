
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

Double_t meanA[93];
Double_t meanB[93];
Double_t meanD[93];

Double_t rmsA[93];
Double_t rmsB[93];
Double_t rmsD[93];

TDatime da[93];
Double_t x[93];

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
             if (line.compare(mask) ==0) 
               {  
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
                 j--;
                 break;
               }
              else{ 
             double drift= (tmp_constant-1.)*100;
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

n++;
meanA[n-1] = (A->GetRMS());
meanB[n-1] =(B->GetRMS());
meanD[n-1] = (D->GetRMS());

rmsA[n-1] = A->GetRMSError();
rmsB[n-1] = B->GetRMSError();
rmsD[n-1] = D->GetRMSError();

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



TGraphErrors* grA = new TGraphErrors(93,x,meanA,nullptr,rmsA); 
TGraphErrors* grB = new TGraphErrors(93,x,meanB,nullptr,rmsB); 
TGraphErrors* grD = new TGraphErrors(93,x,meanD,nullptr,rmsD);




    




Double_t _fx1[432] = {
   1.49557e+09,
   1.49557e+09,
   1.495596e+09,
   1.495596e+09,
   1.495694e+09,
   1.495694e+09,
   1.495887e+09,
   1.495887e+09,
   1.495926e+09,
   1.495926e+09,
   1.495954e+09,
   1.495954e+09,
   1.496032e+09,
   1.496032e+09,
   1.496183e+09,
   1.496183e+09,
   1.496481e+09,
   1.496481e+09,
   1.496649e+09,
   1.496649e+09,
   1.496686e+09,
   1.496686e+09,
   1.496721e+09,
   1.496721e+09,
   1.497255e+09,
   1.497255e+09,
   1.497294e+09,
   1.497294e+09,
   1.497303e+09,
   1.497303e+09,
   1.49744e+09,
   1.49744e+09,
   1.497479e+09,
   1.497479e+09,
   1.49751e+09,
   1.49751e+09,
   1.497602e+09,
   1.497602e+09,
   1.497621e+09,
   1.497621e+09,
   1.497652e+09,
   1.497652e+09,
   1.497691e+09,
   1.497691e+09,
   1.497756e+09,
   1.497756e+09,
   1.49779e+09,
   1.49779e+09,
   1.497942e+09,
   1.497942e+09,
   1.49803e+09,
   1.49803e+09,
   1.498123e+09,
   1.498123e+09,
   1.498184e+09,
   1.498184e+09,
   1.498232e+09,
   1.498232e+09,
   1.498282e+09,
   1.498282e+09,
   1.498333e+09,
   1.498333e+09,
   1.498347e+09,
   1.498347e+09,
   1.498421e+09,
   1.498421e+09,
   1.498466e+09,
   1.498466e+09,
   1.498554e+09,
   1.498554e+09,
   1.49861e+09,
   1.49861e+09,
   1.498644e+09,
   1.498644e+09,
   1.498701e+09,
   1.498701e+09,
   1.498728e+09,
   1.498728e+09,
   1.498743e+09,
   1.498743e+09,
   1.498773e+09,
   1.498773e+09,
   1.499764e+09,
   1.499764e+09,
   1.499837e+09,
   1.499837e+09,
   1.50002e+09,
   1.50002e+09,
   1.500046e+09,
   1.500046e+09,
   1.500111e+09,
   1.500111e+09,
   1.500129e+09,
   1.500129e+09,
   1.500171e+09,
   1.500171e+09,
   1.500219e+09,
   1.500219e+09,
   1.500362e+09,
   1.500362e+09,
   1.500413e+09,
   1.500413e+09,
   1.50043e+09,
   1.50043e+09,
   1.500449e+09,
   1.500449e+09,
   1.500459e+09,
   1.500459e+09,
   1.500503e+09,
   1.500503e+09,
   1.500522e+09,
   1.500522e+09,
   1.500597e+09,
   1.500597e+09,
   1.500692e+09,
   1.500692e+09,
   1.500737e+09,
   1.500737e+09,
   1.500858e+09,
   1.500858e+09,
   1.501182e+09,
   1.501182e+09,
   1.501265e+09,
   1.501265e+09,
   1.501282e+09,
   1.501282e+09,
   1.5013e+09,
   1.5013e+09,
   1.501326e+09,
   1.501326e+09,
   1.501336e+09,
   1.501336e+09,
   1.501416e+09,
   1.501416e+09,
   1.501488e+09,
   1.501488e+09,
   1.501521e+09,
   1.501521e+09,
   1.501582e+09,
   1.501582e+09,
   1.50168e+09,
   1.50168e+09,
   1.501762e+09,
   1.501762e+09,
   1.501821e+09,
   1.501821e+09,
   1.501886e+09,
   1.501886e+09,
   1.501915e+09,
   1.501915e+09,
   1.501977e+09,
   1.501977e+09,
   1.502046e+09,
   1.502046e+09,
   1.502098e+09,
   1.502098e+09,
   1.502166e+09,
   1.502166e+09,
   1.502182e+09,
   1.502182e+09,
   1.502209e+09,
   1.502209e+09,
   1.502294e+09,
   1.502294e+09,
   1.50234e+09,
   1.50234e+09,
   1.502491e+09,
   1.502491e+09,
   1.502675e+09,
   1.502675e+09,
   1.502725e+09,
   1.502725e+09,
   1.502763e+09,
   1.502763e+09,
   1.502789e+09,
   1.502789e+09,
   1.502854e+09,
   1.502854e+09,
   1.502884e+09,
   1.502884e+09,
   1.502919e+09,
   1.502919e+09,
   1.502947e+09,
   1.502947e+09,
   1.502963e+09,
   1.502963e+09,
   1.502991e+09,
   1.502991e+09,
   1.503039e+09,
   1.503039e+09,
   1.50306e+09,
   1.50306e+09,
   1.503166e+09,
   1.503166e+09,
   1.503189e+09,
   1.503189e+09,
   1.503243e+09,
   1.503243e+09,
   1.503355e+09,
   1.503355e+09,
   1.503395e+09,
   1.503395e+09,
   1.50343e+09,
   1.50343e+09,
   1.503507e+09,
   1.503507e+09,
   1.503534e+09,
   1.503534e+09,
   1.503766e+09,
   1.503766e+09,
   1.503831e+09,
   1.503831e+09,
   1.503897e+09,
   1.503897e+09,
   1.503912e+09,
   1.503912e+09,
   1.503972e+09,
   1.503972e+09,
   1.50403e+09,
   1.50403e+09,
   1.504057e+09,
   1.504057e+09,
   1.504118e+09,
   1.504118e+09,
   1.504194e+09,
   1.504194e+09,
   1.504223e+09,
   1.504223e+09,
   1.50426e+09,
   1.50426e+09,
   1.504302e+09,
   1.504302e+09,
   1.50433e+09,
   1.50433e+09,
   1.504364e+09,
   1.504364e+09,
   1.504412e+09,
   1.504412e+09,
   1.504509e+09,
   1.504509e+09,
   1.504552e+09,
   1.504552e+09,
   1.504584e+09,
   1.504584e+09,
   1.504617e+09,
   1.504617e+09,
   1.504644e+09,
   1.504644e+09,
   1.504683e+09,
   1.504683e+09,
   1.504767e+09,
   1.504767e+09,
   1.504818e+09,
   1.504818e+09,
   1.504833e+09,
   1.504833e+09,
   1.504886e+09,
   1.504886e+09,
   1.504903e+09,
   1.504903e+09,
   1.504915e+09,
   1.504915e+09,
   1.504967e+09,
   1.504967e+09,
   1.505008e+09,
   1.505008e+09,
   1.505061e+09,
   1.505061e+09,
   1.505136e+09,
   1.505136e+09,
   1.505215e+09,
   1.505215e+09,
   1.505236e+09,
   1.505236e+09,
   1.505265e+09,
   1.505265e+09,
   1.505287e+09,
   1.505287e+09,
   1.506173e+09,
   1.506173e+09,
   1.506208e+09,
   1.506208e+09,
   1.506235e+09,
   1.506235e+09,
   1.506287e+09,
   1.506287e+09,
   1.506342e+09,
   1.506342e+09,
   1.506405e+09,
   1.506405e+09,
   1.506463e+09,
   1.506463e+09,
   1.50651e+09,
   1.50651e+09,
   1.506585e+09,
   1.506585e+09,
   1.506641e+09,
   1.506641e+09,
   1.506717e+09,
   1.506717e+09,
   1.506787e+09,
   1.506787e+09,
   1.506855e+09,
   1.506855e+09,
   1.506884e+09,
   1.506884e+09,
   1.506905e+09,
   1.506905e+09,
   1.506925e+09,
   1.506925e+09,
   1.506999e+09,
   1.506999e+09,
   1.507059e+09,
   1.507059e+09,
   1.507074e+09,
   1.507074e+09,
   1.507142e+09,
   1.507142e+09,
   1.507196e+09,
   1.507196e+09,
   1.50727e+09,
   1.50727e+09,
   1.507295e+09,
   1.507295e+09,
   1.507341e+09,
   1.507341e+09,
   1.507362e+09,
   1.507362e+09,
   1.507464e+09,
   1.507464e+09,
   1.507494e+09,
   1.507494e+09,
   1.507531e+09,
   1.507531e+09,
   1.507612e+09,
   1.507612e+09,
   1.507671e+09,
   1.507671e+09,
   1.507765e+09,
   1.507765e+09,
   1.507848e+09,
   1.507848e+09,
   1.507957e+09,
   1.507957e+09,
   1.508035e+09,
   1.508035e+09,
   1.508074e+09,
   1.508074e+09,
   1.508148e+09,
   1.508148e+09,
   1.508202e+09,
   1.508202e+09,
   1.508266e+09,
   1.508266e+09,
   1.508308e+09,
   1.508308e+09,
   1.508347e+09,
   1.508347e+09,
   1.5084e+09,
   1.5084e+09,
   1.508421e+09,
   1.508421e+09,
   1.508482e+09,
   1.508482e+09,
   1.508534e+09,
   1.508534e+09,
   1.508548e+09,
   1.508548e+09,
   1.508608e+09,
   1.508608e+09,
   1.508662e+09,
   1.508662e+09,
   1.508735e+09,
   1.508735e+09,
   1.508753e+09,
   1.508753e+09,
   1.508841e+09,
   1.508841e+09,
   1.508907e+09,
   1.508907e+09,
   1.50898e+09,
   1.50898e+09,
   1.509054e+09,
   1.509054e+09,
   1.509085e+09,
   1.509085e+09,
   1.50914e+09,
   1.50914e+09,
   1.509201e+09,
   1.509201e+09,
   1.509233e+09,
   1.509233e+09,
   1.509321e+09,
   1.509321e+09,
   1.509354e+09,
   1.509354e+09,
   1.509395e+09,
   1.509395e+09,
   1.509423e+09,
   1.509423e+09,
   1.509458e+09,
   1.509458e+09,
   1.509559e+09,
   1.509559e+09,
   1.509609e+09,
   1.509609e+09,
   1.509685e+09,
   1.509685e+09,
   1.50979e+09,
   1.50979e+09,
   1.509895e+09,
   1.509895e+09,
   1.510013e+09,
   1.510013e+09,
   1.510214e+09,
   1.510214e+09,
   1.510311e+09,
   1.510311e+09,
   1.511282e+09,
   1.511282e+09,
   1.51133e+09,
   1.51133e+09,
   1.511437e+09,
   1.511437e+09,
   1.511558e+09,
   1.511558e+09,
   1.511638e+09,
   1.511638e+09,
   1.511650e+09,
   1.511660e+09,
   1.511670e+09,
   1.51168e+09};
   Double_t _fy1[432] = {
   0,
   0.19,
   0.19,
   0.317,
   0.317,
   1.128,
   1.128,
   8.661,
   8.661,
   15.216,
   15.216,
   18.862,
   18.862,
   66.469,
   66.469,
   119.962,
   119.962,
   153.018,
   153.018,
   241.583,
   241.583,
   244.524,
   244.524,
   295.308,
   295.308,
   295.396,
   295.396,
   338.545,
   338.545,
   353.203,
   353.203,
   489.399,
   489.399,
   560.983,
   560.983,
   675.427,
   675.427,
   888.473,
   888.473,
   959.715,
   959.715,
   1103.318,
   1103.318,
   1262.499,
   1262.499,
   1524.046,
   1524.046,
   1620.104,
   1620.104,
   1993.014,
   1993.014,
   2511.43,
   2511.43,
   2934.993,
   2934.993,
   2981.667,
   2981.667,
   3195.797,
   3195.797,
   3539.14,
   3539.14,
   3631.864,
   3631.864,
   3640.961,
   3640.961,
   4158.213,
   4158.213,
   4543.639,
   4543.639,
   4936.717,
   4936.717,
   5368.155,
   5368.155,
   5456.459,
   5456.459,
   5862.515,
   5862.515,
   6120.457,
   6120.457,
   6134.038,
   6134.038,
   6322.297,
   6322.297,
   6324.739,
   6324.739,
   6347.775,
   6347.775,
   6399.841,
   6399.841,
   6404.354,
   6404.354,
   6728.477,
   6728.477,
   6785.301,
   6785.301,
   6906.845,
   6906.845,
   7126.5,
   7126.5,
   7307.68,
   7307.68,
   7624.515,
   7624.515,
   7708.545,
   7708.545,
   7804.444,
   7804.444,
   7820.209,
   7820.209,
   7864.362,
   7864.362,
   7880.308,
   7880.308,
   8390.835,
   8390.835,
   8866.472,
   8866.472,
   9005.405,
   9005.405,
   9114.653,
   9114.653,
   9114.708,
   9114.708,
   9114.83,
   9114.83,
   9117.173,
   9117.173,
   9130.674,
   9130.674,
   9185.449,
   9185.449,
   9197.031,
   9197.031,
   9720.493,
   9720.493,
   10236.46,
   10236.46,
   10370.28,
   10370.28,
   10888.36,
   10888.36,
   11315.73,
   11315.73,
   11425.07,
   11425.07,
   11863.55,
   11863.55,
   12353.81,
   12353.81,
   12421.49,
   12421.49,
   12950.91,
   12950.91,
   13449.26,
   13449.26,
   13810.3,
   13810.3,
   14367.31,
   14367.31,
   14416.39,
   14416.39,
   14521.2,
   14521.2,
   15014.05,
   15014.05,
   15440.74,
   15440.74,
   15446.34,
   15446.34,
   15469.35,
   15469.35,
   15573.19,
   15573.19,
   15696.52,
   15696.52,
   15767.94,
   15767.94,
   15915.94,
   15915.94,
   16082.03,
   16082.03,
   16259.45,
   16259.45,
   16330.49,
   16330.49,
   16393.76,
   16393.76,
   16496.96,
   16496.96,
   16792.77,
   16792.77,
   16867.62,
   16867.62,
   17082.02,
   17082.02,
   17189.76,
   17189.76,
   17485.32,
   17485.32,
   17739.35,
   17739.35,
   17807.93,
   17807.93,
   17942.49,
   17942.49,
   18091.53,
   18091.53,
   18113.61,
   18113.61,
   18204.94,
   18204.94,
   18283.53,
   18283.53,
   18588.74,
   18588.74,
   18622.12,
   18622.12,
   18931.54,
   18931.54,
   19228.94,
   19228.94,
   19239.92,
   19239.92,
   19554.78,
   19554.78,
   19677.7,
   19677.7,
   19695.92,
   19695.92,
   19862.51,
   19862.51,
   19984.34,
   19984.34,
   20096.27,
   20096.27,
   20243.53,
   20243.53,
   20497.9,
   20497.9,
   20664.86,
   20664.86,
   20766.34,
   20766.34,
   20866.31,
   20866.31,
   21024.35,
   21024.35,
   21159.67,
   21159.67,
   21403.35,
   21403.35,
   21721.95,
   21721.95,
   22041.59,
   22041.59,
   22115.22,
   22115.22,
   22454.27,
   22454.27,
   22470.3,
   22470.3,
   22485.02,
   22485.02,
   22751.42,
   22751.42,
   22834.33,
   22834.33,
   23149.18,
   23149.18,
   23487.07,
   23487.07,
   23813.62,
   23813.62,
   23923.96,
   23923.96,
   24088.13,
   24088.13,
   24103.13,
   24103.13,
   24106.5,
   24106.5,
   24123.14,
   24123.14,
   24257.83,
   24257.83,
   24599.35,
   24599.35,
   25007.39,
   25007.39,
   25434.2,
   25434.2,
   25847.65,
   25847.65,
   26126.38,
   26126.38,
   26530.17,
   26530.17,
   26994.8,
   26994.8,
   27568.64,
   27568.64,
   27940.99,
   27940.99,
   28378.87,
   28378.87,
   28534.46,
   28534.46,
   28702.3,
   28702.3,
   28838.41,
   28838.41,
   29232.16,
   29232.16,
   29623.28,
   29623.28,
   29696.94,
   29696.94,
   30141.47,
   30141.47,
   30579.94,
   30579.94,
   31049.43,
   31049.43,
   31292,
   31292,
   31567.49,
   31567.49,
   31731.48,
   31731.48,
   32186.25,
   32186.25,
   32492.44,
   32492.44,
   32837.5,
   32837.5,
   33345.05,
   33345.05,
   33815.93,
   33815.93,
   34411.25,
   34411.25,
   34411.25,
   34411.25,
   34921.22,
   34921.22,
   35525.87,
   35525.87,
   35866.56,
   35866.56,
   36399.59,
   36399.59,
   36904.62,
   36904.62,
   37426.91,
   37426.91,
   37852.36,
   37852.36,
   38250.64,
   38250.64,
   38782.06,
   38782.06,
   38831.42,
   38831.42,
   39363.12,
   39363.12,
   39842.44,
   39842.44,
   39912.54,
   39912.54,
   40469.67,
   40469.67,
   40945.76,
   40945.76,
   41471.19,
   41471.19,
   41609.58,
   41609.58,
   42123.77,
   42123.77,
   42724.84,
   42724.84,
   43286.48,
   43286.48,
   43288.44,
   43288.44,
   43557.19,
   43557.19,
   43608.46,
   43608.46,
   44090.02,
   44090.02,
   44350.21,
   44350.21,
   44966.8,
   44966.8,
   45184.79,
   45184.79,
   45429.82,
   45429.82,
   45595.28,
   45595.28,
   45648.35,
   45648.35,
   45699.71,
   45699.71,
   46198.79,
   46198.79,
   46771.5,
   46771.5,
   47495.19,
   47495.19,
   48220.45,
   48220.45,
   48987.28,
   48987.28,
   49407.82,
   49407.82,
   50112.41,
   50112.41,
   50123.92,
   50123.92,
   50161.9,
   50161.9,
   50194.61,
   50194.61,
   50248.79,
   50248.79,
   50281.55,
   50281.55,
   50293.83,
   50293.83,
   0};
   TGraph *graph = new TGraph(432,_fx1,_fy1);
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
   graph->GetXaxis()->SetLimits(x[0],x[92]);
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
grA->GetYaxis()->SetTitle("RMS of Drift [%]");
/*grA->GetXaxis()->SetLabelOffset(999);
grA->GetXaxis()->SetLabelSize(0);
grA->GetXaxis()->SetTickSize(0);*/
grA->GetXaxis()->SetLimits(x[0],x[92]);
grA->GetXaxis()->SetTimeDisplay(1);
grA->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%Y}%F1970-01-01 00:00:00s0");
grA->GetXaxis()->SetNdivisions(509);
grA->GetXaxis()->SetTimeDisplay(1);
grA->GetXaxis()->SetTimeOffset(0,"gmt");
grA->GetXaxis()->SetLabelOffset(0.013);
grA->GetXaxis()->SetLabelSize(0.02);
grA->GetXaxis()->SetTitleSize(0.035);
grA->GetXaxis()->SetTitleOffset(1.36);
cout<<x[92];
grA->SetLineColor(kRed);
grB->SetLineColor(kBlue);
grD->SetLineColor(kGreen);
grA->Draw("AY+PE");
grA->GetYaxis()->SetRangeUser(0.,14.);
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
    printInfo.DrawLatex(0.15, 0.75, "#scale[0.9]{Laser Combined Method 2017}");




  
   
   

   


   
    
    
 } 

