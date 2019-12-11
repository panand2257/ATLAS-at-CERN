
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


TH1F* A1 = new TH1F("BC1", "Channel BC_{1}",50,-0.5,8.);

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

void cha()
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
   ifstream f_Laser;
   ifstream f_Masked;
    

    //--------- new update
    f_Laser.open("run.txt");
    f_Masked.open("masked_channels_313662_28nov2016.txt");
    
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
              //cout<<line<<"\n";
             if (line.compare(mask) ==0) 
               {  
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -100. ;
                 break;
               }
              else{ 
             double drift= (1-tmp_constant)*100;
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Laser.close();
    f_Masked.close();
    
    // Fill histograms

    
    for(int m = 1; m < 65; m++)
    {

           A1 -> Fill(Las_constant[1][m-1][7]);
           A1 -> Fill(Las_constant[2][m-1][7]);     
       
    }//loop over modules

    
    auto myCanvas = new TCanvas("myCanvas", "",484,86,800,750);
    
A1 -> SetLineColor(1);
A1-> SetLineWidth(3);
A1-> GetXaxis() -> SetTitleOffset(1.2);
A1-> GetYaxis() -> SetTitleOffset(1.2);
A1-> GetXaxis() ->  SetTitleFont(62);    //for setting title font
A1-> GetYaxis() ->  SetTitleFont(62);
   
A1-> GetXaxis() -> SetTitleColor(1);
A1-> GetXaxis() ->  SetLabelFont(62);    //for setting Label font
A1-> GetYaxis() ->  SetLabelFont(62);
A1-> GetXaxis() ->  SetLabelSize(0.04);     
A1-> GetYaxis() ->  SetLabelSize(0.04);
A1-> GetXaxis() ->  SetTitleSize(0.04);
A1-> GetYaxis() ->  SetTitleSize(0.04);
A1-> GetYaxis() -> SetTitleColor(1);
   
A1-> GetYaxis() -> SetTitle("Entries");    

A1-> GetXaxis() -> SetTitle("Drift (%)");    
    





    TLegend *legendTopLeft = new TLegend(0.09,0.9,0.95,0.95, "");
    legendTopLeft-> SetNColumns(2);
    legendTopLeft->SetFillStyle(0);
    legendTopLeft->SetFillColor(0);
    legendTopLeft->SetBorderSize(0);
    A1 ->Draw();
    myCanvas -> Draw();
   


    
    
    
 } 

