
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




//Reading files

double Las_constant[4][64][48];
double Cesium_constant[4][64][48];
double Cs_Las[4][64][48];
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

void phi()
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
   ifstream f_Laser1;
   ifstream f_Laser2;
   ifstream f_Cesium1;
   ifstream f_Cesium2;
   
   ifstream f_Masked;
    

    //--------- new update

 /* f_Laser1.open("LaserCombined_267534.txt");                                                                          //IOV1
    f_Laser2.open("LaserCombined_272493.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_263962_11jun2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_270000_17jul2015.txt");
*/

/*   f_Laser1.open("LaserCombined_272493.txt");                                                                          //IOV2
    f_Laser2.open("LaserCombined_277320.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_270000_17jul2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_277321_27aug2015.txt");
*/

 f_Laser1.open("LaserCombined_277320.txt");                                                                        //IOV3
    f_Laser2.open("LaserCombined_284682.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_277321_27aug2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_284600_03nov2015.txt");

   



      //f_Masked.open("masked_channels_267534.txt");      //IOV1
      //f_Masked.open("masked_channels_272493.txt");      //IOV2
     f_Masked.open("masked_channels_277320.txt");      //IOV3

    
         while(!f_Laser1.eof())
           {
            f_Laser1 >> tmp_part;
	    f_Laser1 >> tmp_channel;
            mask = tmp_part + " " + std::to_string(tmp_channel);
            separatePartMod(tmp_part, tmp_module);

            f_Laser1 >> tmp_zero;
            f_Laser1 >> tmp_constant;
            f_Laser1 >> tmp_HV;
            
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
                 break;
               }
              else{ 
             double drift1= (tmp_constant-1.)*100;
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift1;
                 }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Laser1.close();
f_Masked.close();


    //f_Masked.open("masked_channels_272493.txt");     //IOV1
   //f_Masked.open("masked_channels_277320.txt");     //IOV2
    f_Masked.open("masked_channels_284682.txt");     //IOV3
 


while(!f_Laser2.eof())
           {
            f_Laser2 >> tmp_part;
	    f_Laser2 >> tmp_channel;
            mask = tmp_part + " " + std::to_string(tmp_channel);
            separatePartMod(tmp_part, tmp_module);

            f_Laser2 >> tmp_zero;
            f_Laser2 >> tmp_constant;
            f_Laser2 >> tmp_HV;
            
            int myPart_ctr = -1;
            if (tmp_part == "EBA") myPart_ctr = Part::partEBA;
            if (tmp_part == "LBA") myPart_ctr = Part::partLBA;
            if (tmp_part == "LBC") myPart_ctr = Part::partLBC;
            if (tmp_part == "EBC") myPart_ctr = Part::partEBC;
            double a = Las_constant[myPart_ctr][tmp_module-1][tmp_channel];
           
           while(std::getline(f_Masked,line))
            {

             if (line.compare(mask) ==0) 
               {  
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
                 break;
               }
              else{ 
             double drift2= (tmp_constant-1.)*100;
             
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift2 - a;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Laser2.close();
f_Masked.close();


    //f_Masked.open("masked_channels_267534.txt");      //IOV1
    //f_Masked.open("masked_channels_272493.txt");      //IOV2
   f_Masked.open("masked_channels_277320.txt");      //IOV3


 while(!f_Cesium1.eof())
           {
            f_Cesium1 >> tmp_part;
	    f_Cesium1 >> tmp_channel;
            mask = tmp_part + " " + std::to_string(tmp_channel);
            separatePartMod(tmp_part, tmp_module);

            f_Cesium1 >> tmp_zero;
            f_Cesium1 >> tmp_constant;
            f_Cesium1 >> tmp_HV;
            
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
                 
                 Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
                 break;
               }
              else{ 
             double drift3= (tmp_constant-1.)*100;
           Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift3;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Cesium1.close();
f_Masked.close();


    //f_Masked.open("masked_channels_272493.txt");     //IOV1
    //f_Masked.open("masked_channels_277320.txt");     //IOV2
    f_Masked.open("masked_channels_284682.txt");     //IOV3
 
 while(!f_Cesium2.eof())
           {
            f_Cesium2 >> tmp_part;
	    f_Cesium2 >> tmp_channel;
            mask = tmp_part + " " + std::to_string(tmp_channel);
            separatePartMod(tmp_part, tmp_module);

            f_Cesium2 >> tmp_zero;
            f_Cesium2 >> tmp_constant;
            f_Cesium2 >> tmp_HV;
            
            int myPart_ctr = -1;
            if (tmp_part == "EBA") myPart_ctr = Part::partEBA;
            if (tmp_part == "LBA") myPart_ctr = Part::partLBA;
            if (tmp_part == "LBC") myPart_ctr = Part::partLBC;
            if (tmp_part == "EBC") myPart_ctr = Part::partEBC;
            double c = Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel];
            double d = Las_constant[myPart_ctr][tmp_module-1][tmp_channel];

          while(std::getline(f_Masked,line))
            {
              //cout<<line<<"\n";
             if (line.compare(mask) ==0) 
               {  
                 
                 Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
                 Cs_Las[myPart_ctr][tmp_module-1][tmp_channel] = -1000.;
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000.;
                 break;
               }
              else{ 
             double drift4= (tmp_constant-1.)*100;
           Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift4 - c;
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = (Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] -d)   ;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Cesium2.close();


f_Masked.close();
    
    // Fill histograms
    
    
TH1F* A = new TH1F("A", "Channel A",50,-100.,100.);
TH1F* B = new TH1F("B", "Channel B",50,-100.,100.);
TH1F* D = new TH1F("D", "Channel D",50,-100.,100.);

double meanA[64];
double meanB[64];
double meanD[64];
double rmsA[64];
double rmsB[64];
double rmsD[64];
double module[64];

for(int i=0;i<64;i++)
{
module[i]=i+1;
}


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

 meanA[m-1] = A->GetMean();
 meanB[m-1] = B->GetMean();
 meanD[m-1] = D->GetMean();

 rmsA[m-1] = A->GetMeanError();
 rmsB[m-1] = B->GetMeanError();
 rmsD[m-1] = D->GetMeanError();    

A->Reset("ICES");
B->Reset("ICES");
D->Reset("ICES");




}//loop over modules
    
    
    
 auto myCanvas =  new TCanvas("myCanvas", "",484,86,800,750);
 





TGraphErrors* grA = new TGraphErrors(64,module,meanA,nullptr,rmsA); 


TGraphErrors* grB = new TGraphErrors(64,module,meanB,nullptr,rmsB); 
TGraphErrors* grD = new TGraphErrors(64,module,meanD,nullptr,rmsD); 
grA -> SetMarkerStyle(21);
grA ->SetMarkerColor(kRed);
grA-> SetMarkerSize(1.5);
grA->SetLineColor(kRed);

grB -> SetMarkerStyle(22);
grB ->SetMarkerColor(kBlue);
grB-> SetMarkerSize(1.5);
grB->SetLineColor(kBlue);

grD -> SetMarkerStyle(20);
grD ->SetMarkerColor(kGreen);
grD-> SetMarkerSize(1.5);
grD->SetLineColor(kGreen);

grA ->Draw("APE");

//grA->SetTitle("; #phi(Module number); Mean of (Cs - Las) Drift (%) [IOV1]");
//grA->SetTitle("; #phi(Module number); Mean of (Cs - Las) Drift (%) [IOV2]");
grA->SetTitle("; #phi(Module number); Mean of (Cs - Las) Drift (%) [IOV3]");

grA->GetYaxis()->SetRangeUser(-2.,2.);
grB ->Draw("PE");
grD ->Draw("PE");




   
    TLegend *leg = new TLegend(0.7401669,0.7092683,0.9398093,0.9092683,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetHeader("Channels","C"); // option "C" allows to center the header
   leg->AddEntry(grA,"Channel A","p");
   leg->AddEntry(grB,"Channel BC","p");
   leg->AddEntry(grD,"Channel D","p");
   
  leg->Draw();
myCanvas -> Draw();



/*A[i] -> SetLineColor(1);
A[i]-> SetLineWidth(3);
A[i]-> GetXaxis() -> SetTitleOffset(1.2);
A[i]-> GetYaxis() -> SetTitleOffset(1.2);
A[i]-> GetXaxis() ->  SetTitleFont(62);    //for setting title font
A[i]-> GetYaxis() ->  SetTitleFont(62);
   
A[i]-> GetXaxis() -> SetTitleColor(1);
A[i]-> GetXaxis() ->  SetLabelFont(62);    //for setting Label font
A[i]-> GetYaxis() ->  SetLabelFont(62);
A[i]-> GetXaxis() ->  SetLabelSize(0.04);     
A[i]-> GetYaxis() ->  SetLabelSize(0.04);
A[i]-> GetXaxis() ->  SetTitleSize(0.04);
A[i]-> GetYaxis() ->  SetTitleSize(0.04);
A[i]-> GetYaxis() -> SetTitleColor(1);
   
A[i]-> GetYaxis() -> SetTitle("Entries");    

A[i]-> GetXaxis() -> SetTitle("Drift (%)");    
    





    TLegend *legendTopLeft = new TLegend(0.09,0.9,0.95,0.95, "");
    legendTopLeft-> SetNColumns(2);
    legendTopLeft->SetFillStyle(0);
    legendTopLeft->SetFillColor(0);
    legendTopLeft->SetBorderSize(0);
    A[i] ->Draw();*/
 TLatex printInfo;
    printInfo.SetNDC();
    printInfo.SetTextFont(72);
    printInfo.DrawLatex(0.15, 0.85, "#scale[0.5]{ATLAS}");
    printInfo.SetTextFont(42);
    printInfo.DrawLatex(0.25, 0.85, "#scale[0.5]{Internal}");
    printInfo.DrawLatex(0.15, 0.80, "#scale[0.5]{Tile Calorimeter}");
    printInfo.DrawLatex(0.15, 0.75, "#scale[0.5]{2015 Data #sqrt{s} = 13 TeV}");  


  
    
    
 } 

