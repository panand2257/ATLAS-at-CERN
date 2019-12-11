
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

void pmteta()
{   //-------- read the data
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
//f_Laser.open("LaserCombined_314511.txt");               //start
//f_Masked.open("masked_channels_314511.txt");

//f_Laser.open("LaserCombined_325677.txt");               //middle
//f_Masked.open("masked_channels_325677.txt");

f_Laser.open("LaserCombined_341673.txt");              //end
f_Masked.open("masked_channels_341673.txt");
    int j=0;
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
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -100.-100*j ;
                 j--;
                 break;
               }
              else{ 
             double drift= (1-tmp_constant)*100;
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = tmp_constant;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Laser.close();
    f_Masked.close();
    // Fill histograms
    
    
       TH1F* A[30];  
         A[0] =   new TH1F("A1", "Channel A_{1}",50,-5.,5.);
         A[1] =   new TH1F("A2", "Channel A_{2}",50,-5.,5.);
         A[2] =   new TH1F("A3", "Channel A_{3}",50,-5.,5.);
         A[3] =   new TH1F("A4", "Channel A_{4}",50,-5.,5.);
         A[4] =   new TH1F("A5", "Channel A_{5}",50,-5.,5.);
         A[5] =   new TH1F("A6", "Channel A_{6}",50,-5.,5.);
         A[6] =   new TH1F("A7", "Channel A_{7}",50,-5.,5.);
         A[7] =   new TH1F("A8", "Channel A_{8}",50,-5.,5.);
         A[8] =   new TH1F("A9", "Channel A_{9}",50,-5.,5.);
         A[9] =   new TH1F("A10", "Channel A_{10}",50,-5.,5.);
         A[10] =   new TH1F("A12", "Channel A_{12}",50,-5.,5.);
         A[11] =   new TH1F("A12", "Channel A_{13}",50,-5.,5.);
         A[12] =   new TH1F("A13", "Channel A_{14}",50,-5.,5.);
         A[13] =   new TH1F("A14", "Channel A_{15}",50,-5.,5.); 
         A[14] =   new TH1F("A10", "Channel A_{16}",50,-5.,5.);
         A[15] =   new TH1F("A-1", "Channel A_{-1}",50,-5.,5.);
         A[16] =   new TH1F("A-2", "Channel A_{-2}",50,-5.,5.);
         A[17] =   new TH1F("A-3", "Channel A_{-3}",50,-5.,5.);
         A[18] =   new TH1F("A-4", "Channel A_{-4}",50,-5.,5.);
         A[19] =   new TH1F("A-5", "Channel A_{-5}",50,-5.,5.);
         A[20] =   new TH1F("A-6", "Channel A_{-6}",50,-5.,5.);
         A[21] =   new TH1F("A-7", "Channel A_{-7}",50,-5.,5.);
         A[22] =   new TH1F("A-8", "Channel A_{-8}",50,-5.,5.);
         A[23] =   new TH1F("A-9", "Channel A_{-9}",50,-5.,5.);
         A[24] =   new TH1F("A-10", "Channel A_{-10}",50,-5.,5.);
         A[25] =   new TH1F("A-12", "Channel A_{-12}",50,-5.,5.);
         A[26] =   new TH1F("A-13", "Channel A_{-13}",50,-5.,5.);
         A[27] =   new TH1F("A-14", "Channel A_{-14}",50,-5.,5.);
         A[28] =   new TH1F("A-15", "Channel A_{-15}",50,-5.,5.);
         A[29] =   new TH1F("A-16", "Channel A_{-16}",50,-5.,5.);

          TH1F* B[28];  
         B[0] =   new TH1F("B1", "Channel B_{1}",50,-5.,5.);
         B[1] =   new TH1F("B2", "Channel B_{2}",50,-5.,5.);
         B[2] =   new TH1F("B3", "Channel B_{3}",50,-5.,5.);
         B[3] =   new TH1F("B4", "Channel B_{4}",50,-5.,5.);
         B[4] =   new TH1F("B5", "Channel B_{5}",50,-5.,5.);
         B[5] =   new TH1F("B6", "Channel B_{6}",50,-5.,5.);
         B[6] =   new TH1F("B7", "Channel B_{7}",50,-5.,5.);
         B[7] =   new TH1F("B8", "Channel B_{8}",50,-5.,5.);
         B[8] =   new TH1F("B9", "Channel B_{9}",50,-5.,5.);
         B[9] =   new TH1F("B11", "Channel B_{11}",50,-5.,5.);
         B[10] =   new TH1F("B12", "Channel B_{12}",50,-5.,5.);
         B[11] =   new TH1F("B13", "Channel B_{13}",50,-5.,5.);
         B[12] =   new TH1F("B14", "Channel B_{14}",50,-5.,5.);
         B[13] =   new TH1F("B15", "Channel B_{15}",50,-5.,5.); 
         B[14] =   new TH1F("B-1", "Channel B_{-1}",50,-5.,5.);
         B[15] =   new TH1F("B-2", "Channel B_{-2}",50,-5.,5.);
         B[16] =   new TH1F("B-3", "Channel A_{-3}",50,-5.,5.);
         B[17] =   new TH1F("B-4", "Channel B_{-4}",50,-5.,5.);
         B[18] =   new TH1F("B-5", "Channel B_{-5}",50,-5.,5.);
         B[19] =   new TH1F("B-6", "Channel B_{-6}",50,-5.,5.);
         B[20] =   new TH1F("B-7", "Channel B_{-7}",50,-5.,5.);
         B[21] =   new TH1F("B-8", "Channel B_{-8}",50,-5.,5.);
         B[22] =   new TH1F("B-9", "Channel B_{-9}",50,-5.,5.);
         B[23] =   new TH1F("B-11", "Channel B_{-11}",50,-5.,5.);
         B[24] =   new TH1F("B-12", "Channel B_{-12}",50,-5.,5.);
         B[25] =   new TH1F("B-13", "Channel B_{-13}",50,-5.,5.);
         B[26] =   new TH1F("B-14", "Channel B_{-14}",50,-5.,5.);
         B[27] =   new TH1F("B-15", "Channel B_{-15}",50,-5.,5.);
         
         TH1F* D[14];  
         D[0] =   new TH1F("D0", "Channel D_{0}",50,-5.,5.);
         D[1] =   new TH1F("D1", "Channel D_{1}",50,-5.,5.);
         D[2] =   new TH1F("D2", "Channel D_{2}",50,-5.,5.);
         D[3] =   new TH1F("D3", "Channel D_{3}",50,-5.,5.);
         D[4] =   new TH1F("D4", "Channel D_{4}",50,-5.,5.);
         D[5] =   new TH1F("D5", "Channel D_{5}",50,-5.,5.);
         D[6] =   new TH1F("D6", "Channel D_{6}",50,-5.,5.);
         D[7] =   new TH1F("D-1", "Channel D_{-1}",50,-5.,5.);
         D[8] =   new TH1F("D-2", "Channel D_{-2}",50,-5.,5.);
         D[9] =   new TH1F("D-3", "Channel D_{-3}",50,-5.,5.);
         D[10] =   new TH1F("D-4", "Channel D_{-4}",50,-5.,5.);
         D[11] =   new TH1F("D-5", "Channel D_{-5}",50,-5.,5.);
         D[12] =   new TH1F("D-6", "Channel D_{-6}",50,-5.,5.);
         D[13] =   new TH1F("D-6", "Channel D_{-6}",50,-5.,5.);








    for(int m = 1; m < 65; m++)
    {



          // Left Modules
           A[0] -> Fill(Las_constant[1][m-1][1]-Las_constant[1][m-1][4]);        //A1
           A[1] -> Fill(Las_constant[1][m-1][5]-Las_constant[1][m-1][8]);        //A2
           A[2] -> Fill(Las_constant[1][m-1][9]-Las_constant[1][m-1][10]);        //A3
           A[3] -> Fill(Las_constant[1][m-1][15]-Las_constant[1][m-1][18]);       //A4
           A[4] -> Fill(Las_constant[1][m-1][19]-Las_constant[1][m-1][20]);       //A5
           A[5] -> Fill(Las_constant[1][m-1][23]-Las_constant[1][m-1][26]);       //A6
           A[6] -> Fill(Las_constant[1][m-1][29]-Las_constant[1][m-1][32]);       //A7
           A[7] -> Fill(Las_constant[1][m-1][35]-Las_constant[1][m-1][38]);       //A8
           A[8] -> Fill(Las_constant[1][m-1][37]-Las_constant[1][m-1][36]);       //A9
           A[9] -> Fill(Las_constant[1][m-1][45]-Las_constant[1][m-1][46]);       //A10
           A[10] -> Fill(Las_constant[0][m-1][7]-Las_constant[0][m-1][6]);         //A12
           A[11] -> Fill(Las_constant[0][m-1][11]-Las_constant[0][m-1][10]);         //A13
           A[12] -> Fill(Las_constant[0][m-1][21]-Las_constant[0][m-1][20]);        //A14
           A[13] -> Fill(Las_constant[0][m-1][32]-Las_constant[0][m-1][31]);        //A15
           A[14] -> Fill(Las_constant[0][m-1][40]-Las_constant[0][m-1][41]);        //A16




           A[15] -> Fill(Las_constant[2][m-1][1]-Las_constant[2][m-1][4]);
            

           A[16] -> Fill(Las_constant[2][m-1][5]-Las_constant[2][m-1][8]);
            

           A[17] -> Fill(Las_constant[2][m-1][9]-Las_constant[2][m-1][10]);
            

           A[18] -> Fill(Las_constant[2][m-1][15]-Las_constant[2][m-1][18]);
            

           A[19] -> Fill(Las_constant[2][m-1][19]-Las_constant[2][m-1][20]);
             

           A[20] -> Fill(Las_constant[2][m-1][23]-Las_constant[2][m-1][26]);
            

           A[21] -> Fill(Las_constant[2][m-1][29]-Las_constant[2][m-1][32]);
            

           A[22] -> Fill(Las_constant[2][m-1][35]-Las_constant[2][m-1][38]);
            

           A[23] -> Fill(Las_constant[2][m-1][37]-Las_constant[2][m-1][36]);
            

           A[24] -> Fill(Las_constant[2][m-1][45]-Las_constant[2][m-1][46]);
            

           A[25] -> Fill(Las_constant[3][m-1][7]-Las_constant[3][m-1][6]);


           A[26] -> Fill(Las_constant[3][m-1][11]-Las_constant[3][m-1][10]);
    

           A[27] -> Fill(Las_constant[3][m-1][21]-Las_constant[3][m-1][20]);


           A[28] -> Fill(Las_constant[3][m-1][32]-Las_constant[3][m-1][31]);
 

           A[29] -> Fill(Las_constant[3][m-1][40]-Las_constant[3][m-1][41]);
 
          


          // Left Modules
           B[0] -> Fill(Las_constant[1][m-1][3]-Las_constant[1][m-1][2]);         //B1
           B[14] -> Fill(Las_constant[2][m-1][3]-Las_constant[2][m-1][2]);
            
           B[1] -> Fill(Las_constant[1][m-1][7]-Las_constant[1][m-1][6]);         //B2
           B[15] -> Fill(Las_constant[2][m-1][7]-Las_constant[2][m-1][6]);
            
           B[2] -> Fill(Las_constant[1][m-1][11]-Las_constant[1][m-1][12]);        //B3
           B[16] -> Fill(Las_constant[2][m-1][11]-Las_constant[2][m-1][12]);
            
           B[3] -> Fill(Las_constant[1][m-1][17]-Las_constant[1][m-1][16]);        //B4
           B[17] -> Fill(Las_constant[2][m-1][17]-Las_constant[2][m-1][16]);
            
           B[4] -> Fill(Las_constant[1][m-1][21]-Las_constant[1][m-1][22]);        //B5
           B[18] -> Fill(Las_constant[2][m-1][21]-Las_constant[2][m-1][22]);
           
           B[5] -> Fill(Las_constant[1][m-1][27]-Las_constant[1][m-1][28]);       //B6
           B[19] -> Fill(Las_constant[2][m-1][27]-Las_constant[2][m-1][28]);
            
           B[6] -> Fill(Las_constant[1][m-1][33]-Las_constant[1][m-1][34]);        //B7
           B[20] -> Fill(Las_constant[2][m-1][33]-Las_constant[2][m-1][34]);
            
           B[7] -> Fill(Las_constant[1][m-1][39]-Las_constant[1][m-1][40]);       //B8
           B[21] -> Fill(Las_constant[2][m-1][39]-Las_constant[2][m-1][40]);
            
           B[8] -> Fill(Las_constant[1][m-1][47]-Las_constant[1][m-1][42]);        //B9
           B[22] -> Fill(Las_constant[2][m-1][47]-Las_constant[2][m-1][42]);
            
           B[9] -> Fill(Las_constant[0][m-1][9]-Las_constant[0][m-1][8]);           //B11
           B[23] -> Fill(Las_constant[3][m-1][9]-Las_constant[3][m-1][8]);
            
           B[10] -> Fill(Las_constant[0][m-1][15]-Las_constant[0][m-1][14]);         //B12
           B[24] -> Fill(Las_constant[3][m-1][15]-Las_constant[3][m-1][14]);

           B[11] -> Fill(Las_constant[0][m-1][23]-Las_constant[0][m-1][22]);        //B13
           B[25] -> Fill(Las_constant[3][m-1][23]-Las_constant[3][m-1][22]);
    
           B[12] -> Fill(Las_constant[0][m-1][35]-Las_constant[0][m-1][30]);        //B14
           B[26] -> Fill(Las_constant[3][m-1][35]-Las_constant[3][m-1][30]);

           B[13] -> Fill(Las_constant[0][m-1][36]-Las_constant[0][m-1][39]);        //B15
           B[27] -> Fill(Las_constant[3][m-1][36]-Las_constant[3][m-1][39]);
 


          
 
           D[0] -> Fill(Las_constant[1][m-1][0]-Las_constant[1][m-1][0]);         //D0
           D[7] -> Fill(Las_constant[2][m-1][0]-Las_constant[2][m-1][0]);        
            
           D[1] -> Fill(Las_constant[1][m-1][13]-Las_constant[1][m-1][14]);         //D1
           D[8] -> Fill(Las_constant[2][m-1][13]-Las_constant[2][m-1][14]);
            
           D[2] -> Fill(Las_constant[1][m-1][25]-Las_constant[1][m-1][24]);        //D2
           D[9] -> Fill(Las_constant[2][m-1][25]-Las_constant[2][m-1][24]);
            
           D[3] -> Fill(Las_constant[1][m-1][41]-Las_constant[1][m-1][44]);        //D3
           D[10] -> Fill(Las_constant[2][m-1][41]-Las_constant[2][m-1][44]);
            
           D[4] -> Fill(Las_constant[0][m-1][3]-Las_constant[0][m-1][2]);        //D4
           D[11] -> Fill(Las_constant[3][m-1][3]-Las_constant[3][m-1][2]);
           
           D[5] -> Fill(Las_constant[0][m-1][17]-Las_constant[0][m-1][16]);       //D5
           D[12] -> Fill(Las_constant[3][m-1][17]-Las_constant[3][m-1][16]);
            
           D[6] -> Fill(Las_constant[0][m-1][37]-Las_constant[0][m-1][38]);        //D6
           D[13] -> Fill(Las_constant[3][m-1][37]-Las_constant[3][m-1][38]);
            
         




}//loop over modules
    
    
    
 auto myCanvas =  new TCanvas("myCanvas", "",484,86,800,750);
 
double erreta[30];
double erretaD[13];

double meanA[30];
double etaA[30];
double rmsA[30];

double meanB[28];
double etaB[28];
double rmsB[28];

double meanD[14];
double etaD[14] = {0.,0.2,0.4,0.6,0.8,1.0,1.245,0.,-0.2,-0.4,-0.6,-0.8,-1.0,-1.245};
double rmsD[14];

for(int i=0;i<30;i++){
 

meanA[i] = (A[i]->GetRMS());
rmsA[i] = A[i]->GetRMSError();
erreta[i] = 0.05;
if(i<10)
etaA[i] = (0.05 + 0.1*i);

if(i>9 && i<15)
etaA[i] = (0.15 +0.1*i);

if(i>14)
etaA[i] = -1*etaA[i-15];

}

for(int i=0;i<28;i++){
 

meanB[i] = (B[i]->GetRMS());
rmsB[i] = B[i]->GetRMSError();

if(i<9)
etaB[i] = (0.05 + 0.1*i);

if(i>8 && i<14)
etaB[i] = (0.15 +0.1*i);

if(i>13)
etaB[i] = -1*etaB[i-14];

}


for(int i=0;i<14;i++){
 

meanD[i] = (D[i]->GetRMS());
rmsD[i] = D[i]->GetRMSError();
erretaD[i] = 0.1;
}



TGraphErrors* grA = new TGraphErrors(30,etaA,meanA,erreta,rmsA); 


TGraphErrors* grB = new TGraphErrors(28,etaB,meanB,erreta,rmsB); 
TGraphErrors* grD = new TGraphErrors(14,etaD,meanD,erretaD,rmsD); 
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
//grA->SetTitle("; #eta; RMS of #DeltaPMT [Run 314511]");          //start
//grA->SetTitle("; #eta; RMS of #DeltaPMT [Run 325677]");          //middle
grA->SetTitle("; #eta; RMS of #DeltaPMT [Run 341673]");          //end
//grA->GetYaxis()->SetRangeUser(0.,0.07);
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
    printInfo.DrawLatex(0.15, 0.75, "#scale[0.5]{2017 Data #sqrt{s} = 13 TeV}");
   


  
    
    
 } 

