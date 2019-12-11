
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

void scat1()
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

 /*  f_Laser1.open("LaserCombined_267534.txt");                                                                          //IOV1
    f_Laser2.open("LaserCombined_272493.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_263962_11jun2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_270000_17jul2015.txt");
*/

 /*  f_Laser1.open("LaserCombined_272493.txt");                                                                          //IOV2
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

                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Cesium2.close();


f_Masked.close();
    // Fill histograms
    
    
TH2F* A9 = new TH2F("Channel A_{9}", "IOV3",50,-5.,2.,50,-5.,2.);
TH2F* A10 = new TH2F("Channel A_{10}", "",50,-3.,2.,50,-3.,2.);
TH2F* A12 = new TH2F("Channel A_{12}", "",50,-3.,2.,50,-3.,2.);
TH2F* A13 = new TH2F("Channel A_{13}", "",50,-3.,2.,50,-3.,2.);
TH2F* A14 = new TH2F("Channel A_{14}", "",50,-3.,2.,50,-3.,2.);


    for(int m = 1; m < 65; m++)
    {



         
            
           A9 -> Fill(Las_constant[1][m-1][36],Cesium_constant[1][m-1][36]);       //A9
           A9 -> Fill(Las_constant[1][m-1][37],Cesium_constant[1][m-1][37]);
            
           A10 -> Fill(Las_constant[1][m-1][45],Cesium_constant[1][m-1][45]);       //A10
           A10 -> Fill(Las_constant[1][m-1][46],Cesium_constant[1][m-1][46]);
            
           A12 -> Fill(Las_constant[0][m-1][6],Cesium_constant[0][m-1][6]);         //A12
           A12 -> Fill(Las_constant[0][m-1][7],Cesium_constant[0][m-1][7]);

           A13 -> Fill(Las_constant[0][m-1][10],Cesium_constant[0][m-1][10]);         //A13
           A13 -> Fill(Las_constant[0][m-1][11],Cesium_constant[0][m-1][11]);
    
           A14 -> Fill(Las_constant[0][m-1][20],Cesium_constant[0][m-1][20]);        //A14
           A14 -> Fill(Las_constant[0][m-1][21],Cesium_constant[0][m-1][21]);

       
           // Right Modules

            
           A9 -> Fill(Las_constant[2][m-1][36],Cesium_constant[2][m-1][36]);       //A-9
           A9 -> Fill(Las_constant[2][m-1][37],Cesium_constant[2][m-1][37]);
            
           A10 -> Fill(Las_constant[2][m-1][45],Cesium_constant[2][m-1][45]);       //A-10
           A10 -> Fill(Las_constant[2][m-1][46],Cesium_constant[2][m-1][46]);
            
           A12 -> Fill(Las_constant[3][m-1][6],Cesium_constant[3][m-1][6]);         //A-12
           A12 -> Fill(Las_constant[3][m-1][7],Cesium_constant[3][m-1][7]);

           A13 -> Fill(Las_constant[3][m-1][10],Cesium_constant[3][m-1][10]);         //A-13
           A13 -> Fill(Las_constant[3][m-1][11],Cesium_constant[3][m-1][11]);
    
           A14 -> Fill(Las_constant[3][m-1][20],Cesium_constant[3][m-1][20]);        //A-14
           A14 -> Fill(Las_constant[3][m-1][21],Cesium_constant[3][m-1][21]);

         




        
}//loop over modules
    
    
    
 auto myCanvas =  new TCanvas("myCanvas", "",484,86,800,750);
 
A9->GetXaxis()->SetTitle("Laser Drift (%)");
A9->GetYaxis()->SetTitle("Cesium Drift (%)");
A9->SetMarkerColor(kRed);
A9->SetMarkerStyle(31);
A9->SetMarkerSize(1.5);
A9->Draw("SCAT");

A10->SetMarkerColor(kBlue);
A10->SetMarkerStyle(31);
A10->SetMarkerSize(1.5);
A10->Draw("SCATSame");

A12->SetMarkerColor(kGreen);
A12->SetMarkerStyle(31);
A12->SetMarkerSize(1.5);
A12->Draw("SCATSame");

A13->SetMarkerColor(kBlack);
A13->SetMarkerStyle(31);
A13->SetMarkerSize(1.5);
A13->Draw("SCATSame");

A14->SetMarkerColor(kOrange);
A14->SetMarkerStyle(31);
A14->SetMarkerSize(1.5);
A14->Draw("SCATSame");

TF1* f = new TF1("f","x",-10.,10.);
f->Draw("Same");
f->SetLineColor(kBlack);





   
    TLegend *leg = new TLegend(0.7401669,0.7092683,0.9398093,0.9092683,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);

   leg->SetHeader("Channels","C"); // option "C" allows to center the header
   leg->AddEntry(A9,"Channel A_{9}","p");
   leg->AddEntry(A10,"Channel A_{10}","p");
   leg->AddEntry(A12,"Channel A_{12}","p");
   leg->AddEntry(A13,"Channel A_{13}","p");
   leg->AddEntry(A14,"Channel A_{14}","p");
  

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

