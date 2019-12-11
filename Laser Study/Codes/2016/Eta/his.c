
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
#include "TPaveStats.h"

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

void his()
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
//f_Laser.open("LaserCombined_288060.txt");           //start
//f_Masked.open("masked_channels_288060.txt");

//f_Laser.open("LaserCombined_300585.txt");           //middle
//f_Masked.open("masked_channels_300585.txt");

f_Laser.open("LaserCombined_313662.txt");           //end
f_Masked.open("masked_channels_313662.txt");
    
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
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
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
    f_Laser.close();
    f_Masked.close();
    
    // Fill histoAms
    
    
       TH1F* A =   new TH1F("A", "Channel A",50,-5.,5.);
        

          TH1F* B=   new TH1F("BC", "Channel BC",50,-10.,25.);
        
         TH1F* D =   new TH1F("D", "Channel D",50,-5.,25.);









    for(int m = 1; m < 65; m++)
    {



          // Left Modules
           A -> Fill(Las_constant[1][m-1][1]);        //A1
           A -> Fill(Las_constant[1][m-1][4]);
            
           A -> Fill(Las_constant[1][m-1][5]);        //A2
           A -> Fill(Las_constant[1][m-1][8]);
            
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

           A-> Fill(Las_constant[2][m-1][1]);         //A-1
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
           B-> Fill(Las_constant[1][m-1][34]);
            
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

 
           D-> Fill(Las_constant[1][m-1][0]);         //D0
           D-> Fill(Las_constant[1][m-1][0]);
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






}//loop over modules
    
    
    
 auto myCanvas =  new TCanvas("myCanvas", "",484,86,800,750);
 

A -> SetMarkerStyle(21);
A ->SetMarkerColor(kRed);
A-> SetMarkerSize(1.5);

A->SetLineColor(kRed);
B->SetLineColor(kBlue);
D->SetLineColor(kGreen);

B -> SetMarkerStyle(22);
B ->SetMarkerColor(kBlue);
B-> SetMarkerSize(1.5);

D -> SetMarkerStyle(20);
D ->SetMarkerColor(kGreen);
D-> SetMarkerSize(1.5);






B ->Draw();
//B->SetTitle(";Drift (%) [Run 288060];Number of Entries");          //start
//B->SetTitle(";Drift (%) [Run 300585];Number of Entries");          //middle
B->SetTitle(";Drift (%) [Run 313662];Number of Entries");          //end
//A->GetYaxis()->SetRangeUser(0.4,1.4);
D ->Draw("Sames");
A ->Draw("Sames");



   
    TLegend *leg = new TLegend(0.7401669,0.7092683,0.9398093,0.9092683,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetHeader("Channels","C"); // option "C" allows to center the header
   leg->AddEntry(A,"Channel A","l");
   leg->AddEntry(B,"Channel BC","l");
   leg->AddEntry(D,"Channel D","l");
   
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
    printInfo.DrawLatex(0.15, 0.75, "#scale[0.5]{2016 Data #sqrt{s} = 13 TeV}");

  


  
    
    
 } 

