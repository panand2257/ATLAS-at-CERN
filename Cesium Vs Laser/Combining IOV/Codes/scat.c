
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

void scat()
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
    

/*     f_Laser1.open("LaserCombined_267534.txt");                                                                          //IOV1
    f_Laser2.open("LaserCombined_272493.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_263962_11jun2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_270000_17jul2015.txt");
*/

  f_Laser1.open("LaserCombined_272493.txt");                                                                          //IOV2
    f_Laser2.open("LaserCombined_284682.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_270000_17jul2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_277321_27aug2015.txt");



   



    //  f_Masked.open("masked_channels_267534.txt");      //IOV1
     f_Masked.open("masked_channels_272493.txt");      //IOV2


    
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


   // f_Masked.open("masked_channels_272493.txt");     //IOV1
   f_Masked.open("masked_channels_284682.txt");     //IOV2

 


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


   // f_Masked.open("masked_channels_267534.txt");      //IOV1
   f_Masked.open("masked_channels_272493.txt");      //IOV2



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


   // f_Masked.open("masked_channels_272493.txt");     //IOV1
   f_Masked.open("masked_channels_284682.txt");     //IOV2

 
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
    
    
TH2F* A = new TH2F("Channel A", "IOV2",50,-6.,4.,50,-6.,4.);
TH2F* B = new TH2F("Channel BC", "IOV2",50,-6.,4.,50,-6.,4.);
TH2F* D = new TH2F("Channel D", "IOV2",50,-6.,4.,50,-6.,4.);


    for(int m = 1; m < 65; m++)
    {



          A -> Fill(Las_constant[1][m-1][1],Cesium_constant[1][m-1][1]);        //A1
           A -> Fill(Las_constant[1][m-1][4],Cesium_constant[1][m-1][4]);
            
           A-> Fill(Las_constant[1][m-1][5],Cesium_constant[1][m-1][5]);        //A2
           A-> Fill(Las_constant[1][m-1][8],Cesium_constant[1][m-1][8]);
            
           A -> Fill(Las_constant[1][m-1][9],Cesium_constant[1][m-1][9]);        //A3
           A -> Fill(Las_constant[1][m-1][10],Cesium_constant[1][m-1][10]);
            
           A -> Fill(Las_constant[1][m-1][15],Cesium_constant[1][m-1][15]);       //A4
           A -> Fill(Las_constant[1][m-1][18],Cesium_constant[1][m-1][18]);
            
           A -> Fill(Las_constant[1][m-1][19],Cesium_constant[1][m-1][19]);       //A5
           A -> Fill(Las_constant[1][m-1][20],Cesium_constant[1][m-1][20]);
             
           A -> Fill(Las_constant[1][m-1][23],Cesium_constant[1][m-1][23]);       //A6
           A -> Fill(Las_constant[1][m-1][26],Cesium_constant[1][m-1][26]);
           
           A -> Fill(Las_constant[1][m-1][29],Cesium_constant[1][m-1][29]);       //A7
           A -> Fill(Las_constant[1][m-1][32],Cesium_constant[1][m-1][32]);
            
           A -> Fill(Las_constant[1][m-1][35],Cesium_constant[1][m-1][35]);       //A8
           A -> Fill(Las_constant[1][m-1][38],Cesium_constant[1][m-1][38]);
            
           A -> Fill(Las_constant[1][m-1][36],Cesium_constant[1][m-1][36]);       //A9
           A -> Fill(Las_constant[1][m-1][37],Cesium_constant[1][m-1][37]);
            
           A -> Fill(Las_constant[1][m-1][45],Cesium_constant[1][m-1][45]);       //A10
           A -> Fill(Las_constant[1][m-1][46],Cesium_constant[1][m-1][46]);
            
           A -> Fill(Las_constant[0][m-1][6],Cesium_constant[0][m-1][6]);         //A12
           A -> Fill(Las_constant[0][m-1][7],Cesium_constant[0][m-1][7]);

           A -> Fill(Las_constant[0][m-1][10],Cesium_constant[0][m-1][10]);         //A13
           A -> Fill(Las_constant[0][m-1][11],Cesium_constant[0][m-1][11]);
    
           A -> Fill(Las_constant[0][m-1][20],Cesium_constant[0][m-1][20]);        //A14
           A -> Fill(Las_constant[0][m-1][21],Cesium_constant[0][m-1][21]);

           A -> Fill(Las_constant[0][m-1][31],Cesium_constant[0][m-1][31]);        //A15
           A -> Fill(Las_constant[0][m-1][32],Cesium_constant[0][m-1][32]);
 
           A -> Fill(Las_constant[0][m-1][40],Cesium_constant[0][m-1][40]);        //A16
           A -> Fill(Las_constant[0][m-1][41],Cesium_constant[0][m-1][41]);
 
           // Right Modules

           A -> Fill(Las_constant[2][m-1][1],Cesium_constant[2][m-1][1]);         //A-1
           A -> Fill(Las_constant[2][m-1][4],Cesium_constant[2][m-1][4]);
            
           A -> Fill(Las_constant[2][m-1][5],Cesium_constant[2][m-1][5]);         //A-2
           A -> Fill(Las_constant[2][m-1][8],Cesium_constant[2][m-1][8]);
            
           A -> Fill(Las_constant[2][m-1][9],Cesium_constant[2][m-1][9]);        //A-3
           A -> Fill(Las_constant[2][m-1][10],Cesium_constant[2][m-1][10]);
            
           A -> Fill(Las_constant[2][m-1][15],Cesium_constant[2][m-1][15]);        //A-4
           A -> Fill(Las_constant[2][m-1][18],Cesium_constant[2][m-1][18]);
            
           A -> Fill(Las_constant[2][m-1][19],Cesium_constant[2][m-1][19]);       //A-5
           A -> Fill(Las_constant[2][m-1][20],Cesium_constant[2][m-1][20]);
            
           A -> Fill(Las_constant[2][m-1][23],Cesium_constant[2][m-1][23]);       //A-6
           A -> Fill(Las_constant[2][m-1][26],Cesium_constant[2][m-1][26]);
            
           A -> Fill(Las_constant[2][m-1][29],Cesium_constant[2][m-1][29]);       //A-7
           A -> Fill(Las_constant[2][m-1][32],Cesium_constant[2][m-1][32]);
            
           A -> Fill(Las_constant[2][m-1][35],Cesium_constant[2][m-1][35]);      //A-8
           A -> Fill(Las_constant[2][m-1][38],Cesium_constant[2][m-1][38]);
            
           A -> Fill(Las_constant[2][m-1][36],Cesium_constant[2][m-1][36]);       //A-9
           A -> Fill(Las_constant[2][m-1][37],Cesium_constant[2][m-1][37]);
            
           A -> Fill(Las_constant[2][m-1][45],Cesium_constant[2][m-1][45]);       //A-10
           A -> Fill(Las_constant[2][m-1][46],Cesium_constant[2][m-1][46]);
            
           A -> Fill(Las_constant[3][m-1][6],Cesium_constant[3][m-1][6]);         //A-12
           A -> Fill(Las_constant[3][m-1][7],Cesium_constant[3][m-1][7]);

           A -> Fill(Las_constant[3][m-1][10],Cesium_constant[3][m-1][10]);         //A-13
           A -> Fill(Las_constant[3][m-1][11],Cesium_constant[3][m-1][11]);
    
           A -> Fill(Las_constant[3][m-1][20],Cesium_constant[3][m-1][20]);        //A-14
           A -> Fill(Las_constant[3][m-1][21],Cesium_constant[3][m-1][21]);

           A -> Fill(Las_constant[3][m-1][31],Cesium_constant[3][m-1][31]);        //A-15
           A -> Fill(Las_constant[3][m-1][32],Cesium_constant[3][m-1][32]);
 
           A -> Fill(Las_constant[3][m-1][40],Cesium_constant[3][m-1][40]);        //A-16
           A -> Fill(Las_constant[3][m-1][41],Cesium_constant[3][m-1][41]);






          // Left Modules
           B -> Fill(Las_constant[1][m-1][2],Cesium_constant[1][m-1][2]);         //B1
           B -> Fill(Las_constant[1][m-1][3],Cesium_constant[1][m-1][3]);
            
           B -> Fill(Las_constant[1][m-1][6],Cesium_constant[1][m-1][6]);         //B2
           B -> Fill(Las_constant[1][m-1][7],Cesium_constant[1][m-1][7]);
            
           B -> Fill(Las_constant[1][m-1][11],Cesium_constant[1][m-1][11]);        //B3
           B -> Fill(Las_constant[1][m-1][12],Cesium_constant[1][m-1][12]);
            
           B -> Fill(Las_constant[1][m-1][16],Cesium_constant[1][m-1][16]);        //B4
           B -> Fill(Las_constant[1][m-1][17],Cesium_constant[1][m-1][17]);
            
           B -> Fill(Las_constant[1][m-1][21],Cesium_constant[1][m-1][21]);        //B5
           B -> Fill(Las_constant[1][m-1][22],Cesium_constant[1][m-1][22]);
           
           B -> Fill(Las_constant[1][m-1][27],Cesium_constant[1][m-1][27]);       //B6
           B -> Fill(Las_constant[1][m-1][28],Cesium_constant[1][m-1][28]);
            
           B -> Fill(Las_constant[1][m-1][33],Cesium_constant[1][m-1][33]);        //B7
           B -> Fill(Las_constant[1][m-1][34],Cesium_constant[1][m-1][34]);
            
           B -> Fill(Las_constant[1][m-1][39],Cesium_constant[1][m-1][39]);       //B8
           B -> Fill(Las_constant[1][m-1][40],Cesium_constant[1][m-1][40]);
            
           B -> Fill(Las_constant[1][m-1][42],Cesium_constant[1][m-1][42]);        //B9
           B -> Fill(Las_constant[1][m-1][47],Cesium_constant[1][m-1][47]);
            
           B -> Fill(Las_constant[0][m-1][8],Cesium_constant[0][m-1][8]);           //B11
           B -> Fill(Las_constant[0][m-1][9],Cesium_constant[0][m-1][9]);
            
           B -> Fill(Las_constant[0][m-1][14],Cesium_constant[0][m-1][14]);         //B12
           B -> Fill(Las_constant[0][m-1][15],Cesium_constant[0][m-1][15]);

           B -> Fill(Las_constant[0][m-1][22],Cesium_constant[0][m-1][22]);        //B13
           B -> Fill(Las_constant[0][m-1][23],Cesium_constant[0][m-1][23]);
    
           B -> Fill(Las_constant[0][m-1][30],Cesium_constant[0][m-1][30]);        //B14
           B -> Fill(Las_constant[0][m-1][35],Cesium_constant[0][m-1][35]);

           B -> Fill(Las_constant[0][m-1][36],Cesium_constant[0][m-1][36]);        //B15
           B -> Fill(Las_constant[0][m-1][39],Cesium_constant[0][m-1][39]);
 


           // Right Modules


           B -> Fill(Las_constant[2][m-1][2],Cesium_constant[2][m-1][2]);        //B-1
           B -> Fill(Las_constant[2][m-1][3],Cesium_constant[2][m-1][3]);
 
           B -> Fill(Las_constant[2][m-1][6],Cesium_constant[2][m-1][6]);        //B-2
           B -> Fill(Las_constant[2][m-1][7],Cesium_constant[2][m-1][7]);
            
           B -> Fill(Las_constant[2][m-1][11],Cesium_constant[2][m-1][11]);        //B-3
           B -> Fill(Las_constant[2][m-1][12],Cesium_constant[2][m-1][12]);
            
           B -> Fill(Las_constant[2][m-1][16],Cesium_constant[2][m-1][16]);      //B-4
           B -> Fill(Las_constant[2][m-1][17],Cesium_constant[2][m-1][17]);
            
           B -> Fill(Las_constant[2][m-1][21],Cesium_constant[2][m-1][21]);      //B-5
           B -> Fill(Las_constant[2][m-1][22],Cesium_constant[2][m-1][22]);
            
           B -> Fill(Las_constant[2][m-1][27],Cesium_constant[2][m-1][27]);        //B-6
           B -> Fill(Las_constant[2][m-1][28],Cesium_constant[2][m-1][28]);
            
           B -> Fill(Las_constant[2][m-1][33],Cesium_constant[2][m-1][33]);      //B-7
           B -> Fill(Las_constant[2][m-1][34],Cesium_constant[2][m-1][34]);
            
           B -> Fill(Las_constant[2][m-1][39],Cesium_constant[2][m-1][39]);      //B-8
           B -> Fill(Las_constant[2][m-1][40],Cesium_constant[2][m-1][40]);
            
           B -> Fill(Las_constant[2][m-1][42],Cesium_constant[2][m-1][42]);      //B-9
           B -> Fill(Las_constant[2][m-1][47],Cesium_constant[2][m-1][47]);
            
           B -> Fill(Las_constant[3][m-1][8],Cesium_constant[3][m-1][8]);     //B-11
           B -> Fill(Las_constant[3][m-1][9],Cesium_constant[3][m-1][9]);
            
           B -> Fill(Las_constant[3][m-1][14],Cesium_constant[3][m-1][14]);      //B-12
           B -> Fill(Las_constant[3][m-1][15],Cesium_constant[3][m-1][15]);
            
           B -> Fill(Las_constant[3][m-1][22],Cesium_constant[3][m-1][22]);         //B-13
           B -> Fill(Las_constant[3][m-1][23],Cesium_constant[3][m-1][23]);

           B -> Fill(Las_constant[3][m-1][30],Cesium_constant[3][m-1][30]);         //B-14
           B -> Fill(Las_constant[3][m-1][35],Cesium_constant[3][m-1][35]);
    
           B -> Fill(Las_constant[3][m-1][36],Cesium_constant[3][m-1][36]);        //B-15
           B -> Fill(Las_constant[3][m-1][39],Cesium_constant[3][m-1][39]);

 
           D -> Fill(Las_constant[1][m-1][0],Cesium_constant[1][m-1][0]);         //D0
           D -> Fill(Las_constant[1][m-1][0],Cesium_constant[1][m-1][0]);
           D -> Fill(Las_constant[2][m-1][0],Cesium_constant[2][m-1][0]);         //D0
           D -> Fill(Las_constant[2][m-1][0],Cesium_constant[2][m-1][0]);        
            
           D -> Fill(Las_constant[1][m-1][13],Cesium_constant[1][m-1][13]);         //D1
           D -> Fill(Las_constant[1][m-1][14],Cesium_constant[1][m-1][14]);
            
           D -> Fill(Las_constant[1][m-1][24],Cesium_constant[1][m-1][24]);        //D2
           D -> Fill(Las_constant[1][m-1][25],Cesium_constant[1][m-1][25]);
            
           D -> Fill(Las_constant[1][m-1][41],Cesium_constant[1][m-1][41]);        //D3
           D -> Fill(Las_constant[1][m-1][44],Cesium_constant[1][m-1][44]);
            
           D -> Fill(Las_constant[0][m-1][2],Cesium_constant[0][m-1][2]);        //D4
           D -> Fill(Las_constant[0][m-1][3],Cesium_constant[0][m-1][3]);
           
           D -> Fill(Las_constant[0][m-1][16],Cesium_constant[0][m-1][16]);       //D5
           D -> Fill(Las_constant[0][m-1][17],Cesium_constant[0][m-1][17]);
            
           D -> Fill(Las_constant[0][m-1][37],Cesium_constant[0][m-1][37]);        //D6
           D -> Fill(Las_constant[0][m-1][38],Cesium_constant[0][m-1][38]);
            
           D -> Fill(Las_constant[2][m-1][13],Cesium_constant[2][m-1][13]);       //D-1
           D -> Fill(Las_constant[2][m-1][14],Cesium_constant[2][m-1][14]);
            
           D -> Fill(Las_constant[2][m-1][24],Cesium_constant[2][m-1][24]);        //D-2
           D -> Fill(Las_constant[2][m-1][25],Cesium_constant[2][m-1][25]);
            
           D -> Fill(Las_constant[2][m-1][41],Cesium_constant[2][m-1][41]);           //D-3
           D -> Fill(Las_constant[2][m-1][44],Cesium_constant[2][m-1][44]);
            
           D -> Fill(Las_constant[3][m-1][2],Cesium_constant[3][m-1][2]);         //D-4
           D -> Fill(Las_constant[3][m-1][3],Cesium_constant[3][m-1][3]);

           D -> Fill(Las_constant[3][m-1][16],Cesium_constant[3][m-1][16]);        //D-5
           D -> Fill(Las_constant[3][m-1][17],Cesium_constant[3][m-1][17]);
    
           D -> Fill(Las_constant[3][m-1][37],Cesium_constant[3][m-1][37]);        //D-6
           D -> Fill(Las_constant[3][m-1][38],Cesium_constant[3][m-1][38]);



}//loop over modules
    
    
    
 auto myCanvas =  new TCanvas("myCanvas", "",484,86,800,750); 
A->GetXaxis()->SetTitle("Laser Drift (%)");
A->GetYaxis()->SetTitle("Cesium Drift (%)");
A->SetMarkerColor(kRed);
A->SetMarkerStyle(31);
A->SetMarkerSize(1.0);
A->Draw("SCAT");

B->GetXaxis()->SetTitle("Laser Drift (%)");
B->GetYaxis()->SetTitle("Cesium Drift (%)");
B->SetMarkerColor(kBlue);
B->SetMarkerStyle(31);
B->SetMarkerSize(1.0);
B->Draw("SCATSame");

D->GetXaxis()->SetTitle("Laser Drift (%)");
D->GetYaxis()->SetTitle("Cesium Drift (%)");
D->SetMarkerColor(kGreen);
D->SetMarkerStyle(31);
D->SetMarkerSize(1.0);
D->Draw("SCATSame");

TF1* f = new TF1("f","x",-10.,10.);
f->Draw("Same");
f->SetLineColor(kBlack);

/*TH1D * AX = A->ProjectionX();
TH1D * AY = A->ProjectionY();

TH1D * BX = B->ProjectionX();
TH1D * BY = B->ProjectionY();

TH1D * DX = D->ProjectionX();
TH1D * DY = D->ProjectionY();

DX->GetXaxis()->SetTitle("Drift (%)");
DX->SetTitle("IOV3");
DX->GetYaxis()->SetTitle("Number of entries");
DX->SetLineColor(kRed);
DY->SetLineColor(kBlue);

DX->Draw();
DY->Draw("Same");*/

   
   TLegend *leg = new TLegend(0.7401669,0.7092683,0.9398093,0.9092683,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
leg->SetTextSize(0.03);
   //leg->SetHeader("Channels","C"); // option "C" allows to center the header
   leg->AddEntry(A,"Channel A","p");
   leg->AddEntry(B,"Channel BC","p");
   leg->AddEntry(D,"Channel D","p");
  

  leg->Draw();

/*TLegend *leg = new TLegend(0.7401669,0.7092683,0.9398093,0.9092683,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.3);
   leg->SetHeader("Channel D","C"); // option "C" allows to center the header
   leg->AddEntry(DX,"Laser","l");
   leg->AddEntry(DY,"Cesium","l");*/

  

  leg->Draw();

myCanvas -> Draw();

TLatex printInfo;
    printInfo.SetNDC();
    printInfo.SetTextFont(72);
    printInfo.DrawLatex(0.15, 0.85, "#scale[0.9]{ATLAS}");
    printInfo.SetTextFont(42);
    printInfo.DrawLatex(0.25, 0.85, "#scale[0.9]{Internal}");
    printInfo.DrawLatex(0.15, 0.80, "#scale[0.9]{Tile Calorimeter}");
    printInfo.DrawLatex(0.15, 0.75, "#scale[0.9]{Laser Combined Method 2015}");

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

    
    
 } 

