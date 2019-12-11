
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

void eta()
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
      f_Laser1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/LaserData/LaserCombined_267534.txt");                //IOV1
    f_Laser2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/LaserData/LaserCombined_271584.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_263962_11jun2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_270000_17jul2015.txt");


  /*  f_Laser1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/LaserData/LaserCombined_271584.txt");                //IOV2
    f_Laser2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/LaserData/LaserCombined_277320.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_270000_17jul2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_277321_27aug2015.txt");
*/

   /* f_Laser1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/LaserData/LaserCombined_277320.txt");                //IOV3
    f_Laser2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/LaserData/LaserCombined_284682.txt");
    f_Cesium1.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_277321_27aug2015.txt");
    f_Cesium2.open("/home/pratyush/Documents/Project/CERN/Cesium Vs Laser/CesiumData/cs_284600_03nov2015.txt");*/

    f_Masked.open("masked_channels_287886_12dec2015.txt");
    
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
              //cout<<line<<"\n";
             if (line.compare(mask) ==0) 
               {  
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
                 break;
               }
              else{ 
             double drift1= (tmp_constant-1)*100;
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift1;
                 }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Laser1.close();
f_Masked.close();


    f_Masked.open("masked_channels_287886_12dec2015.txt");
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
              //cout<<line<<"\n";
             if (line.compare(mask) ==0) 
               {  
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -1000. ;
                 break;
               }
              else{ 
             double drift2= (tmp_constant-1)*100;
             
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift2 - a;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Laser2.close();
f_Masked.close();


    f_Masked.open("masked_channels_287886_12dec2015.txt");
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
             double drift3= (tmp_constant-1)*100;
           Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift3;
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Cesium1.close();
f_Masked.close();


    f_Masked.open("masked_channels_287886_12dec2015.txt");
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
             double drift4= (tmp_constant-1)*100;
           Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] = drift4 - c;
           Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = abs(Cesium_constant[myPart_ctr][tmp_module-1][tmp_channel] -d);
                  }
        }
     f_Masked.clear();
     f_Masked.seekg(0, ios::beg);
}
    f_Cesium2.close();


f_Masked.close();
    
    // Fill histograms
    
    
       TH1F* A[30];  
         A[0] =   new TH1F("A1", "Channel A_{1}",50,-100.,100.);
         A[1] =   new TH1F("A2", "Channel A_{2}",50,-100.,100.);
         A[2] =   new TH1F("A3", "Channel A_{3}",50,-100.,100.);
         A[3] =   new TH1F("A4", "Channel A_{4}",50,-100.,100.);
         A[4] =   new TH1F("A5", "Channel A_{5}",50,-100.,100.);
         A[5] =   new TH1F("A6", "Channel A_{6}",50,-100.,100.);
         A[6] =   new TH1F("A7", "Channel A_{7}",50,-100.,100.);
         A[7] =   new TH1F("A8", "Channel A_{8}",50,-100.,100.);
         A[8] =   new TH1F("A9", "Channel A_{9}",50,-100.,100.);
         A[9] =   new TH1F("A10", "Channel A_{10}",50,-100.,100.);
         A[10] =   new TH1F("A12", "Channel A_{12}",50,-100.,100.);
         A[11] =   new TH1F("A12", "Channel A_{13}",50,-100.,100.);
         A[12] =   new TH1F("A13", "Channel A_{14}",50,-100.,100.);
         A[13] =   new TH1F("A14", "Channel A_{15}",50,-100.,100.); 
         A[14] =   new TH1F("A10", "Channel A_{16}",50,-100.,100.);
         A[15] =   new TH1F("A-1", "Channel A_{-1}",50,-100.,100.);
         A[16] =   new TH1F("A-2", "Channel A_{-2}",50,-100.,100.);
         A[17] =   new TH1F("A-3", "Channel A_{-3}",50,-100.,100.);
         A[18] =   new TH1F("A-4", "Channel A_{-4}",50,-100.,100.);
         A[19] =   new TH1F("A-5", "Channel A_{-5}",50,-100.,100.);
         A[20] =   new TH1F("A-6", "Channel A_{-6}",50,-100.,100.);
         A[21] =   new TH1F("A-7", "Channel A_{-7}",50,-100.,100.);
         A[22] =   new TH1F("A-8", "Channel A_{-8}",50,-100.,100.);
         A[23] =   new TH1F("A-9", "Channel A_{-9}",50,-100.,100.);
         A[24] =   new TH1F("A-10", "Channel A_{-10}",50,-100.,100.);
         A[25] =   new TH1F("A-12", "Channel A_{-12}",50,-100.,100.);
         A[26] =   new TH1F("A-13", "Channel A_{-13}",50,-100.,100.);
         A[27] =   new TH1F("A-14", "Channel A_{-14}",50,-100.,100.);
         A[28] =   new TH1F("A-15", "Channel A_{-15}",50,-100.,100.);
         A[29] =   new TH1F("A-16", "Channel A_{-16}",50,-100.,100.);

          TH1F* B[28];  
         B[0] =   new TH1F("B1", "Channel B_{1}",50,-100.,100.);
         B[1] =   new TH1F("B2", "Channel B_{2}",50,-100.,100.);
         B[2] =   new TH1F("B3", "Channel B_{3}",50,-100.,100.);
         B[3] =   new TH1F("B4", "Channel B_{4}",50,-100.,100.);
         B[4] =   new TH1F("B5", "Channel B_{5}",50,-100.,100.);
         B[5] =   new TH1F("B6", "Channel B_{6}",50,-100.,100.);
         B[6] =   new TH1F("B7", "Channel B_{7}",50,-100.,100.);
         B[7] =   new TH1F("B8", "Channel B_{8}",50,-100.,100.);
         B[8] =   new TH1F("B9", "Channel B_{9}",50,-100.,100.);
         B[9] =   new TH1F("B11", "Channel B_{11}",50,-100.,100.);
         B[10] =   new TH1F("B12", "Channel B_{12}",50,-100.,100.);
         B[11] =   new TH1F("B13", "Channel B_{13}",50,-100.,100.);
         B[12] =   new TH1F("B14", "Channel B_{14}",50,-100.,100.);
         B[13] =   new TH1F("B15", "Channel B_{15}",50,-100.,100.); 
         B[14] =   new TH1F("B-1", "Channel B_{-1}",50,-100.,100.);
         B[15] =   new TH1F("B-2", "Channel B_{-2}",50,-100.,100.);
         B[16] =   new TH1F("B-3", "Channel A_{-3}",50,-100.,100.);
         B[17] =   new TH1F("B-4", "Channel B_{-4}",50,-100.,100.);
         B[18] =   new TH1F("B-5", "Channel B_{-5}",50,-100.,100.);
         B[19] =   new TH1F("B-6", "Channel B_{-6}",50,-100.,100.);
         B[20] =   new TH1F("B-7", "Channel B_{-7}",50,-100.,100.);
         B[21] =   new TH1F("B-8", "Channel B_{-8}",50,-100.,100.);
         B[22] =   new TH1F("B-9", "Channel B_{-9}",50,-100.,100.);
         B[23] =   new TH1F("B-11", "Channel B_{-11}",50,-100.,100.);
         B[24] =   new TH1F("B-12", "Channel B_{-12}",50,-100.,100.);
         B[25] =   new TH1F("B-13", "Channel B_{-13}",50,-100.,100.);
         B[26] =   new TH1F("B-14", "Channel B_{-14}",50,-100.,100.);
         B[27] =   new TH1F("B-15", "Channel B_{-15}",50,-100.,100.);
         
         TH1F* D[13];  
         D[0] =   new TH1F("D0", "Channel D_{0}",50,-100.,100.);
         D[1] =   new TH1F("D1", "Channel D_{1}",50,-100.,100.);
         D[2] =   new TH1F("D2", "Channel D_{2}",50,-100.,100.);
         D[3] =   new TH1F("D3", "Channel D_{3}",50,-100.,100.);
         D[4] =   new TH1F("D4", "Channel D_{4}",50,-100.,100.);
         D[5] =   new TH1F("D5", "Channel D_{5}",50,-100.,100.);
         D[6] =   new TH1F("D6", "Channel D_{6}",50,-100.,100.);
         D[7] =   new TH1F("D-1", "Channel D_{-1}",50,-100.,100.);
         D[8] =   new TH1F("D-2", "Channel D_{-2}",50,-100.,100.);
         D[9] =   new TH1F("D-3", "Channel D_{-3}",50,-100.,100.);
         D[10] =   new TH1F("D-4", "Channel D_{-4}",50,-100.,100.);
         D[11] =   new TH1F("D-5", "Channel D_{-5}",50,-100.,100.);
         D[12] =   new TH1F("D-6", "Channel D_{-6}",50,-100.,100.);









    for(int m = 1; m < 65; m++)
    {



          // Left Modules
           A[0] -> Fill(Las_constant[1][m-1][1]);        //A1
           A[0] -> Fill(Las_constant[1][m-1][4]);
            
           A[1] -> Fill(Las_constant[1][m-1][5]);        //A2
           A[1] -> Fill(Las_constant[1][m-1][8]);
            
           A[2] -> Fill(Las_constant[1][m-1][9]);        //A3
           A[2] -> Fill(Las_constant[1][m-1][10]);
            
           A[3] -> Fill(Las_constant[1][m-1][15]);       //A4
           A[3] -> Fill(Las_constant[1][m-1][18]);
            
           A[4] -> Fill(Las_constant[1][m-1][19]);       //A5
           A[4] -> Fill(Las_constant[1][m-1][20]);
             
           A[5] -> Fill(Las_constant[1][m-1][23]);       //A6
           A[5] -> Fill(Las_constant[1][m-1][26]);
           
           A[6] -> Fill(Las_constant[1][m-1][29]);       //A7
           A[6] -> Fill(Las_constant[1][m-1][32]);
            
           A[7] -> Fill(Las_constant[1][m-1][35]);       //A8
           A[7] -> Fill(Las_constant[1][m-1][38]);
            
           A[8] -> Fill(Las_constant[1][m-1][36]);       //A9
           A[8] -> Fill(Las_constant[1][m-1][37]);
            
           A[9] -> Fill(Las_constant[1][m-1][45]);       //A10
           A[9] -> Fill(Las_constant[1][m-1][46]);
            
           A[10] -> Fill(Las_constant[0][m-1][6]);         //A12
           A[10] -> Fill(Las_constant[0][m-1][7]);

           A[11] -> Fill(Las_constant[0][m-1][10]);         //A13
           A[11] -> Fill(Las_constant[0][m-1][11]);
    
           A[12] -> Fill(Las_constant[0][m-1][20]);        //A14
           A[12] -> Fill(Las_constant[0][m-1][21]);

           A[13] -> Fill(Las_constant[0][m-1][31]);        //A15
           A[13] -> Fill(Las_constant[0][m-1][32]);
 
           A[14] -> Fill(Las_constant[0][m-1][40]);        //A16
           A[14] -> Fill(Las_constant[0][m-1][41]);
 
           // Right Modules

           A[15] -> Fill(Las_constant[2][m-1][1]);         //A-1
           A[15] -> Fill(Las_constant[2][m-1][4]);
            
           A[16] -> Fill(Las_constant[2][m-1][5]);         //A-2
           A[16] -> Fill(Las_constant[2][m-1][8]);
            
           A[17] -> Fill(Las_constant[2][m-1][9]);        //A-3
           A[17] -> Fill(Las_constant[2][m-1][10]);
            
           A[18] -> Fill(Las_constant[2][m-1][15]);        //A-4
           A[18] -> Fill(Las_constant[2][m-1][18]);
            
           A[19] -> Fill(Las_constant[2][m-1][19]);       //A-5
           A[19] -> Fill(Las_constant[2][m-1][20]);
            
           A[20] -> Fill(Las_constant[2][m-1][23]);       //A-6
           A[20] -> Fill(Las_constant[2][m-1][26]);
            
           A[21] -> Fill(Las_constant[2][m-1][29]);       //A-7
           A[21] -> Fill(Las_constant[2][m-1][32]);
            
           A[22] -> Fill(Las_constant[2][m-1][35]);      //A-8
           A[22] -> Fill(Las_constant[2][m-1][38]);
            
           A[23] -> Fill(Las_constant[2][m-1][36]);       //A-9
           A[23] -> Fill(Las_constant[2][m-1][37]);
            
           A[24] -> Fill(Las_constant[2][m-1][45]);       //A-10
           A[24] -> Fill(Las_constant[2][m-1][46]);
            
           A[25] -> Fill(Las_constant[3][m-1][6]);         //A-12
           A[25] -> Fill(Las_constant[3][m-1][7]);

           A[26] -> Fill(Las_constant[3][m-1][10]);         //A-13
           A[26] -> Fill(Las_constant[3][m-1][11]);
    
           A[27] -> Fill(Las_constant[3][m-1][20]);        //A-14
           A[27] -> Fill(Las_constant[3][m-1][21]);

           A[28] -> Fill(Las_constant[3][m-1][31]);        //A-15
           A[28] -> Fill(Las_constant[3][m-1][32]);
 
           A[29] -> Fill(Las_constant[3][m-1][40]);        //A-16
           A[29] -> Fill(Las_constant[3][m-1][41]);






          // Left Modules
           B[0] -> Fill(Las_constant[1][m-1][2]);         //B1
           B[0] -> Fill(Las_constant[1][m-1][3]);
            
           B[1] -> Fill(Las_constant[1][m-1][6]);         //B2
           B[1] -> Fill(Las_constant[1][m-1][7]);
            
           B[2] -> Fill(Las_constant[1][m-1][11]);        //B3
           B[2] -> Fill(Las_constant[1][m-1][12]);
            
           B[3] -> Fill(Las_constant[1][m-1][16]);        //B4
           B[3] -> Fill(Las_constant[1][m-1][17]);
            
           B[4] -> Fill(Las_constant[1][m-1][21]);        //B5
           B[4] -> Fill(Las_constant[1][m-1][22]);
           
           B[5] -> Fill(Las_constant[1][m-1][27]);       //B6
           B[5] -> Fill(Las_constant[1][m-1][28]);
            
           B[6] -> Fill(Las_constant[1][m-1][33]);        //B7
           B[6] -> Fill(Las_constant[1][m-1][34]);
            
           B[7] -> Fill(Las_constant[1][m-1][39]);       //B8
           B[7] -> Fill(Las_constant[1][m-1][40]);
            
           B[8] -> Fill(Las_constant[1][m-1][42]);        //B9
           B[8] -> Fill(Las_constant[1][m-1][47]);
            
           B[9] -> Fill(Las_constant[0][m-1][8]);           //B11
           B[9] -> Fill(Las_constant[0][m-1][9]);
            
           B[10] -> Fill(Las_constant[0][m-1][14]);         //B12
           B[10] -> Fill(Las_constant[0][m-1][15]);

           B[11] -> Fill(Las_constant[0][m-1][22]);        //B13
           B[11] -> Fill(Las_constant[0][m-1][23]);
    
           B[12] -> Fill(Las_constant[0][m-1][30]);        //B14
           B[12] -> Fill(Las_constant[0][m-1][35]);

           B[13] -> Fill(Las_constant[0][m-1][36]);        //B15
           B[13] -> Fill(Las_constant[0][m-1][39]);
 


           // Right Modules


           B[14] -> Fill(Las_constant[2][m-1][2]);        //B-1
           B[14] -> Fill(Las_constant[2][m-1][3]);
 
           B[15] -> Fill(Las_constant[2][m-1][6]);        //B-2
           B[15] -> Fill(Las_constant[2][m-1][7]);
            
           B[16] -> Fill(Las_constant[2][m-1][11]);        //B-3
           B[16] -> Fill(Las_constant[2][m-1][12]);
            
           B[17] -> Fill(Las_constant[2][m-1][16]);      //B-4
           B[17] -> Fill(Las_constant[2][m-1][17]);
            
           B[18] -> Fill(Las_constant[2][m-1][21]);      //B-5
           B[18] -> Fill(Las_constant[2][m-1][22]);
            
           B[19] -> Fill(Las_constant[2][m-1][27]);        //B-6
           B[19] -> Fill(Las_constant[2][m-1][28]);
            
           B[20] -> Fill(Las_constant[2][m-1][33]);      //B-7
           B[20] -> Fill(Las_constant[2][m-1][34]);
            
           B[21] -> Fill(Las_constant[2][m-1][39]);      //B-8
           B[21] -> Fill(Las_constant[2][m-1][40]);
            
           B[22] -> Fill(Las_constant[2][m-1][42]);      //B-9
           B[22] -> Fill(Las_constant[2][m-1][47]);
            
           B[23] -> Fill(Las_constant[3][m-1][8]);     //B-11
           B[23] -> Fill(Las_constant[3][m-1][9]);
            
           B[24] -> Fill(Las_constant[3][m-1][14]);      //B-12
           B[24] -> Fill(Las_constant[3][m-1][15]);
            
           B[25] -> Fill(Las_constant[3][m-1][22]);         //B-13
           B[25] -> Fill(Las_constant[3][m-1][23]);

           B[26] -> Fill(Las_constant[3][m-1][30]);         //B-14
           B[26] -> Fill(Las_constant[3][m-1][35]);
    
           B[27] -> Fill(Las_constant[3][m-1][36]);        //B-15
           B[27] -> Fill(Las_constant[3][m-1][39]);

 
           D[0] -> Fill(Las_constant[1][m-1][0]);         //D0
           D[0] -> Fill(Las_constant[1][m-1][0]);
           D[0] -> Fill(Las_constant[2][m-1][0]);         //D0
           D[0] -> Fill(Las_constant[2][m-1][0]);        
            
           D[1] -> Fill(Las_constant[1][m-1][13]);         //D1
           D[1] -> Fill(Las_constant[1][m-1][14]);
            
           D[2] -> Fill(Las_constant[1][m-1][24]);        //D2
           D[2] -> Fill(Las_constant[1][m-1][25]);
            
           D[3] -> Fill(Las_constant[1][m-1][41]);        //D3
           D[3] -> Fill(Las_constant[1][m-1][44]);
            
           D[4] -> Fill(Las_constant[0][m-1][2]);        //D4
           D[4] -> Fill(Las_constant[0][m-1][3]);
           
           D[5] -> Fill(Las_constant[0][m-1][16]);       //D5
           D[5] -> Fill(Las_constant[0][m-1][17]);
            
           D[6] -> Fill(Las_constant[0][m-1][37]);        //D6
           D[6] -> Fill(Las_constant[0][m-1][38]);
            
           D[7] -> Fill(Las_constant[2][m-1][13]);       //D-1
           D[7] -> Fill(Las_constant[2][m-1][14]);
            
           D[8] -> Fill(Las_constant[2][m-1][24]);        //D-2
           D[8] -> Fill(Las_constant[2][m-1][25]);
            
           D[9] -> Fill(Las_constant[2][m-1][41]);           //D-3
           D[9] -> Fill(Las_constant[2][m-1][44]);
            
           D[10] -> Fill(Las_constant[3][m-1][2]);         //D-4
           D[10] -> Fill(Las_constant[3][m-1][3]);

           D[11] -> Fill(Las_constant[3][m-1][16]);        //D-5
           D[11] -> Fill(Las_constant[3][m-1][17]);
    
           D[12] -> Fill(Las_constant[3][m-1][37]);        //D-6
           D[12] -> Fill(Las_constant[3][m-1][38]);






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

double meanD[13];
double etaD[13] = {0.,0.2,0.4,0.6,0.8,1.0,1.245,-0.2,-0.4,-0.6,-0.8,-1.0,-1.245};
double rmsD[30];

for(int i=0;i<30;i++){
 

meanA[i] = (A[i] ->GetMean());
rmsA[i] = A[i] ->GetMeanError();
erreta[i] = 0.05;
if(i<10)
etaA[i] = (0.05 + 0.1*i);

if(i>9 && i<15)
etaA[i] = (0.15 +0.1*i);

if(i>14)
etaA[i] = -1*etaA[i-15];

}

for(int i=0;i<28;i++){
 

meanB[i] = (B[i] ->GetMean());
rmsB[i] = B[i] ->GetMeanError();

if(i<9)
etaB[i] = (0.05 + 0.1*i);

if(i>8 && i<14)
etaB[i] = (0.15 +0.1*i);

if(i>13)
etaB[i] = -1*etaB[i-14];

}


for(int i=0;i<13;i++){
 

meanD[i] = (D[i] ->GetMean());
rmsD[i] = D[i] ->GetMeanError();
erretaD[i] = 0.1;
}




TGraphErrors* grA = new TGraphErrors(30,etaA,meanA,erreta,rmsA); 


TGraphErrors* grB = new TGraphErrors(28,etaB,meanB,erreta,rmsB); 
TGraphErrors* grD = new TGraphErrors(13,etaD,meanD,erretaD,rmsD); 
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
grA->SetTitle("; #eta; |Cs - Laser| Difference in Drift (%) [IOV1]");
//grA->GetYaxis()->SetRangeUser(0.5,1.3);
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
   


  
    
    
 } 

