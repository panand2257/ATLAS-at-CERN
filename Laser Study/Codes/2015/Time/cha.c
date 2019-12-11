
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


TH1F* A = new TH1F("A", "Channel A",50,-0.5,8.);
TH1F* B = new TH1F("B", "Channel B",50,-0.5,8.);
TH1F* D = new TH1F("D", "Channel D",50,-0.5,8.);

Double_t meanA[58];
Double_t meanB[58];
Double_t meanD[58];

Double_t rmsA[58];
Double_t rmsB[58];
Double_t rmsD[58];

TDatime da[58];
Double_t x[58];

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
    
     ifstream f_Laser;
     ifstream f_Masked;
        
         f_Laser.open(f);
         f_Masked.open("masked_channels_287886_12dec2015.txt");
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
                 
                 Las_constant[myPart_ctr][tmp_module-1][tmp_channel] = -100.-j ;
                 j--;
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
f_Masked.close();    
f_Laser.close();
    
    // Fill histograms
    
    
    for(int m = 1; m < 65; m++)
    {
        
          A -> Fill(Las_constant[1][m-1][1]);        //A1
           A -> Fill(Las_constant[1][m-1][4]);
            
          
           // Right Modules

           A -> Fill(Las_constant[2][m-1][1]);         //A-1
           A -> Fill(Las_constant[2][m-1][4]);
            
           





          // Left Modules
           B -> Fill(Las_constant[1][m-1][2]);         //B1
           B -> Fill(Las_constant[1][m-1][3]);
         
           // Right Modules


           B -> Fill(Las_constant[2][m-1][2]);        //B-1
           B -> Fill(Las_constant[2][m-1][3]);
 
         
 
            
           D -> Fill(Las_constant[1][m-1][13]);         //D1
           D -> Fill(Las_constant[1][m-1][14]);
         
            
           D -> Fill(Las_constant[2][m-1][13]);       //D-1
           D -> Fill(Las_constant[2][m-1][14]);
         

       
 //loop over modules
    

}
//cout<<A->GetMean()<<endl;
n++;
cout<<date<<f<<A->GetMean()<<"\n";
cout<<"\n";
meanA[n-1] = -1*(A->GetMean());
meanB[n-1] =-1*(B->GetMean());
meanD[n-1] = -1*(D->GetMean());

rmsA[n-1] = A->GetStdDev();
rmsB[n-1] = B->GetStdDev();
rmsD[n-1] = D->GetStdDev();

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











TGraphErrors* grA = new TGraphErrors(58,x,meanA,nullptr,rmsA); 
TGraphErrors* grB = new TGraphErrors(58,x,meanB,nullptr,rmsB); 
TGraphErrors* grD = new TGraphErrors(58,x,meanD,nullptr,rmsD);




    




Double_t _fx1[259] = {
   1.433315e+09,
   1.433317e+09,
   1.433317e+09,
   1.433351e+09,
   1.433351e+09,
   1.433439e+09,
   1.433439e+09,
   1.433539e+09,
   1.433539e+09,
   1.433655e+09,
   1.433655e+09,
   1.433697e+09,
   1.433697e+09,
   1.433879e+09,
   1.433879e+09,
   1.433895e+09,
   1.433895e+09,
   1.433917e+09,
   1.433917e+09,
   1.433959e+09,
   1.433959e+09,
   1.433977e+09,
   1.433977e+09,
   1.434118e+09,
   1.434118e+09,
   1.434204e+09,
   1.434204e+09,
   1.434255e+09,
   1.434255e+09,
   1.43603e+09,
   1.43603e+09,
   1.436065e+09,
   1.436065e+09,
   1.436128e+09,
   1.436128e+09,
   1.436229e+09,
   1.436229e+09,
   1.436296e+09,
   1.436296e+09,
   1.436354e+09,
   1.436354e+09,
   1.436541e+09,
   1.436541e+09,
   1.436599e+09,
   1.436599e+09,
   1.436638e+09,
   1.436638e+09,
   1.43666e+09,
   1.43666e+09,
   1.436759e+09,
   1.436759e+09,
   1.436826e+09,
   1.436826e+09,
   1.436911e+09,
   1.436911e+09,
   1.436977e+09,
   1.436977e+09,
   1.436994e+09,
   1.436994e+09,
   1.437349e+09,
   1.437349e+09,
   1.43736e+09,
   1.43736e+09,
   1.439428e+09,
   1.439428e+09,
   1.439473e+09,
   1.439473e+09,
   1.439536e+09,
   1.439536e+09,
   1.439565e+09,
   1.439565e+09,
   1.439621e+09,
   1.439621e+09,
   1.439639e+09,
   1.439639e+09,
   1.439656e+09,
   1.439656e+09,
   1.439688e+09,
   1.439688e+09,
   1.439738e+09,
   1.439738e+09,
   1.439744e+09,
   1.439744e+09,
   1.439844e+09,
   1.439844e+09,
   1.439857e+09,
   1.439857e+09,
   1.439947e+09,
   1.439947e+09,
   1.440123e+09,
   1.440124e+09,
   1.44018e+09,
   1.44018e+09,
   1.440244e+09,
   1.440244e+09,
   1.44029e+09,
   1.44029e+09,
   1.440332e+09,
   1.440332e+09,
   1.440352e+09,
   1.440352e+09,
   1.440441e+09,
   1.440441e+09,
   1.440494e+09,
   1.440494e+09,
   1.440503e+09,
   1.440503e+09,
   1.441561e+09,
   1.441561e+09,
   1.441578e+09,
   1.441578e+09,
   1.441693e+09,
   1.441693e+09,
   1.441772e+09,
   1.441772e+09,
   1.441829e+09,
   1.441829e+09,
   1.441845e+09,
   1.441845e+09,
   1.44198e+09,
   1.44198e+09,
   1.442104e+09,
   1.442104e+09,
   1.442155e+09,
   1.442155e+09,
   1.442194e+09,
   1.442194e+09,
   1.442241e+09,
   1.442241e+09,
   1.442397e+09,
   1.442397e+09,
   1.44244e+09,
   1.44244e+09,
   1.442512e+09,
   1.442512e+09,
   1.44258e+09,
   1.44258e+09,
   1.442628e+09,
   1.442628e+09,
   1.442691e+09,
   1.442691e+09,
   1.442768e+09,
   1.442768e+09,
   1.442788e+09,
   1.442788e+09,
   1.442854e+09,
   1.443048e+09,
   1.443048e+09,
   1.443146e+09,
   1.443146e+09,
   1.44317e+09,
   1.44317e+09,
   1.443249e+09,
   1.443249e+09,
   1.443299e+09,
   1.443299e+09,
   1.443319e+09,
   1.443319e+09,
   1.443388e+09,
   1.443388e+09,
   1.443435e+09,
   1.443435e+09,
   1.443449e+09,
   1.443449e+09,
   1.443492e+09,
   1.443492e+09,
   1.443525e+09,
   1.443525e+09,
   1.443641e+09,
   1.443641e+09,
   1.443765e+09,
   1.443765e+09,
   1.443786e+09,
   1.443786e+09,
   1.443892e+09,
   1.443892e+09,
   1.443965e+09,
   1.443965e+09,
   1.444054e+09,
   1.444054e+09,
   1.444067e+09,
   1.444067e+09,
   1.444094e+09,
   1.444094e+09,
   1.44413e+09,
   1.44413e+09,
   1.444149e+09,
   1.444149e+09,
   1.444379e+09,
   1.444379e+09,
   1.444399e+09,
   1.444399e+09,
   1.444428e+09,
   1.444428e+09,
   1.444564e+09,
   1.444564e+09,
   1.44467e+09,
   1.44467e+09,
   1.44474e+09,
   1.44474e+09,
   1.444864e+09,
   1.444864e+09,
   1.444894e+09,
   1.444894e+09,
   1.444931e+09,
   1.444931e+09,
   1.445019e+09,
   1.445019e+09,
   1.445053e+09,
   1.445053e+09,
   1.445102e+09,
   1.445102e+09,
   1.445127e+09,
   1.445127e+09,
   1.445176e+09,
   1.445176e+09,
   1.445309e+09,
   1.445309e+09,
   1.445327e+09,
   1.445327e+09,
   1.445386e+09,
   1.445386e+09,
   1.445452e+09,
   1.445452e+09,
   1.445551e+09,
   1.44555e+09,
   1.445615e+09,
   1.445615e+09,
   1.445696e+09,
   1.445696e+09,
   1.445825e+09,
   1.445825e+09,
   1.445872e+09,
   1.445872e+09,
   1.445987e+09,
   1.445987e+09,
   1.44608e+09,
   1.44608e+09,
   1.446182e+09,
   1.446182e+09,
   1.446221e+09,
   1.446221e+09,
   1.446311e+09,
   1.446311e+09,
   1.446407e+09,
   1.446407e+09,
   1.446442e+09,
   1.446442e+09,
   1.446482e+09,
   1.446482e+09,
   1.446592e+09,
   1.446592e+09,
   1.447e+09,
   1.448e+09,
   1.449e+09,
1.4491e+09,
1.4495e+09,
1.44985e+09,
1.4499e+09
   };
   Double_t _fy1[259] = {
   0,
   0,
   0.064,
   0.064,
   0.128,
   0.128,
   0.238,
   0.238,
   0.671,
   0.671,
   1.003,
   1.003,
   1.391,
   1.391,
   1.391,
   1.391,
   1.391,
   1.391,
   1.392,
   1.392,
   1.394,
   1.394,
   1.4,
   1.4,
   1.406,
   1.406,
   4.975,
   4.975,
   8.233,
   8.233,
   8.287,
   8.287,
   11.798,
   11.798,
   13.264,
   13.264,
   16.056,
   16.056,
   22.432,
   22.432,
   30.524,
   30.524,
   43.706,
   43.706,
   44.435,
   44.435,
   45.019,
   45.019,
   61.255,
   61.255,
   81.897,
   81.897,
   100.743,
   100.743,
   100.752,
   100.752,
   101.486,
   101.486,
   108.81,
   108.81,
   110.907,
   110.907,
   114.141,
   114.141,
   114.681,
   114.681,
   119.967,
   119.967,
   122.527,
   122.527,
   125.491,
   125.491,
   127.028,
   127.028,
   129.409,
   129.409,
   131.08,
   131.08,
   138.486,
   138.486,
   138.595,
   138.595,
   157.895,
   157.895,
   157.995,
   157.995,
   163.545,
   163.545,
   174.599,
   174.599,
   186.673,
   186.673,
   213.904,
   213.904,
   214.907,
   214.907,
   216.961,
   216.961,
   224.676,
   224.676,
   225.724,
   225.724,
   225.79,
   225.79,
   225.79,
   225.79,
   225.81,
   225.81,
   226.872,
   226.872,
   236.963,
   236.963,
   265.198,
   265.198,
   287.59,
   287.59,
   299.151,
   299.151,
   301.815,
   301.815,
   368.724,
   368.724,
   377.31,
   377.31,
   399.605,
   399.605,
   434.581,
   434.581,
   493.413,
   493.413,
   494.64,
   494.64,
   569.801,
   569.801,
   654.953,
   654.953,
   663.804,
   663.804,
   718.612,
   718.612,
   752.786,
   752.786,
   754.76,
   754.76,
   803.514,
   895.391,
   895.391,
   1001.511,
   1001.511,
   1005.185,
   1005.185,
   1111.167,
   1111.167,
   1121.312,
   1121.312,
   1127.099,
   1127.099,
   1211.179,
   1211.179,
   1280.161,
   1280.161,
   1294.814,
   1294.814,
   1324.47,
   1324.47,
   1354.79,
   1354.79,
   1527.499,
   1527.499,
   1588.084,
   1588.084,
   1599.239,
   1599.239,
   1746.982,
   1746.982,
   1890.188,
   1890.188,
   1929.55,
   1929.55,
   1944.559,
   1944.559,
   1995.723,
   1995.723,
   1997.294,
   1997.294,
   2018.075,
   2018.075,
   2248.247,
   2248.247,
   2277.907,
   2277.907,
   2292.887,
   2292.887,
   2490.245,
   2490.245,
   2656.005,
   2656.005,
   2656.005,
   2656.005,
   2656.006,
   2656.006,
   2656.017,
   2656.017,
   2656.055,
   2656.055,
   2656.17,
   2656.17,
   2656.273,
   2656.273,
   2656.553,
   2656.553,
   2656.641,
   2656.641,
   2656.767,
   2656.767,
   2656.769,
   2656.769,
   2670.089,
   2670.089,
   2700.724,
   2700.724,
   2806.437,
   2806.437,
   2819.803,
   2819.803,
   2943.308,
   2943.308,
   3008.102,
   3008.102,
   3041.103,
   3041.103,
   3057.004,
   3057.004,
   3356.63,
   3356.63,
   3517.996,
   3517.996,
   3563.759,
   3563.759,
   3575.413,
   3575.413,
   3797.165,
   3797.165,
   4079.722,
   4079.722,
   4148.013,
   4148.013,
   4187.28,
   4187.28,
   4370.512,
   4370.512,
   4370.512,
      4370.512,
     4370.512,
4370.512,
4370.512,
4370.512,
4370.512,
     0.};

   

   TGraph *graph = new TGraph(259,_fx1,_fy1);
  // graph->SetName("");
   graph->SetTitle("Luminosity");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cccccc");
   graph->SetFillColor(ci);

   ci = TColor::GetColor("#cccccc");
   graph->SetLineColor(ci);
   graph->Draw("AYFB1");
   graph->GetXaxis()->SetLimits(x[0],x[57]);
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
grA->GetXaxis()->SetLimits(x[0],x[57]);
grA->GetXaxis()->SetTimeDisplay(1);
grA->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%Y}%F1970-01-01 00:00:00s0");
grA->GetXaxis()->SetNdivisions(509);
grA->GetXaxis()->SetTimeDisplay(1);
grA->GetXaxis()->SetTimeOffset(0,"gmt");
grA->GetXaxis()->SetLabelOffset(0.013);
grA->GetXaxis()->SetLabelSize(0.02);
grA->GetXaxis()->SetTitleSize(0.035);
grA->GetXaxis()->SetTitleOffset(1.36);


grA->Draw("AY+PE");

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
   leg->SetHeader("Channels","C"); // option "C" allows to center the header
   leg->AddEntry(grA,"Channel A","p");
   leg->AddEntry(grB,"Channel BC","p");
   leg->AddEntry(grD,"Channel D","p");
   
  leg->Draw();
  


  
   
   

   


   
    
    
 } 

