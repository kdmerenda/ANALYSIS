#include "ComparaisonStudy.h"
#include <evt/Event.h>
#include <fevt/FEvent.h>
#include <fevt/Eye.h>
#include <fevt/TelescopeSimData.h>
#include <fevt/Telescope.h>
#include <fevt/PixelSimData.h>
#include <fevt/Pixel.h>
#include <fevt/FdConstants.h>
#include <det/Detector.h>
#include <fdet/FDetector.h>
#include <fdet/Eye.h>
#include <fdet/Telescope.h>
#include <fdet/Camera.h>
#include <fdet/Pixel.h>
#include <fdet/Channel.h>
#include <fdet/Mirror.h>
#include <fdet/Filter.h>
#include <fdet/Corrector.h>

#include <utl/Point.h>
#include <utl/UTMPoint.h>
#include <utl/Particle.h>
#include <utl/Vector.h>
#include <utl/TimeStamp.h>
#include <utl/Particle.h>
#include <utl/ErrorLogger.h>
#include <utl/ReferenceEllipsoid.h>
#include <utl/RandomEngine.h>
#include <fwk/LocalCoordinateSystem.h>
#include <fwk/CoordinateSystemRegistry.h>
#include <fwk/CentralConfig.h>
#include <fwk/SVNGlobalRevision.h>
#include <fwk/RandomEngineRegistry.h>
#include <utl/Math.h>

#include <utl/MathConstants.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/AxialVector.h>
#include <utl/Vector.h>
#include <utl/PhysicalConstants.h>
#include <utl/AugerUnits.h>
#include <utl/UTCDateTime.h>
#include <CLHEP/Random/Randomize.h>

#include <TMap.h>
#include <TMath.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TVector3.h>
#include <vector>
#include <AugerEvent.h>
#include <EyeEvent.hh>
#include <EyePixelList.hh>
#include <FadcData.hh>
#include <TPaletteAxis.h>

using namespace std;
using namespace utl;
using namespace fwk;
using namespace evt;
using namespace fevt;
using namespace fdet;
using namespace det;
using CLHEP::RandFlat;

namespace ComparaisonStudyNS {

// Initialize your static data members.

// Define the methods preferibly in the same order as listed in the header file.

ComparaisonStudy::ComparaisonStudy() 
{
  // Initialize your class.
  // Prefer initialization lists over = to 
  // initialize your members.
}

ComparaisonStudy::~ComparaisonStudy()
{
}

VModule::ResultFlag 
ComparaisonStudy::Init()
{
  INFO("Init()");
  //here we want to start by comparing only the first page of the traces. later we will do a full comparaison when the simulation will be doing followers properly.   
  fEventCounter = 0;
  fSimEventCounter = 0;
  fNBinTrace = 300; //num bins in full trace
  fMaxBinTrace = 2999; //max num of bins in full trace
  fMinBinTrace = 0;
  fNBinPage = int(1000./((fMaxBinTrace-fMinBinTrace)/fNBinTrace)); //num bin per page
  int BinSize = int(((fMaxBinTrace-fMinBinTrace)/fNBinTrace));
  TFile *fTMP = new TFile("traces/traces.root","recreate");
  fTMP->Close();

  outputPlots = new TFile("traces/FullTraces.root","recreate");

  for(int i = 0; i < fNRows; i++){  
    //same title for both data and simulation
    TString hPixelRowTitle("Row ");hPixelRowTitle+= (i+1);hPixelRowTitle+= ";Time Bin (100 ns);Column Number;ADC Counts / ";hPixelRowTitle+= BinSize; " (100 ns) Time Bins";
    
    for(int j=0; j < fNumFiles; j++){
      TString hPixelRowNameSim("hSimPixelRow");hPixelRowNameSim+= (i+1);
      hPixelRowNameSim+="_";hPixelRowNameSim+=j;
      hPixelRow[j][i] = new TH2F(hPixelRowNameSim,hPixelRowTitle,fNBinTrace,fMinBinTrace,fMaxBinTrace,fNColumns,0,fNColumns);
    }
  }

  const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);
  //  detTelGlobal = detTelTMP;
  for(int i = 0; i < fNTels; i++){
    for(int j = 0; j < fNPixels; j++){
      for(int k=0; k < fNumFiles; k++){      
	TString hnameSim("hSim");hnameSim+=k;hnameSim+="_tel";
	hnameSim+=i+1;hnameSim+="_r";hnameSim+=detTelGlobal.GetPixel(j+1).GetRow();hnameSim+="_c";hnameSim+=detTelGlobal.GetPixel(j+1).GetColumn();
	fhRawPixel[k][i][j] = new TH1F(hnameSim,hnameSim,fNBinTrace,fMinBinTrace,fMaxBinTrace);
      }

    }
  }
  return eSuccess;
}

VModule::ResultFlag 
ComparaisonStudy::Run(evt::Event& event)
{
  const evt::Header& header = event.GetHeader();  
  
  AugerEvent& rawEvent = event.GetRawEvent();

  cout << header.GetTime().GetGPSNanoSecond() << endl;
  //simluated the events at three frame long for now. 
  if(fEventCounter%3 == 0){
    INFO("New SIM Event!");
    fSimEventCounter++; //sim event counter
    fEventCounter = 0; //page counter
  }
  fEventCounter++;
  fData = false;
  
  for (AugerEvent::EyeIterator eyeIter = rawEvent.EyesBegin();
       eyeIter != rawEvent.EyesEnd(); ++eyeIter) {
    
    TEyeEvent& eyeEvent = *eyeIter;
    TEyePixelList* eyePixelList = eyeEvent.GetPixelList();
    TEyeFADCData* eyeFADCData = eyeEvent.GetFADCData();
    const unsigned int numpixels = eyePixelList->GetNumPixels();
    TEyeEventHeader *eyeHeader = eyeEvent.GetEventHeader();
    const int eyeId = eyeEvent.GetEventHeader()->GetEyeNo();
    const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
    FEvent& fdEvent = event.GetFEvent();
    fevt::Eye& eye = fdEvent.GetEye(eyeId);
    const fdet::Eye& detEye = detFD.GetEye(eyeId);


    cout << fEventCounter <<" " << fSimEventCounter<< " " << eyeHeader->GetEventNo() <<  " "  <<  eyePixelList->GetEyeNo() <<  " " << eyePixelList->GetNumPixels() << endl;
    
    TFile *fTMP = new TFile("traces/traces.root","update");

    for (unsigned int iPixel = 0; iPixel < eyePixelList->GetNumPixels(); ++iPixel) {
      
      FdUtil::Fd::FdPixelNumber pixelNumber = eyePixelList->GetPixel(iPixel);
      const unsigned int mirrorId = FdUtil::Fd::GetEyeMirrorNo(pixelNumber);
      fevt::Telescope& tel = eye.GetTelescope(mirrorId);
      //cout << " detEye " << eyeId<<endl;
      const fdet::Telescope& detTel = detEye.GetTelescope(mirrorId);
      const unsigned int channelId = FdUtil::Fd::GetMirrorPixelNo(pixelNumber);
      //     const unsigned int pixelId = FdUtil::Fd::GetPixelNo(pixelNumber);
      
      const TFADCData* fadcData = eyeFADCData->GetFADCData(pixelNumber);

      
      const fdet::Channel& thisChannel = detTel.GetChannel(channelId);
      const unsigned int pixelId = thisChannel.GetPixelId();


      if (!fadcData || thisChannel.IsVirtual() || pixelId > 440) {
	continue;
      }
      //cout << mirrorId << " " << channelId << " " << pixelId << endl;
      
      
      const unsigned int startBin = fadcData->GetTraceStartBin();
      const unsigned int endBin = fadcData->GetTraceEndBin();
      
      TFADCData::FADCDataWord fadcword = fadcData->GetFADCTrace();
      
      
      baseline=0.;
      baselineRMS=0.;
      baselinetail=0.;
      
      for (unsigned int pos = startBin; pos <startBin+100; ++pos) baseline += int(FADCDataWordGetData(&fadcword[pos]))/100.;
      
      for (unsigned int pos = endBin-100; pos < endBin; ++pos) baselinetail += int(FADCDataWordGetData(&fadcword[pos]))/100.;
      
      TString hName("h");
      hName+="Sim";
      hName+=fSimEventCounter;
      hName+="Evt";
      hName+=fEventCounter;
      hName+="_m";hName+=mirrorId;hName+="_r";hName+=detTel.GetPixel(pixelId).GetRow();hName+="_c";hName+=detTel.GetPixel(pixelId).GetColumn();

      
      TH1F* trace = new TH1F(hName,hName,fNBinPage,startBin,endBin);
      
      for (unsigned int pos = startBin; pos <= endBin; ++pos) {
	charge = FADCDataWordGetData(&fadcword[pos]);
 	if (pixelId <= 440 && abs(baseline-baselinetail) > 10 && baseline >baselinetail){
	  trace->SetBinContent(int(pos/((fMaxBinTrace-fMinBinTrace)/fNBinTrace)),charge-baselinetail);	  
	}else{
	  trace->SetBinContent(int(pos/((fMaxBinTrace-fMinBinTrace)/fNBinTrace)),charge-baseline);	  
	}
      }//trace loop

      GlueTrace(fSimEventCounter-1,trace, mirrorId, pixelId);

    }//pixel loop
    fTMP->Write();
    fTMP->Close();
    
  }//mirror or eye loop?

  //reset page counter
  //everything below needs to happen only once fhRawPixel is fully filled!
  // if(!fEventCounter == 3){
  //   return eSuccess;
  // }
  if(fSimEventCounter == fNumFiles && fEventCounter == 3){
    

    //FIT THE TRACES
    const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
    Detector::GetInstance().Update(times);
    const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
    const fdet::Eye& detEye = detFD.GetEye(4);
    const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);
    //    double longitudes[fNumFiles] = {-62,-62.5,-63,-63.5,-64,-64.5,-65,-65.5,-66,-66.5,-67};
    double longitudes[fNumFiles] = {-65,-65.5,-66};
    
    double Integrals[fNumFiles];
    double IntegralsFromFit[fNumFiles];
    double TraceCounter[fNumFiles];
    double FirstPixel[fNumFiles];
    double FirstPixelTime[fNumFiles];
    double FirstPixelRow[fNumFiles];
    double FirstPixelColumn[fNumFiles];
    double PeakBrightPixel[fNumFiles];
    double PeakBrightPixelValue[fNumFiles];
    double InteBrightPixel[fNumFiles];
    double InteBrightPixelValue[fNumFiles];
    TF1* MINFIT;    
    vector<TF1*> allFits[fNumFiles];//this will save the fits done below, readjust for nTELS
    vector<double> allFitsRows[fNumFiles];//this will save the fits done below, readjust for nTELS
    vector<double> allFitsColumns[fNumFiles];//this will save the fits done below, readjust for nTELS
    int skewed = 0;

    for(int k = 0; k < fNumFiles; k++){
      double MINTIME = 99999.;
      double MAXAMPLITUDE = 0.;
      double MAXINTEGRAL = 0.;
      for(int i = 0; i < fNTels; i++){
	Integrals[k]=0;
	IntegralsFromFit[k]=0;
	TraceCounter[k]=0;
	for(int j = 0; j < fNPixels; j++){

	  if(fhRawPixel[k][i][j]->GetEntries() == 0){
	    delete fhRawPixel[k][i][j];
	    continue;
	  }

	  //fill the 2D hists per rows
	  for(int iBin = 1; iBin <= fhRawPixel[k][i][j]->GetNbinsX(); iBin++){ 
	    hPixelRow[k][detTelGlobal.GetPixel(j+1).GetRow()-1]->SetBinContent(iBin,detTelGlobal.GetPixel(j+1).GetColumn(),fhRawPixel[k][i][j]->GetBinContent(iBin));
	  }

	  
	  //need a way better fit...
	  TString FitName("g1");FitName+=i;FitName+=j;FitName+=k;
	  TF1* g1;

	  if(!skewed) g1 = new TF1(FitName,"gaus",200,3000);
	  if(skewed)  g1 = new TF1(FitName,asymGaussian,200,3000,4);

	  if(!skewed) g1->SetParameters(fhRawPixel[k][i][j]->GetMaximum(),fhRawPixel[k][i][j]->GetMaximumBin(),fhRawPixel[k][i][j]->GetStdDev());//initial guess
	  if(skewed) g1->SetParameters(fhRawPixel[k][i][j]->GetMaximum(),fhRawPixel[k][i][j]->GetMaximumBin(),fhRawPixel[k][i][j]->GetStdDev(),fhRawPixel[k][i][j]->GetStdDev());//initial guess
	  
	  fhRawPixel[k][i][j]->Fit(g1,"RQ");

	  //clean up the fits
	  if(!skewed && (
	     g1->GetParameter(1) < 280                      // above trigger bin
	     || g1->GetParameter(1) < 2*g1->GetParameter(2) // have a peak
	     || g1->GetParameter(0) < 20                     // positive peak above bg noise
	     || g1->GetParameter(0) <  g1->GetParError(0)
	     || g1->GetParameter(1) <  g1->GetParError(1)
	     || g1->GetParameter(2) <  g1->GetParError(2)
	     || g1->Integral(0,3000) < 0
	     || g1->GetParameter(0) > 8000
	     || g1->GetChisquare()/g1->GetNDF() > 40)
 	     )continue;

	  //ALERT: the 2*sigma stuff should be the amplitude not the mean..	  
	  double integralTMPHIST = fhRawPixel[k][i][j]->Integral();
	  //save some stuff of interest
	  Integrals[k]+=integralTMPHIST;
	  TraceCounter[k]++;
	  
	  //find first pixel 
	  if(g1->GetParameter(1) < MINTIME) {
	    MINTIME = g1->GetParameter(1);
	    FirstPixel[k] = j+1;
	  }
	  //find peak brightness
	  if(g1->GetParameter(0) > MAXAMPLITUDE) {
	    MAXAMPLITUDE = g1->GetParameter(0);
	    PeakBrightPixel[k] = j+1;
	  }
	  //find largest integral
	  double integralTMPFIT = g1->Integral(0,3000);
	  if(integralTMPFIT > MAXINTEGRAL) {
	    MAXINTEGRAL = integralTMPFIT;
	    InteBrightPixel[k] = j+1;
	  }
	  IntegralsFromFit[k]+=integralTMPFIT;

	  allFits[k].push_back(g1);
	  allFitsRows[k].push_back(detTelGlobal.GetPixel(j+1).GetRow());
	  allFitsColumns[k].push_back(detTelGlobal.GetPixel(j+1).GetColumn());
 	}
      }
      if(!int(FirstPixel[k])) continue;//skip if no fits were done.
      
      PeakBrightPixelValue[k] = MAXAMPLITUDE;
      InteBrightPixelValue[k] = MAXINTEGRAL;
      FirstPixelTime[k] = MINTIME;
      FirstPixelRow[k] = detTelGlobal.GetPixel(FirstPixel[k]).GetRow();
      FirstPixelColumn[k] = detTelGlobal.GetPixel(FirstPixel[k]).GetColumn();
      
      cout << endl;
      cout << "Longitude: "<< longitudes[k] << endl;
      cout << "     TOT Integral: " <<  Integrals[k] << " TOT Integral (FITS): " <<  IntegralsFromFit[k] << " Pixel Count: " << TraceCounter[k] << endl;
      cout << "     First Pixel Time: R" << detTelGlobal.GetPixel(FirstPixel[k]).GetRow()<< " C" <<detTelGlobal.GetPixel(FirstPixel[k]).GetColumn() << " @ " <<MINTIME << endl;
      cout << "     Max. Peak Brightness Pixel: R"<<  detTelGlobal.GetPixel(PeakBrightPixel[k]).GetRow()<< " C" <<detTelGlobal.GetPixel(PeakBrightPixel[k]).GetColumn() << " @ " << MAXAMPLITUDE << endl;
      cout << "     Max. Integrated Brightness: R"<<  detTelGlobal.GetPixel(InteBrightPixel[k]).GetRow()<< " C" <<detTelGlobal.GetPixel(InteBrightPixel[k]).GetColumn() << " @ " << MAXINTEGRAL <<endl;
      
    }  

    //plot the values
    TCanvas * cstats = new TCanvas("cstats","cstats",1000,600);

    TGraph gTraceCounter(fNumFiles,longitudes,TraceCounter);
    gTraceCounter.SetTitle("Number of Pixels Passing Cuts; Longitude; Number of Pixels");
    gTraceCounter.SetMarkerStyle(20);
    gTraceCounter.SetMarkerSize(2);
    gTraceCounter.Draw("AP");
    cstats->SaveAs("outputs/TraceCounter.png");

    TGraph gFirstPixelTime(fNumFiles,longitudes,FirstPixelTime);
    gFirstPixelTime.SetTitle("Time Location of First Pixel; Longitude; Time Bin");
    gFirstPixelTime.SetMarkerStyle(20);
    gFirstPixelTime.SetMarkerSize(2);
    gFirstPixelTime.Draw("AP");
    cstats->SaveAs("outputs/FirstPixelTime.png");

    TGraph gFirstPixelRow(fNumFiles,longitudes,FirstPixelRow);
    gFirstPixelRow.SetTitle("Row Location of First Pixel; Longitude; Row Number");
    gFirstPixelRow.SetMarkerStyle(20);
    gFirstPixelRow.SetMarkerSize(2);
    gFirstPixelRow.Draw("AP");
    cstats->SaveAs("outputs/FirstPixelRow.png");

    TGraph gFirstPixelColumn(fNumFiles,longitudes,FirstPixelColumn);
    gFirstPixelColumn.SetTitle("Column Location of First Pixel; Longitude; Column Number");
    gFirstPixelColumn.SetMarkerStyle(20);
    gFirstPixelColumn.SetMarkerSize(2);
    gFirstPixelColumn.Draw("AP");
    cstats->SaveAs("outputs/FirstPixelColumn.png");

    TGraph gIntegralsFromFit(fNumFiles,longitudes,IntegralsFromFit);
    gIntegralsFromFit.SetTitle("Total Integrated Light (FITS); Longitude; Total ADC Counts");
    gIntegralsFromFit.SetMarkerStyle(20);
    gIntegralsFromFit.SetMarkerSize(2);
    gIntegralsFromFit.Draw("AP");
    cstats->SaveAs("outputs/IntegralsFromFit.png");

    TGraph gIntegrals(fNumFiles,longitudes,Integrals);
    gIntegrals.SetTitle("Total Integrated Light (HIST'S); Longitude; Total ADC Counts");
    gIntegrals.SetMarkerStyle(20);
    gIntegrals.SetMarkerSize(2);
    gIntegrals.Draw("AP");
    cstats->SaveAs("outputs/Integrals.png");

    Color_t colorList[11] = {kPink,kRed,kOrange,kSpring+8,kTeal+8,kTeal,kCyan,kAzure-4,kBlue,kViolet+10,kBlack};
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Amplitude vs Time of Quality Fits; Time Bin;Amplitude (ADC)");
    TMultiGraph *mgintegral = new TMultiGraph();
    mgintegral->SetTitle("Integral vs Time of Quality Fits; Time Bin;Amplitude (ADC)");
    for(int k = 0; k < fNumFiles; k++){
      const int arraysize = (int)TraceCounter[k];

      //amplitude plots
      TGraph* gscatter = new TGraph(arraysize);
      for(int j = 0; j < arraysize; j++){
	gscatter->SetPoint(j,allFits[k][j]->GetParameter(1),allFits[k][j]->GetParameter(0));
      }
      gscatter->SetMarkerStyle(20);
      gscatter->SetMarkerSize(1);
      mg->Add(gscatter);
      TString gscattertitle("Amplitude vs Time of Quality Fits");gscattertitle+=longitudes[k];gscattertitle+="; Time Bin; Amplitude (ADC)";
      gscatter->SetTitle(gscattertitle);
      gscatter->Draw("AP");
      TString gscattername("outputs/Scatter");gscattername+=longitudes[k];gscattername+=".png";
      cstats->SaveAs(gscattername);

      //integral plots
      TGraph* gscatterintegral = new TGraph(arraysize);
      for(int j = 0; j < arraysize; j++){
	gscatterintegral->SetPoint(j,allFits[k][j]->GetParameter(1),allFits[k][j]->Integral(0,3000));
      }
      gscatterintegral->SetMarkerStyle(20);
      gscatterintegral->SetMarkerSize(1);
      mgintegral->Add(gscatterintegral);
      TString gscattertitleintegral("Integral vs Time of Quality Fits");gscattertitleintegral+=longitudes[k];gscattertitleintegral+="; Time Bin; Integral (ADC)";
      gscatterintegral->SetTitle(gscattertitleintegral);
      gscatterintegral->Draw("AP");
      TString gscatternameintegral("outputs/ScatterIntegral");gscatternameintegral+=longitudes[k];gscatternameintegral+=".png";
      cstats->SaveAs(gscatternameintegral);
    }
    mg->Draw("AP");
    cstats->SaveAs("outputs/Scatter.png");    
    mgintegral->Draw("AP");
    cstats->SaveAs("outputs/ScatterIntegral.png");    
    
    //overlay of same column as first pixel on same page. 
    TCanvas * cfits = new TCanvas("cfits","cfits",1000,800);
    gStyle->SetOptStat(0);
    for(int k = 0; k < fNumFiles; k++){
      const int arraysize = (int)TraceCounter[k];
      TString  hCanvasName("Column ");hCanvasName+= FirstPixelColumn[k]; hCanvasName+=" : Longitude = ";hCanvasName+=longitudes[k];hCanvasName+=";Time Bins (100ns);ADC Counts";
      TH1F hCanvas("hCanvas",hCanvasName,100,0,2000);
            hCanvas.GetYaxis()->SetRangeUser(0,2000);
      hCanvas.GetYaxis()->SetTitleOffset(1.4);
      hCanvas.Draw();
      TLegend lcstack(.15,.45,.35,.85);
      lcstack.SetTextSize(0.03);
      lcstack.SetTextFont(52);
      //      gStyle->SetPalette(53);
      int colorindex = 0;
      for(int j = 0; j < arraysize; j++){
	if(allFitsColumns[k][j] != FirstPixelColumn[k]) continue;
	allFits[k][j]->SetLineColor(colorList[colorindex]);
	//	allFits[k][j]->SetLineColor(gStyle->GetColorPalette(5*colorindex));
	colorindex++;
	allFits[k][j]->SetLineWidth(3);
	allFits[k][j]->SetRange(200,2000);
	allFits[k][j]->Draw("same");
	TString lcstackName("Row ");lcstackName+=allFitsRows[k][j];
	lcstack.AddEntry(allFits[k][j],lcstackName,"l");
      }
      TString cstackName("traces/traces_");cstackName+=longitudes[k];cstackName+=".png";
      lcstack.Draw();
      cfits->Update();
      cfits->SaveAs(cstackName.Data());  
    }

    //save all in mem
    outputPlots->Write();
    outputPlots->Close();
	  
  }
  return eSuccess;
}

VModule::ResultFlag 
ComparaisonStudy::Finish() 
{
  
  return eSuccess;
}

VModule::ResultFlag
ComparaisonStudy::GlueTrace(int iTagtmp, TH1F* theTrace, int telIdtmp, int pixelIdtmp) 
{
  double PreviousBinContent = 0;

  for(int iBin = 1; iBin <= theTrace->GetNbinsX(); iBin++){
    double BinContent = theTrace->GetBinContent(iBin);
    double normalizedBinContent;
    
    //if (fEventCounter == 1 && iBin > theTrace->GetNbinsX() - 10 ) continue;
    //    if( BinContent < PreviousBinContent - 1000) continue;
    fhRawPixel[iTagtmp][telIdtmp-1][pixelIdtmp-1]->SetBinContent(iBin+int((fEventCounter-1)*1000./((fMaxBinTrace-fMinBinTrace)/fNBinTrace)),BinContent);
    PreviousBinContent = BinContent;
    
    
  }
  return eSuccess;
}

Double_t ComparaisonStudy::asymGaussian(Double_t *x, Double_t *par)
{
    Double_t y=0.,xr=0.;
    if (x[0]<=par[1]) {xr = (x[0]-par[1])/par[2]; y = par[0]*exp(-0.5*xr*xr);}
    if (x[0]>par[1]) {xr = (x[0]-par[1])/(par[2]*par[3]); y = par[0]*exp(-0.5*xr*xr);}
    return y;
}
Double_t ComparaisonStudy::skewedGaussian(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t xr = (xx-par[1])/par[2];
  return (par[0]*exp(-0.5*xr*xr)*(1+TMath::Erf(par[3]*(xx-par[1])/(TMath::Sqrt(2)*par[2]))));
}



}