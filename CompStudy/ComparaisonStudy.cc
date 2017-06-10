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
#include <TTree.h>
#include <TGraph.h>
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
  fDataEventCounter = 0;
  fSimEventCounter = 0;
  outputPlots = new TFile("traces/FullTraces.root","recreate");

  for(int i = 0; i < fNRows; i++){  
    //same title for both data and simulation
    TString hPixelRowTitle("Row ");hPixelRowTitle+= (i+1);hPixelRowTitle+= ";Time Bin (100 ns);Column Number;ADC Counts / 20 Time Bins";
    
    //initialize data first
    TString hPixelRowNameDat("hDataPixelRow");hPixelRowNameDat+= (i+1);
    hPixelRow[0][i] = new TH2F(hPixelRowNameDat,hPixelRowTitle,100,0,999,fNColumns,0,fNColumns);
    //intialize simulation
    for(int j=1; j < fNumFiles; j++){
      TString hPixelRowNameSim("hSimPixelRow");hPixelRowNameSim+= (i+1);
      hPixelRowNameSim+="_";hPixelRowNameSim+=j;
      hPixelRow[j][i] = new TH2F(hPixelRowNameSim,hPixelRowTitle,100,0,999,fNColumns,0,fNColumns);
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

      //and now for data
      TString hnameDat("hDat_tel");hnameDat+=i+1;hnameDat+="_r";hnameDat+=detTelGlobal.GetPixel(j+1).GetRow();hnameDat+="_c";hnameDat+=detTelGlobal.GetPixel(j+1).GetColumn();
      fhRawPixel[0][i][j] = new TH1F(hnameDat,hnameDat,100,0,999);

      //intialize storing for sim
      for(int k=1; k < fNumFiles; k++){      
	TString hnameSim("hSim");hnameSim+=k;hnameSim+="_tel";
	hnameSim+=i+1;hnameSim+="_r";hnameSim+=detTelGlobal.GetPixel(j+1).GetRow();hnameSim+="_c";hnameSim+=detTelGlobal.GetPixel(j+1).GetColumn();
	fhRawPixel[k][i][j] = new TH1F(hnameSim,hnameSim,100,0,999);
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

  fEventCounter++;
  //to deal with sim and data in input file
  if(!(header.GetId() == "eye1_run1_event1")) {
    fData = true;
    if (!(header.GetId() == "eye1_run5254_event6336")) return eSuccess;
    fDataEventCounter++;
    INFO("New Data Event!");
  }else{
    fData = false;
    INFO("New SIM Event!");
    fSimEventCounter++;
  }
  cout << fDataEventCounter << " " << fSimEventCounter << endl;
  //  if(fDataEventCounter > 1 || fSimEventCounter > 1) return eSuccess;
  
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


    cout << eyeEvent.GetEventHeader()->GetEventNo() << " " << eyePixelList->GetNumPixels() << " " <<  eyePixelList->GetEyeNo() << endl;
    
    TString fName("traces.root");
    TFile *fTMP;
    if( fEventCounter ==1) {
      fTMP = new TFile(fName,"recreate");
    }else{
      fTMP = new TFile(fName,"update");
    }

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
      
      for (unsigned int pos = startBin; pos <startBin+280; ++pos) baseline += int(FADCDataWordGetData(&fadcword[pos]))/280.;
      
      
      for (unsigned int pos = startBin; pos <startBin+280; ++pos) baselineRMS += pow(int(FADCDataWordGetData(&fadcword[pos])),2)/280.;
      
      baselineRMS -= pow(baseline,2);
      if (baselineRMS>0) baselineRMS=sqrt(baselineRMS);
      else baselineRMS=0.;


      TString hName("h");
      if(fData) {
	hName+="Data";
	hName+=fDataEventCounter;
      }
      if(!fData){
	hName+="Sim";
	hName+=fSimEventCounter;
      }
      hName+="_m";hName+=mirrorId;hName+="_r";hName+=detTel.GetPixel(pixelId).GetRow();hName+="_c";hName+=detTel.GetPixel(pixelId).GetColumn();
      //      cout << hName.Data() << endl;
      TH1F* trace = new TH1F(hName,hName,100,startBin,endBin);
      
      for (unsigned int pos = startBin; pos <= endBin; ++pos) {
	charge = FADCDataWordGetData(&fadcword[pos]);
	trace->SetBinContent(int(pos/10.),charge-baseline);	  
      }//trace loop

      //      cout << iTag << " " << trace << " " << mirrorId << " " << pixelId << endl;

      //this itag below is needed as we are trying to select one event in a long list of data event and we are setting it as the first element in the array of histograms.
      //it turns out fSimCounter is alread at least 1 when we are here...
      if(fData)iTag = 0;
      if(!fData)iTag = fSimEventCounter;
      GlueTrace(iTag,trace, mirrorId, pixelId);//here glue means save in memory. 

    }//pixel loop
    fTMP->Write();
    fTMP->Close();

    }//mirror or eye loop?

  //everything below needs to happen only once fhRawPixel is fully filled!
  if(!(fDataEventCounter+fSimEventCounter == fNumFiles))return eSuccess;


  //here I am looking  at creating the 2d histos for each rows of each files
  const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);

  
  for(int i = 0; i < fNTels; i++){
    for(int j = 0; j < fNPixels; j++){
      for(int k = 0; k < fNumFiles; k++){
	//know we can treat the full array the same way (data and sim)
     	if(fhRawPixel[k][i][j]->GetEntries() == 0){
	  delete fhRawPixel[k][i][j];
	  continue;
	}
	for(int iBin = 1; iBin <= fhRawPixel[k][i][j]->GetNbinsX(); iBin++){ 
	  hPixelRow[k][detTelGlobal.GetPixel(j+1).GetRow()-1]->SetBinContent(iBin,detTelGlobal.GetPixel(j+1).GetColumn(),fhRawPixel[k][i][j]->GetBinContent(iBin));
	} 
      }
    } 
  }  
  

  //here I make plot with superimposed traces for data and simulaton.
  //  double I0[fNumFiles-1] = {8e04,8.2e04,8.4e04,8.6e04,8.8e04,9e04,9.2e04,9.4e04,9.6e04,9.8e04,1e05,1.02e05,1.04e05,1.06e05,1.08e05,1.1e05,1.12e05,1.14e05,1.16e05,1.18e05,1.2e05,1.22e05,1.24e05,1.26e05,1.28e05,1.3e05,1.32e05};
  //  double I0[fNumFiles-1] = {9.6e04,9.8e04,1e05,1.02e05,1.04e05,1.06e05,1.08e05,1.1e05,1.12e05,1.14e05,1.16e05};
  //I0 used for altitude too
  double I0[fNumFiles-1] = {2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000};
  Color_t colorList[11] = {kPink,kRed,kOrange,kSpring+8,kTeal+8,kTeal,kCyan,kAzure-4,kBlue,kViolet+10,kBlack};

  TCanvas *cstack = new TCanvas("cstack","cstack",1000,800);
  gStyle->SetOptStat(0);
  TF1* fits[fNumFiles][fNRows];//this will save the fits done below, readjust
  for(int k = 1; k < fNumFiles; k++){
    //    TString  hCanvasName("Column 10: Peak Current = ");hCanvasName+=I0[k-1];
    TString  hCanvasName("Column 10: Source Altitude = ");hCanvasName+=I0[k-1];
    hCanvasName+="m;Time Bins (100ns);ADC Counts";
    TH1F hCanvas("hCanvas",hCanvasName,100,0,999);
    hCanvas.GetYaxis()->SetRangeUser(0,800);
    hCanvas.GetYaxis()->SetTitleOffset(1.4);
    hCanvas.Draw();
    TLegend lcstack(.15,.60,.35,.85);
    TLegend lcstack2(.15,.45,.35,.55);
    lcstack.SetTextSize(0.03);
    lcstack.SetTextFont(52);
    lcstack2.SetTextSize(0.03);
    lcstack2.SetTextFont(62);
    lcstack.SetBorderSize(0);
    lcstack2.SetBorderSize(0);
    int traceCounter = 0;
    for(int j = 6; j <= 10; j++){
      TString htmpName("htmp");htmpName+=k;htmpName+=j;
      TH1F *htmp= new TH1F(htmpName,"htmp",100,0,999);
      for(int iBinx = 1; iBinx <= hPixelRow[k][j]->GetNbinsX(); iBinx++){
	htmp->SetBinContent(int(iBinx),hPixelRow[k][j]->GetBinContent(iBinx,10));//choosing column 10
      }
      TF1* g1 = new TF1(htmpName+="fit","gaus");
      if(j==10 || j==9) g1->SetRange(200,500);
      if(j==8) g1->SetRange(300,600);
      if(j==7) g1->SetRange(400,800);
      if(j==6) g1->SetRange(600,999);
      htmp->Fit(g1,"RNQ");
      g1->SetLineColor(colorList[int(2*(j-6))]);
      g1->SetLineWidth(3);
      g1->SetLineStyle(2);
      g1->SetRange(200,999);

      if(g1->GetParameter(0) > g1->GetParError(0)) g1->Draw("same");
      
      TString htmpDatName("htmpDat");htmpDatName+=k;htmpDatName+=j;
      TH1F *htmpDat= new TH1F(htmpDatName,"htmpDat",100,0,999);
      for(int iBinx = 1; iBinx <= hPixelRow[0][j]->GetNbinsX(); iBinx++){
	htmpDat->SetBinContent(int(iBinx),hPixelRow[0][j]->GetBinContent(iBinx,10));
      }
      TF1* g2 = new TF1(htmpDatName+="fit","gaus",200,999);
      htmpDat->Fit(g2,"RNQ");
      g2->SetLineColor(colorList[int(2*(j-6))]);
      g2->SetLineWidth(3);
      g2->Draw("same");
      TString lcstackName("Row ");lcstackName+=j+1;
      lcstack.AddEntry(g2,lcstackName,"l");

      if(traceCounter ==0){
	TF1* g1tmp = (TF1*)g1->Clone();
	TF1* g2tmp = (TF1*)g2->Clone();
	g1tmp->SetLineColor(kBlack);
	lcstack2.AddEntry(g1tmp,"Simulation","l");
	g2tmp->SetLineColor(kBlack);
	lcstack2.AddEntry(g2tmp,"Data","l");
      }
      traceCounter ++;
      if(traceCounter < 6) fits[0][j] = g2;//only need to save data fit for one go across pixels
      fits[k][j] = g1; 
    }

    TString cstackName("traces/traces_");cstackName+=I0[k-1];cstackName+=".png";
    lcstack.Draw();
    lcstack2.Draw();
    cstack->Update();
    cstack->SaveAs(cstackName.Data());  
  }
  

  //below is the analysis on the fits done above...
  TCanvas* cfits = new TCanvas("cfits","cfits",1000,800);
  TMultiGraph *mgFitAmplitude = new TMultiGraph();
  TMultiGraph *mgFitTimeDiff = new TMultiGraph();
  TMultiGraph *mgFitChi2 = new TMultiGraph();
  int plotCounter = 0;
  TLegend lcamplitude(.15,.60,.35,.85);
  TLegend lcamplitude2(.15,.45,.35,.55);
  lcamplitude.SetTextSize(0.03);
  lcamplitude.SetTextFont(52);
  lcamplitude2.SetTextSize(0.03);
  lcamplitude2.SetTextFont(62);
  lcamplitude.SetBorderSize(0);
  lcamplitude2.SetBorderSize(0);
  for(int j = 6; j <= 10; j++){
    //data
    TGraphErrors* gFitDataAmplitude = new TGraphErrors(1);
    double peakMaximumData = fits[0][j]->GetParameter(0);
    double peakMeanData = fits[0][j]->GetParameter(1);
    double peakStdDevData = fits[0][j]->GetParameter(2);
    if (peakMaximumData > 0){
      gFitDataAmplitude->SetPoint(1,7200,peakMaximumData);
      gFitDataAmplitude->SetPointError(1,0,fits[0][j]->GetParError(0));
    }
    gFitDataAmplitude->SetMarkerColor(colorList[int(2*(j-6))]);
    gFitDataAmplitude->SetMarkerSize(2);
    gFitDataAmplitude->SetMarkerStyle(23);
    mgFitAmplitude->Add(gFitDataAmplitude);
    TString lcstackNameDat("Data Guess");
    if(plotCounter==0){
      TGraphErrors* gTMPAmplitude = (TGraphErrors*)gFitDataAmplitude->Clone();
      gTMPAmplitude->SetMarkerColor(kBlack);
      lcamplitude2.AddEntry(gTMPAmplitude,lcstackNameDat,"p");
    }

    
    //simulation
    TGraphErrors* gFitAmplitude = new TGraphErrors(fNumFiles-1);
    TGraph* gFitChi2 = new TGraph(fNumFiles-1);
    TGraphErrors* gFitTimeDiff = new TGraphErrors(fNumFiles-1);
    for(int k = 1; k < fNumFiles; k++){

      //chi2 analysis
      double chiSquareFits = fits[k][j]->GetChisquare()/fits[k][j]->GetNDF();//get reduced chi2

      //amplitude and timing analysis
      double amplitudeFits = fits[k][j]->GetParameter(0);//get amp
      double amplitudeFitsErrors = fits[k][j]->GetParError(0);//get amp
      double meanFits = fits[k][j]->GetParameter(1);//get mean
      double meanFitsErrors = fits[k][j]->GetParError(1);//get mean
      double stddevFits = fits[k][j]->GetParameter(2);//get stddev
      double stddevFitsErrors = fits[k][j]->GetParError(2);//get stddev
      cout << k+1 << " " << j+1 << " " << amplitudeFits << " "  << meanFits << " " << stddevFits << endl;
      cout << k+1 << " " << j+1 << " " << amplitudeFitsErrors << " "  << meanFitsErrors << " " << stddevFitsErrors << endl<< endl;
      if (amplitudeFits > 20 && meanFits > 3*stddevFits && amplitudeFits > amplitudeFitsErrors) {
	gFitChi2->SetPoint(k-1,I0[k-1],chiSquareFits);
	gFitAmplitude->SetPoint(k-1,I0[k-1],amplitudeFits);
	gFitAmplitude->SetPointError(1,0,amplitudeFitsErrors);	
	gFitTimeDiff->SetPoint(k-1,I0[k-1],meanFits-peakMeanData);
	gFitTimeDiff->SetPointError(1,0,0);	
      }
    }
    gFitAmplitude->SetMarkerColor(colorList[int(2*(j-6))]);
    gFitAmplitude->SetMarkerStyle(20);
    gFitTimeDiff->SetMarkerColor(colorList[int(2*(j-6))]);
    gFitTimeDiff->SetMarkerStyle(20);
    gFitChi2->SetMarkerColor(colorList[int(2*(j-6))]);
    gFitChi2->SetMarkerStyle(20);
    mgFitAmplitude->Add(gFitAmplitude);
    mgFitTimeDiff->Add(gFitTimeDiff);
    mgFitChi2->Add(gFitChi2);

    //legend fill and aplicable to all three plots
    TString lcstackName("Simulation Row ");lcstackName+=j+1;
    lcamplitude.AddEntry(gFitAmplitude,lcstackName,"p");
    
    plotCounter++;
  }
  
  //save amplitude plots
  //range for I0.
  // double xRangeMin = 78000;
  // double xRangeMax = 134000;
  //range for alt.
  double xRangeMin = 2000;
  double xRangeMax = 8500;
  mgFitAmplitude->SetTitle("Amplitude Variation for Column 10; Altitude (m); Peak Maximum");
  //  mgFitAmplitude->SetTitle("Amplitude Variation for Column 10; Peak Current (A); Peak Maximum");
  mgFitAmplitude->Draw("AP");
  mgFitAmplitude->GetYaxis()->SetTitleOffset(1.3);
  mgFitAmplitude->GetXaxis()->SetRangeUser(xRangeMin,xRangeMax);
  mgFitAmplitude->GetYaxis()->SetRangeUser(0,800);
  lcamplitude2.Draw();
  lcamplitude.Draw();
  cfits->Update();
  cfits->SaveAs("traces/gAmplitudes.png");

  //save timediff plots
  //  mgFitChi2->SetTitle("Simulation Reduced #chi^{2} for Column 10; Peak Current (A); #chi^{2}");
  mgFitChi2->SetTitle("Simulation Reduced #chi^{2} for Column 10; Altitude (m); #chi^{2}");
  mgFitChi2->Draw("AP");
  mgFitChi2->GetYaxis()->SetTitleOffset(1.3);
  mgFitChi2->GetYaxis()->SetRangeUser(0,20);
  mgFitChi2->GetXaxis()->SetRangeUser(xRangeMin,xRangeMax);
  lcamplitude.Draw();
  cfits->Update();
  cfits->SaveAs("traces/gChi2.png");

  //save timediff plots
  mgFitTimeDiff->SetTitle("Peak Time Difference w/ Data for Column 10; Altitude (m); Time Difference (Time Bins of 100ns)");
  //  mgFitTimeDiff->SetTitle("Peak Time Difference w/ Data for Column 10; Peak Current (A); Time Difference (Time Bins of 100ns)");
  mgFitTimeDiff->Draw("AP");
  mgFitTimeDiff->GetYaxis()->SetTitleOffset(1.3);
  mgFitTimeDiff->GetYaxis()->SetRangeUser(-200,0);
  mgFitTimeDiff->GetXaxis()->SetRangeUser(xRangeMin,xRangeMax);
  // lcamplitude.SetX1NDC(.55);
  // lcamplitude.SetY1NDC(.15);
  // lcamplitude.SetX2NDC(.75);
  // lcamplitude.SetY2NDC(.40);
  lcamplitude.Draw();
  cfits->Update();
  cfits->SaveAs("traces/gTimeDiff.png");

  
  //below one is studying the residuals of the 2 d histogram and for pixels. 
  TH2F *hResidual[fNumFiles][fNRows];  
  for(int iRow = 0; iRow < fNRows; iRow++){
    for(int k = 0; k < fNumFiles; k++){
      TString hResidualName("hResidual");hResidualName+=k;hResidualName+="_";hResidualName+=iRow+1;
      TString hResidualTitle("Residual for Row ");hResidualTitle+=iRow+1;
      //get the data histogram
      
      hResidual[k][iRow] = (TH2F*) hPixelRow[0][iRow]->Clone();
      hResidual[k][iRow]->SetNameTitle(hResidualName, hResidualTitle);
      hResidual[k][iRow]->Add(hPixelRow[k][iRow],-1.);
    }
  }

  //plotting the residual for all selected rows in 2D hist form
  cout << "PLOTTING THE RESIDUALS..." << endl;
  TCanvas* cglues = new TCanvas("cglues","cglues",1000,800);
  for(int i=6; i <= 10; i++){
    for(int k = 0; k < fNumFiles; k++){
      hResidual[k][i]->Draw("colz");
      hResidual[k][i]->GetZaxis()->SetTitleOffset(1.02); 
      hResidual[k][i]->GetZaxis()->SetTitle("Residual ADC/2#mus"); 
      cglues->Update();
      
      TString fNameGlues("glues/");fNameGlues+=hResidual[k][i]->GetName();
      fNameGlues+=".png";
      //cglues->SaveAs(fNameGlues.Data());  
    }
  }

  
  TCanvas* cresidual = new TCanvas("cresidual","cresidual",1000,800);
  TMultiGraph *mg = new TMultiGraph();
  TLegend* lc = new TLegend(.15,.15,.25,.45);
  
  for(int i=6; i <= 10; i++){  //look through rows 6-10
    TGraph* gResidual = new TGraph(fNumFiles-1);
    gResidual->SetMarkerStyle(21);
    gResidual->SetMarkerColor(colorList[i-5]);
    gResidual->SetMarkerSize(2.0);
    for(int k = 1; k < fNumFiles; k++){//skip the 0 value from data/data
      double residualTotal=0;
      for(int iBinx=1;iBinx<=hResidual[k][i]->GetNbinsX(); iBinx++){
	residualTotal += hResidual[k][i]->GetBinContent(iBinx, 10);//look at c10
      }
      residualTotal = residualTotal/hResidual[k][i]->GetNbinsX();   
      //  cout << I0[k-1] << " " << residualTotal << endl;
      gResidual->SetPoint(k-1,I0[k-1],residualTotal);
    }
    
    TString lcName("Row");lcName+=i+1;
    mg->Add(gResidual);
    lc->AddEntry(gResidual,lcName,"p");
    
  }
  //  mg->SetTitle("Residual Average Per Pixel: Column 10; Peak Current (A); Data-Simulation Residual");
  mg->SetTitle("Residual Average Per Pixel: Column 10; Altitude (m); Data-Simulation Residual");
  mg->Draw("AP");
  mg->GetYaxis()->SetTitleOffset(1.3);
  lc->Draw("same");
  cresidual->Update();
  cresidual->SaveAs("traces/gResidual.png");
  
  
  
  outputPlots->Write();
  outputPlots->Close();
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
  // cout << theTrace->GetNbinsX() << " " << fhRawPixel[iTagtmp][telIdtmp-1][pixelIdtmp-1]->GetNbinsX() <<endl;
  for(int iBin = 1; iBin <= theTrace->GetNbinsX(); iBin++){
    double normalizedBinContent = theTrace->GetBinContent(iBin);
    fhRawPixel[iTagtmp][telIdtmp-1][pixelIdtmp-1]->SetBinContent(iBin,normalizedBinContent);
  }
  return eSuccess;
}




}
