#include "DoubleELVESAnalysis.h"
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
#include "TVirtualFFT.h"
#include <TGraph.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TH2D.h>
#include <THStack.h>
#include <TPaletteAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include "TSpectrum.h"
#include <TVector3.h>
#include <vector>
#include <AugerEvent.h>
#include <EyeEvent.hh>
#include <EyePixelList.hh>
#include <FadcData.hh>
#include "TPolyMarker.h"

using namespace std;
using namespace utl;
using namespace fwk;
using namespace evt;
using namespace fevt;
using namespace fdet;
using namespace det;
using CLHEP::RandFlat;

namespace DoubleELVESAnalysis {

// Initialize your static data members.

// Define the methods preferibly in the same order as listed in the header file.

DoubleELVESAnalysis::DoubleELVESAnalysis() 
{
  // Initialize your class.
  // Prefer initialization lists over = to 
  // initialize your members.
}

DoubleELVESAnalysis::~DoubleELVESAnalysis()
{
}

VModule::ResultFlag 
DoubleELVESAnalysis::Init()
{
  INFO("=====================>>>>>>>>>>>>>>>> Init!");
  hdeltaTMedian = new TH1F("deltaTMedian","Median Difference; time bin (100ns); counts", 100, 200, 2000);
  hdeltaTMean = new TH1F("deltaTMean","Mean Time Difference; time bin (100ns); counts", 100, 200, 2000);
  hdeltaTMeanCenter = new TH1F("deltaTMeanCenter","Mean Time Difference for Center Column; time bin (100ns); counts", 100, 200, 2000);

  fBinning = 2; //in us
  fNPages = 10; //make it longer than necessary
  
  fEventCounter = 0;
  fEventTracker = 0;
  fDataEventCounter = 0;
  fRebin = 10*fBinning;
  fNBinTrace = 1000./(fRebin)*fNPages; //num bins in full trace --> 10th of detector
  fMaxBinTrace = 1000*fNPages; //upper bin value in full trace
  fMinBinTrace = 1;//lower bin value in full trace
  fNBinPage = 1000./(fRebin) ; //num bin per page
  cout << fRebin << " " << fNBinTrace << " " << fNBinPage << " " << fMinBinTrace << " " << fMaxBinTrace << endl;

  
  int BinSize = int(((fMaxBinTrace-fMinBinTrace)/fNBinTrace));
  const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);


  outputPlots = new TFile("traces/FullTraces.root","recreate");
  outputPlots->Close();
  //create the template to hold the long traces of all events. 
  for(int i = 0; i < fNTels; i++){
    for(int j = 0; j < fNPixels; j++){
      TString hname("m");hname+=i+1;hname+="_r";hname+=detTelGlobal.GetPixel(j+1).GetRow();hname+="_c";hname+=detTelGlobal.GetPixel(j+1).GetColumn();
      hnameBase[i][j].Append(hname);
      fhRawPixel[i][j] = new TH1F(hname,hname,fNBinTrace,fMinBinTrace,fMaxBinTrace);
    }
  }
  return eSuccess;
}

VModule::ResultFlag 
DoubleELVESAnalysis::Run(evt::Event& event)
{

  AugerEvent& rawEvent = event.GetRawEvent();

  
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


    cout << eyeEvent.GetEventHeader()->GetEventNo() << " " << eyePixelList->GetNumPixels() << " " <<  eyePixelList->GetEyeNo() << " " << eyeHeader->GetTimeStamp()->GetSec() + 1e-9*eyeHeader->GetTimeStamp()->GetNanoSec()<<endl;

    //only look at first page.
    // for (int pageNum = 1; pageNum < 9; pageNum++){
    //   if ( eyeEvent.GetEventHeader()->GetEventNo() == fEventCounter+pageNum) return eSuccess;
    // } fEventCounter = eyeEvent.GetEventHeader()->GetEventNo();

    //look at first two pages only...170us is 25.5 km alt lightning, first need to keep track of things...
    //tag first event in the series of 2-10 pages
    // if(fEventCounter == 0 || (fEventCounter > 0 && eyeEvent.GetEventHeader()->GetEventNo() == fEventTracker+1)){
    //   fEventCounter++;
    // }
    // if (fEventCounter > 0 && eyeEvent.GetEventHeader()->GetEventNo() != fEventTracker+1){
    // }

    // if (fEventCounter > 0 && eyeEvent.GetEventHeader()->GetEventNo() != fEventTracker+1){
    //   fEventTracker = eyeEvent.GetEventHeader()->GetEventNo();
    //   fEventCounter = 0;
    // }

    if( eyeEvent.GetEventHeader()->GetEventNo() == fEventTracker+1 ){
      fEventCounter++;
    }else{
      fDataEventCounter++;
      fEventCounter = 0;
    }
    cout << fEventCounter << " " << fEventTracker  <<  " " << eyeEvent.GetEventHeader()->GetEventNo()<< " " << fDataEventCounter << endl;
    fEventTracker = eyeEvent.GetEventHeader()->GetEventNo();

    if(fEventCounter == 0) {
      // Save previous trace once we looped through one event. Make sure to put this once more in the finish to look at last event. 
      if(fDataEventCounter > 1){
	Analysis();
      }

      // SetUpFullPixelTrace
      INFO("New Raw Event!\n");
      for(int i = 0; i < fNTels; i++){
	for(int j = 0; j < fNPixels; j++){
	  TString hnameUpdate("hRun");
	  hnameUpdate+=eyeEvent.GetEventHeader()->GetRunNo();
	  hnameUpdate+="_Evt";
	  hnameUpdate+=fEventTracker;
	  hnameUpdate+="_";
	  hnameUpdate+=hnameBase[i][j];
	  fhRawPixel[i][j]->SetName(hnameUpdate);
	  //clean histogram for reuse at each event. 
	  for(int iBin = 1; iBin < fhRawPixel[i][j]->GetNbinsX()+1; iBin++){
	    fhRawPixel[i][j]->SetBinContent(iBin,0.0);
	  }
	}
      }
    }
    // if(fEventCounter > 1) {
    //   INFO("Skipping Extra Pages.");
    //   return eSuccess;
    // }

    TString fName("/home/kswiss/Workspace/workoffline/analysis/DoubleELVESAnalysis/traces/");
    fName+="run";
    fName+=eyeEvent.GetEventHeader()->GetRunNo();
    fName+=".event";
    fName+=eyeEvent.GetEventHeader()->GetEventNo();
    fName+=".eye";
    fName+=eyeId;
    fName+=".root";
    TFile fTMP(fName,"recreate");
    
    int mirrorIdGlobal;//each eye event has all pixels lit in one list... multiple mirros as well... This global id is the first mirror lit.
    
    for (unsigned int iPixel = 0; iPixel < eyePixelList->GetNumPixels(); ++iPixel) {
      
      FdUtil::Fd::FdPixelNumber pixelNumber = eyePixelList->GetPixel(iPixel);
      const unsigned int mirrorId = FdUtil::Fd::GetEyeMirrorNo(pixelNumber);
      if(iPixel==0) mirrorIdGlobal = mirrorId;
      const unsigned int channelId = FdUtil::Fd::GetMirrorPixelNo(pixelNumber);
      //     const unsigned int pixelId = FdUtil::Fd::GetPixelNo(pixelNumber);
      
      const TFADCData* fadcData = eyeFADCData->GetFADCData(pixelNumber);
      
      fevt::Telescope& tel = eye.GetTelescope(mirrorId);
      const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
      const fdet::Eye& detEye = detFD.GetEye(eyeId);
      //cout << " detEye " << eyeId<<endl;
      const fdet::Telescope& detTel = detEye.GetTelescope(mirrorId);
      
      const fdet::Channel& thisChannel = detTel.GetChannel(channelId);
      const unsigned int pixelId = thisChannel.GetPixelId();
      
      
      if (!fadcData || thisChannel.IsVirtual()) {
	continue;
      }
      
      const unsigned int startBin = fadcData->GetTraceStartBin();
      const unsigned int endBin = fadcData->GetTraceEndBin();
      
      TFADCData::FADCDataWord fadcword = fadcData->GetFADCTrace();
      
      double baseline = 0.;
      double baselinetail=0.;
      double charge=0.;
      
      for (unsigned int pos = startBin; pos < startBin+100; ++pos) baseline += int(FADCDataWordGetData(&fadcword[pos]))/100.;
      for (unsigned int pos = endBin-100; pos < endBin; ++pos) baselinetail += int(FADCDataWordGetData(&fadcword[pos]))/100.;
      
      TString hName("h"); hName+=fEventCounter;hName+="_m";hName+=mirrorId;hName+="_r";hName+=detTel.GetPixel(pixelId).GetRow();hName+="_c";hName+=detTel.GetPixel(pixelId).GetColumn();
      TH1F* trace = new TH1F(hName,hName,1000,startBin,endBin);
      
      for (unsigned int pos = startBin; pos <= endBin; ++pos) {
	charge = FADCDataWordGetData(&fadcword[pos]);
 	if (pixelId <= 440 && abs(baseline-baselinetail) > 10 && baseline >baselinetail){
	  //	  trace->SetBinContent(int(pos/((fMaxBinTrace-fMinBinTrace)/fNBinTrace)),charge-baselinetail);
	  trace->SetBinContent(pos,charge-baselinetail);	  
	  
	}else{
	  trace->SetBinContent(pos,charge-baseline);	  
	}
      }
      trace->Rebin(fRebin);
      
      GlueTrace(fEventCounter,trace,mirrorId,pixelId);
      //      delete trace;
    }//pixel loop
    fTMP.Write();
    fTMP.Close();
  }//eye loop 

  return eSuccess;
  

}

  //this function is called when fEventCounter = 0, before the fHRawPixel are cleared
VModule::ResultFlag
DoubleELVESAnalysis::Analysis()
{
  //sample info to get col and row from pixelId
  // const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  //  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);

  vector<double> timeDifference;
  vector<double> timeFirstPeak;
  vector<int> columnFirstPeak;
  struct ELVESTraceData {
    double t0;//time of First Peak
    double dt;// time difference between the two peaks
    int colnum;// column number
    int rownum;//row number
    int mirnum;//mirror num
    double chi2; ///chi2
  };
  vector <ELVESTraceData> FirstPeak;
  
  struct by_t0 { 
    bool operator()(ELVESTraceData const &a, ELVESTraceData const &b) { 
      return a.t0 < b.t0;
    }
  };

  int npages = 3;
  int nbins = npages*fNBinPage;//look at three pages
  
  int doubleELVESPixelCounter = 0;
  for(int i = 0; i < fNTels; i++){
    for(int j = 0; j < fNPixels; j++){
      int pixelId = j+1;
      TSpectrum *spectrum = new TSpectrum(100);
      TSpectrum *spectrumbg = new TSpectrum(100);
      Float_t *source = new Float_t[nbins];
      Float_t *dest = new Float_t[nbins];
      TString hName("");hName+=fhRawPixel[i][j]->GetName();
      TString spectrumName(hName);
      TString spectrumbgName(hName);
      spectrum->SetName(spectrumName.Append("_spectrum"));
      //      cout << spectrumName << endl;
      spectrumbg->SetName(spectrumbgName.Append("_spectrumbg"));
      TH1F *d = new TH1F(hName.Append("d"),"",nbins,fMinBinTrace,npages*1000);
      TH1F *hbraw = new TH1F(hName.Append("braw"),"",nbins,fMinBinTrace,npages*1000);
      for (int k = 0; k < nbins; k++) {
	source[k]=fhRawPixel[i][j]->GetBinContent(k+1);
	hbraw->SetBinContent(k+1,fhRawPixel[i][j]->GetBinContent(k+1));
      }
      Int_t nfound = spectrum->SearchHighRes(source, dest, nbins, 1.5, 30, kTRUE, 3, kTRUE, 3);
      TH1 *hb = spectrum->Background(hbraw,2,"same");
      hb->SetName(hName.Append("_background"));
      Int_t nfoundbg = spectrumbg->Search(hb,3,"",0.20);
      //      cout << nfound << " " << nfoundbg << endl;

      Float_t *xpeaks = spectrum->GetPositionX();
      Double_t fPositionX[100];
      Double_t fPositionY[100];
      for (int k = 0; k < nfound; k++) {
	int bin = 1 + Int_t(xpeaks[k] + 0.5);
	fPositionX[k] = fhRawPixel[i][j]->GetBinCenter(bin);
	fPositionY[k] = fhRawPixel[i][j]->GetBinContent(bin);
      }
      TPolyMarker * pm = (TPolyMarker*)fhRawPixel[i][j]->GetListOfFunctions()->FindObject("TPolyMarker");
      if (pm) {
	fhRawPixel[i][j]->GetListOfFunctions()->Remove(pm);
	delete pm;
      }

      pm = new TPolyMarker(nfound, fPositionX, fPositionY);
      fhRawPixel[i][j]->GetListOfFunctions()->Add(pm);
      pm->SetMarkerStyle(23);
      pm->SetMarkerColor(kRed);
      pm->SetMarkerSize(1.3);
      for (int k = 1; k <= nbins; k++) d->SetBinContent(k,dest[k]);
      d->SetLineColor(kRed);
  
      Double_t par[6];
      Double_t chi2tmp;
      if(nfoundbg == 2   &&  fPositionY[0] > 20 &&  fPositionY[1] > 20 && TMath::Abs(fPositionX[1]-fPositionX[0]) < 2700){
      // if(fhRawPixel[i][j]->GetMean()!=0){
	//	cout << "Found "<< nfound<<" candidate peaks to fit\n" << endl;
	//cout << "Found "<< nfoundbg<<" candidate bg peaks to fit\n" << endl;
	//	cout <<"Time Difference: " << fPositionX[1]-fPositionX[0] << endl;
	par[0]=fPositionY[0];par[1]=fPositionX[0];par[2]=50;
	par[3]=fPositionY[1];par[4]=fPositionX[1];par[5]=50;
	TF1 total("total","gaus(0)+gaus(3)",200,npages*1000);
	total.SetParameters(par);
	fhRawPixel[i][j]->Fit(&total,"RQ");
	//	cout << total.GetChisquare()/total.GetNDF() << endl;
	chi2tmp = total.GetChisquare()/total.GetNDF();
	outputPlots = new TFile("traces/FullTraces.root","update");
	fhRawPixel[i][j]->Write();
	hb->Write();
	outputPlots->Close();
      }else{
	// delete d;
	// delete hb;
	continue;         
      }
      
      doubleELVESPixelCounter++;
      //      cout << "============" << par[1] << " " << TMath::Abs(par[4]-par[1]) << " " << doubleELVESPixelCounter << endl;
      ELVESTraceData FirstPeakTmp;
      FirstPeakTmp.t0 = par[1];
      FirstPeakTmp.dt =TMath::Abs(par[4]-par[1]);
      FirstPeakTmp.colnum = detTelGlobal.GetPixel(pixelId).GetColumn();
      FirstPeakTmp.rownum = detTelGlobal.GetPixel(pixelId).GetRow();
      FirstPeakTmp.mirnum = i+1;
      FirstPeakTmp.chi2 = chi2tmp;
      timeFirstPeak.push_back(fPositionX[0]);
      columnFirstPeak.push_back(detTelGlobal.GetPixel(pixelId).GetColumn());
      timeDifference.push_back(TMath::Abs(fPositionX[1]-fPositionX[0]));
      FirstPeak.push_back(FirstPeakTmp);

      
    }//end pixel loop
  }
  //for reference
  // struct ELVESTraceData {
  //   double t0;//time of First Peak
  //   double dt;// time difference between the two peaks
  //   int colnum;// column number
  //   int rownum;//row number
  // };

  
  //figure this out
  //  int mirrorIdGlobal = 1;
  sort(FirstPeak.begin(),FirstPeak.end(),by_t0());
  double averageDeltaT = accumulate(timeDifference.begin(), timeDifference.end(), 0.0)/timeDifference.size();
  double averageMinColT=0;
  int averageCounter=0;

  TH2F * heatmapDeltaT[6];
  TH2F * heatmapChi2[6];
  for(int i=0; i < 6; i++){
    TString hNameheatmapDeltaT("");hNameheatmapDeltaT+=fhRawPixel[0][0]->GetName();
    hNameheatmapDeltaT.Remove(hNameheatmapDeltaT.Length()-9,hNameheatmapDeltaT.Length());
    hNameheatmapDeltaT+="_heatmapDeltaT_m";hNameheatmapDeltaT+=i+1;
    heatmapDeltaT[i] = new TH2F(hNameheatmapDeltaT,"heatmapDeltaT",20,1,20,22,1,22);
    TString hNameheatmapChi2("");hNameheatmapChi2+=fhRawPixel[0][0]->GetName();
    hNameheatmapChi2.Remove(hNameheatmapChi2.Length()-9,hNameheatmapChi2.Length());
    hNameheatmapChi2+="_heatmapChi2_m";hNameheatmapChi2+=i+1;
    heatmapChi2[i] = new TH2F(hNameheatmapChi2,"heatmapChi2",20,1,20,22,1,22);
  }
  if (FirstPeak.size() != 0){
    for(int i=0; i < FirstPeak.size(); i++){
      heatmapDeltaT[FirstPeak[i].mirnum-1]->Fill(FirstPeak[i].colnum,FirstPeak[i].rownum,FirstPeak[i].dt);
      heatmapChi2[FirstPeak[i].mirnum-1]->Fill(FirstPeak[i].colnum,FirstPeak[i].rownum,FirstPeak[i].chi2);
      TPaletteAxis *palette = (TPaletteAxis*)heatmapDeltaT[FirstPeak[i].mirnum-1]->GetListOfFunctions()->FindObject("palette");
      delete palette;
      //calculate average based on center column
      if(FirstPeak[i].colnum == FirstPeak[0].colnum){
	averageMinColT += FirstPeak[i].dt;
	averageCounter++;
      }
    }
    averageMinColT = averageMinColT / averageCounter;
  }
  outputPlots = new TFile("traces/FullTraces.root","update");
  for(int i=0; i < 6; i++){    
    //    heatmapDeltaT[i]->GetZaxis()->SetRangeUser(averageDeltaT-50,averageDeltaT+50);
    if(heatmapDeltaT[i]->GetMean()!=0) heatmapDeltaT[i]->Write();
    if(heatmapChi2[i]->GetMean()!=0) heatmapChi2[i]->Write();
  }
  outputPlots->Close();
  
  cout <<endl<< averageMinColT << " " << averageDeltaT << " " << doubleELVESPixelCounter<< endl<<endl;
  
      //making sure a full-bodied elves is lookedc at and take care of stereo
  if(doubleELVESPixelCounter >= 10){
    // if(EventTimeList.size()>=1){
    // 	for(int i = 0; i < EventTimeList.size(); i++){
    // 	  // if(abs((eyeHeader->GetTimeStamp()-EventTimeList[i]).double()) < 1e-6){
    // 	  //   cout << "ALERT" << endl;
    // 	  // }
    // 	}
    //      }
    // EventTimeList.push_back(eyeHeader->GetTimeStamp());
    hdeltaTMeanCenter->Fill(averageMinColT);
    hdeltaTMean->Fill(averageDeltaT);
    hdeltaTMedian->Fill(CalcMHWScore(timeDifference));
    //      fTMP.Write();
    //      fTMP.Close();
  }else{}// fTMP.Close();   fTMP.Delete();}
  
  
  return eSuccess;  
}


VModule::ResultFlag 
DoubleELVESAnalysis::Finish() 
{
  //analysis and save the last event
  Analysis();
  // outputPlots = new TFile("traces/FullTraces.root","update");
  // for(int i = 0; i < fNTels; i++){
  //   for(int j = 0; j < fNPixels; j++){
  //     if(fhRawPixel[i][j]->GetMean()!=0) fhRawPixel[i][j]->Write();
  //   }
  // }
  // outputPlots->Close();

  TCanvas* cT = new TCanvas("cT","cT",800,600);
  hdeltaTMedian->Draw();
  cT->SaveAs("DeltaTMedian.png");
  hdeltaTMean->Draw();
  cT->SaveAs("DeltaTMean.png");
  hdeltaTMeanCenter->Draw();
  cT->SaveAs("DeltaTMeanCenter.png");
  
  outputPlots->Close();
  return eSuccess;
  
}


VModule::ResultFlag
DoubleELVESAnalysis::GlueTrace(int iTagtmp, TH1F* theTrace, int telIdtmp, int pixelIdtmp) 
{
  for(int iBin = 1; iBin < theTrace->GetNbinsX()+1; iBin++){
    double BinContent = theTrace->GetBinContent(iBin);
    fhRawPixel[telIdtmp-1][pixelIdtmp-1]->SetBinContent(iBin+int((fEventCounter)*1000./((fMaxBinTrace-fMinBinTrace)/fNBinTrace)),BinContent);  
    
  }
  return eSuccess;
}


double DoubleELVESAnalysis::CalcMHWScore(vector<double> scores)
{
    double median;
    size_t size = scores.size();
    
    sort(scores.begin(), scores.end());
    
    if (size  % 2 == 0)
      {
	median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
      }
    else
      {
	median = scores[size / 2];
      }
    
    return median;
}
  
}


// For special applications.
void AugerOfflineUser()
{
}
