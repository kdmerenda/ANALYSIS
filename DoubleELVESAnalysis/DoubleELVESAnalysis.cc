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
  hdeltaTMedian = new TH1F("deltaTMedian","Median Difference; time bin (100ns); counts", 50, 200, 700);
  hdeltaTMean = new TH1F("deltaTMean","Mean Time Difference; time bin (100ns); counts", 50, 200, 700);
  hdeltaTMeanCenter = new TH1F("deltaTMeanCenter","Mean Time Difference for Center Column; time bin (100ns); counts", 50, 00, 700);

  fEventCounter = 0;
  fNBinTrace = 300; //num bins in full trace
  fMaxBinTrace = 2999; //max num of bins in full trace
  fMinBinTrace = 0;
  fNBinPage = int(1000./((fMaxBinTrace-fMinBinTrace)/fNBinTrace)); //num bin per page
  int BinSize = int(((fMaxBinTrace-fMinBinTrace)/fNBinTrace));
  const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);
  for(int i = 0; i < fNTels; i++){
    for(int j = 0; j < fNPixels; j++){
	TString hnameSim("hSim");hnameSim+="_tel";
	hnameSim+=i+1;hnameSim+="_r";hnameSim+=detTelGlobal.GetPixel(j+1).GetRow();hnameSim+="_c";hnameSim+=detTelGlobal.GetPixel(j+1).GetColumn();
	fhRawPixel[i][j] = new TH1F(hnameSim,hnameSim,fNBinTrace,fMinBinTrace,fMaxBinTrace);
    }
  }
  return eSuccess;
}

VModule::ResultFlag 
DoubleELVESAnalysis::Run(evt::Event& event)
{

  AugerEvent& rawEvent = event.GetRawEvent();
  INFO("New Raw Event!");
  
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


    cout << eyeEvent.GetEventHeader()->GetEventNo() << " " << eyePixelList->GetNumPixels() << " " <<  eyePixelList->GetEyeNo() << endl;

    //only look at first page.
    for (int pageNum = 1; pageNum < 9; pageNum++){
      if ( eyeEvent.GetEventHeader()->GetEventNo() == fEventCounter+pageNum) return eSuccess;
    } fEventCounter = eyeEvent.GetEventHeader()->GetEventNo();

    
    TString fName("/home/kswiss/Workspace/workoffline/analysis/DoubleELVESAnalysis/traces/");
    fName+="run";
    fName+=eyeEvent.GetEventHeader()->GetRunNo();
    fName+=".event";
    fName+=eyeEvent.GetEventHeader()->GetEventNo();
    fName+=".eye";
    fName+=eyeId;
    fName+=".root";
    TFile fTMP(fName,"recreate");

    
    int doubleELVESPixelCounter = 0;
    vector<double> timeDifference;
    vector<double> timeFirstPeak;
    vector<int> columnFirstPeak;
    struct TripleTriple {
      double A;
      double B;
      int C;
      int D;
    };
    
    struct by_A { 
      bool operator()(TripleTriple const &a, TripleTriple const &b) { 
	return a.A < b.A;
      }
    };
    
    vector<TripleTriple> FirstPeak;
    int mirrorIdGlobal;
    for (unsigned int iPixel = 0; iPixel < eyePixelList->GetNumPixels(); ++iPixel) {
      
      FdUtil::Fd::FdPixelNumber pixelNumber = eyePixelList->GetPixel(iPixel);
      const unsigned int mirrorId = FdUtil::Fd::GetEyeMirrorNo(pixelNumber);
      mirrorIdGlobal = mirrorId;
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
      //  cout << mirrorId << " " << channelId << " " << pixelId << endl;
      
      
      const unsigned int startBin = fadcData->GetTraceStartBin();
      const unsigned int endBin = fadcData->GetTraceEndBin();
      
      TFADCData::FADCDataWord fadcword = fadcData->GetFADCTrace();
      
      double baseline = 0.;
      double baselineRMS = 0.;
      double charge=0.;
            
             
      // PEDESTAL CALCULATION : the first 280 bins
      // for (unsigned int pos = startBin; pos <startBin+280; ++pos) baseline += int(FADCDataWordGetData(&fadcword[pos]))/280.;
      
      // for (unsigned int pos = startBin; pos <startBin+280; ++pos) baselineRMS += pow(int(FADCDataWordGetData(&fadcword[pos])),2)/280.;

      for (unsigned int pos = startBin; pos < startBin+100; ++pos) baseline += int(FADCDataWordGetData(&fadcword[pos]))/100.;
      
      for (unsigned int pos = startBin; pos <startBin+100; ++pos) baselineRMS += pow(int(FADCDataWordGetData(&fadcword[pos])),2)/100.;
      
      baselineRMS -= pow(baseline,2);
      if (baselineRMS>0) baselineRMS=sqrt(baselineRMS);
      else baselineRMS=0.;

      TString hName("h"); hName+=fEventCounter;hName+="_m";hName+=mirrorId;hName+="_r";hName+=detTel.GetPixel(pixelId).GetRow();hName+="_c";hName+=detTel.GetPixel(pixelId).GetColumn();
      TH1F* trace = new TH1F(hName,hName,1000,startBin,endBin);

      for (unsigned int pos = startBin; pos <= endBin; ++pos) {
	charge = FADCDataWordGetData(&fadcword[pos]);
       	if (pixelId <= 440) trace->Fill(pos,charge-baseline);
      }

      // if(trace->GetMaximum() < 50) {
      //  	delete trace;
      //  	continue;
      // }

      // TH1 *hm=0;
      // TVirtualFFT::SetTransform(0);
      // hm = trace->FFT(hm, "MAG");
      // hm->SetNameTitle(hName.Append("_FFT"),"Magnitude of the 1st transform");
      // hm->Draw();
      
      // TF1 p1("p1","pol1",0,200);
      // TF1 g1("g1","gaus",200,600);
      // TF1 g2("g2","gaus",600,1000);
      // // p1.SetLineColor(kBlack);
      // g1.SetLineColor(kRed);
      // g2.SetLineColor(kBlue);
      // // trace->Fit(&p1,"R");
      // trace->Fit(&g1,"RQ");
      // trace->Fit(&g2,"RQ+");
      // // g1.GetParameters(&par[0]);
      // // g2.GetParameters(&par[3]);
      
      // //normalize      
      // if( //Difference between the mean of the gaussian (dt)
      // 	 TMath::Abs( g1.GetParameter(1) - g2.GetParameter(1)) < 150 ||
      // 	 //Amplitude of the fit at the mean
      // 	 g1.Eval(g1.GetParameter(1)) < 50 || g2.Eval(g2.GetParameter(1)) < 50 ||
      // 	 //The difference in amplitudes at the mean
      // 	 TMath::Abs(g1.GetParameter(0) - g2.GetParameter(0)) > 150 ||
      // 	 //The quality of the fits
      // 	 g2.GetChisquare()/g2.GetNDF() > 10 || g1.GetChisquare()/g1.GetNDF() > 10 ||
      // 	 //The standard deviation of the fits (width of peaks)
      // 	  g1.GetParameter(2) > 200 || g2.GetParameter(2) > 200 ||
      // 	 //significance of the peak? nope need p-value
      // 	 // g1.GetParameter(0) < 2*g1.GetParameter(2) || g2.GetParameter(0) < 2*g2.GetParameter(2) ||
      // 	 //make sure the paramters are greater than their error.
      // 	 g1.GetParameter(0) <  g1.GetParError(0)    || g1.GetParameter(1) <  g1.GetParError(1)|| g1.GetParameter(2) <  g1.GetParError(2) ||
      // 	 g2.GetParameter(0) <  g2.GetParError(0)    || g2.GetParameter(1) <  g2.GetParError(1)|| g2.GetParameter(2) <  g2.GetParError(2) 
      // 	  ){
      // 	delete trace;
      // 	continue;
      // }
      
      // doubleELVESPixelCounter++;
      // TripleTriple FirstPeakTmp;
      // FirstPeakTmp.A = g1.GetParameter(1);
      // FirstPeakTmp.C = TMath::Abs(g1.GetParameter(1) - g2.GetParameter(1));
      // FirstPeakTmp.B = detTel.GetPixel(pixelId).GetColumn();
      // timeFirstPeak.push_back(g1.GetParameter(1));
      // columnFirstPeak.push_back(detTel.GetPixel(pixelId).GetColumn());
      // timeDifference.push_back(TMath::Abs(g1.GetParameter(1) - g2.GetParameter(1)));
      // FirstPeak.push_back(FirstPeakTmp);
      //||
      //	 p1.GetChisquare()/p1.GetNDF() > 2 )
      // Double_t par[6];
      // TF1 *g1    = new TF1("g1","gaus",200,600);
      // TF1 *g2    = new TF1("g2","gaus",600,1000);
      // TF1 *total = new TF1("total","gaus(0)+gaus(3)",200,1000);
      // total->SetLineColor(2);
      // trace->Fit(g1,"R");
      // trace->Fit(g2,"R+");
      // g1->GetParameters(&par[0]);
      // g2->GetParameters(&par[3]);
      // total->SetParameters(par);
      // trace->Fit(total,"R+");
      //delete g1, g2, total;

      Int_t nbins = 1000;
      TSpectrum *spectrum = new TSpectrum(100);
      Float_t * source = new Float_t[nbins];
      Float_t * dest   = new Float_t[nbins];
      TH1F *d = new TH1F(hName.Append("d"),"",nbins,0,999);
      for (int i = 0; i < nbins; i++) source[i]=trace->GetBinContent(i + 1);
      spectrum->SetName(hName.Append("_spectrum"));
      //      Int_t nfound = spectrum->Search(trace,50,"",0.10);
      Int_t nfound = spectrum->SearchHighRes(source, dest, nbins, 50, 20, kTRUE, 2, kTRUE, 3);
      printf("Found %d candidate peaks to fit\n",nfound);
      TH1 *hb = spectrum->Background(trace,20,"same");
      hb->SetName(hName.Append("_background"));
      Float_t *xpeaks = spectrum->GetPositionX();
      Double_t fPositionX[100];
      Double_t fPositionY[100];
      for (int i = 0; i < nfound; i++) {
	int a=xpeaks[i];
	int bin = 1 + Int_t(a + 0.5);
	fPositionX[i] = trace->GetBinCenter(bin);
	fPositionY[i] = trace->GetBinContent(bin);
      }
      TPolyMarker * pm = (TPolyMarker*)trace->GetListOfFunctions()->FindObject("TPolyMarker");
      if (pm) {
	trace->GetListOfFunctions()->Remove(pm);
	delete pm;
      }
      pm = new TPolyMarker(nfound, fPositionX, fPositionY);
      trace->GetListOfFunctions()->Add(pm);
      pm->SetMarkerStyle(23);
      pm->SetMarkerColor(kRed);
      pm->SetMarkerSize(1.3);
      for (int i = 0; i < nbins; i++) d->SetBinContent(i + 1,dest[i]);
      d->SetLineColor(kRed);

      Double_t par[6];
      if(nfound = 2  &&  fPositionY[0] > 20 &&  fPositionY[1] > 20){
	cout <<"Time Difference: " << fPositionX[1]-fPositionX[0] << endl;
	par[0]=fPositionY[0];par[1]=fPositionX[0];par[2]=100;
	par[3]=fPositionY[1];par[4]=fPositionX[1];par[5]=100;
	TF1 total("total","gaus(0)+gaus(3)",200,1000);
	total.SetParameters(par);
	trace->Fit(&total,"RQ");
      }else{
	delete trace;
	delete d;
	delete hb;
	continue;         
      }


      doubleELVESPixelCounter++;
      TripleTriple FirstPeakTmp;
      FirstPeakTmp.A = par[1];
      FirstPeakTmp.B =TMath::Abs(par[4]-par[1]);
      FirstPeakTmp.C = detTel.GetPixel(pixelId).GetColumn();
      FirstPeakTmp.D = detTel.GetPixel(pixelId).GetRow();
      timeFirstPeak.push_back(fPositionX[0]);
      columnFirstPeak.push_back(detTel.GetPixel(pixelId).GetColumn());
      timeDifference.push_back(fPositionX[1]-fPositionX[0]);
      FirstPeak.push_back(FirstPeakTmp);


    }//end pixel loop

    sort(FirstPeak.begin(),FirstPeak.end(),by_A());
    double averageDeltaT = accumulate(timeDifference.begin(), timeDifference.end(), 0.0)/timeDifference.size();
    //A: mean g1
    //C: Abs (mean g1 - mean g2)
    //B: Pixel Num Column
    double averageMinColT=0;
    int averageCounter=0;
    TString heatmapDeltaTName("heatmapDeltaT_m");heatmapDeltaTName+=mirrorIdGlobal;
    TH2F heatmapDeltaT(heatmapDeltaTName,"heatmapDeltaT",20,1,20,22,1,22);
    if (FirstPeak.size() != 0){
      for(int i=0; i < FirstPeak.size(); i++){
	cout << FirstPeak[i].A  << endl;
 	cout << FirstPeak[i].B  << endl;
	cout << FirstPeak[i].C  << endl;
	cout << FirstPeak[i].D  << endl;
	cout << averageDeltaT << endl << endl;
	heatmapDeltaT.Fill(FirstPeak[i].C,FirstPeak[i].D,FirstPeak[i].B);
	//calculate average based on center column
	if(FirstPeak[i].C == FirstPeak[0].C){
	  averageMinColT += FirstPeak[i].B;
	  averageCounter++;
	}
      }
      averageMinColT = averageMinColT / averageCounter;
    }

    
    //making sure a full event is looked at. 
    if(doubleELVESPixelCounter >= 5){
      hdeltaTMeanCenter->Fill(averageMinColT);
      hdeltaTMean->Fill(averageDeltaT);
      hdeltaTMedian->Fill(CalcMHWScore(timeDifference));
      fTMP.Write();
      fTMP.Close();
    }else{ fTMP.Close();   fTMP.Delete();}
    
  }
  
  return eSuccess;
}

VModule::ResultFlag 
DoubleELVESAnalysis::Finish() 
{
  TCanvas* cT = new TCanvas("cT","cT",800,600);
  hdeltaTMedian->Draw();
  cT->SaveAs("DeltaTMedian.png");
  hdeltaTMean->Draw();
  cT->SaveAs("DeltaTMean.png");
  hdeltaTMeanCenter->Draw();
  cT->SaveAs("DeltaTMeanCenter.png");
  
  return eSuccess;
  
}


VModule::ResultFlag
DoubleELVESAnalysis::GlueTrace(int iTagtmp, TH1F* theTrace, int telIdtmp, int pixelIdtmp) 
{
  double PreviousBinContent = 0;

  for(int iBin = 1; iBin <= theTrace->GetNbinsX(); iBin++){
    double BinContent = theTrace->GetBinContent(iBin);
    double normalizedBinContent;
    
    //if (fEventCounter == 1 && iBin > theTrace->GetNbinsX() - 10 ) continue;
    //    if( BinContent < PreviousBinContent - 1000) continue;
    fhRawPixel[telIdtmp-1][pixelIdtmp-1]->SetBinContent(iBin+int((fEventCounter-1)*1000./((fMaxBinTrace-fMinBinTrace)/fNBinTrace)),BinContent);
    PreviousBinContent = BinContent;
    
    
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
