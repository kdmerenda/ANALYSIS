#include "NetProd.h"
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
#include <fevt/EyeHeader.h>
#include <fdet/Telescope.h>
#include <fdet/Camera.h>
#include <fdet/Pixel.h>
#include <fdet/Channel.h>
#include <fdet/Mirror.h>
#include <fdet/Filter.h>
#include <fdet/Corrector.h>

#include <EyeEventFile.hh>

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
#include <TCanvas.h>
#include <TCutG.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <vector>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

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


namespace NetProd {


NetProd::NetProd() 
{
}

NetProd::~NetProd()
{
}

VModule::ResultFlag 
NetProd::Init()
{
  exec("rm -f /media/kswiss/ExtraDrive1/ELVESNet/merged/TruthMerged.list");

  fData=false;
  //how many pages will be glues together and outputed?
  fNPages = 3;
  //Time binning in microseconds
  fBinning = 2;
  /*
    1000 time bins = 0.1 us \
    100  time bins = 1   us  } us*10*timebins= 1000
    50   time bins = 2   us /
   */
  //trackers of the number of events belonging to the same one. 
  fOutFileName = new TString("");
  fEventCounter = 0;
  fEventTracker = 0;
  fDataEventCounter = 0;
  //automatic binning based on num pages
  fNBinTrace = 1000./(10*fBinning)*fNPages; //num bins in full trace --> 10th of detector
  fMaxBinTrace = 1000*fNPages; //upper bin value in full trace
  fMinBinTrace = 1;//lower bin value in full trace
  fNBinPage = 1000./(fBinning*10) ; //num bin per page

  //merging all mirrors
  for(int j = 0; j < fNTels*fNPixels; j++){
    TString hname("Merge");hname+="_";hname+=j+1;
    fhRawPixelMerge[j] = new TH1F(hname,hname,fNBinTrace,fMinBinTrace,fMaxBinTrace);
  }
  
  ostringstream info;
  info << "\nOuput format of data: \n"
    "        Number of Pages = " << fNPages << "\n"
    "Number of Bins in Trace = " << fNBinTrace << "\n"
    " Number of Bins in Page = " << fNBinPage << "\n"
    "          Bin Size (us) = " << fBinning << "\n";
  INFO(info);

  TFile *fTMP = new TFile("traces/traces.root","recreate");
  fTMP->Close();
  
  return eSuccess;
}

VModule::ResultFlag 
NetProd::Run(evt::Event& event)
{

  AugerEvent& rawEvent = event.GetRawEvent();
  const evt::Header& header = event.GetHeader();


  double EventTimeCurrent = header.GetTime().GetGPSSecond() + 1e-9*header.GetTime().GetGPSNanoSecond();
  bool toadd;

    
  
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


    if(eyeEvent.GetEventHeader()->GetEventNo() != fEventTracker+1) {
      fEventCounter = 0;
    }

    ostringstream info;
    std::string fileName;
    if(fEventCounter == 0){
      TString stmp("");stmp+= header.GetId();
      stmp.Remove(0,5);//remove the eye from the ID to make hname
      fFirstEventId = stmp;
      info << "NEW EVENT!!! Starting work on page #" << fEventCounter+1 << " of " << fFirstEventId << endl;
      fDataEventCounter++; //count the amount of full events we are outputting. 

      //get filename out of file list
      ifstream fileList;
      fileList.open("sloppy.list");
      int fileListCounter = 0;
      while (std::getline(fileList,fileName)){
	fileListCounter++;
	if(fileListCounter == fDataEventCounter) {
	  cout << fileName << endl;
	  break;
	}
      }
      fOutFileName->Clear();
      fOutFileName->Append(fileName.data());
      fileList.close();
      for(int i = 0; i < fNTels*fNPixels; i++){
	for(int j = 1; j <= fhRawPixelMerge[i]->GetNbinsX()+1; j++){
	  fhRawPixelMerge[i]->SetBinContent(j,0);
	}
      }
      

    }else{
      info << "RUNNING >>>>>>>> Page #" << fEventCounter+1 << endl;
    }
    INFO(info);
    
    fEventTracker = eyeEvent.GetEventHeader()->GetEventNo();
    fEventCounter++;

    for (unsigned int iPixel = 0; iPixel < eyePixelList->GetNumPixels(); ++iPixel) {

      FdUtil::Fd::FdPixelNumber pixelNumber = eyePixelList->GetPixel(iPixel);
      const unsigned int mirrorId = FdUtil::Fd::GetEyeMirrorNo(pixelNumber);
      const unsigned int channelId = FdUtil::Fd::GetMirrorPixelNo(pixelNumber);
  
      const TFADCData* fadcData = eyeFADCData->GetFADCData(pixelNumber);

      fevt::Telescope& tel = eye.GetTelescope(mirrorId);
      const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
      const fdet::Eye& detEye = detFD.GetEye(eyeId);
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

      TString hName("h"); hName+=fEventTracker;hName+="_";hName+=fEventCounter;hName+="_m";hName+=mirrorId;hName+="_r";hName+=detTel.GetPixel(pixelId).GetRow();hName+="_c";hName+=detTel.GetPixel(pixelId).GetColumn();
      TH1F* trace = new TH1F(hName,hName,1000,startBin,endBin);

      //need to reindex the pixel id. 
      int globalPixelId = pixelId+(fNTels-mirrorId)*440;
      
      if(pixelId > 440) continue;
      for (unsigned int pos = startBin; pos <= endBin; ++pos) {
	charge = FADCDataWordGetData(&fadcword[pos]);
 	if (abs(baseline-baselinetail) > 10 && baseline >baselinetail){
	  trace->SetBinContent(pos,charge-baselinetail);	  
	}else{
	  trace->SetBinContent(pos,charge-baseline);	  
	}
      }
      trace->Rebin(int(fBinning*10));
  
      //gluing traces for the merged event. 
      for(int iBin = 1; iBin <= trace->GetNbinsX(); iBin++){
	fhRawPixelMerge[globalPixelId-1]->SetBinContent(iBin+int((fEventCounter-1)*trace->GetNbinsX()),trace->GetBinContent(iBin));
      }
      delete trace;
 
    }//pixel loop
  }//eye loop 

  //once we looped through three pages save them
  if(fEventCounter != fNPages) return eSuccess;

  //obtain the truth location of this event. 
  TString com2exe(fOutFileName->Data());com2exe.Remove(com2exe.Length()-9,9);
  com2exe.Prepend("grep -i '");com2exe.Append("' /media/kswiss/ExtraDrive1/ELVESNet/raw/Truth.list");  
  cout << exec(com2exe.Data()) << endl;
  std::istringstream iss(exec(com2exe.Data()));
  float TruthLat, TruthLon;
  string tmpdump;
  iss >> tmpdump >> TruthLat >> TruthLon;
  cout << TruthLat << " "  << TruthLon << endl;

  // TFile *fTMP = new TFile("traces/traces.root","update");
  // for(int j = 0; j < fNTels*fNPixels; j++)
  //   if(fhRawPixelMerge[j]->GetEntries()!=0)fhRawPixelMerge[j]->Write();
  // fTMP->Close();

  //to output the 6 telescopes, change the foolliwng variable to fNTels
  int wantedfNTels = 3;//fNTels;
  TH2F theELVES("theELVES","",wantedfNTels*fNColumns,1,wantedfNTels*fNColumns,fNRows,1,fNRows);
  const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);
  
  for(int i = 0; i < wantedfNTels; i++){
    for(int j = 0; j < fNPixels; j++){
      int pixelId = j+1;
      int COL = detTelGlobal.GetPixel(pixelId).GetColumn() + i*20;
      int ROW = detTelGlobal.GetPixel(pixelId).GetRow();
      theELVES.Fill(COL,ROW,fhRawPixelMerge[j+i*fNPixels]->Integral());
    }
  }

  TCanvas theC("theC","theC",1000,500);
  gStyle->SetOptStat(0);
  theELVES.Draw("col");
  fOutFileName->Replace(35,3,"merged");
  fOutFileName->Remove(fOutFileName->Length()-5,5);fOutFileName->Append(".png");
  theC.SaveAs(fOutFileName->Data());


  //dump the histogram values in text file
  fstream outFile;
  fOutFileName->Remove(fOutFileName->Length()-4,4);fOutFileName->Append(".dat");
  cout << fOutFileName->Data() << endl;
  outFile.open(fOutFileName->Data(),std::ios::out);
  outFile << wantedfNTels*fNPixels << endl;
  outFile << fNBinTrace << endl;
  for(int i = 0; i < wantedfNTels*fNPixels; i++){
    for(int j = 1; j <= fhRawPixelMerge[i]->GetNbinsX(); j++){
      outFile << fhRawPixelMerge[i]->GetBinContent(j) << endl;
    }
  }
  outFile.close();

  fstream outFile2;
  outFile2.open("/media/kswiss/ExtraDrive1/ELVESNet/merged/TruthMerged.list",std::ios::out | std::ios::app);
  outFile2 << fOutFileName->Data() << " " << TruthLat << " " << TruthLon << endl;
  outFile2.close();
  return eSuccess;
}

  VModule::ResultFlag 
NetProd::Finish() 
{
  return eSuccess;
}
  
std::string NetProd::exec(const char* cmd) {
  std::array<char, 128> buffer;
  std::string result;
  std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
  if (!pipe) throw std::runtime_error("popen() failed!");
  while (!feof(pipe.get())) {
    if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
      result += buffer.data();
  }
    return result;
}

  
}



// For special applications.
void AugerOfflineUser()
{
}
