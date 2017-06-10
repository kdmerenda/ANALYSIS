/**
   \file
   implementation of FdElvesFinder

   \author:  Roberto Mussa
   \version $Id: FdProcessElves.cc  $
   \date 16 Jan 2017
 
   Modified version of FdSDPFinder to reconstruct elves
 
 - v0 : basic
 - v1 : extended readout 
 - v1.1: corrected trace gluing
 - 2017-02-01: BoltFitter debugging:
     - site coordinates and backwall info passed as parameters to t0linfitter
 - 2017-02-03: PlotTimeFit5D_Output added, to plot the results of the fit to the bolt direction
 -             Fit3D because we fit only Time0, Lat, Long; 
 -             it will soon become Fit5D to account also for Hbolt (lightning altitude) 
 -                 and Hd (light emission altitude) : but we need to well understand the contributions to the chisquare
 - 2017-02-07: FdSingleElvesTime_CameraViewPlot added : we want to see things in the two perspectives: FDview at basic level
 - 2017-02-12: FdSingleElvesPull_CameraViewPlot added : Probably it is better to keep dummy pixels in a separate dethod
 - 2017-02-14: FdSingleElvesPixel_CameraViewPlot added:  just dummy pixels
 - 2017-03-
 - 2017-03-
 - 2017-03-
 - 2017-03-
 - 2017-04-
 - 2017-04-
 
*/

//static const char CVSId[] =
//"$Id: FdSDPFinder.cc 20863 2012-04-09 20:16:57Z paul $";

#include <iostream>

#include <evt/Event.h>


#include <fevt/FEvent.h>
#include <fevt/Eye.h>
#include <fevt/EyeHeader.h>
#include <fevt/EyeRecData.h>
#include <fevt/Pixel.h>
#include <fevt/Telescope.h>
#include <fevt/PixelRecData.h>
#include <fevt/TelescopeTriggerData.h>

#include <fwk/CentralConfig.h>
#include <fwk/SVNGlobalRevision.h>
#include <fwk/LocalCoordinateSystem.h>

#include <utl/Point.h>
#include <utl/AxialVector.h>
#include <utl/ErrorLogger.h>
#include <utl/MathConstants.h>
#include <utl/Reader.h>

#include <det/Detector.h>

#include <fdet/FDetector.h>
#include <fdet/Eye.h>
#include <fdet/Pixel.h>
#include <fdet/Channel.h>
#include <fdet/Telescope.h>
#include <fdet/Camera.h>

#include <TMinuit.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TImage.h>
#include <TLegend.h>

#include <TH2F.h>
#include <TH1.h>
#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPaletteAxis.h>
#include <TMath.h>
#include <TF1.h>
#include <TVectorF.h>
#include <TVector3.h>
#include <TFitResult.h>
#include <TProfile.h>

#include <AugerEvent.h>
#include <EyeEvent.hh>
#include <EyePixelList.hh>
#include <FadcData.hh>

#include <atm/CloudResult.h>


#include "FdProcessElves.h"

using namespace std;
//using namespace FdElvesFinder;
using namespace evt;
using namespace fevt;
using namespace fwk;
using namespace utl;
using namespace det;
using namespace fdet;
using namespace atm;


VModule::ResultFlag FdProcessElves::Init(){
    

   // fstream fIn("XR_e1_2014.d", ios::in);
    fstream fIn("testfile.d", ios::in);
    
    int nl=0;
    while (!fIn.eof()) {
        if (NMAXPAGES==9)fIn >>  GPSsec[nl] >> GPSns[nl] >> evno[nl] >> nwd[nl] >> bayno[nl] >> dts[0][nl] >> dts[1][nl] >> dts[2][nl] >> dts[3][nl] >> dts[4][nl] >> dts[5][nl] >> dts[6][nl] >> dts[7][nl] ;
        if (NMAXPAGES==3)fIn >>  GPSsec[nl] >> GPSns[nl] >> evno[nl] >> nwd[nl] >> bayno[nl] >> dts[0][nl] >> dts[1][nl] ;

        nl++;
    }
    NpsPages = nl;
    
    cout << "PRESCAN file read " << NpsPages-1 << " PAGES  "<< endl;
    
    fIn.close();
    
    
    float BW[4]  = {-30.00,60.03,-171.85,-116.68}; 	// with respect to the Est, in degrees
    float H[4] = {1.45,1.45,1.45,1.7}; // altitude asl
    
    float LG[4] = {-69.449673,-69.012203,-69.210845,-69.599633};
    float LT[4] = {-35.495759,-35.291974,-34.935916,-35.114138};
    int CS[4] = {4,6,42,3};
    
    TLEAltitude=85.;
    

    TCanvas *ctest = new TCanvas("ctest"," test ", 1400, 1000);
    ctest -> Divide(2,2);
    
    TH2F *testLL  = new TH2F("testLL","test; Long; Lat ",100,-78.244,-54.756,100,-45,-25);
    
    
    for (int i=0; i<4; i++){
        backwall[i]= BW[i];
        hsite[i] = H[i];
        longsite[i] =LG[i];
        latsite[i]=LT[i];
        sitecolor[i]=CS[i];
        
        
        
        
        Site1Location.SetMag(Rearth+hsite[i]);
        Site1Location.SetPhi(longsite[i]*degree);
        Site1Location.SetTheta((90.-latsite[i])*degree);
        
        East1 = PoloNord.Cross(Site1Location);
        
        East1.SetMag(1.);
        

        
        ctest->cd(i+1);
        testLL->Draw();
        
        for (int j=1; j<6; j+=2){
            
            for (int itc = 0 ; itc<2;itc++){
                for (int ipp=0; ipp<NPIXELS; ipp++){
                    
                    int itp = itc*NPIXELS+ipp;
                    
                    
                    double GeomCorr = CalculateGeomCorr(itp,i,j) ;
                    PixelEdgesLL[itp]->Draw("L");
                    
                    
                }
            }
        }

        
        
        
    }
    
    
    cout << " Saving TEST.gif"<< endl;
    ctest->SaveAs("TEST.gif");

    
    
    
    gROOT->SetStyle("Plain");
    int         palette[MaxColors];
    
    Int_t FI =    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, MaxColors);
    for (int i=0;i<MaxColors;i++) palette[i] = FI+i;
    
    TColor::SetPalette(MaxColors,palette);
    gStyle->SetNumberContours(MaxColors);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);

    Branch topB =
     CentralConfig::GetInstance()->GetTopBranch("FdProcessElves");

  if (!topB) {
    ERROR("Could not find branch FdProcessElves");
    return eFailure;
  }

  topB.GetChild("minPixels").GetData(fMinPixels);
  topB.GetChild("maxPixels").GetData(fMaxPixels);

  // info output
  ostringstream info;
  info << " Version: "
       << GetVersionInfo(VModule::eRevisionNumber) << "\n"
          " Parameters:\n"
          "       min no. pixels: " << fMinPixels << "\n"
          "       max no. pixels: " << fMaxPixels << "\n";

  INFO(info);
    //fOut = new TFile("Output.root","recreate");
    
//    ResetAllVectors();
    fPixels = new TObjArray;
    fPixels2 = new TObjArray;
    fPixels3 = new TObjArray;
    fPixels4 = new TObjArray;
    fPixels5 = new TObjArray;
    fPixels6 = new TObjArray;
    fPixels7 = new TObjArray;
    fPixels8 = new TObjArray;

    fSites = new TObjArray;
    fSites2 = new TObjArray;

    for (int np = 0 ; np < NMAXPAGES ; np++ ){
        for (int itt=0; itt<6;itt++) TelHit[np][itt]=false;
        
        for (int npp = 0; npp < NPIXELS; npp++)pixHasData[np][npp]=false;
    }
    for (int itt=0; itt<3; itt++) HitTel[itt]=0;

    for (int itc = 0 ; itc<3;itc++){
      for (int ipp=0; ipp< NPIXELS; ipp++){
          int ir=ipp%NROWS;
          int row = ipp%NROWS+1;
          int col = int(ipp/NROWS)+1;

	stringstream ss;
	ss << "QvsTime" << itc <<  "_R" << row << "_C" << col ;
    stringstream ssc;
    ssc << "QCvsTime" << itc << "_R" << row << "_C" << col ;
          stringstream sscc;
          sscc << "QCCvsTime" << itc << "_R" << row << "_C" << col ;
          stringstream ssccc;
          ssccc << "QCCCvsTime" << itc << "_R" << row << "_C" << col ;
    stringstream ssc2;
    ssc2 << "dQCvsTime" << itc <<  "_R" << row << "_C" << col ;

	//cout << " Initializing histo " << ss.str().c_str() << endl;                                                    
	// the histogram has 3*NBINS; currently we integrate on 2 microsec, 20 bins                                      
	// we allocate a max of 3 TELESCOPES for the elve
    // CAREFUL: NMAXPAGES = 1 in 2013, 3 in 2014-6, ? in 2017
        //  int nslots = NMAXPAGES*size/nSlotInterval;
	      hRawPixels[itc][ipp]=new TH1F(ss.str().c_str(),"Raw Signal vs Time, in one pixel",nslots,0.,nslots);

          hCalPixels[itc][ipp]=new TH1F(ssc.str().c_str(),"Cal Signal vs Time, in one pixel",nslots,0.,nslots);
          hGCorrPixels[itc][ipp]=new TH1F(sscc.str().c_str(),"Signal Geometry Corrected  vs Time, in one pixel",nslots,0.,nslots);
          hTCorrPixels[itc][ipp]=new TH1F(ssccc.str().c_str(),"Signal Corrected for Geometry and Atmosphere vs Time, in one pixel",nslots,0.,nslots);

          hCalResiduals[itc][ipp]=new TH1F(ssc2.str().c_str(),"Cal Signal vs Time, in one pixel",nslots,0.,nslots);

          stringstream ss2;
          ss2 << "PeakTime_Tel"<< itc << "_Row_" << ir ;
          stringstream ss3;
          ss3 << "Peak Time vs Column, for Row "<< ir ;
          
          if(ipp<NROWS)hPixelPeakTimes[itc][ir]=new TH1D(ss2.str().c_str(),ss3.str().c_str(),NCOLS,0.,NCOLS);
         
      }}
    int ncx=60;
    for (int ir=0; ir<NROWS; ir++){
        stringstream ss2;
        ss2 << "TimeMax_vs_Col"<< ir+1  ;
        stringstream ss3;
        ss3 << "Peak Time vs Column, for Row "<< ir+1 ;
        
        TimeMax_vs_Col[ir]=new TH1F(ss2.str().c_str(),ss3.str().c_str(),ncx,0.5,ncx+0.5);
        
        
        stringstream ss2a;
        ss2a << "NphMax_vs_Col"<< ir+1  ;
        stringstream ss3a;
        ss3a << "NphMax vs Column, for Row "<< ir+1 ;
        
        NphMax_vs_Col[ir]=new TH1F(ss2a.str().c_str(),ss3a.str().c_str(),ncx,0.5,ncx+0.5);
        
        stringstream ss2b;
        ss2b << "SigmaT_vs_Col"<< ir+1  ;
        stringstream ss3b;
        ss3b << "Sigma Time Raising Edge vs Column, for Row "<< ir+1 ;
        
        
        SigmaT_vs_Col[ir]=new TH1F(ss2b.str().c_str(),ss3b.str().c_str(),ncx,0.5,ncx+0.5);
        
        
        stringstream ss2c;
        ss2c << "Skew_vs_Col"<< ir+1  ;
        stringstream ss3c;
        ss3c << "Skewness of the Peak vs Column, for Row "<< ir+1 ;
        
        Skew_vs_Col[ir]=new TH1F(ss2c.str().c_str(),ss3c.str().c_str(),ncx,0.5,ncx+0.5);
        

    }

    
    hRawTmax = new TH1F("hRawTmax"," Highest Pulse Peak Time, Raw; iCOL*22+iROW ",NPIXELS*3,0,NPIXELS*3);
    hRawQmax = new TH1F("hRawQmax"," Highest Pulse Peak Nphotons, Raw; iCOL*22+iROW ",NPIXELS*3,0,NPIXELS*3);
    hRawTmax2= new TH1F("hRawTmax2"," 2nd highest Pulse Peak Time, Raw ",NPIXELS*3,0,NPIXELS*3);
    hRawQmax2= new TH1F("hRawQmax2"," 2nd highest Pulse Peak Nphotons, Raw ",NPIXELS*3,0,NPIXELS*3);

    hFineTmax= new TH1F("hFineTmax"," Highest Pulse Peak Time, Fine ; R+22*C; Time ",NPIXELS*3,0,NPIXELS*3);
    hFineQmax= new TH1F("hFineQmax"," Highest Pulse Peak Nphotons, Fine ; R+22*C; Nph",NPIXELS*3,0,NPIXELS*3);
    hFineQtot= new TH1F("hFineQtot"," Highest Pulse Tot Nphotons, Fine ; R+22*C; Nph",NPIXELS*3,0,NPIXELS*3);

    hFineQCmax= new TH1F("hFineQCmax"," Highest Pulse Peak Nphotons, Geometry corrected ; R+22*C; Nph/kmq",NPIXELS*3,0,NPIXELS*3);
    hFineQCCmax= new TH1F("hFineQCCmax"," Highest Pulse Peak Nphotons, Geom+AtmoCorr ; R+22*C; Nph/kmq",NPIXELS*3,0,NPIXELS*3);
    hFineQCtot= new TH1F("hFineQCtot"," Highest Pulse Tot Nphotons, Geometry corrected ; R+22*C; Nph",NPIXELS*3,0,NPIXELS*3);
    hFineQCCtot= new TH1F("hFineQCCtot"," Highest Pulse Tot Nphotons, Geom+AtmoCorr  ; R+22*C; Nph",NPIXELS*3,0,NPIXELS*3);
    
    
    hFineTmax2= new TH1F("hFineTmax2"," Highest 2nd xPulse Peak Time, Fine ; R+22*C; Time ",NPIXELS*3,0,NPIXELS*3);
    hFineQmax2= new TH1F("hFineQmax2"," Highest Pulse Peak Nphotons, Fine ; R+22*C; Nph",NPIXELS*3,0,NPIXELS*3);
    hFineSigmaT= new TH1F("hFineSigmaT"," Highest Pulse Left Sigma, Fine ; R+22*C; #sigma_{T}",NPIXELS*3,0,NPIXELS*3);
    hFineSkew= new TH1F("hFineSkew"," Highest Pulse Skewness, Fine ; R+22*C; Skew",NPIXELS*3,0,NPIXELS*3);

    hFIT3Tmax= new TH1F("hFIT3Tmax"," Highest Pulse Peak Time, Fine; R+22*C; Time ",NPIXELS*3,0,NPIXELS*3);
    hFIT3Qmax= new TH1F("hFIT3Qmax"," Highest Pulse Peak Photons, Fine; R+22*C; Nph ",NPIXELS*3,0,NPIXELS*3);
    hFIT3SigmaT= new TH1F("hFIT3SigmaT"," Highest Pulse Peak Left Sigma , Fine; R+22*C; delta Time ",NPIXELS*3,0,NPIXELS*3);
    hFIT3Skew= new TH1F("hFIT3Skew"," Highest Pulse Skewness , Fine; R+22*C; skewness ",NPIXELS*3,0,NPIXELS*3);

    hFitPullsN= new TH1F("hFitPullsN"," Negative Pull in Time Fit; R+22*C; Pull ",NPIXELS*3,0,NPIXELS*3);
    hFitPullsP= new TH1F("hFitPullsP"," Positive Pull in Time Fit; R+22*C; Pull ",NPIXELS*3,0,NPIXELS*3);

    double MinQC=-0.015e15;
    double MaxQC= 1.00e15;

    hQC_vs_Dist[0]= new TH1F("QC_vs_Dist0","Corrected Nph/Area vs distance from bolt; ArcR(km); dNph/dArea",10,0., 500.);
    hQC_vs_Dist[1]= new TH1F("QC_vs_Dist1","Corrected Nph/Area vs distance from bolt; ArcR(km); dNph/dArea",10,0., 500.);
    hQC_vs_Dist[2]= new TH1F("QC_vs_Dist2","Corrected Nph/Area vs distance from bolt; ArcR(km); dNph/dArea",10,0., 500.);
    
    hQC_vs_SinThBE[0]= new TH1F("QC_vs_SinThBE0","Corrected Nph/Area vs sin(theta,dipole) from bolt; sin#theta_{dipole}; dNph/dArea",10,0., 1.7);
    hQC_vs_SinThBE[1]= new TH1F("QC_vs_SinThBE1","Corrected Nph/Area vs sin(theta,dipole) from bolt; sin#theta_{dipole}; dNph/dArea",10,0., 1.7);
    hQC_vs_SinThBE[2]= new TH1F("QC_vs_SinThBE2","Corrected Nph/Area vs sin(theta,dipole) from bolt; sin#theta_{dipole}; dNph/dArea",10,0., 1.7);
    
    for (int i=0; i<NDP1*NDP2; i++){
        int ii=i%NDP1;
        int jj=i/NDP1;
        double delta = (ii+1)*DDELTA_KM; double phi0 = jj*(360./NDP2);
        stringstream ss3;ss3 << "pQCC_vs_ElliDist" << i << endl;
        stringstream ss4;ss4 << "hQCC_Spread" << i << endl;

        pQCC_vs_ElliDist[i] = new TProfile(ss3.str().c_str(),ss3.str().c_str(),NDBbins,0., 400.,0.,3.e15);
        hQCC_Spread[i] = new TH1F(ss4.str().c_str(),ss4.str().c_str(),40,0.,2.);

        
    }

    
    IsThisFirstEvent=true;
  return eSuccess;
}


/** Ignore eye if less than fMinPixels (datacard file\code <minPixels>\endcode)
 */
VModule::ResultFlag FdProcessElves::Run(evt::Event& event){
    
    if (!event.HasFEvent()) return eSuccess;
    // the pedestal array has to be reset only at the beginning of the elve
    //for (int imir=0; imir< 6; imir++){
    //	for (int ipixx =0; ipixx<440; ipixx ++){ rawPed[imir][ipixx]=0.;}}
    
    AugerEvent& rawEvent = event.GetRawEvent();
    
    for (AugerEvent::EyeIterator eyeIter = rawEvent.EyesBegin();
         eyeIter != rawEvent.EyesEnd(); ++eyeIter) {
        
        TEyeEvent& eyeEvent = *eyeIter;
        //TEyeEventHeader* eyeHeader =  eyeEvent.GetEventHeader();
        
        TEyePixelList* eyePixelList = eyeEvent.GetPixelList();
        
        const unsigned int numpixels = eyePixelList->GetNumPixels();
        
        cout << "NUM PIXELS " << numpixels << endl;
        
        const int eyeId = eyeEvent.GetEventHeader()->GetEyeNo();
        const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
        
        
        FEvent& fdEvent = event.GetFEvent();
        fevt::Eye& eye = fdEvent.GetEye(eyeId);

        site1 = eyeId-1;
        Site1Location.SetMag(Rearth+hsite[site1]);
        Site1Location.SetPhi(longsite[site1]*degree);
        Site1Location.SetTheta((90.-latsite[site1])*degree);
        
        East1 = PoloNord.Cross(Site1Location);
        
        East1.SetMag(1.);
        
        

        
        TimeStamp	timeMy = eye.GetHeader().GetTimeStamp();
        
        cout << "GPS " << timeMy <<  " NUM PIXELS " << numpixels << " EVTYPE  ";
        cout << eye.GetHeader().GetEventType() << " EVCLASS " << eye.GetHeader().GetEventClass() << endl;
      
        fTelescopeCount = 1;
        CalculateDelaysBetweenBays(eyeEvent,eye);
        // here we check the timing: for now on we just check the GPSsecond :
        // as we know that in some cases we have 2 elves in the same second,
        // we will need to at least check also the msec
        // IF this is NOT an elve follower , we have to reset all histos and variables
        // ||( timeMy.GetGPSSecond() == ThisGPSsec && deltaNanoTime>1000000))
       /* if ( NewEventFound(eye))  {
            
            if (!IsThisFirstEvent){
                    cout << "Writing Output File ... " << endl;
                    ProcessHistos();
                    WriteOutputFile(ThisGPSsec,EventNanoTime);
                
                } else IsThisFirstEvent = false;
                ResetAllForNewEvent(timeMy.GetGPSSecond(),timeMy.GetGPSNanoSecond());
            }
            else nPages++;
        
            ProcessRawData(eyeEvent,eye);
        */
	//if(fTelescopeCount == 1)DrawAmplitudeColor_LongLatMap(eyeEvent,eye);
        
        
        
        if (!ProcessingEvent){

            
            
            if (fTelescopeCount == NbaysEventFoundinPreScan(eye)){
                 ProcessingEvent=true;
                 fL7TelescopeCount = fTelescopeCount;
                 nPages=0;
            }
            
        }
        
        if (ProcessingEvent){
            nPages++;
            cout << ThisGPSsec << " Processing Raw data ... " << nPages << " " << EffNpages << endl;
            ProcessRawData(eyeEvent,eye);
            
            
            if (nPages == EffNpages){
                cout << "Writing Output File ... " << endl;
                ProcessHistos();
                WriteOutputFile(ThisGPSsec,EventNanoTime);
                ResetAllForNewEvent();
                
            }
        }
       
    }
    
    return eSuccess;
}

VModule::ResultFlag FdProcessElves::Finish(){
    
    LastReset();
    
  return eSuccess;
}



/*
 * \brief Find an inital guess of the SDP.
 *
 * Find the cross products of the pixels with a pulse
 * Pick the one that accomodates most trackbed pixels
 *
 */
bool FdProcessElves::NewEventFound(fevt::Eye& eye){
    TimeStamp	timeMy = eye.GetHeader().GetTimeStamp();
    bool isnewevent=false;
    int deltaNanoTime = (int) timeMy.GetGPSNanoSecond() - (int) EventNanoTime;
    if (timeMy.GetGPSSecond() != ThisGPSsec) deltaNanoTime = 999999999;
    if (deltaNanoTime>1000000){
        isnewevent = true;
        cout << timeMy.GetGPSSecond() << " NewEventFound " << (int) timeMy.GetGPSNanoSecond()
        << " - " << EventNanoTime << " = "  << deltaNanoTime <<  endl;
    }
    return isnewevent;
}

int FdProcessElves::NbaysEventFoundinPreScan(fevt::Eye& eye){
    TimeStamp	timeMy = eye.GetHeader().GetTimeStamp();
    bool isfound = false;
    EffNpages=0;
    int nbays=0;
    int ips=0;
    //cout << "Searching NbaysEvent in PRESCAN " << endl;
    
    while (ips < NpsPages  && !isfound){
        cout << ips << "..";
        int NextEvno = 0;
        int ThisEvno = evno[ips];
        
        int deltaGPS = timeMy.GetGPSSecond() - GPSsec[ips];
        
        //cout << "Searching PRESCAN " << ips << " "  << GPSsec[ips] << " " << deltaGPS << " " << rawMinGPSns - GPSns[ips] <<endl;
        if (deltaGPS == 0){
                int dGPS = rawMinGPSns - GPSns[ips];
            
        if (dGPS ==0){
            isfound = true;
            EffNpages=1;
            PreScanIndex = ips;
            ThisGPSsec = GPSsec[ips];
            EventNanoTime = GPSns[ips];


        if (ips+1<NpsPages) NextEvno = evno[ips+1];
        if (NextEvno > ThisEvno){
            //ProcessingEvent = true;
            //    One Bay event found
            ThisBay1=bayno[ips]; ThisBay2 = 0; nbays = 1;
            LeftMostBay = ThisBay1;
            EffNpages=1; int timesum=0;
            for (int ipg=0; ipg<NMAXPAGES-1; ipg++){
                timesum+=dts[ipg][ips];
                
                if (timesum<(105+ipg*100)*1000)EffNpages=ipg+2;
                
            }
            //if (dts[0][ips]+dts[1][ips]<205000)EffNpages=3;
            //else  if (dts[0][ips]<105000)EffNpages=2;
            cout << "1BAYEventFound " << GPSsec[ips] << " " << GPSns[ips] << " ";
            cout << evno[ips]  << " " << EffNpages<< " " << bayno[ips] << " ";
            cout << dts[0][ips] << " " << dts[1][ips] << endl ;
        } else if (NextEvno == ThisEvno) {
            //ProcessingEvent = true;
            //    Two Bays event found
            ThisBay1=bayno[ips]; ThisBay2 = bayno[ips+1]; nbays = 2;
            LeftMostBay = ThisBay1;
            if (ThisBay2 > LeftMostBay)LeftMostBay = ThisBay2;
  /*
            if (dts[0][ips]+dts[1][ips]<205000)EffNpages=3;
            else  if (dts[0][ips]<105000)EffNpages=2;
            
            if (dts[0][ips+1]+dts[1][ips+1]<205000 && EffNpages<3)EffNpages=3;
            else  if (dts[0][ips+1]<105000  && EffNpages<2)EffNpages=2;
    */
            
            EffNpages=1; int timesum=0;
            for (int ipg=0; ipg<NMAXPAGES-1; ipg++){
                timesum+=dts[ipg][ips];
                if (timesum<(105+ipg*100)*1000)EffNpages=ipg+2;
            }
            
            int tmpEffNpages=1; timesum=0;
            for (int ipg=0; ipg<NMAXPAGES-1; ipg++){
                timesum+=dts[ipg][ips+1];
                if (timesum<(105+ipg*100)*1000)tmpEffNpages=ipg+2;
            }
            if (tmpEffNpages>EffNpages)EffNpages=tmpEffNpages;

            cout << "2BAYEventFound " << GPSsec[ips] << " " << GPSns[ips] << " ";
            cout << evno[ips]  << " " << EffNpages<< " " << bayno[ips] << " ";
            cout << dts[0][ips] << " " << dts[1][ips] << endl;
            cout << "2BAYEventFound " << GPSsec[ips+1] << " " << GPSns[ips+1] << " ";
            cout << evno[ips+1]  << " " << EffNpages<< " " << bayno[ips+1] << " ";
            cout << dts[0][ips+1] << " " << dts[1][ips+1] << endl;
        }
        }}
               ips++;
    }
    if (ips == NpsPages) cout << " None found in PRESCAN " << endl;
    
        /*
        int NextGPS = 0;
        int ThisGPS = GPSsec[ips];
        
        int deltaGPS = timeMy.GetGPSSecond() - GPSsec[ips];
        int deltanGPS = timeMy.GetGPSNanoSecond() - GPSns[ips];
        
        if (deltaGPS == 0 && deltanGPS ==0){
            for (int jps=ips+1; jps<Npspages; jps++)
                if (GPSsec[jps]-ThisGPS == 0)

          if (ips+1<NpsPages) NextGPS = GPSsec[ips+1];
            if (NextGPS > ThisGPS){
                    ThisBay1=bayno[ips]; ThisBay2 = 0; nbays = 1;
                    if (dts[0][ips]+dts[1][ips]<205000)EffNpages=3;
                    else  if (dts[0][ips]<105000)EffNpages=2;
                
                 }  else {
                    ThisBay1=bayno[ips]; ThisBay2 = bayno[ips+1]; nbays = 2;
                     
                     if (dts[0][ips]+dts[1][ips]<205000)EffNpages=3;
                     else  if (dts[0][ips]<105000)EffNpages=2;
                     
                     if (dts[0][ips+1]+dts[1][ips+1]<205000 && EffNpages<3)EffNpages=3;
                     else  if (dts[0][ips+1]<105000  && EffNpages<2)EffNpages=2;
                 }
    
        }
    */
    
    
    return nbays;
}


void FdProcessElves::CalculateDelaysBetweenBays(TEyeEvent& eyeEvent, fevt::Eye& eye){

    
    TimeStamp	timeMy = eye.GetHeader().GetTimeStamp();

    TEyeEventHeader* eyeHeader =  eyeEvent.GetEventHeader();
    
    rawMinGPSns = 1000000000;
    for (UInt_t i = 0; i < FdUtil::Fd::kEYE_NMIRRORS; ++i){TelGPSns[i]=1000000000;TelHasElveTrigger[i]=false;}
    
    for (UInt_t i = 0; i < FdUtil::Fd::kEYE_NMIRRORS; ++i) {
        const UInt_t j = i + FdUtil::Fd::kEYE_FIRST_MIRROR;
        if (eyeHeader->IsMirrorDataPresent(j)) {
            const FdRoot::TTimeStamp mirrorTimeStamp =
            *eyeHeader->GetRawMirrorTimeStamp(j);
            cout  << timeMy.GetGPSSecond() << " RawGPSTime -- mirror "
            << j << "  " << mirrorTimeStamp << " Label " << eyeHeader->GetMirrorEventLabel(j)<< endl;
            if (eyeHeader->GetMirrorEventLabel(j) ==7)TelHasElveTrigger[j]=true;
            if (eyeHeader->GetMirrorEventLabel(j) ==8)TelHasElveTrigger[j]=true;
            TelGPSns[j-1]=eyeHeader->GetRawMirrorTimeStamp(j)->GetNanoSec(); cout << " DELAY " << TelGPSns[j-1] << endl;
            if (TelGPSns[j-1]<rawMinGPSns)rawMinGPSns=TelGPSns[j-1];
        }}
    
    
    for (UInt_t i = 0; i < FdUtil::Fd::kEYE_NMIRRORS; ++i){
        if(TelGPSns[i]<1000000000){TelGPSns[i]-=rawMinGPSns; cout << "TelescopeDelay mirror " << i+1 << " : " << TelGPSns[i] << endl;}
        if(TelGPSns[i]<1000000000 && TelHasElveTrigger[i])fTelescopeCount++;
    }
    //        if (telcnt<2)continue;
    
    
    const int eyeId = eyeEvent.GetEventHeader()->GetEyeNo();
    const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
    const fdet::Eye& detEye = detFD.GetEye(eyeId);
   
    for (fdet::Eye::TelescopeIterator telIter = detEye.TelescopesBegin();
         telIter != detEye.TelescopesEnd(); ++telIter) {
        
        const unsigned int mirrorId = telIter->GetId();
        
        if (eyeHeader->IsMirrorDataPresent(mirrorId)) {
            
            if (!eye.HasTelescope(mirrorId))
                eye.MakeTelescope(mirrorId);
            
            fevt::Telescope& tel = eye.GetTelescope(mirrorId);
            
            if (!tel.HasTriggerData()){ cout << "MAKE TRIGGER DATA " << endl ;
                tel.MakeTriggerData(); }
            fevt::TelescopeTriggerData& telTrigger = tel.GetTriggerData();
            
            cout << "TLTLABEL mirror " << (eyeId-1)*6+mirrorId << " "  << eyeHeader->GetMirrorEventLabel(mirrorId) << " GPStime= " << timeMy.GetGPSSecond()
            <<  "." << eyeHeader->GetRawMirrorTimeStamp(mirrorId)->GetNanoSec() << endl;
            //				}
        }
    }
    


}

void FdProcessElves::WriteOutputFile(long unsigned int gpssec, long unsigned int gpsnsec){

    fOut = new TFile("myOutput.root","recreate");
    cout << gpssec << "." << gpsnsec << " Preparing to put all histos in new file " << endl ;
    
    WriteAllHistos();

    fOut->Close();

    

    stringstream ss;
    ss<<"EVT"<<gpssec<<"-"<<gpsnsec;    string tit = ss.str();       tit += ".root";
    cout << " new file " << tit << endl ;

    char str1[110];
    
    strcpy (str1,"mv myOutput.root ");
    strcat (str1,tit.c_str());
    
    cout << str1 << endl;
    Int_t sys_err;
    
    sys_err = system(str1);
    
}

void FdProcessElves::WriteAllHistos(){

    cout << "Writing TMax, NphMax, SigmaT, Skew histos"<< endl;
    for (int ir=0; ir<NROWS; ir++){
        TimeMax_vs_Col[ir]->Write();
        NphMax_vs_Col[ir]->Write();
        SigmaT_vs_Col[ir]->Write();
        Skew_vs_Col[ir]->Write();
    }

    cout << "Writing Pulsefitting  histos"<< endl;

    
    hChi2BoltFinder->Write();
    hPhiMinBoltFinder->Write();
    hTMinBoltFinder->Write();

    cout << "Writing All trace  histos"<< endl;

    for (int itc = 0 ; itc<fL7TelescopeCount;itc++){
        for (int ipp=0; ipp<NPIXELS; ipp++){
            
            int itp = itc*NPIXELS+ipp;
            hRawPixels[itc][ipp]->Write();
            hCalPixels[itc][ipp]->Write();
            hGCorrPixels[itc][ipp]->Write();
            hCalResiduals[itc][ipp]->Write();
        }
    }

    cout << "Writing Raw and Fine Tmax histos"<< endl;

    hRawTmax->Write();
    hRawQmax->Write();
    hRawTmax2->Write();
    hRawQmax2->Write();
    hFineTmax->Write();
    hFineQmax->Write();
    hFineQCmax->Write();
    hFineQCCmax->Write();
    hFineQtot->Write();
    hFineQCtot->Write();
    hFineQCCtot->Write();
    hFineTmax2->Write();
    hFineQmax2->Write();
    hFineSigmaT->Write();
    hFineSkew->Write();
    hFIT3Tmax->Write();
    hFIT3Qmax->Write();
    hFIT3SigmaT->Write();
    hFIT3Skew->Write();
    hFitPullsN->Write();
    hFitPullsP->Write();
    hCountElves->Write();
    
    pQCC_vs_DistBolt->Write();
    hQCC_Asym->Write();
    
    for(int i=0; i<NDP1*NDP2; i++){
        pQCC_vs_ElliDist[i]->Write();
        hQCC_Spread[i]->Write();
    }
    
    hQCCAveSpread ->Write();

  
}

void FdProcessElves::ResetAllCounters(){
    nPages = 1;
    QtotTrace = 0;
    minphi = 360.;
    maxphi = 0;
    change_phi = false;
    min_charge =999999999.;
    max_charge = -9999999.;
    LeftMostBay = 0;


}

void FdProcessElves::ResetAllVectors(){
    
    cout << "Deleting all fPixels .." << endl;
    
    fPixels->Delete();
    fPixels2->Delete();
    fPixels3->Delete();
    fPixels4->Delete();
    fPixels5->Delete();
    fPixels6->Delete();
    fPixels7->Delete();
    fPixels8->Delete();
    fSites->Delete();
    fSites2->Delete();
    
    if (fB100km) delete fB100km;
    if (fB200km) delete fB200km;
    
    
    cout << "Deleting fBE1_100km ... ";
    if (fBE1_100km) delete fBE1_100km;
    cout << "Deleting fBE1_200km ... " << endl;
    if (fBE1_200km) delete fBE1_200km;
    
    
    
    cout << "Deleting all TelHits .." << endl;
    
    for (int np = 0 ; np < NMAXPAGES ; np++ ){
        for (int itt=0; itt<6;itt++) TelHit[np][itt]=false;
        for (int itt=0; itt<3; itt++) HitTel[itt]=0;
        for (int npp = 0; npp < 440; npp++)pixHasData[np][npp]=false;
    }

    cout << "Deleting all Qmx,Tmx,rawPed,caltrace..." << endl;

    for (int itc = 0 ; itc<3;itc++){
        
        for (int ipp=0; ipp<NPIXELS; ipp++){
            
            int itp = itc*NPIXELS+ipp;
            
            Qmx[itp] = 0.;
            Tmx[itp] = 0.;
            R_EFD[itp] = 0.;
            R_BE[itp] = 0.;
            Dist_BE[itp] = 0.;
            SinTh_BE[itp] = 0.;
            dT_EFD[itp] = 0.;
            rawPed[itc][ipp] = 0.;
            rawPed2[itc][ipp] = 0.;
            QtotPixel[itc][ipp]=0;
            LongPixelVert[itc]=0;
            LatPixelVert[itc]=0;
            
            for (int ns=0; ns<nslots; ns++){
                //rawTrace[itc][ipp][ns]=0;
                calTrace[itc][ipp][ns]=0;
                EmTime_ns[itc][ipp][ns]=0;
            }
        }
    }


}

void FdProcessElves::ResetAllHistos(){

    for (int ir=0; ir<NROWS; ir++){
        delete TimeMax_vs_Col[ir];
        delete NphMax_vs_Col[ir];
        delete SigmaT_vs_Col[ir];
        delete Skew_vs_Col[ir];
    }

    hChi2BoltFinder->Reset();
    hPhiMinBoltFinder->Reset();
    hTMinBoltFinder->Reset();
    for (int itc = 0 ; itc<3;itc++){
        
        for (int ir=0; ir<NROWS; ir++)hPixelPeakTimes[itc][ir]->Reset();
        for (int ipp=0; ipp<NPIXELS; ipp++){
            int itp = itc*NPIXELS+ipp;
            hRawPixels[itc][ipp]->Reset();
            hCalPixels[itc][ipp]->Reset();
            hGCorrPixels[itc][ipp]->Reset();
            
            hCalResiduals[itc][ipp]->Reset();
        }
    }

    hRawTmax->Reset();
    hRawQmax->Reset();
    hRawTmax2->Reset();
    hRawQmax2->Reset();
    
    hFineTmax->Reset();
    hFineQmax->Reset();
    hFineQCmax->Reset();
    hFineQCCmax->Reset();
    
    hFineQtot->Reset();
    hFineQCtot->Reset();
    hFineQCCtot->Reset();
    
    hFineTmax2->Reset();
    hFineQmax2->Reset();
    
    hFineSigmaT->Reset();
    hFineSkew->Reset();

    hFIT3Tmax->Reset();
    hFIT3SigmaT->Reset();
    hFIT3Skew->Reset();
    hFIT3Qmax->Reset();
    hFitPullsN->Reset();
    hFitPullsP->Reset();
    
    pQCC_vs_DistBolt->Reset();
    hQCC_Asym->Reset();
    
    for(int i=0; i<NDP1*NDP2; i++){
        pQCC_vs_ElliDist[i]->Reset();
        hQCC_Spread[i]->Reset();
    }
    hQCCAveSpread->Reset();
    

    
//    hCountElves->Reset();


}


void FdProcessElves::LastReset(){
    
    cout << " Starting LastReset " << endl;
    ResetAllHistos();
    ResetAllVectors();
    ResetAllCounters();
    
    
}

void FdProcessElves::ResetAllForNewEvent(){
    //long unsigned int gpssec, long unsigned int gpsnsec){
 
 
 ResetAllHistos();
 ResetAllVectors();
 ResetAllCounters();
 
 ProcessingEvent = false ;
 EffNpages=0;
 PreScanIndex = -1;

 
 }
/*
 
void FdProcessElves::ResetAllForNewEvent(long unsigned int gpssec, long unsigned int gpsnsec, int iopt=0){


    ResetAllHistos();
    ResetAllVectors();
    ResetAllCounters();
    
    
    ThisGPSsec = gpssec;
    EventNanoTime = gpsnsec;

    
    fPixels->Delete();
    fPixels2->Delete();
    nPages = 1;
    QtotTrace = 0;
    minphi = 360.;
    maxphi = 0;
    change_phi = false;
    
    for (int np = 0 ; np < 3 ; np++ ){
        for (int itt=0; itt<6;itt++) TelHit[np][itt]=false;
        for (int itt=0; itt<3; itt++) HitTel[itt]=0;
        for (int npp = 0; npp < 440; npp++)pixHasData[np][npp]=false;
        
    }
    min_charge =999999999.;
    max_charge = -9999999.;
    
    
    
    stringstream ss;
    ss<<"EVT"<<ThisGPSsec<<"-"<<EventNanoTime;
    string tit = ss.str();
    tit += ".root";
    
    cout << " new file " << tit << endl ;
    
    ThisGPSsec = gpssec;
    EventNanoTime = gpsnsec;

    if (iopt==2) return;
    
    fOut = new TFile("myOutput.root","recreate");


    for (int ir=0; ir<NROWS; ir++){
        TimeMax_vs_Col[ir]->Write();
        delete TimeMax_vs_Col[ir];
        NphMax_vs_Col[ir]->Write();
        delete NphMax_vs_Col[ir];
        SigmaT_vs_Col[ir]->Write();
        delete SigmaT_vs_Col[ir];
        Skew_vs_Col[ir]->Write();
        delete Skew_vs_Col[ir];

    }
    
    hChi2BoltFinder->Write();
    hChi2BoltFinder->Reset();

    hPhiMinBoltFinder->Write();
    hPhiMinBoltFinder->Reset();
    hTMinBoltFinder->Write();
    hTMinBoltFinder->Reset();
    
    for (int itc = 0 ; itc<3;itc++){
        
        for (int ir=0; ir<NROWS; ir++)hPixelPeakTimes[itc][ir]->Reset();
            
        
        for (int ipp=0; ipp<NPIXELS; ipp++){
            
            int itp = itc*NPIXELS+ipp;
            hRawPixels[itc][ipp]->Write();
            hCalPixels[itc][ipp]->Write();

            hCalResiduals[itc][ipp]->Write();

            //cout << " Resetting hRawPixels" << ipp << endl;
            hRawPixels[itc][ipp]->Reset();
            hCalPixels[itc][ipp]->Reset();
            hCalResiduals[itc][ipp]->Reset();

            //hRawTmax2->Reset();
            //hRawQmax2->Reset();
            //hFineTmax->Reset();
            //hFineQmax->Reset();

           // if (itp%10 == 0) cout << itp << endl;
            
            Qmx[itp] = 0.;
            Tmx[itp] = 0.;
            rawPed[itc][ipp] = 0.;
            rawPed2[itc][ipp] = 0.;
            QtotPixel[itc][ipp]=0;
            for (int ns=0; ns<nslots; ns++){
                //rawTrace[itc][ipp][ns]=0;
                calTrace[itc][ipp][ns]=0;
                EmTime_ns[itc][ipp][ns]=0;
            }
            //for (int ip=0; ip<3; ip++){
            //pixPed[itc][ip][ipp] = -100.;
            //}
            //for (int isl=0; isl<150; isl++) {
            //	 pixel_slot[itc][ipp][isl]=0.;
            //}
        }
    }
    
    hRawTmax->Write();    hRawTmax->Reset();

    hRawQmax->Write();    hRawQmax->Reset();
 
    hFineTmax->Write();    hFineTmax->Reset();
    
    hFineQmax->Write();    hFineQmax->Reset();

    hFineSigmaT->Write();    hFineSigmaT->Reset();
    
    hFineSkew->Write();    hFineSkew->Reset();
    hFIT3Tmax->Write();    hFIT3Tmax->Reset();

    
    fOut->Close();
    char str1[110];
   
    strcpy (str1,"mv myOutput.root ");
    strcat (str1,tit.c_str());
    
    cout << str1 << endl;
    Int_t sys_err;

    sys_err = system(str1);
    
}
*/

void FdProcessElves::ProcessRawData(TEyeEvent& eyeEvent, fevt::Eye& eye){

    
    GlueTraces(eyeEvent,eye);
    
    
}

void FdProcessElves::ProcessHistos(){

    //if( nPages < NMAXPAGES ) return;
    NpixCount1=0;
    NpixCount2=0;
    
    // pixel loop
    for (int it=0; it<fL7TelescopeCount; it++){
        
        cout << "TELESCOPE " << it+1 << endl;
        for (int ip=0; ip<NPIXELS ; ip++){
            // for each pixel , fits pulse with asym gaussian and gives PeakAmpli, PeakTime, SigmaRight, SigmaLeft
            // with proper errors
            FdElvesPulseFinder(it,ip);
            // for each pixel, after 1st fit, check residuals and tries to refit the residual with a 2nd asym gauss
            //FdElvesExtraPulseFinder(it,ip);
        }
    }
    
    hCountElves->Fill(NpixCount1,NpixCount2);
    
    cout << ThisGPSsec << " DOUBLESearch: Npixels with 1,2 pulses: " << NpixCount1 << " " <<  NpixCount2 << endl;
    
    bool IsDoubleElves=false;
    if (NpixCount2>NP2MIN)IsDoubleElves=true;
  
    
    if (!IsDoubleElves) {
    // preliminary estimate of bolt location, assuming Hd=80 km and Hbolt=0km
        TVector3 RawBoltLocation=FdSingleElvesBoltFinder();

        double RawBoltLatitude = 90.-RawBoltLocation.Theta()/degree;
        if (RawBoltLatitude>89.) {cout << " Raw Bolt Location FAILED : event REJECTED " << ThisGPSsec<< endl; return;}
        double RawBoltLongitude = RawBoltLocation.Phi()/degree;

        double RawBoltDistance = Rearth*(RawBoltLocation.Angle(Site1Location));
        cout << " RAWBOLT Latitude Longitude Distance" << RawBoltLatitude << " " << RawBoltLongitude << " "  << RawBoltDistance << endl;
    
        double RawBoltT0 = 28.-RawBoltDistance/SpeedOfLight;
    
    
        // the fitter loops also on Hd , therefore  also geometry corrections have to be postponed
        //TVector3 FineBoltLocation=
    
        // Summary PLot of all traces
        if (FdAllTraces_CameraViewPlot()>0){
    
        // Prepares the dummy markers for showing the elves Camera View
        FdSingleElvesPixel_CameraViewIni();

        // Fits the lighning location starting from the RawBoltLocation
            if (FdSingleElvesBoltFitter(RawBoltLocation,RawBoltT0 )>10){
        
                CameraView_4TPLOTS();
                
                FdElvesGeometryCorrection();
        
                FdElvesAtmosphericCorrection();
        
        

                CameraView_4QPLOTS();
 
                LongLatView_QPLOTS();
            }
        }
        //LightEmission_vs_EllipticDistancePlot();
        
        //FdSingleElvesQCCtot_LongLatViewPlot();

        //LightEmission_vs_DistancePlot();
        //FdSingleElvesQCCtot_LongLatViewPlot();
        //IsoRDistance_LongLatViewDraw();

        
        //LightEmission_vs_EllipticDistancePlot();
        //FdSingleElvesQCCtot_LongLatViewPlot();
        //IsoEDistance_LongLatViewDraw();
        
    } else {

        TVector3 RawBoltLocation=FdDoubleElvesBoltFinder();
        
        double RawBoltLatitude = 90.-RawBoltLocation.Theta()/degree;
        if (RawBoltLatitude>89.) {cout << " Raw Bolt Location FAILED : event REJECTED " << ThisGPSsec<< endl; return;}
        double RawBoltLongitude = RawBoltLocation.Phi()/degree;
        
        double RawBoltDistance = Rearth*(RawBoltLocation.Angle(Site1Location));
        cout << " RAWBOLT Latitude Longitude Distance" << RawBoltLatitude << " " << RawBoltLongitude << " "  << RawBoltDistance << endl;
        
        double RawBoltT0 = 28.-RawBoltDistance/SpeedOfLight;
        

        
        // Summary PLot of all traces
        if (FdAllTraces_CameraViewPlot()>0){
        
        // Prepares the dummy markers for showing the elves Camera View
        FdSingleElvesPixel_CameraViewIni();

        
        // Fits the lighning location starting from the RawBoltLocation
            if (FdSingleElvesBoltFitter(RawBoltLocation,RawBoltT0 )>10){
        
                CameraView_4TPLOTS();
        
                FdElvesGeometryCorrection();
        
                FdElvesAtmosphericCorrection();
        
        

                CameraView_4QPLOTS();
        
                LongLatView_QPLOTS();
            }
        //LightEmission_vs_EllipticDistancePlot();

        //FdSingleElvesQCCtot_LongLatViewPlot();
        }
    
    }
//    FdSingleElvesSigmaT_CameraViewPlot();
    
    // after this we have Long, Lat, H_emission, H_bolt, T_bolt
    
    // Following parts will probably go to UserModule
    //FdElvesCalculateIonosphereLightEmissionDensity();
    //  - FdElvesGeometryCorrection;       // depends on Hd
    //  - FdElvesAtmosphericCorrection;    // depends on local weather
    //FdElvesAngularDistribution();
    
    
}



void  FdProcessElves::CameraView_4TPLOTS(){

    cFDViewPlot = new TCanvas("cFDViewPlot","FD View ", 1400, 1000);
    cFDViewPlot -> Divide(2,2);

    //  shows FD Camera View Color Plots of : Time, SigmaT_left, Photons at Diaphragm,  SigmaT_right
    FdSingleElvesTime_CameraViewPlot(1);
    //FdSingleElvesPull_CameraViewPlot(2);
    FdSingleElvesSigmaT_CameraViewPlot(3);
    FdSingleElvesSigmaR_CameraViewPlot(4);

    stringstream ss3;
    ss3<< "EVTS" << ThisGPSsec << "." << EventNanoTime << "_4TPLOTS_FDVIEW"  << ".png";

    cFDViewPlot->SaveAs(ss3.str().c_str());

    delete cFDViewPlot;

}



void  FdProcessElves::CameraView_4QPLOTS(){
    
    cFDViewPlot2 = new TCanvas("cFDViewPlot2","FD View ", 1400, 1000);
    cFDViewPlot2 -> Divide(2,2);
    
    FdSingleElvesQmax_CameraViewPlot(1);
    FdSingleElvesQCtot_CameraViewPlot(3);
    FdSingleElvesQCCtot_CameraViewPlot(4);
    stringstream ss4;
    ss4<< "EVTS" << ThisGPSsec << "." << EventNanoTime << "_4QPLOTS_FDVIEW"  << ".png";
    
    cFDViewPlot2->SaveAs(ss4.str().c_str());
    
    delete cFDViewPlot2;
    
}

void  FdProcessElves::LongLatView_QPLOTS(){
    
    TCanvas * cQC_vs_RDist = new TCanvas("cQC_vs_RDist","Corrected charge vs Radial Distance",200,10,1600,800);
    cout << "Drawing cQC_vs_RDist " << endl;
    cQC_vs_RDist->Divide(2,1);
    cQC_vs_RDist->cd(1);                          LightEmission_vs_DistancePlot();
    cQC_vs_RDist->cd(2);                          FdSingleElvesQCCtot_LongLatViewPlot(cQC_vs_RDist);
    //IsoRDistance_LongLatViewDraw();
    
    
    stringstream ss4;
    ss4<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_QPLOTS_LLVIEW_RDIST.png"  ;
    
    cQC_vs_RDist->SaveAs(ss4.str().c_str());
    
    delete cQC_vs_RDist;
 

    
     TCanvas * cQC_vs_EDist = new TCanvas("cQC_vs_EDist","Corrected charge vs Elliptic Distance",200,10,1600,800);
    cout << "Drawing cQC_vs_EDist " << endl;

    cQC_vs_EDist->Divide(2,1);
    cQC_vs_EDist->cd(1);                   int imin = LightEmission_vs_EllipticDistancePlot();
    cQC_vs_EDist->cd(2);                          FdSingleElvesQCCtot_LongLatViewPlot2(cQC_vs_EDist);

    stringstream ss4b;
    ss4b<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_QPLOTS_LLVIEW_EDIST.png"  ;
    
    cQC_vs_EDist->SaveAs(ss4b.str().c_str());
    
    delete cQC_vs_EDist;
    
    FdElvesBestEllipseSearch(imin);

    
}



void  FdProcessElves::FdElvesGeometryCorrection(){
// Calculate Geometry correction : projects into Ionosphere surface, calculate Area

    cout << ThisGPSsec << " LEFTMOSTBAY,site  " << LeftMostBay << " " << site1 << endl;
    for (int itc = 0 ; itc<fL7TelescopeCount;itc++){
        for (int ipp=0; ipp<NPIXELS; ipp++){
        
            int itp = itc*NPIXELS+ipp;
            
            
            
            int lmbay = (LeftMostBay-1)%6;
            double GeomCorr = CalculateGeomCorr(itp,site1,lmbay) ;
            for (int its=0; its<50*NMAXPAGES; its++){
                float Qt = hCalPixels[itc][itp]->GetBinContent(its+1);
                hGCorrPixels[itc][itp]->SetBinContent(its+1,Qt*GeomCorr);
            }
            float  Qtot2s = GeomCorr*hFineQtot->GetBinContent(itp+1);
            hFineQCtot->SetBinContent(itp+1,Qtot2s);


            
        }
    }
}


void  FdProcessElves::FdElvesAtmosphericCorrection(){
    // Calculate Atmospheric correction : uses total Optical Depth of the Atmosphere
    
    for (int itc = 0 ; itc<fL7TelescopeCount;itc++){
        for (int ipp=0; ipp<NPIXELS; ipp++){
            
            int itp = itc*NPIXELS+ipp;
            
            double AtmoCorr = CalculateAtmoCorr(ipp) ;
            for (int its=0; its<50*NMAXPAGES; its++){
                float Qt = hGCorrPixels[itc][itp]->GetBinContent(its+1);
                hTCorrPixels[itc][itp]->SetBinContent(its+1,Qt*AtmoCorr);
            }
            float Qtot2s = AtmoCorr*hFineQCtot->GetBinContent(itp+1);
            hFineQCCtot->SetBinContent(itp+1,Qtot2s);
            
            
            
        }
    }
    
}


double FdProcessElves::CalculateAtmoCorr(int npix){

    //double Airmass_as= 1./(cos(Theta_s)+0.050572*pow(96.07995-Theta_s/degree,-1.6364));
    double VMOD = 0.54; // just from a simple calculation of Vertical Total Optical Depth
    double VAOD = 0.08;
    double VTOD = VAOD+VMOD;
    // then I assume that airmass accounts only for the sectheta term
    double Theta_o = Theta_pixel[npix];
    
    double Airmass_oa = 1./(cos(Theta_o)+0.50572*pow(96.07995-Theta_o/degree,-1.6364));
    
    
    if (cos(Theta_o)<0.) Airmass_oa=200.;
    
    //Weigh_oa *= exp(-Airmass_as*vtod);//this is negligible
    //cout << Theta_o << " " <<  arg_oa << " AM O-A " << Airmass_oa << endl;
    double atmocorr =  exp(Airmass_oa*VTOD);// /pow(MirrorArea,2);

    return atmocorr;
}


void  FdProcessElves::LightEmission_vs_AnglePlot(){

    TGraph *QC_vs_SinThBE[3]; // we define THREE TGraph  to account for fully contained , partially constinaed
    
    QC_vs_SinThBE[0]  =new TGraph();
    QC_vs_SinThBE[0]->SetMarkerSize(1.0);
    QC_vs_SinThBE[0]->SetMarkerColor(3);
    QC_vs_SinThBE[0]->SetMarkerStyle(24);
    
    QC_vs_SinThBE[1]  =new TGraph();
    QC_vs_SinThBE[1]->SetMarkerSize(1.0);
    QC_vs_SinThBE[1]->SetMarkerColor(2);
    QC_vs_SinThBE[1]->SetMarkerStyle(24);
    
    QC_vs_SinThBE[2]  =new TGraph();
    QC_vs_SinThBE[2]->SetMarkerSize(1.0);
    QC_vs_SinThBE[2]->SetMarkerColor(1);
    QC_vs_SinThBE[2]->SetMarkerStyle(24);
    
    
    for (int rp = 0; rp<NROWS; rp++) {
        //cout << " Drawing QC for row " << rp << endl;
        int rowcolor = rp+45 + 15*(rp%3);
        QC_vs_SinThBE_RP[rp] = new TGraph();
        QC_vs_SinThBE_RP[rp]->SetMarkerSize(1.0);
        //QC_vs_DistBoltRP[rp]->SetMarkerColor(44*(rp%2+1)-2*rp+15);
        QC_vs_SinThBE_RP[rp]->SetMarkerColor(rowcolor);
        QC_vs_SinThBE_RP[rp]->SetMarkerStyle(20+rp%4);
        stringstream ssp;
        ssp<< "RowN" << rp+1;
        QC_vs_SinThBE_RP[rp]->SetName(ssp.str().c_str());
        
    }
   
    
    for (int cp = 0; cp<NCOLS*fL7TelescopeCount; cp++) {
        
        QC_vs_SinThBE_CP[cp] = new TGraph();
        float ms = 0.5+0.1*abs(cp-10*fL7TelescopeCount+0.5);
        QC_vs_SinThBE_CP[cp]->SetMarkerSize(ms);
        QC_vs_SinThBE_CP[cp]->SetMarkerColor(1);
        QC_vs_SinThBE_CP[cp]->SetMarkerStyle(32);
        if (cp<10*fL7TelescopeCount)QC_vs_SinThBE_CP[cp]->SetMarkerStyle(26);
        stringstream ssp;
        ssp<< "ColN" << cp+1;
        QC_vs_SinThBE_CP[cp]->SetName(ssp.str().c_str());
        
        
    }

    
    
    int telTmin = HitTel[0];
    
    //  if (fL7TelescopeCount == 2 && HitTel[0] == HitTel[1]+1)xc[1]=1;
    //  else if (fL7TelescopeCount == 2 && HitTel[1] == HitTel[0]+1){xc[0]=1; telTmin = HitTel[1];}
    
    double NphMax[4] = {};
    
    int MINROW = 1; int MAXROW = 22;

    int np[3]={}, nrp[NROWS]={}, ncp[NCOLS]={};
    
    //for (int npix=0 ; npix < NPIXELS; npix++){
    for (int itc = 0 ; itc<fL7TelescopeCount;itc++){
        for (int ipp=0; ipp<NPIXELS; ipp++){
            
            int npix = itc*NPIXELS+ipp;
            int rp = npix%NROWS;
            int cp = int(npix/NROWS);
            
            double SinThDipole=SinTh_BE[npix];
            
            double Nph = hFineQCCtot->GetBinContent(npix+1);
            for (int i=0; i<4; i++){
                if (Nph>NphMax[i]){
                    for (int j=3;j>i;j--)NphMax[j]=NphMax[j-1];
                    
                    NphMax[i]=Nph;
                    break;
                }
            }
            
            
            if (rp>=MINROW-1 && rp<MAXROW-1 && cp<NCOLS*fL7TelescopeCount&& cp>=0 && Nph > 0.){
                //QC_vs_Dist[0]->Fill(IonosphereBoltPixelDist,Nph);
                
                
                QC_vs_SinThBE[0]->SetPoint(np[0]++,SinThDipole,Nph);//*Geom_corr[npix]);
                //pQCC_vs_SinThBE->Fill(SinThDipole,Nph);
                
                QC_vs_SinThBE_RP[rp]->SetPoint(nrp[rp]++,SinThDipole,Nph);
                QC_vs_SinThBE_CP[cp]->SetPoint(ncp[cp]++,SinThDipole,Nph);
                
                //if(Nph>Qttmax)Qttmax=Nph;
                
            }
        }
        
    }
    
    
    TCanvas * cQC_vs_SinThBE = new TCanvas("cQC_vs_SinThBE","Corrected charge vs SinThBE",200,10,800,800);
    // cQC_vs_Dist->Divide(3,1);
    //QC_vs_Dist[2]->Reset();
    cout << "Drawing cQC_vs_SinThBE " << endl;
    
    cQC_vs_SinThBE->cd();
    //for (int ii=0; ii<3; ii++){cQC_vs_Dist->cd(ii+1);QC_vs_Dist[ii]->Draw("zcol");}
    hQC_vs_SinThBE[0]->SetMaximum(1.2*NphMax[3]); hQC_vs_SinThBE[0]->Draw();
    //QC_vs_SinThBE[0]->Draw("P");
    
    
    TLegend *leg1b = new TLegend(0.75,0.15,0.9,0.9);
    leg1b->SetTextSize(0.025);
    
    for (int rp=0; rp<NROWS; rp++){
        QC_vs_SinThBE_RP[rp]->Draw("P");
        stringstream ssp;
        ssp<< "Row" << rp+1;
        stringstream sspn;
        sspn<< "RowN" << rp+1;
        leg1b->AddEntry(sspn.str().c_str(),ssp.str().c_str(),"p");
    }
    leg1b->Draw("P");
    
    TLegend *leg2b = new TLegend(0.6,0.15,0.75,0.9);
    leg2b->SetTextSize(0.025/fL7TelescopeCount);
    
    for (int cp=0; cp<NCOLS*fL7TelescopeCount; cp++){
        QC_vs_SinThBE_CP[cp]->Draw("P");
        stringstream ssp;
        ssp<< "Col" << cp+1;
        stringstream sspn;
        sspn<< "ColN" << cp+1;
        leg2b->AddEntry(sspn.str().c_str(),ssp.str().c_str(),"p");
    }
    leg2b->Draw("P");
    
    stringstream ss3b;
    ss3b<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_NPHCC_vs_SinTh"  << ".png";
    
    cQC_vs_SinThBE->SaveAs(ss3b.str().c_str());
    
    delete cQC_vs_SinThBE;
    


}

void  FdProcessElves::LightEmission_vs_DistancePlot(){

    // initialize QC_vs_DistBolt[3]
    TGraph *QC_vs_DistBolt[3]; // we define THREE TGraph  to account for fully contained , partially constinaed
    
    QC_vs_DistBolt[0]  =new TGraph();
    QC_vs_DistBolt[0]->SetMarkerSize(1.0);
    QC_vs_DistBolt[0]->SetMarkerColor(3);
    QC_vs_DistBolt[0]->SetMarkerStyle(24);
    
    QC_vs_DistBolt[1]  =new TGraph();
    QC_vs_DistBolt[1]->SetMarkerSize(1.0);
    QC_vs_DistBolt[1]->SetMarkerColor(2);
    QC_vs_DistBolt[1]->SetMarkerStyle(24);
    
    QC_vs_DistBolt[2]  =new TGraph();
    QC_vs_DistBolt[2]->SetMarkerSize(1.0);
    QC_vs_DistBolt[2]->SetMarkerColor(1);
    QC_vs_DistBolt[2]->SetMarkerStyle(24);

    // initialize QC_vs_DistBoltRP
    int np[3]={}, nrp[NROWS]={}, ncp[NCOLS]={};
    for (int rp = 0; rp<NROWS; rp++) {
        //cout << " Drawing QC for row " << rp << endl;
        int rowcolor = rp+45 + 15*(rp%3);
        QC_vs_DistBoltRP[rp] = new TGraph();
        QC_vs_DistBoltRP[rp]->SetMarkerSize(1.0);
        //QC_vs_DistBoltRP[rp]->SetMarkerColor(44*(rp%2+1)-2*rp+15);
        QC_vs_DistBoltRP[rp]->SetMarkerColor(rowcolor);
        QC_vs_DistBoltRP[rp]->SetMarkerStyle(20+rp%4);
        stringstream ssp;
        ssp<< "Row" << rp+1;
        QC_vs_DistBoltRP[rp]->SetName(ssp.str().c_str());
        
    }
    
    // initialize QC_vs_DistBoltCP
    for (int cp = 0; cp<NCOLS*fL7TelescopeCount; cp++) {
        
        QC_vs_DistBoltCP[cp] = new TGraph();
        float ms = 0.5+0.1*abs(cp-10*fL7TelescopeCount+0.5);
        QC_vs_DistBoltCP[cp]->SetMarkerSize(ms);
        QC_vs_DistBoltCP[cp]->SetMarkerColor(1);
        QC_vs_DistBoltCP[cp]->SetMarkerStyle(32);
        if (cp<10*fL7TelescopeCount)QC_vs_DistBoltCP[cp]->SetMarkerStyle(26);
        stringstream ssp;
        ssp<< "Col" << cp+1;
        QC_vs_DistBoltCP[cp]->SetName(ssp.str().c_str());

        
    }
 
    int MINROW = 1; int MAXROW = 22;
    
    int telTmin = HitTel[0];
    
  //  if (fL7TelescopeCount == 2 && HitTel[0] == HitTel[1]+1)xc[1]=1;
  //  else if (fL7TelescopeCount == 2 && HitTel[1] == HitTel[0]+1){xc[0]=1; telTmin = HitTel[1];}
    
    double NphMax[4] = {};
    
    //for (int npix=0 ; npix < NPIXELS; npix++){
    for (int itc = 0 ; itc<fL7TelescopeCount;itc++){
        for (int ipp=0; ipp<NPIXELS; ipp++){
                
            int npix = itc*NPIXELS+ipp;
            int rp = npix%NROWS;
            int cp = int(npix/NROWS);

        //    dT_BFD[npix] = PathFromBolt_to_FD(npix,BestFitParameter)/SpeedOfLight;
            double IonosphereBoltPixelDist=R_BE[npix];
            double SinThDipole=SinTh_BE[npix];

            double Nph = hFineQCCtot->GetBinContent(npix+1);
            for (int i=0; i<4; i++){
                if (Nph>NphMax[i]){
                    for (int j=3;j>i;j--)NphMax[j]=NphMax[j-1];
                    
                    NphMax[i]=Nph;
                    break;
                }
            }
            
            //cout << "QCCHECK_R"<<rp+1<<"_C"<<cp+1<<" "<< IonosphereBoltPixelDist << " "<< Nph <<" "<<endl;
            cout << IonosphereBoltPixelDist << " " << Nph << " NPHMAX " << NphMax[0] << " " << NphMax[1] << " "<< NphMax[2] << " "<< NphMax[3] << endl;
          
            if (rp>=MINROW-1 && rp<MAXROW-1 && cp<20*fL7TelescopeCount&& cp>=0 && Nph > 0.){
                //QC_vs_Dist[0]->Fill(IonosphereBoltPixelDist,Nph);
                
                pQCC_vs_DistBolt->Fill(IonosphereBoltPixelDist,Nph);

                QC_vs_DistBoltRP[rp]->SetPoint(nrp[rp]++,IonosphereBoltPixelDist,Nph);
                QC_vs_DistBoltCP[cp]->SetPoint(ncp[cp]++,IonosphereBoltPixelDist,Nph);
            
            }
        }
    
    }
    
    //for (int ii=0; ii<3; ii++){cQC_vs_Dist->cd(ii+1);QC_vs_Dist[ii]->Draw("zcol");}
    hQC_vs_Dist[0]->SetMaximum(1.2*NphMax[3]); hQC_vs_Dist[0]->Draw();
    //QC_vs_DistBolt[0]->Draw("P");
    
    
    TLegend *leg1 = new TLegend(0.75,0.15,0.9,0.9);
    leg1->SetTextSize(0.025);

    for (int rp=0; rp<NROWS; rp++){
        QC_vs_DistBoltRP[rp]->Draw("P");
        stringstream ssp;
        ssp<< "Row" << rp+1;
        leg1->AddEntry(ssp.str().c_str(),ssp.str().c_str(),"p");
    }
    leg1->Draw("P");

    TLegend *leg2 = new TLegend(0.6,0.15,0.75,0.9);
    leg2->SetTextSize(0.025/fL7TelescopeCount);
    
    for (int cp=0; cp<NCOLS*fL7TelescopeCount; cp++){
        QC_vs_DistBoltCP[cp]->Draw("P");
        stringstream ssp;
        ssp<< "Col" << cp+1;
        leg2->AddEntry(ssp.str().c_str(),ssp.str().c_str(),"p");
    }
    leg2->Draw("P");

 

   hQC_vs_Dist[0]->Reset();

    
}


void  FdProcessElves::LightEmission_vs_DistanceFit(){



cout << "Fitting cQC_vs_Dist " << endl;

double NPHmaxRaw = 0.;
double RawRingR  = 0.;

for (int k=0; k<NDBbins; k++){
    double NPH_R = pQCC_vs_DistBolt->GetBinContent(k+1);
    
    if (NPH_R>NPHmaxRaw){NPHmaxRaw=NPH_R;RawRingR = k*400./NDBbins;}
    
}




if (NPHmaxRaw>0.){
    
    TCanvas * cQC_vs_Dist2 = new TCanvas("cQC_vs_Dist2","Corrected charge vs Distance",200,10,1600,800);
    
    cQC_vs_Dist2->Divide(2,1);
    
    //pQCC_vs_DistBolt->Fit("gaus");
    
    int ntotbins = nPages*size/nSlotInterval;
    fitFcn = new TF1("fitFcn",P2exp,0,500,2);
    /*
     fitFcn -> SetParameters(NPHmaxRaw,RawRingR , 0.2*RawRingR ,1.);
     fitFcn -> SetParLimits(1,0.,250.);
     fitFcn -> SetParLimits(2,0.02*RawRingR,RawRingR);
     fitFcn -> SetParLimits(3,0.2,1.5);
     */
    
    float B=RawRingR;
    float A=NPHmaxRaw/(exp(-2.)*B*B );
    fitFcn -> SetParameters(A,B);
    pQCC_vs_DistBolt->Fit("fitFcn","R");
    //hCalPixels[it][ip]->Fit("fitFcn","R");
    float RingRadius  = fitFcn->GetParameter(1) ;
    
    float fitNPHmax = pow(RingRadius,2)*exp(-2)*fitFcn->GetParameter(0) ;
    //float RingSigma   = fitFcn->GetParameter(2) ;
    //float RingSkew    = fitFcn->GetParameter(3)  ;
    
    
    // TFitResultPtr r = pQCC_vs_DistBolt->Fit("gaus","S");
    // This is the fit to get the column with the first pixel in a row
    /*Double_t chi2   = r->Chi2(); // to retrieve the fit chi2
     Double_t fitNPHmax   = r->Value(0); // retrieve the value for the parameter 0
     Double_t RingRadius    = r->Value(1); // retrieve the value for the parameter 1
     Double_t RingSigma   = r->Value(2); // retrieve the value for the parameter 2
     */
    cQC_vs_Dist2->cd(1);hQCC_Asym->Draw();
    
    for (int k=0; k<NDBbins; k++){
        double sigmaQcc = pQCC_vs_DistBolt->GetBinError(k+1);
        double NPH_R = pQCC_vs_DistBolt->GetBinContent(k+1);
        
        if (NPH_R/fitNPHmax>0.01) hQCC_Asym->Fill(sigmaQcc/NPH_R);
        
    }
    
    cQC_vs_Dist2->cd(2); pQCC_vs_DistBolt->Draw();
    
    
    
    
    double RelAziAsymmetry = hQCC_Asym->GetMean();
    // double RingFWHM = 1.18*(RingSigma*(1.+RingSkew));
    
    cout << ThisGPSsec << " NPH_vs_DistBolt Parameters Ring Radius,nPHmax,spread " << RingRadius  << " " << fitNPHmax << " " << RelAziAsymmetry <<  endl;
    
    
    
    stringstream ss3b;
    ss3b<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_NPHCC_vs_DIST2"  << ".png";
    
    cQC_vs_Dist2->SaveAs(ss3b.str().c_str());
    
    
    
}

}

int  FdProcessElves::LightEmission_vs_EllipticDistancePlot(){

    
   
    int MINROW = 1; int MAXROW = 22;

    for (int itc = 0 ; itc<fL7TelescopeCount;itc++){
        for (int ipp=0; ipp<NPIXELS; ipp++){
            
            int npix = itc*NPIXELS+ipp;
            int rp = npix%NROWS;
            int cp = int(npix/NROWS);
            
            //    dT_BFD[npix] = PathFromBolt_to_FD(npix,BestFitParameter)/SpeedOfLight;
            float Nph = hFineQCCtot->GetBinContent(npix+1);
            
            
            if (rp>=MINROW-1 && rp<MAXROW-1 && cp<NCOLS*fL7TelescopeCount&& cp>=0 && Nph > 0.){
                
                // loop on ellipses, changing one focus with respect to the bolt location
                for (int i=0; i<NDP1*NDP2; i++){
                    int ii=i%NDP1;
                    int jj=i/NDP1;
                    double delta = (ii+1)*DDELTA_KM; double phi0 = jj*(360./NDP2)*degree;
                    float ElliDist1 = CalculateEllipticDistance(npix,delta,phi0);
                    pQCC_vs_ElliDist[i]->Fill(ElliDist1,Nph);
                    
                }
                
                
            }
            
            // int mx=22-cp;
            
            //hqtot[rp]->Fill(mx,Nph);
            //hdist[rp]->Fill(mx,BoltPixelDist[npix]);
        }// ends loop on NPIXELS
        
    }// ends loop on TELESCOPES
    
    int imin=-1; double MinRelSpread=100000.;
  
    
    // Loop on all values of delta and phi0 to find the minimum spread
    
    for(int i=0;i<NDP1*NDP2; i++){
        
        double NphMax[4]={};
        
        for (int k=0; k<NDBbins; k++){
            double NPH_R = pQCC_vs_ElliDist[i]->GetBinContent(k+1);
 
            for (int i=0; i<4; i++){
                if (NPH_R>NphMax[i]){
                    for (int j=3;j>i;j--)NphMax[j]=NphMax[j-1];
                    
                    NphMax[i]=NPH_R;
                    break;
                }
            }

            
        }
    
        for (int k=0; k<NDBbins; k++){
            double NPH_R = pQCC_vs_ElliDist[i]->GetBinContent(k+1);
    
    
            if (NPH_R<NphMax[3]*1.2 && NPH_R > 0  ){
        
        
                double sigmaNPH = pQCC_vs_ElliDist[i]->GetBinError(k+1);
                double NPH_R = pQCC_vs_ElliDist[i]->GetBinContent(k+1);
            
                if (sigmaNPH>0){hQCC_Spread[i]->Fill(sigmaNPH/NPH_R);
                    cout<< " " << sigmaNPH/NPH_R ;}
            }
        }
        cout << endl;
        float RelSpread = hQCC_Spread[i]->GetMean();
        
        if (hQCC_Spread[i]->GetEntries()<4)RelSpread=0;
        
        int ii=i%NDP1;
        int jj=i/NDP1;
        double delta = (ii+1)*DDELTA_KM; double phi0 = jj*(360./NDP2);
        cout << ThisGPSsec << " AVESPREAD d(km),phi(deg),RelSpread " << delta << " " << phi0 << " " << RelSpread << endl;
        
         hQCCAveSpread->SetBinContent(ii+1,jj+1,RelSpread);
        if (RelSpread<MinRelSpread && RelSpread>0.){imin = i; MinRelSpread = RelSpread;}
        
    }
    
    //  if a minimum was found we plot the final results
    
    
    
    int nrp[NROWS]={}; int ncp[NCOLS]={};
   
    if (imin>=0){
        cout << ThisGPSsec << " MINIMUM SPREAD IN BIN " << imin << " " << MinRelSpread << endl;
        for (int rp = 0; rp<NROWS; rp++) {
            //cout << " Drawing QC for row " << rp << endl;
            int rowcolor = rp+45 + 15*(rp%3);
            QC_vs_ElliDistRP[rp] = new TGraph();
            QC_vs_ElliDistRP[rp]->SetMarkerSize(1.0);
            //QC_vs_ElliDistRP[rp]->SetMarkerColor(44*(rp%2+1)-2*rp+15);
            QC_vs_ElliDistRP[rp]->SetMarkerColor(rowcolor);
            QC_vs_ElliDistRP[rp]->SetMarkerStyle(20+rp%4);
            stringstream ssp;
            ssp<< "Row" << rp+1;
            QC_vs_ElliDistRP[rp]->SetName(ssp.str().c_str());
            
        }
        
        
        
        for (int cp = 0; cp<NCOLS*fL7TelescopeCount; cp++) {
            
            QC_vs_ElliDistCP[cp] = new TGraph();
            float ms = 0.5+0.1*abs(cp-10*fL7TelescopeCount+0.5);
            QC_vs_ElliDistCP[cp]->SetMarkerSize(ms);
            QC_vs_ElliDistCP[cp]->SetMarkerColor(1);
            QC_vs_ElliDistCP[cp]->SetMarkerStyle(32);
            if (cp<10*fL7TelescopeCount)QC_vs_ElliDistCP[cp]->SetMarkerStyle(26);
            stringstream ssp;
            ssp<< "Col" << cp+1;
            QC_vs_ElliDistCP[cp]->SetName(ssp.str().c_str());
            
            
        }
        
        
        
        int ii=imin%NDP1;
        int jj=imin/NDP1;
        
        double delta = (ii+1)*DDELTA_KM; double phi0 = jj*(360./NDP2)*degree;
        DrawSkewedDipole(delta,phi0);
        
        
        for (int npix=0; npix<fTelescopeCount*NPIXELS; npix++){
            int nnpix = npix%NPIXELS;
            int bay=npix/NPIXELS;
            int rp = npix%NROWS;
            int cp = int(npix/NROWS);

            float ElliDist1 = CalculateEllipticDistance(npix,delta,phi0);
            double Nph = hFineQCCtot->GetBinContent(npix+1);
        
            QC_vs_ElliDistRP[rp]->SetPoint(nrp[rp]++,ElliDist1,Nph);
            QC_vs_ElliDistCP[cp]->SetPoint(ncp[cp]++,ElliDist1,Nph);
            
            

        }
        
        //double delta = 1000.; double phi0 = 90.*degree;
        //DrawSkewedDipole(delta,phi0);

    }
    
    
        pQCC_vs_ElliDist[imin]->SetMarkerStyle(21);
        pQCC_vs_ElliDist[imin]->SetMarkerColor(2);
        //pQCC_vs_DistBolt->SetMarkerStyle(25);
        //pQCC_vs_DistBolt->SetMarkerColor(1);
        pQCC_vs_ElliDist[imin]->Draw("E");
        //pQCC_vs_DistBolt->Draw("SAME");
    for (int rp=0; rp<NROWS; rp++)                   if(nrp[rp]>0)   QC_vs_ElliDistRP[rp]->Draw("P");
    for (int cp=0; cp<NCOLS*fTelescopeCount;cp++)    if(ncp[cp]>0)   QC_vs_ElliDistCP[cp]->Draw("P");
    
    
    
    return imin;
}
    

void FdProcessElves::FdElvesBestEllipseSearch(int imin){

    TCanvas * cBestElli = new TCanvas("cBestElli","Corrected charge vs Elliptic Distance",200,10,1600,800);
    
    cBestElli->Divide(3,1);
    
    cBestElli->cd(1);          hQCCAveSpread->Draw("ZCOL");
    
    cBestElli->cd(2);          hQCC_Spread[imin]->Draw();
    
    
    stringstream ss3;
    ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_NPHCC_BestEllipse.png";
    
    cBestElli->SaveAs(ss3.str().c_str());

}


void FdProcessElves::FdSingleElvesPixel_CameraViewIni(){
    
    /*
      FdSingleElvesPixel_CameraViewIni
       Just prepares dummy pixels for two bays
     
    */
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    for (int npix=0; npix<2*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        Theta_pixel[nnpix] = XPixel.Theta();
        //double timepix = hFIT3Tmax->GetBinContent(npix+1);
        
        TMarker * mark = new TMarker(-azipixel, elepixel ,4); //empty dot
        mark->SetMarkerColor(2);    mark->SetMarkerSize(1.5);
        fPixels->Add(mark);
        
    }
    
    
    
    
    /*
     stringstream ss3;
     ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_TPEAK_FDVIEW"  << ".png";
     
     cTpeakPlot->SaveAs(ss3.str().c_str());
     
     delete cTpeakPlot;
     */
    
    
    
}

void FdProcessElves::FdSingleElvesPixel_LongLatViewIni(double Hd=85.){
    
    /*
     FdSingleElvesPixel_CameraViewIni
     Just prepares dummy pixels for two bays
     
     */
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    double Rsite = Rearth+hsite[site1];
    
    for (int npix=0; npix<fTelescopeCount*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ((backwall[site1]+15.+30.*(LeftMostBay-bay-1) )*degree);
        //    XPixel.RotateZ((LeftBayCenter-30.*bay)*degree);
        //XPixel.RotateZ((backwall[site1]+15.+30.*tel)*degree);
      
        
        East1 = PoloNord.Cross(Site1Location);
        
        double phi    = XPixel.Phi()/degree;
        double theta  = XPixel.Theta()/degree;

        TVector3 AziView1 = East1;
        AziView1.Rotate(phi*degree,Site1Location);
        TVector3 MaxCirc1 = Site1Location.Cross(AziView1);
        TVector3 PixelHdVert = Site1Location;
        double EFD = RootOf(1.,2.*Rsite*cos(theta*degree),Rsite*Rsite-pow(Rearth+Hd,2),1);
        double cta = (EFD*EFD-pow(Rearth+Hd,2)-Rsite*Rsite)/(2*Rsite*(Rearth+Hd));
        if (cta>=1 || cta <=-1)cta = 1;
        double alpha_EFD = acos(cta);
        // The vertical at FD site is rotated on the max circle to the vertical at light emission
        // and Latitude+Longitude are calculated
        PixelHdVert.Rotate(-alpha_EFD,MaxCirc1);
        if (PixelHdVert.Angle(Site1Location)>30.*degree) PixelHdVert = - PixelHdVert ;
        
        double LatPixelHdV = 90.-PixelHdVert.Theta()/degree;
        double LongPixelHdV = PixelHdVert.Phi()/degree;

        //PixelDirs2->SetPoint(npix,XPixel2.Phi()/degree,90.-XPixel2.Theta()/degree);
        double oo2[7], ph2[7];
        oo2[0] = oo+dOmega/2.; ph2[0] = ph-dPhi/3.;
        oo2[1] = oo+dOmega/2.; ph2[1] = ph+dPhi/3.;
        oo2[2] = oo;                     ph2[2] = ph+2.*dPhi/3.;
        oo2[3] = oo-dOmega/2.; ph2[3] = ph+dPhi/3.;
        oo2[4] = oo-dOmega/2.; ph2[4] = ph-dPhi/3.;
        oo2[5] = oo;                     ph2[5] = ph-2.*dPhi/3.;
        oo2[6] = oo2[0];               ph2[6] = ph2[0];
//        PixelEdges[npix] = new TGraph();
        PixelEdgesLL[npix] = new TGraph();
        
        for (int iV = 0; iV<7 ; iV++){
            
            double mcosOO2 = -cos(oo2[iV]);
            double z2 = mcosOO2 * cos(ph2[iV]);
            double x2 = mcosOO2 * sin(ph2[iV]);
            double y2 = sin(oo2[iV]);
            
            TVector3 *XPixEdge = new TVector3(-x2,-y2,-z2);
            XPixEdge->RotateY(73.5*degree);
            XPixEdge->RotateZ((backwall[site1]+15.+30.*(LeftMostBay-bay))*degree);
            
            PixelEdges[npix]->SetPoint(iV,XPixEdge->Phi()/degree,90.-XPixEdge->Theta()/degree);
            TVector3 AziView2 = East1;
            AziView2.Rotate(XPixEdge->Phi(),Site1Location);
            TVector3 MaxCirc2 = Site1Location.Cross(AziView2);
            TVector3 PixEdgeHdVert = Site1Location;
            // this is the distance between the FD and the point of light emission,
            // at Hd altitude (default 90 km)
            EFD = RootOf(1.,2.*Rsite*cos(XPixEdge->Theta()),Rsite*Rsite-pow(Rearth+Hd,2),1);
            cta = (EFD*EFD-pow(Rearth+Hd,2)-Rsite*Rsite)/(2*Rsite*(Rearth+Hd));
            if (cta>=1 || cta <=-1)cta = 1;
            alpha_EFD = acos(cta);
            
            PixEdgeHdVert.Rotate(-alpha_EFD,MaxCirc2);
            if (PixEdgeHdVert.Angle(Site1Location)>30.*degree) PixEdgeHdVert = - PixEdgeHdVert ;
            
            double LatPixEdgeHdV = 90.-PixEdgeHdVert.Theta()/degree;
            double LongPixEdgeHdV = PixEdgeHdVert.Phi()/degree;
            PixelEdgesLL[npix]->SetPoint(iV,LongPixEdgeHdV,LatPixEdgeHdV);
            
            cout << " : " << LongPixelHdV << " " << LatPixelHdV;

        }// end loop over Hexagon

        cout << endl;
        
    }
    
 
    
    
}



int FdProcessElves::FdAllTraces_CameraViewPlot(){
    /*
     FdAllTraces_CameraViewPlot() 
        input: hFIT3Qmax, filled in SingleElvesBoltFinder to normalize histos
               ThisGPSsec, EventNanoTime ti write the name of the PNG file
     
        output: AllTracesFDV, created to produce Canvas cFDAllTracesPlot
    
     */

    gStyle->SetOptStat(0);
    double Qtmax = 0.; int minRow=0;
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        double Qmaxpix = hFIT3Qmax->GetBinContent(npix+1);
        if (Qmaxpix>Qtmax){Qtmax=Qmaxpix; minRow=npix%NROWS+1;}
    }
    
    if (minRow == 0)return minRow;
    
    minRow -=4; if (minRow<1)minRow=1;
    
    if (minRow>13)minRow=13;

    // in case of two bay event, the column number goes from 1 to 40
    //int xc[3]={0,0,0};
    
    int telTmin = HitTel[0];
    
    //if (fL7TelescopeCount == 2 && HitTel[0] == HitTel[1]+1)xc[1]=1;
    //else if (fL7TelescopeCount == 2 && HitTel[1] == HitTel[0]+1){xc[0]=1; telTmin = HitTel[1];}
    
    
    stringstream ss2;
    
    if (fL7TelescopeCount ==1){
    ss2<< "All traces , Camera View, Eye: " << site1+1 << " Bay: " << HitTel[0] << " Event: " << ThisGPSsec << "." << EventNanoTime << " minRow = "  << minRow << " Nph,max = " << Qtmax ;
        ss2<< " ; Column; Row" ;}
    if (fL7TelescopeCount ==2){
        ss2<< "All traces , Camera View, Eye: " << site1+1 << " Bays: " << HitTel[0] << "," << HitTel[1] << " Event: " << ThisGPSsec << "." << EventNanoTime << " minRow = "  << minRow << " Nph,max = " << Qtmax ;
        ss2<< " ; Column; Row" ;}



    TH2F *AllTracesFDV = new TH2F("AllTracesFDV",ss2.str().c_str(),2250,0.,45.,10*nslots,0.,10*nslots);

    
    for (int itel=0;itel<fL7TelescopeCount;itel++){
    for (int ir=minRow;ir<minRow+10;ir++){
        for (int ic=1; ic<=NCOLS; ic++){
            for (int its=0; its<50*nPages; its++){
                int itp=(ir-1)+NROWS*(ic-1);
                float Qt = hCalPixels[itel][itp]->GetBinContent(its+1);
                int ixp = 1+ itel*50*NCOLS + 50*ic -int(Qt*75./Qtmax);
                if (ir%2 ==0)ixp -= 25;
                int iyp = 1+ its + (ir-minRow)*50*NMAXPAGES;
                float cr = (ic);
                if (ir%2 ==1)cr += 40.;
                AllTracesFDV->SetBinContent(ixp,iyp,cr);
        }
       }
      }
    }
    
    TCanvas *cFDAllTracesPlot = new TCanvas("cFDAllTracesPlot","FD View ", 3000, 2000);
    cFDAllTracesPlot ->cd();
    
   
    AllTracesFDV->Draw("zcol");
    
    stringstream ss3;
    ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_ALLTRACES_FDVIEW"  << ".png";
    
    cFDAllTracesPlot->SaveAs(ss3.str().c_str());
    
    delete cFDAllTracesPlot;
    delete AllTracesFDV;

    return minRow;

}

void FdProcessElves::FdSingleElvesTime_CameraViewPlot(int iv){
    const int NMAXCOLS = 60;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        double timepix = hFIT3Tmax->GetBinContent(npix+1);
        
        int MaxTime = nPages*100;
        if (MaxTime>400) MaxTime =400;
        
        if(timepix>0 && timepix<MaxTime) {//if(bay==1)maxphi = 50.;
            int timepixcolor = (int) ((timepix-28.)*NMAXCOLS/(MaxTime))+COLOROFFSET;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(timepixcolor);
            mark2->SetMarkerSize(1.5);
            
            fPixels2->Add(mark2);
            
        }
        /*if(npix>2*NPIXELS && npix < 2*NPIXELS+256){
            int timepixcolor = nnpix;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(timepixcolor);
            mark2->SetMarkerSize(1.5);
            fPixels2->Add(mark2);
        }*/
    }
    
    for (int i=0; i<50; i++) {
        int pixcolor = (int) (i*NMAXCOLS/50.)+COLOROFFSET;
        double azipixel = -10.+i*1.;
        double elepixel = 35.;
        TMarker * mark2 = new TMarker(azipixel, elepixel ,8);
        mark2->SetMarkerColor(pixcolor);
        mark2->SetMarkerSize(1.5);
        
        fPixels2->Add(mark2);
    }
    
    //TCanvas *cTpeakPlot = new TCanvas("cTpeakPlot","FD View Time Peak", 1400, 1000);
    cFDViewPlot->cd(iv);
    
    
    TH2D * TpeakPlot = new TH2D("Tpeak","Tpeak",90,minphi,maxphi,120,0.,40.);
    
    TpeakPlot->SetStats(kFALSE);
    TpeakPlot->GetXaxis()->SetTitle("-#phi [deg]");
    TpeakPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    TpeakPlot->GetZaxis()->SetTitle("T_{peak}[#mu s]");
    
//    for (int i=0; i<120 ; i++)TpeakPlot->Fill(79,i/3.+0.5,i*2.5);
    
    TpeakPlot->Draw();
    
    
    
    
    TIter next(fPixels);
    
    while( TMarker * mark = (TMarker *) next() ) {
        cFDViewPlot->cd(iv);
        mark->Draw("SAME");
    }
    
    
    TIter next2(fPixels2);
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cFDViewPlot->cd(iv);
        mark2->Draw("SAME");
    }
    
    
    /*
     stringstream ss3;
     ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_TPEAK_FDVIEW"  << ".png";
     
     cTpeakPlot->SaveAs(ss3.str().c_str());
     
     delete cTpeakPlot;
     */
    
    
    
}



void FdProcessElves::FdSingleElvesPull_CameraViewPlot(int iv ){
    const int NMAXCOLS = 50;
    const double PULLRANGE = 10.;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        double pullPpix = hFitPullsP->GetBinContent(npix+1); if (pullPpix>PULLRANGE)pullPpix=PULLRANGE;
        double pullNpix = -hFitPullsN->GetBinContent(npix+1); if (pullNpix<-PULLRANGE)pullNpix=-PULLRANGE;
        double pullpix = pullNpix; if (pullpix==0.)pullpix=pullPpix;
         
        if(pullpix!=0.) {//if(bay==1)maxphi = 50.;
            int pullpixcolor = (int) ((pullpix+PULLRANGE)*NMAXCOLS/(2.*PULLRANGE))+COLOROFFSET;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(pullpixcolor);
            mark2->SetMarkerSize(1.5);
            
            fPixels3->Add(mark2);
            
        }
        
        
    }

    
    for (int i=0; i<50; i++) {
        int pixcolor = (int) (i*NMAXCOLS/50.)+COLOROFFSET;
        double azipixel = -10.+i*1.;
        double elepixel = 35.;
        TMarker * mark2 = new TMarker(azipixel, elepixel ,8);
        mark2->SetMarkerColor(pixcolor);
        mark2->SetMarkerSize(1.5);
        
        fPixels3->Add(mark2);
    }

    
    //TCanvas *cPullPlot = new TCanvas("cPullPlot","FD View Pull", 1400, 1000);
    
    
    TH2D * PullPlot = new TH2D("PullPlot","Pulls ",90,minphi,maxphi,120,0.,40.);
    
    PullPlot->SetStats(kFALSE);
    PullPlot->GetXaxis()->SetTitle("-#phi [deg]");
    PullPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    PullPlot->GetZaxis()->SetTitle("Pull");
    
//    for (int i=0; i<120 ; i++)PullPlot->Fill(49,i/3.+0.5,i*2.5);
    
    //cFDViewPlot->cd(iv);
    PullPlot->Draw();
    
    
    delete PullPlot;
    
    TIter next(fPixels);
    
    while( TMarker * mark = (TMarker *) next() ) {
        //cFDViewPlot->cd(iv);
        mark->Draw("SAME");
    }
    
    
    TIter next2(fPixels3);
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        //cFDViewPlot->cd(iv);
        mark2->Draw("SAME");
    }
    
    
    
    
}


void FdProcessElves::FdSingleElvesSigmaT_CameraViewPlot(int iv ){
    const int NMAXCOLS = 50;
    const double SIGMALRANGE = 100.;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        double sigmaTpix = hFIT3SigmaT->GetBinContent(npix+1); if (sigmaTpix>SIGMALRANGE)sigmaTpix=SIGMALRANGE;
        
        if(sigmaTpix>0.) {//if(bay==1)maxphi = 50.;
            int sigmaTpixcolor = (int) (sigmaTpix)*NMAXCOLS/(SIGMALRANGE)+COLOROFFSET;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(sigmaTpixcolor);
            mark2->SetMarkerSize(1.5);
            
            fPixels4->Add(mark2);
            
        }
        
        
    }

    
    for (int i=0; i<50; i++) {
        int pixcolor = (int) (i*NMAXCOLS/50.)+COLOROFFSET;
        double azipixel = -10.+i*1.;
        double elepixel = 35.;
        TMarker * mark2 = new TMarker(azipixel, elepixel ,8);
        mark2->SetMarkerColor(pixcolor);
        mark2->SetMarkerSize(1.5);
        
        fPixels4->Add(mark2);
    }

    
    //TCanvas *cPullPlot = new TCanvas("cPullPlot","FD View Pull", 1400, 1000);
    
    
    TH2D * SigmaLPlot = new TH2D("SigmaLPlot","Pulse LEFT Sigma ",90,minphi,maxphi,120,0.,40.);
    
    SigmaLPlot->SetStats(kFALSE);
    SigmaLPlot->GetXaxis()->SetTitle("-#phi [deg]");
    SigmaLPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    SigmaLPlot->GetZaxis()->SetTitle("Pulse LEFT Sigma (#mu s)");
    
    //    for (int i=0; i<120 ; i++)PullPlot->Fill(49,i/3.+0.5,i*2.5);
    
    cFDViewPlot->cd(iv);    SigmaLPlot->Draw();
    
   // delete SigmaLPlot;
    
    
    TIter next(fPixels);
    
    while( TMarker * mark = (TMarker *) next() ) {
        cFDViewPlot->cd(iv);
        mark->Draw("SAME");
    }
    
    
    TIter next2(fPixels4);
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cFDViewPlot->cd(iv);
        mark2->Draw("SAME");
    }
    
    
    
    
}


void FdProcessElves::FdSingleElvesSigmaR_CameraViewPlot(int iv ){
    const int NMAXCOLS = 50;
    const double SIGMARRANGE = 200.;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        double sigmaRpix = hFIT3SigmaT->GetBinContent(npix+1)*hFIT3Skew->GetBinContent(npix+1); if (sigmaRpix>SIGMARRANGE)sigmaRpix=SIGMARRANGE;
        
        if(sigmaRpix>0.) {//if(bay==1)maxphi = 50.;
            int sigmaRpixcolor = (int) (sigmaRpix)*NMAXCOLS/(SIGMARRANGE)+COLOROFFSET;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(sigmaRpixcolor);
            mark2->SetMarkerSize(1.5);
            
            fPixels6->Add(mark2);
            
        }
        
        
    }
    
    
    
    for (int i=0; i<50; i++) {
        int pixcolor = (int) (i*NMAXCOLS/50.)+COLOROFFSET;
        double azipixel = -10.+i*1.;
        double elepixel = 35.;
        TMarker * mark2 = new TMarker(azipixel, elepixel ,8);
        mark2->SetMarkerColor(pixcolor);
        mark2->SetMarkerSize(1.5);
        
        fPixels6->Add(mark2);
    }


    
    //TCanvas *cPullPlot = new TCanvas("cPullPlot","FD View Pull", 1400, 1000);
    
    
    TH2D * SigmaRPlot = new TH2D("SigmaRPlot","Pulse RIGHT Sigma ",90,minphi,maxphi,120,0.,40.);
    
    SigmaRPlot->SetStats(kFALSE);
    SigmaRPlot->GetXaxis()->SetTitle("-#phi [deg]");
    SigmaRPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    SigmaRPlot->GetZaxis()->SetTitle("Pulse RIGHT Sigma (#mu s)");
    
    //    for (int i=0; i<120 ; i++)PullPlot->Fill(49,i/3.+0.5,i*2.5);
    
    cFDViewPlot->cd(iv);    SigmaRPlot->Draw();
    
    
    
    
    TIter next(fPixels);
    
    while( TMarker * mark = (TMarker *) next() ) {
        cFDViewPlot->cd(iv);
        mark->Draw("SAME");
    }
    
    
    TIter next2(fPixels6);
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cFDViewPlot->cd(iv);
        mark2->Draw("SAME");
    }
    
    
    
    
}



void FdProcessElves::FdSingleElvesQmax_CameraViewPlot(int iv ){
    
    cout << "Plotting QMAX Camera View .... " ;
    
    const int NMAXCOLS = 50;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    
    double QMAXRANGE = 0.;
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        double Qmaxpix = hFIT3Qmax->GetBinContent(npix+1);
        if (Qmaxpix>QMAXRANGE)QMAXRANGE=Qmaxpix;
    }
    
    for (int npix=0; npix<2*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        double Qmaxpix = hFIT3Qmax->GetBinContent(npix+1);
        
        if(Qmaxpix>0.) {//if(bay==1)maxphi = 50.;
            int Qmaxpixcolor = (int) (Qmaxpix)*NMAXCOLS/(QMAXRANGE)+COLOROFFSET;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(Qmaxpixcolor);
            mark2->SetMarkerSize(1.5);
            
            fPixels5->Add(mark2);
            
        }
        
        
    }
    
    
    for (int i=0; i<50; i++) {
        int pixcolor = (int) (i*NMAXCOLS/50.)+COLOROFFSET;
        double azipixel = -10.+i*1.;
        double elepixel = 35.;
        TMarker * mark2 = new TMarker(azipixel, elepixel ,8);
        mark2->SetMarkerColor(pixcolor);
        mark2->SetMarkerSize(1.5);
        
        fPixels5->Add(mark2);
    }

    
    //TCanvas *cPullPlot = new TCanvas("cPullPlot","FD View Pull", 1400, 1000);
    
    
    TH2D * QmaxPlot = new TH2D("QmaxPlot","Pulse Qmax(Photons at Diaphragm) ",90,minphi,maxphi,120,0.,40.);
    
    QmaxPlot->SetStats(kFALSE);
    QmaxPlot->GetXaxis()->SetTitle("-#phi [deg]");
    QmaxPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    QmaxPlot->GetZaxis()->SetTitle("Pulse Qmax(Photons at Diaphragm)");
    
    //    for (int i=0; i<120 ; i++)PullPlot->Fill(49,i/3.+0.5,i*2.5);
    
    cFDViewPlot2->cd(iv);    QmaxPlot->Draw();
    
    
    
    
    TIter next(fPixels);
    
    while( TMarker * mark = (TMarker *) next() ) {
        cFDViewPlot2->cd(iv);
        mark->Draw("SAME");
    }
    
    
    TIter next2(fPixels5);
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cFDViewPlot2->cd(iv);
        mark2->Draw("SAME");
    }
    
    
    cout << " .... Qmax completed " << endl;
    
}

void FdProcessElves::FdSingleElvesQCtot_CameraViewPlot(int iv ){
    
     cout << "Plotting QCtot Camera View .... " ;
    
    const int NMAXCOLS = 50;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    
    double QMAXRANGE = 0.;
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        double Qmaxpix = hFineQCtot->GetBinContent(npix+1);
        if (Qmaxpix>QMAXRANGE)QMAXRANGE=Qmaxpix;
    }
    
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        double Qmaxpix = hFineQCtot->GetBinContent(npix+1);
        
        if(Qmaxpix>0.) {//if(bay==1)maxphi = 50.;
            int Qmaxpixcolor = (int) (Qmaxpix*NMAXCOLS/QMAXRANGE)+COLOROFFSET;
            
            //cout << " QMAXPIX " << npix << " " << Qmaxpix << " " << Qmaxpixcolor << endl;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(Qmaxpixcolor);
            mark2->SetMarkerSize(1.5);
            
            fPixels7->Add(mark2);
            
        }
        
        

    }
    
    
    
    for (int i=0; i<50; i++) {
        int pixcolor = (int) (i*NMAXCOLS/50.)+COLOROFFSET;
        double azipixel = -10.+i*1.;
        double elepixel = 35.;
        TMarker * mark2 = new TMarker(azipixel, elepixel ,8);
        mark2->SetMarkerColor(pixcolor);
        mark2->SetMarkerSize(1.5);
        
        fPixels7->Add(mark2);
    }

    
    //TCanvas *cPullPlot = new TCanvas("cPullPlot","FD View Pull", 1400, 1000);
    
    
    TH2D * QCtotPlot = new TH2D("QCtotPlot","Pulse QCtot(Light Emission Surface Density in Ionosphere) ",90,minphi,maxphi,120,0.,40.);
    
    QCtotPlot->SetStats(kFALSE);
    QCtotPlot->GetXaxis()->SetTitle("-#phi [deg]");
    QCtotPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    QCtotPlot->GetZaxis()->SetTitle("Light Emission Surface Density(Photons per unit Area), NOT corrected for atmo");
    
    //    for (int i=0; i<120 ; i++)PullPlot->Fill(49,i/3.+0.5,i*2.5);
    
    cFDViewPlot2->cd(iv);    QCtotPlot->Draw();
    
    
    
    
    TIter next(fPixels);
    
    while( TMarker * mark = (TMarker *) next() ) {
        cFDViewPlot2->cd(iv);
        mark->Draw("SAME");
    }
    
    
    TIter next2(fPixels7);
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cFDViewPlot2->cd(iv);
        mark2->Draw("SAME");
    }
    
    
    cout << "...... QCtot  .... completed " << endl;
    
}


void FdProcessElves::FdSingleElvesQCCtot_CameraViewPlot(int iv ){
    
     cout << "Plotting QCCtot Camera View .... " ;
    
    const int NMAXCOLS = 50;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    double minphi=-20.;
    double maxphi=50.;
    
    double NphMax[4] = {};
    
    
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        double Nph = hFineQCCtot->GetBinContent(npix+1);
        for (int i=0; i<4; i++){
            if (Nph>NphMax[i]){
                for (int j=3;j>i;j--)NphMax[j]=NphMax[j-1];
                
                NphMax[i]=Nph;
                break;
            }
        }
    }
   double QMAXRANGE = NphMax[3];
  // cout << QMAXRANGE <<  " NPHMAX[4] " << NphMax[0] << " " << NphMax[1] << " "<< NphMax[2] << " "<< NphMax[3] << endl;
    
   // double QMAXRANGE = 6.e14;
    
    for (int npix=0; npix<2*NPIXELS; npix++){
        int nnpix = npix%NPIXELS;
        int bay=npix/NPIXELS;
        int row = nnpix%NROWS+1;
        int col = int(nnpix/NROWS)+1;
        const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
        const double ph = (11.66667 - row) * dPhi;
        
        const double mcosOO = -cos(oo);
        const double z = mcosOO * cos(ph);
        const double u = mcosOO * sin(ph);
        const double y = sin(oo);
        
        TVector3 XPixel(-u,-y,-z);
        XPixel.RotateY(73.5*degree);
        XPixel.RotateZ(-30.*bay*degree);
        
        
        double azipixel = XPixel.Phi()/degree;
        double elepixel = 90.-XPixel.Theta()/degree;
        double Qmaxpix = hFineQCCtot->GetBinContent(npix+1);
        
        if(Qmaxpix>0.) {//if(bay==1)maxphi = 50.;
            int Qmaxpixcolor = (int) (Qmaxpix*NMAXCOLS/QMAXRANGE)+COLOROFFSET;
            
            //cout << " QMAXPIX " << npix << " " << Qmaxpix << " " << Qmaxpixcolor << endl;
            TMarker * mark2 = new TMarker(-azipixel, elepixel ,8);
            mark2->SetMarkerColor(Qmaxpixcolor);
            mark2->SetMarkerSize(1.5);
            
            fPixels8->Add(mark2);
            
        }
       
        
        
    }
    
    
    for (int i=0; i<50; i++) {
        int pixcolor = (int) (i*NMAXCOLS/50.)+COLOROFFSET;
        double azipixel = -10.+i*1.;
        double elepixel = 35.;
        TMarker * mark2 = new TMarker(azipixel, elepixel ,8);
        mark2->SetMarkerColor(pixcolor);
        mark2->SetMarkerSize(1.5);
        
        fPixels8->Add(mark2);
    }

    
    
    TH2D * QCCtotPlot = new TH2D("QCCtotPlot","Pulse QCCtot(Light Emission Surface Density in Ionosphere) ",90,minphi,maxphi,120,0.,40.);
    
    QCCtotPlot->SetStats(kFALSE);
    QCCtotPlot->GetXaxis()->SetTitle("-#phi [deg]");
    QCCtotPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    QCCtotPlot->GetZaxis()->SetTitle("Light Emission Surface Density(Photons per unit Area), corrected for atmo");
    
    //    for (int i=0; i<120 ; i++)PullPlot->Fill(49,i/3.+0.5,i*2.5);
    
    cFDViewPlot2->cd(iv);    QCCtotPlot->Draw();
    
    
    
    
    TIter next(fPixels);
    
    while( TMarker * mark = (TMarker *) next() ) {
        cFDViewPlot2->cd(iv);
        mark->Draw("SAME");
    }
    
    
    TIter next2(fPixels8);
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cFDViewPlot2->cd(iv);
        mark2->Draw("SAME");
    }
    
    
    
    cout << "...... QCCtot  .... completed " << endl;

    
    
}



void FdProcessElves::FdSingleElvesQCCtot_LongLatViewPlot(TCanvas* cLLmap){
    
    cout << "Plotting QCCtot LongLat View .... " ;
    
    const int NMAXCOLS = 50;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    
    //stringstream ss3f;
    //ss3f << " Surface Density of light emission ; Long(degrees); Lat(degrees)";
    
    //    TH2F *PixelsLL  = new TH2F("PixelsLL",ss3f.str().c_str(),100,-72.122,-59.878,100,-40,-30);
    //TH2F *PixelsLL  = new TH2F("PixelsLL",ss3f.str().c_str(),100,-78.244,-54.756,100,-45,-25);
    cLLmap->SetGrid();
    PixelsLL->SetStats(kFALSE);
    PixelsLL->Draw("");
    TImage *img = TImage::Open("Argentina_map_LLsquare.gif");
    if (!img) {
        printf("Could not create an image... exit\n");
        return;
    }else cout << "Drawing Argentina map ..." << endl;
    
    img->SetConstRatio(kFALSE);
    //img->Draw("+");
    
    
    
  //  cLLmap->ls();
//    cLLmap->Update();

    
 //   gPad->Modified();
    
    
    
    
    
    
    
    for (int i=0; i<4; i++){
        
        //cout << " QMAXPIX " << npix << " " << Qmaxpix << " " << Qmaxpixcolor << endl;
        TMarker * mark2 = new TMarker(longsite[i], latsite[i],8);
        // cout << "SITEMARKER LONG,LAT " << longsite[i] << " " << latsite[i] << endl;
        mark2->SetMarkerColor(sitecolor[i]);
        mark2->SetMarkerSize(1.5);
        mark2->SetMarkerStyle(30);
        
        fSites->Add(mark2);
        
    }
    
    TMarker * mark2 = new TMarker(longsite[site1], latsite[site1],8);
    mark2->SetMarkerColor(sitecolor[site1]);
    mark2->SetMarkerSize(1.8);
    mark2->SetMarkerStyle(29);
    //cout << "THISSITEMARKER LONG,LAT " << longsite[site1] << " " << latsite[site1] << endl;
    
    
    fSites->Add(mark2);
    
    TIter next(fSites);
    cout << "Drawing FD sites ..." << endl;
    while( TMarker * mark2 = (TMarker *) next() ) {
        cLLmap->cd(2);
        mark2->Draw("SAME");
    }
    
    
    
    TIter next2(fSites2);
    cout << "Drawing bolt locations ..." << endl;
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cLLmap->cd(2);
        mark2->Draw("SAME");
    }
    
    
    double NphMax[4] = {};
    cout << " Finding Nphmax 1,2,3,4...."<< endl;
    
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        double Nph = hFineQCCtot->GetBinContent(npix+1);
        for (int i=0; i<4; i++){
            if (Nph>NphMax[i]){
                for (int j=3;j>i;j--)NphMax[j]=NphMax[j-1];
                
                NphMax[i]=Nph;
                break;
            }
        }
    }
    double QMAXRANGE = NphMax[3];
     //cout << QMAXRANGE <<  " NPHMAX[4] " << NphMax[0] << " " << NphMax[1] << " "<< NphMax[2] << " "<< NphMax[3] << endl;
    
    // double QMAXRANGE = 6.e14;
    
    for (int npix=0; npix<fTelescopeCount*NPIXELS; npix++){
             double Qmaxpix = hFineQCCtot->GetBinContent(npix+1);
        
        if(Qmaxpix>0.) {//if(bay==1)maxphi = 50.;
            int Qmaxpixcolor = (int) (Qmaxpix*NMAXCOLS/QMAXRANGE)+COLOROFFSET;
            
            if(PixelEdgesLL[npix]->GetN() > 0){
                PixelEdgesLL[npix]->SetFillColor(Qmaxpixcolor);
                
                PixelEdgesLL[npix]->Draw("F");
                
                PixelEdgesLL[npix]->Draw("L");
            }
            
        }else PixelEdgesLL[npix]->Draw("L");
        
        
        
    }
    if (fB100km->GetN() > 10) {
         fB100km->SetLineColor(2);
         fB100km->Draw("L");
         fB200km->SetLineColor(4);
         fB200km->Draw("L");
    }
    /*
     fBE1_100km->SetLineColor(6);
     fBE1_100km->Draw("L");
     fBE1_200km->SetLineColor(1);
     fBE1_200km->Draw("L");
     */
    /*
     stringstream ss4x;
     ss4x<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_NPHCC_LLVIEW"  << ".png";
     
     
     cLLmap->SaveAs(ss4x.str().c_str());
     */
    
    //delete img;
    
}

void FdProcessElves::FdSingleElvesQCCtot_LongLatViewPlot2(TCanvas *cLLmap){
    
    cout << "Plotting QCCtot LongLat View .... " ;
    
    const int NMAXCOLS = 50;
    const int COLOROFFSET = 50;
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;

    
  cLLmap->SetGrid();
   PixelsLL->SetStats(kFALSE);
    PixelsLL->Draw("");
    cLLmap->cd(2);
    
    TImage *img = TImage::Open("Argentina_map_LLsquare.gif");
    if (!img) {
        printf("Could not create an image... exit\n");
        return;
   }else cout << "Drawing Argentina map ..." << endl;
    
    img->SetConstRatio(kFALSE);
    //img->Draw("+");
    
    
    
    //gPad->Modified();
    
    
    
    
    for (int i=0; i<4; i++){
        
        //cout << " QMAXPIX " << npix << " " << Qmaxpix << " " << Qmaxpixcolor << endl;
        TMarker * mark2 = new TMarker(longsite[i], latsite[i],8);
       // cout << "SITEMARKER LONG,LAT " << longsite[i] << " " << latsite[i] << endl;
        mark2->SetMarkerColor(sitecolor[i]);
        mark2->SetMarkerSize(1.5);
        mark2->SetMarkerStyle(30);
        
        fSites->Add(mark2);
        
    }
    
    TMarker * mark2 = new TMarker(longsite[site1], latsite[site1],8);
    mark2->SetMarkerColor(sitecolor[site1]);
    mark2->SetMarkerSize(1.8);
    mark2->SetMarkerStyle(29);
    //cout << "THISSITEMARKER LONG,LAT " << longsite[site1] << " " << latsite[site1] << endl;
    
    
    fSites->Add(mark2);
    
    TIter next(fSites);
    cout << "Drawing FD sites ..." << endl;
    while( TMarker * mark2 = (TMarker *) next() ) {
        cLLmap->cd(2);
        mark2->Draw("SAME");
    }
   

    
    TIter next2(fSites2);
    cout << "Drawing bolt locations ..." << endl;
    
    while( TMarker * mark2 = (TMarker *) next2() ) {
        cLLmap->cd(2);
        mark2->Draw("SAME");
    }

    
    double NphMax[4] = {};
    cout << " Finding Nphmax E1,2,3,4...."<< endl;
    
    for (int npix=0; npix<fL7TelescopeCount*NPIXELS; npix++){
        double Nph = hFineQCCtot->GetBinContent(npix+1);
        for (int i=0; i<4; i++){
            if (Nph>NphMax[i]){
                for (int j=3;j>i;j--)NphMax[j]=NphMax[j-1];
                
                NphMax[i]=Nph;
                break;
            }
        }
    }
    double QMAXRANGE = NphMax[3];
    // cout << QMAXRANGE <<  " NPHMAX[4] " << NphMax[0] << " " << NphMax[1] << " "<< NphMax[2] << " "<< NphMax[3] << endl;
    
    // double QMAXRANGE = 6.e14;
    
    for (int npix=0; npix<fTelescopeCount*NPIXELS; npix++){

        double Qmaxpix = hFineQCCtot->GetBinContent(npix+1);
        
        if(Qmaxpix>0.) {//if(bay==1)maxphi = 50.;
            int Qmaxpixcolor = (int) (Qmaxpix*NMAXCOLS/QMAXRANGE)+COLOROFFSET;
            
            PixelEdgesLL[npix]->SetFillColor(Qmaxpixcolor);
            
            PixelEdgesLL[npix]->Draw("F");
 		         
            PixelEdgesLL[npix]->Draw("L");
            
        }else PixelEdgesLL[npix]->Draw("L");
        
        
        
    }
/*
    fB100km->SetLineColor(2);
    fB100km->Draw("L");
    fB200km->SetLineColor(4);
    fB200km->Draw("L");
 */
    if (fBE1_100km->GetN() > 10 &&  fBE1_200km->GetN() > 10 ) {

     fBE1_100km->SetLineColor(6);
     fBE1_100km->Draw("L");
     fBE1_200km->SetLineColor(1);
     fBE1_200km->Draw("L");
    }
    
    
    //cLLmap->ls();
   /*
    stringstream ss4x;
    ss4x<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_NPHCC_LLVIEW"  << ".png";
    
    
    cLLmap->SaveAs(ss4x.str().c_str());
 */
    //delete img;
    
}

void FdProcessElves::GlueTraces(TEyeEvent& eyeEvent, fevt::Eye& eye){

    TEyePixelList* eyePixelList = eyeEvent.GetPixelList();
    TEyeFADCData* eyeFADCData = eyeEvent.GetFADCData();
    const unsigned int numpixels = eyePixelList->GetNumPixels();
    TEyeEventHeader *eyeHeader = eyeEvent.GetEventHeader();
    
    const int eyeId = eyeEvent.GetEventHeader()->GetEyeNo();
    

    
    cout << "GlueTraces: NUM PIXELS " << numpixels << " Pedestal calculation , GPS =" << ThisGPSsec << endl;

    
    
    // INSERZIONE ROUTINES
    
    int countVirtuals=0;
    int countErrors=0;
    int numTels = 0 ;
    
    
    // calculate pedestals from raw data and fill raw traces
    for (unsigned int iPixel = 0; iPixel < numpixels; ++iPixel) {
        
        FdUtil::Fd::FdPixelNumber pixelNumber = eyePixelList->GetPixel(iPixel);
        
        const unsigned int mirrorId = FdUtil::Fd::GetEyeMirrorNo(pixelNumber);
        const unsigned int mirrorChannelId = FdUtil::Fd::GetMirrorPixelNo(pixelNumber);
        
        // reject all pixels that belong to mirrors without ELVE (7) or ELVE FOLLOWER (8) trigger
        ///cout << " WRONG EVENT LABEL ? " << endl;
        if (WrongEventLabel(eyeHeader->GetMirrorEventLabel(mirrorId)))continue;

        if (!eye.HasTelescope(mirrorId))eye.MakeTelescope(mirrorId);
        
        //cout << " TelHit[" << nPages-1 << "," << mirrorId-1 << "]" << endl;
        
        if (!TelHit[nPages-1][mirrorId-1]){TelHit[nPages-1][mirrorId-1]=true;
            //cout << " TelHit TRUE for mirror "<< mirrorId << " for PAGE " << nPages << endl;

            if (nPages==1)HitTel[numTels]=mirrorId;
            //cout << " HitTel[ " << numTels << "]=" << mirrorId << endl;

            numTels++;}
        
        fevt::Telescope& tel = eye.GetTelescope(mirrorId);
        const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
        const fdet::Eye& detEye = detFD.GetEye(eyeId);
        //cout << " detEye " << eyeId<<endl;
        const fdet::Telescope& detTel = detEye.GetTelescope(mirrorId);
       
        const fdet::Channel& thisChannel = detTel.GetChannel(mirrorChannelId);
        
        // virtual channels must be treated differently
        if (thisChannel.IsVirtual()) {
            ++countVirtuals;
            continue;
        }
        
        
        // Correct for channel/pixel mismatch
        
        
        
        // const fdet::Channel& detChannel = detTel.GetChannel(mirrorChannelId);
        unsigned int pixelId = thisChannel.GetPixelId();
        if (mirrorChannelId != pixelId) {
            ostringstream info;
            info << "remapped pixel during calibration: channel=" << mirrorChannelId
            << " -> pixel=" << pixelId;
            INFO(info);
        }
        
        if (!tel.HasPixel(pixelId))
            tel.MakePixel(pixelId);
        fevt::Pixel& pixel = tel.GetPixel(pixelId);
        
 

        
        const fdet::Pixel& detPixel = detFD.GetPixel(pixel);
        
        const int telId = detPixel.GetTelescopeId();
  //      const int eyeId = detPixel.GetEyeId();
        int pixrow = detPixel.GetRow();
        int pixcol = detPixel.GetColumn();
       
        
        //const int pixelId = detPixel.GetId();
        
        // Get the calibration constant for this channel
        const int channelId = detPixel.GetChannelId();
        const fdet::Pixel& channel = detFD.GetEye(eyeId).GetTelescope(telId).GetPixel(channelId);
        //const TabulatedFunction& calibration = channel.GetEndToEndCalibration();
        
        
        const double bestCalibConst = channel.GetEndToEndCalibrationAtReferenceWavelength();
        
        const TFADCData* fadcData = eyeFADCData->GetFADCData(pixelNumber);
        
        // Protect yourself against the stupid non-existent traces...
        if (!fadcData) {
            ++countErrors;
            continue;
        }
        
        
        TFADCData::FADCDataWord fadcword = fadcData->GetFADCTrace();
        
        //const double fadcBinSize = detCamera.GetFADCBinSize();
        //const unsigned int fadcBinLength = (unsigned int)(detCamera.GetFADCTraceLength());
        
        //  copy the Raw trace into an utl::Trace
        
        //utl::TraceI rawFADCTrace(fadcBinLength, fadcBinSize);
        
        const unsigned int startBin = fadcData->GetTraceStartBin();
        const unsigned int endBin = fadcData->GetTraceEndBin();
        
        double baseline = 0.;
        double baselineRMS = 0.;
        double charge=0.;
        
        
        // kpos is the starting value of the rebinned data: is 0, 50, or 100
        int kpos = size/nSlotInterval*(nPages-1);
        int kpos0 = kpos;
        
        //cout << " Starting RAWPED_calc from " << kpos << endl;
        
        // PEDESTAL CALCULATION : the first 280 bins
        for (unsigned int pos = startBin; pos <startBin+280; ++pos) baseline += int(FADCDataWordGetData(&fadcword[pos]))/280.;
        
        for (unsigned int pos = startBin; pos <startBin+280; ++pos) baselineRMS += pow(int(FADCDataWordGetData(&fadcword[pos])),2)/280.;
        
        baselineRMS -= pow(baseline,2);
        
        if (baselineRMS>0) baselineRMS=sqrt(baselineRMS);
        else baselineRMS=0.;
        
        
        int delay = TelGPSns[mirrorId-1]/100;    // This should be ZERO for the first telescope
        int skipNbins = int(delay/nSlotInterval);
        int i0 = 0;
        
        
        int telIndex = (LeftMostBay - mirrorId)%6;
        
        //cout << "RAWPED " << rawPed[telIndex][pixelId-1] << " for " << pixelId-1 <<endl;
        
        //cout << "STEP 1.." << endl;
        if (rawPed[telIndex][pixelId-1] == 0. && baseline >0)skipNbins+=(nPages-1)*size/nSlotInterval ;
        else if (nPages == 3 && rawPed2[telIndex][pixelId-1] == 0. && baseline >0) {skipNbins+=(nPages-2)*size/nSlotInterval ; i0=size/nSlotInterval; }
        
        //cout << "STEP 2..." << endl;

        // the rawPed is set to baseline ... but ...
        if (rawPed[telIndex][pixelId-1] == 0. )rawPed[telIndex][pixelId-1] = baseline ;
        if (nPages == 2 ) rawPed2[telIndex][pixelId-1]  = baseline ;

        //cout << "STEP 3...." << endl;
        

        int partial_bin = delay%nSlotInterval;  // this regards the last bin, which is only partially filled for the side bay
        double Qpartial=0;
        // The first bins of the second rawTrace, if page=1, are filled with pedestal
        if (nPages ==1 && skipNbins>0) {
            for(int ii=i0;ii<skipNbins;ii++){//rawTrace[telIndex][pixelId-1][ii]=nSlotInterval*baseline;
                hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*ii,nSlotInterval*baseline );   }
            if(partial_bin>0){//rawTrace[telIndex][pixelId-1][skipNbins]=baseline*partial_bin;
                hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*skipNbins,baseline*partial_bin );
            }
        }
        
        //cout << "STEP 4....." << endl;

        if (nPages ==2 && skipNbins>size/nSlotInterval) {
            for(int ii=i0;ii<skipNbins;ii++){//rawTrace[telIndex][pixelId-1][ii]=nSlotInterval*baseline;
                hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*ii,nSlotInterval*baseline );   }
            if(partial_bin>0){//rawTrace[telIndex][pixelId-1][skipNbins]=baseline*partial_bin;
                hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*skipNbins,baseline*partial_bin );
            }
        }
        
        //cout << "STEP 5....." << endl;

        if (nPages ==3 && skipNbins>2*size/nSlotInterval) {
            for(int ii=i0;ii<skipNbins;ii++){//rawTrace[telIndex][pixelId-1][ii]=nSlotInterval*baseline;
                hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*ii,nSlotInterval*baseline );   }
            if(partial_bin>0){//rawTrace[telIndex][pixelId-1][skipNbins]=baseline*partial_bin;
                hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*skipNbins,baseline*partial_bin );
            }
        }
        
       // cout << " DELAY CHECK " << delay << " skipbins " << skipNbins << " numTels "  <<  numTels << " endbin " << endBin << " nPages " << nPages << endl;
        
        // This is the main LOOP for filling rawTrace[TEL][PIX][BIN]
        //cout << "Main LOOP for filling RAW TRACE  " << endl;
        for (unsigned int pos = startBin; pos <= endBin; ++pos) {
            //if (pos%100 == 0)cout << pos << endl;
            // here we get the DATA at address pos;
            if (pos>=startBin)charge += FADCDataWordGetData(&fadcword[pos]);
            
            
            //  this is the core of the rebinning : we save integrated charge only every nSlotInterval=20 bins
            if (pos == endBin && partial_bin>0 ) {
                
                int kposx = int((delay+pos+size*(nPages-1))/nSlotInterval);
                ////cout << " PARTIAL_BINNING " << partial_bin << " KPOSX2b " << kposx << " Q= " << charge << endl;
                
                if (pixelId-1 < 440 && kposx < nslots)hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*kposx,charge );
                //if (pixelId-1 < 440 && kposx < nslots)rawTrace[telIndex][pixelId-1][kposx] = charge;
                charge=0.;
                
            }
            
            if( (pos+delay)%nSlotInterval < nSlotInterval-1 ) continue;
            if( (pos+delay)%nSlotInterval == nSlotInterval-1 ) {
                
                int kposx = int((delay+pos+size*(nPages-1))/nSlotInterval);
                // This is the main charge filling
                //if (partial_bin>0) cout << " KPOSX2 " << kposx << " POS= " << pos << " Q= " << charge << endl;
                //else cout << " KPOSX1 " << kposx << " POS= " << pos << " Q= " << charge << endl;
                
                //cout << "NUMtels " << numTels ;
                
                if (pixelId-1 < 440 && kposx < nslots)hRawPixels[telIndex][pixelId-1]->Fill(0.5+1.*kposx,charge );
                //cout << " pixelId " << pixelId << endl;
                charge=0.;
            }
        }// END of MAIN LOOP for FILLING rawtrace
        
      
        kpos = size/nSlotInterval*(nPages-1);

        //cout << "Main LOOP for filling CAL TRACE  " << kpos << " " << nPages << endl;

        for (unsigned int ppos = 0; ppos < size/nSlotInterval; ++ppos) {
            
            //double Qraw = rawTrace[telIndex][pixelId-1][ppos+kpos];
            double Qraw = hRawPixels[telIndex][pixelId-1]->GetBinContent(ppos+kpos+1);

            
            float Qpos = bestCalibConst*( Qraw-nSlotInterval*rawPed[telIndex][pixelId-1] );
            
            float fudge = 2.5;
            float dQpos = fudge*sqrt(bestCalibConst* Qraw);
            
            

            if (bestCalibConst>900.)Qpos=0;
            //if (ppos%20 == 0)cout << nPages << "-" << ppos+kpos  << " " << bestCalibConst*Qraw << " " << telIndex << " " << pixelId-1<< " " << Qpos << " " << dQpos << endl;
           
            //cout << " Filling CalPixels " ;
            hCalPixels[telIndex][pixelId-1]->SetBinContent(ppos+kpos+1,Qpos);
            hCalPixels[telIndex][pixelId-1]->SetBinError(ppos+kpos+1,dQpos);
            //cout <<  ppos+kpos+1 << endl;
            
           
             //cout << " Filling CalTrace " ;
            calTrace[telIndex][pixelId-1][ppos+kpos]= Qpos;
            
            //cout <<  ppos+kpos << endl;
            QtotPixel[telIndex][pixelId-1]+= Qpos;
            
            
        }
    }

}


bool FdProcessElves::WrongEventLabel(int EventLabel){
    
    if (EventLabel<7 || EventLabel>8) return true;
    else return false;
    
}

void FdProcessElves::FdElvesPulseFinder(int it, int ip){
   // cout << "Entering FdElvesPulseFinder Telescope=" << it << "Pixel=" << ip << endl;

    bool Found1pulse = FindPulse(it,  ip);
    
    if (Found1pulse )ParametrizePixelPulse( it, ip);

    int dTpulses = FindSecondPulse(it, ip);
    
    bool Found2pulses = false;
    if (abs(dTpulses)>0)Found2pulses = true;
    
    bool IsSecondHigh = false;
    if (dTpulses<0)IsSecondHigh = true;
    
  //  if (Found1pulse || Found2pulses )ParametrizePixelPulse( it, ip);
  //  if (Found2pulses )ParametrizePixelPulse2( it, ip);
    
    if (Found2pulses && !IsSecondHigh)ParametrizePixelPulse2( it, ip, 0);
    if (Found2pulses && IsSecondHigh)ParametrizePixelPulse2( it, ip, 1);
    
    if (Found2pulses)NpixCount2++;
    else if (Found1pulse)NpixCount1++;

}


bool FdProcessElves::FindPulse(int it, int ip){
/* 
 *
 *  FindPulse searches for the pulse maximum across all the trace and save the
 *  timeslot when it occurs. It can be mis guided by double elves of type "2nd high"
 *
 *  it: telescope index: 0 in single bay events, 0,1 in double bay events
 *                   0,1,2 is for triple bay events , very rare
 *  ip: pixel index : goes from 0 to NPIX-1
 *
 *  Qpeak and Tpeak are stored in 1D histos (hRawQmax and hRawTmax)
 *
 *  INPUT hCalPixels[tel][pix]
 *  OUTPUT hRawQmax, hRawTmax
 */
    bool found = false;
    float Qpeak = QMIN; // below a certain value (5 sigmas?) we cannot call it a pulse
    float Tpeak = 0.;
    int itp = it*NPIXELS+ip;

    //cout << " hCalPixels[it][ip] has " << nPages*size/nSlotInterval << "slots" <<endl;
    for (unsigned int ppos = 0; ppos < nPages*size/nSlotInterval; ++ppos) {
        //cout << "PPOS=" << ppos << endl;
        double Qpos=hCalPixels[it][ip]->GetBinContent(ppos+1);
        //cout << itp << "-QPOS=" << Qpos << endl;
        
        if (Qpos>Qpeak){
            Qpeak=Qpos;
            Tpeak=(1.*ppos+0.5)*DTSLOT;
        }
    }
    //cout << " Qpeak " << Qpeak << " Tpeak " << Tpeak << endl;
    if (Qpeak>QMIN) {
        //cout << "Filling Raw Qmax and Tmax " << itp << endl;
        hRawQmax->SetBinContent(itp+1,Qpeak);
        hRawTmax->SetBinContent(itp+1,Tpeak);
        found = true;
    }
    //cout << it  << " " <<  ip  << " " << Qpeak << endl;
    return found;
}


int FdProcessElves::FindSecondPulse(int it, int ip){
    /*
     *
     *  FindPulse searches for the pulse maximum across all the trace and save the
     *  timeslot when it occurs. It can be mis guided by double elves of type "2nd high"
     *
     *  it: telescope index: 0 in single bay events, 0,1 in double bay events
     *                   0,1,2 is for triple bay events , very rare
     *  ip: pixel index : goes from 0 to NPIX-1
     *
     *  Qpeak and Tpeak are stored in 1D histos (hRawQmax and hRawTmax)
     *
     *  INPUT hCalPixels[tel][pix]
     *  OUTPUT hRawQmax, hRawTmax
     */
    bool found = false; int dTpulses = 0;
    
    float Qpeak = QMIN; // below a certain value (5 sigmas?) we cannot call it a pulse
    float Tpeak = 0.;
    
    int itp = it*NPIXELS+ip;
    float Tpeak0=hRawTmax->GetBinContent(itp+1);
    float Qpeak0=hRawQmax->GetBinContent(itp+1);

    
   // cout << " hCalResiduals[it][ip] has " << nPages*size/nSlotInterval << "slots" <<endl;
    for (unsigned int ppos = 0; ppos < nPages*size/nSlotInterval; ++ppos) {
        //cout << "PPOS=" << ppos ;

        double Qpos=hCalResiduals[it][ip]->GetBinContent(ppos+1);
        //cout << " Pixel " << itp << " QPOS=" << Qpos << endl;
        
        if (Qpos>Qpeak){
            Qpeak=Qpos;
            Tpeak=(1.*ppos+0.5)*DTSLOT;
        }
    }
    //cout << " 2ndPulse: Qpeak " << Qpeak << " " << Qpeak0<< " Tpeak " << Tpeak<< " " << Tpeak0<< endl;
    if (Qpeak>QMIN && Qpeak0>QMIN) {
        cout << ThisGPSsec << " DOUBLE ELVES: Filling Raw Qmax and Tmax " << itp;
        if (Tpeak>Tpeak0+DTMIN){
            cout << " T1,Q1= " << Tpeak0 << " " << Qpeak0;
            cout << " dT,Q2= " << Tpeak-Tpeak0 << " " << Qpeak << " CaseA" << endl;
            
            hRawQmax2->SetBinContent(itp+1,Qpeak);
            hRawTmax2->SetBinContent(itp+1,Tpeak);
            
            dTpulses = Tpeak-Tpeak0;
            found = true;
        } else if (Tpeak0>Tpeak+DTMIN) {
            cout << " T1,Q1= " << Tpeak << " " << Qpeak;
            cout << " dT,Q2= " << Tpeak0-Tpeak << " " << Qpeak0 << " CaseB" << endl;

            hRawQmax2->SetBinContent(itp+1,Qpeak0);
            hRawTmax2->SetBinContent(itp+1,Tpeak0);
            hRawQmax->SetBinContent(itp+1,Qpeak);
            hRawTmax->SetBinContent(itp+1,Tpeak);
            
            dTpulses = Tpeak-Tpeak0;
            found = true;
            
        }
    }
    //cout << it  << " " <<  ip  << " " << Qpeak << endl;
    return dTpulses;
}




void FdProcessElves::ParametrizePixelPulse(int it,int ip){
    /*
     *  ParametrizePixelPulse fits the pulse with an asymmetric gaussian
     *
     *  it: telescope index: 0 in single bay events, 0,1 in double bay events
     *                   0,1,2 is for triple bay events , very rare
     *  ip: pixel index : goes from 0 to NPIX-1
     *
     *  Qpeak and Tpeak are stored in 1D histos (hRawQmax and hRawTmax)
     *  ToDo: Undershoot correction for large pulses
     *
     *  INPUT: hRawQmax,hRawTmax,hCalPixels[it][ip]
     *  INTERMEDIATE: htrace
     *  OUTPUT: hpixels[itp],hresiduals[itp],hpixelsUC[itp]
     *  * hpixels is FITTED , then copied to hCalPixels[it][ip]
     *  * hresiduals is ready for fitting the DOUBLE ELVE, and also copied to hCalResiduals[it][ip]
     *  Fit results and errors are stored in :
     *  OUTPUT2:        hFineTmax,hFineQmax,hFineSigmaT,hFineSkew
     *
     */
    
    int row = ip%NROWS+1;
    int col = int(ip/NROWS)+1;
    int itp = it*NPIXELS+ip;
    
    Qmx[itp] = hRawQmax->GetBinContent(itp+1);
    Tmx[itp] = hRawTmax->GetBinContent(itp+1);
    // TimePeak_vs_Col[row-1]->SetPoint(col-1,mx,Tmx[ip]*dTbin);
    //if(row%2==0)TimePeak_vs_Col[row-1]->SetPoint(col-1,col,Tmx[ip]*dTbin);
    //if(row%2==1)TimePeak_vs_Col[row-1]->SetPoint(col-1,col+0.5,Tmx[ip]*dTbin);
    int ntotbins = nPages*size/nSlotInterval;
    fitFcn = new TF1("fitFcn",asymGaussian,0,ntotbins*DTSLOT,4);
    //fitFcn = new TF1("fitFcn",asymGaussianU,0,150*dTbin,5);
    
    //cout << " PIXEL: " << ip << " Qmax: " << Qmx[itp] << " at T " << Tmx[itp] << endl;
    fitFcn -> SetParameters(Qmx[itp],Tmx[itp],4*DTSLOT,1.);
    
    if (Qmx[itp]<QMIN)return;
    
    htrace->Reset();hpull->Reset();
    for (int i=0; i<nPages*size/nSlotInterval; i++){
        double Qcal = hCalPixels[it][ip]->GetBinContent(i+1);
        double dQcal = hCalPixels[it][ip]->GetBinError(i+1);
        
        htrace->SetBinContent(i+1,Qcal);
        htrace->SetBinError(i+1,dQcal);
        // htrace->SetBinError(i+1,sqrt(20./200.));
    }
    TH1F *hpixels= (TH1F *)htrace->Clone();
    TH1F *hresiduals= (TH1F *)htrace->Clone();
    TH1F *hpixelsUC= (TH1F *)htrace->Clone();
    
    stringstream shi;
    //ss2<<inputFile1<<" dt:"<<it*2;
    shi<<"htrax" << it << "_Tel" << HitTel[it] << "_R"<<row<<"_C"<<col;
    string tith = shi.str();
    hpixels->SetName(tith.c_str());
    
    stringstream shiUC;
    //ss2<<inputFile1<<" dt:"<<it*2;
    shiUC<<"htraxUC" << it << "_Tel" << HitTel[it] <<   "_R"<<row<<"_C"<<col;
    string tithUC = shiUC.str();
    hpixelsUC->SetName(tithUC.c_str());
    
    stringstream shi2;
    //ss2<<inputFile1<<" dt:"<<it*2;
    shi2<<"hres" << it << "_Tel" << HitTel[it] << "_R"<<row<<"_C"<<col;
    string tith2 = shi2.str();
    hresiduals->SetName(tith2.c_str());
    
    //CorrectUndershoot(ip);
    cout << " Fitting Pixel Pulse " << itp << " Tel" << HitTel[it] << " R"<<row<<" C"<<col<< endl;
    hpixels->Fit("fitFcn","R");
    //hCalPixels[it][ip]->Fit("fitFcn","R");
    float peakq = fitFcn->GetParameter(0) ;  float sigmapeakq= fitFcn->GetParError(0);
    float peakt = fitFcn->GetParameter(1) ;  float sigmapeakt= fitFcn->GetParError(1);
    float sigmaL  = fitFcn->GetParameter(2) ;  float sigmasigmaL = fitFcn->GetParError(2);
    float skew    = fitFcn->GetParameter(3)  ;  float sigmaskew = fitFcn->GetParError(3);
    float TmaxInt = 1500.;
    if (peakt+2*sigmaL*skew<TmaxInt ) TmaxInt = peakt+2*sigmaL*skew;
    float Qtot2s  = fitFcn->Integral(0.,TmaxInt);   // Qtot = int(Q) from 0 to Tpeak+2sigma
    
    float d1 = sigmapeakq/peakq;
    float d2 = sigmasigmaL/sigmaL;
    float d3 = sigmaskew/(1+skew);
    float dQtot2s  = Qtot2s*sqrt(pow(d1,2)+pow(d2,2)+pow(d3,2));
    
    
    if(peakt>25.&& peakt <nPages*size/nSlotInterval*DTSLOT) {
        hFineTmax->SetBinContent(itp+1,peakt);
        hFineTmax->SetBinError(itp+1,sigmapeakt);
        hFineQmax->SetBinContent(itp+1,peakq);
        hFineQmax->SetBinError(itp+1,sigmapeakq);
        hFineSigmaT->SetBinContent(itp+1,sigmaL);
        hFineSigmaT->SetBinError(itp+1,sigmasigmaL);
        hFineSkew->SetBinContent(itp+1,skew);
        hFineSkew->SetBinError(itp+1,sigmaskew);
        hFineQtot->SetBinContent(itp+1,Qtot2s);
        hFineQtot->SetBinError(itp+1,dQtot2s);
        
    }
    //hCalPixels[it][ip]=(TH1F *)hpixels->Clone();
    
    for (int i=1; i<=nPages*size/nSlotInterval; i++){
        double res = hpixels->GetBinContent(i) - fitFcn->Eval(hpixels->GetBinCenter(i));
        //double res =    hCalPixels[it][ip]->GetBinContent(i) - fitFcn->Eval(hCalPixels[it][ip]->GetBinCenter(i));
        //hresiduals
        hCalResiduals[it][ip]->SetBinContent(i,res);
        double sigma_guess = sqrt(20./200.);
        double pul=res/sigma_guess ;
        hpull->Fill(pul);
    }
    
    //hCalResiduals[it][ip]=(TH1F *)hresiduals->Clone();
    
    delete hpixels ;
    delete  hresiduals ;
    /*
     
     hpulls[itp]=(TH1F *)hpull->Clone();
     stringstream shi3;
     //ss2<<inputFile1<<" dt:"<<it*2;
     shi3<<"hpull_"<<row<<"_"<<col;
     string tith3 = shi3.str();
     hpulls[itp]->SetName(tith3.c_str());
     
     */
    
    /*
     
     PeakAmpli->SetBinContent(col,row,abs(fitFcn->GetParameter(0)));
     
     float peakt = fitFcn->GetParameter(1) ;  Sigma_Time[ip]= fitFcn->GetParError(1);
     if(peakt>25.&& peakt <300.) PeakTime->SetBinContent(col,row,peakt );
     
     float pulsewid = abs(fitFcn->GetParameter(2));
     if(pulsewid<100.)PulseSigma->SetBinContent(col,row,pulsewid);
     
     float skewp = abs(fitFcn->GetParameter(3));
     if (skewp<10.)PulseSkew->SetBinContent(col,row,skewp);
     
     //float USratio=abs(fitFcn->GetParameter(4)/fitFcn->GetParameter(0));
     float USratio=0.05;   // FAKE: just to give it a value
     if (USratio<0.1)PulseUS->SetBinContent(col,row,USratio);
     
     // cout << " FIT** MAX Value : " << fitFcn->GetParameter(0) << " at T " << fitFcn->GetParameter(1) ;
     // cout << " sigma_left " << abs(fitFcn->GetParameter(2)) <<  " sigma_right " ;
     // cout <<  abs(fitFcn->GetParameter(2)*fitFcn->GetParameter(3)) << " INTEGRAL " ;
     // cout << fitFcn->Integral(0.,300.) << endl;
     
     QtotC[ip]=fitFcn->Integral(0.,300.);
     double QtotC400 = fitFcn->Integral(0.,400.);
     double QtotC500 = fitFcn->Integral(0.,500.);
     int ipcompl = (row-1)+NROWS*(18-col-1);
     if (row%2==0) ipcompl = (row-1)+NROWS*(19-col-1);
     cout << " Extra_" << row << "_" << col << " " << QtotC400/QtotC[ip] << " " << QtotC500/QtotC[ip] << endl;
     if (col>9 && col <18 && row%2>0) cout << " Asymmetry_row" << row << "_col_" << col << "-" << 18-col << " " << QtotC[ip]/QtotC[ipcompl] << endl;
     if (col>9 && col <19 && row%2==0) cout << " Asymmetry_row" << row << "_col_" << col << "-" << 19-col << " " << QtotC[ip]/QtotC[ipcompl] << endl;
     
     
     double Q=Qtot[ip];
     double Qrel=(Q/Qttmax);
     double QdtThres = 0.2*Qtmx[ip];
     // search for the pulse start (constant fraction = 20%)
     for (int it=Tmx[ip]-20; it<Tmx[ip]-1;it++){
     if (it>10 && it<145 && !startfound[ip]){
     double thisQdt = Qdt2b[2*NPIXELS*(it-1)+ip];
     double thatQdt = Qdt2b[2*NPIXELS*(it-2)+ip];
     if (thisQdt>QdtThres && thatQdt<QdtThres) {timthr1 = (it-1)*dTbin+dTbin*(thisQdt-QdtThres)/(thisQdt-thatQdt);
     cout << "PIXEL " << ip << " MAX " << Qtmx[ip] << " START " << timthr1 << endl;
     startfound[ip]=true;
     Nstartfound[row-1]++;
     TimeStart_vs_Col[row-1]->SetPoint(Nstartfound[row-1],mx,timthr1);
     }
     }}
     // search for the pulse stop (zero crossing)
     // the second threshold is for the falling edge: (again 20% of Qpeak)
     QdtThres = 0.2*Qtmx[ip];
     
     for (int it=Tmx[ip]; it<Tmx[ip]+40;it++){
     if (it>11 && it<148 && !stopfound[ip]){
     double thisQdt = Qdt2b[2*NPIXELS*(it-1)+ip];
     double thatQdt = Qdt2b[2*NPIXELS*(it-2)+ip];
     if (thisQdt<QdtThres && thatQdt>QdtThres) {timthr2 = (it-1)*dTbin+dTbin*(thisQdt-QdtThres)/(thisQdt-thatQdt);
     cout << "PIXEL " << ip << " MAX " << Qtmx[ip] << " STOP" << timthr2 << endl;
     stopfound[ip]=true;
     Nstopfound[row-1]++;
     TimeStop_vs_Col[row-1]->SetPoint(Nstopfound[row-1],mx,timthr2);
     }
     
     }
     }
     */
}
// END of PULSE search

void FdProcessElves::ParametrizePixelPulse2(int it,int ip, int ic ){
    /*
     *  ParametrizePixelPulse2 fits the second pulse with an asymmetric gaussian
     *
     *  it: telescope index: 0 in single bay events, 0,1 in double bay events
     *                   0,1,2 is for triple bay events , very rare
     *  ip: pixel index : goes from 0 to NPIX-1
     *
     *  Qpeak and Tpeak are stored in 1D histos (hRawQmax and hRawTmax)
     *  ToDo: Undershoot correction for large pulses
     *
     *  INPUT: hRawQmax2,hRawTmax2,hCalResiduals[it][ip]
     *  INTERMEDIATE: htrace
     *  OUTPUT: hpixels[itp],hresiduals[itp],hpixelsUC[itp]
     *  * hpixels is FITTED , then copied to hCalPixels[it][ip]
     *  * hresiduals is ready for fitting the DOUBLE ELVE, and also copied to hCalResiduals[it][ip]
     *  Fit results and errors are stored in :
     *  OUTPUT2:        hFineTmax,hFineQmax,hFineSigmaT,hFineSkew
     *
     */
    
    int row = ip%NROWS+1;
    int col = int(ip/NROWS)+1;
    int itp = it*NPIXELS+ip;
    
    float Qmx2 = hRawQmax2->GetBinContent(itp+1);
    float Tmx2 = hRawTmax2->GetBinContent(itp+1);
    
    if (ic > 0) {
        Qmx2 = hRawQmax->GetBinContent(itp+1);
        Tmx2 = hRawTmax->GetBinContent(itp+1);
    }
    // TimePeak_vs_Col[row-1]->SetPoint(col-1,mx,Tmx[ip]*dTbin);
    //if(row%2==0)TimePeak_vs_Col[row-1]->SetPoint(col-1,col,Tmx[ip]*dTbin);
    //if(row%2==1)TimePeak_vs_Col[row-1]->SetPoint(col-1,col+0.5,Tmx[ip]*dTbin);
    if (Qmx2<QMIN)return;
    
    int ntotbins = nPages*size/nSlotInterval;
    fitFcn = new TF1("fitFcn",asymGaussian,0,ntotbins*DTSLOT,4);
    //fitFcn2 = new TF1("fitFcn",asymGaussian2,0,ntotbins*DTSLOT,8);

    //fitFcn = new TF1("fitFcn",asymGaussianU,0,150*dTbin,5);
    
    cout << " PIXEL: " << ip << " Qmax 2ndPeak: " << Qmx2 << " at T " << Tmx2 << endl;
    fitFcn -> SetParameters(Qmx2,Tmx2,4*DTSLOT,1.);
    
    
    htrace->Reset();hpull->Reset();
    for (int i=0; i<nPages*size/nSlotInterval; i++){
        double Qcal = hCalResiduals[it][ip]->GetBinContent(i+1);
        double dQcal = hCalResiduals[it][ip]->GetBinError(i+1);
        
        htrace->SetBinContent(i+1,Qcal);
        htrace->SetBinError(i+1,dQcal);
        // htrace->SetBinError(i+1,sqrt(20./200.));
    }
    TH1F *hpixels= (TH1F *)htrace->Clone();
    TH1F *hresiduals= (TH1F *)htrace->Clone();
    TH1F *hpixelsUC= (TH1F *)htrace->Clone();
    
    stringstream shi;
    //ss2<<inputFile1<<" dt:"<<it*2;
    shi<<"htrax" << it << "_Tel" << HitTel[it] << "_R"<<row<<"_C"<<col;
    string tith = shi.str();
    hpixels->SetName(tith.c_str());
    
    stringstream shiUC;
    //ss2<<inputFile1<<" dt:"<<it*2;
    shiUC<<"htraxUC" << it << "_Tel" << HitTel[it] <<   "_R"<<row<<"_C"<<col;
    string tithUC = shiUC.str();
    hpixelsUC->SetName(tithUC.c_str());
    
    stringstream shi2;
    //ss2<<inputFile1<<" dt:"<<it*2;
    shi2<<"hres" << it << "_Tel" << HitTel[it] << "_R"<<row<<"_C"<<col;
    string tith2 = shi2.str();
    hresiduals->SetName(tith2.c_str());
    
    //CorrectUndershoot(ip);
    cout << " Fitting 2ndPeak" << itp << endl;
    hpixels->Fit("fitFcn","R");
    //hCalPixels[it][ip]->Fit("fitFcn","R");
    float peakq = fitFcn->GetParameter(0) ;  float sigmapeakq= fitFcn->GetParError(0);
    float peakt = fitFcn->GetParameter(1) ;  float sigmapeakt= fitFcn->GetParError(1);
    float sigmaL  = fitFcn->GetParameter(2) ;  float sigmasigmaL = fitFcn->GetParError(2);
    float skew    = fitFcn->GetParameter(3)  ;  float sigmaskew = fitFcn->GetParError(3);
    
    
    if(peakt>25.&& peakt <nPages*size/nSlotInterval*DTSLOT) {
        cout << " FIT 2nd PEAK for pixel " << itp+1 << " " << peakt << " " << peakq << endl;
        hFineTmax2->SetBinContent(itp+1,peakt);
        hFineTmax2->SetBinError(itp+1,sigmapeakt);
        hFineQmax2->SetBinContent(itp+1,peakq);
        hFineQmax2->SetBinError(itp+1,sigmapeakq);
        /*hFineSigmaT->SetBinContent(itp+1,sigmaL);
        hFineSigmaT->SetBinError(itp+1,sigmasigmaL);
        hFineSkew->SetBinContent(itp+1,skew);
        hFineSkew->SetBinError(itp+1,sigmaskew);
        */
    }
    
    
    
    delete hpixels ;
    delete  hresiduals ;
}
// END of 2nd PULSE search


TVector3 FdProcessElves::FdSingleElvesBoltFinder(){
/*
* FdSingleElvesBoltFinder: preliminary estimate of bolt location, assuming Hd=80 km and Hbolt=0km
*  input: hFineTmax, split in rows
*/
    
    // in case of two bay event, the column number goes from 1 to 40
    /*int xc[3]={0,0,0};
    
    int ncx = NCOLS*3 ;
    
    int telTmin = HitTel[0];
    
    if (fL7TelescopeCount == 2 && HitTel[0] == HitTel[1]+1)xc[1]=1;
    else if (fL7TelescopeCount == 2 && HitTel[1] == HitTel[0]+1){xc[0]=1; telTmin = HitTel[1];}
    */
    // TODO1: what if the two Telescopes are not adjacent?
    // TODO2: what if fTelescopeCount == 3? (never seen sofar)
    
    int ncx = NCOLS*3 ;
    
    int telTmin = HitTel[0];
    

//  TimeMax_vs_Col[row-1]  is going to be fitted
// all other histos are used to clean the histos to be fit from spurious values
// Cuts on NphMax, SigmaT, Skew
    
    
    for (int ir=0; ir<NROWS; ir++){
        stringstream ss2;
        ss2 << "TimeMax_vs_Col"<< ir+1  ;
        stringstream ss3;
        ss3 << "Peak Time vs Column, for Row "<< ir+1 ;
        
        TimeMax_vs_Col[ir]=new TH1F(ss2.str().c_str(),ss3.str().c_str(),ncx,0.5,ncx+0.5);
        
        stringstream ss2a;
        ss2a << "NphMax_vs_Col"<< ir+1  ;
        stringstream ss3a;
        ss3a << "NphMax vs Column, for Row "<< ir+1 ;

        NphMax_vs_Col[ir]=new TH1F(ss2a.str().c_str(),ss3a.str().c_str(),ncx,0.5,ncx+0.5);

        stringstream ss2b;
        ss2b << "SigmaT_vs_Col"<< ir+1  ;
        stringstream ss3b;
        ss3b << "Sigma Time Raising Edge vs Column, for Row "<< ir+1 ;
        

        SigmaT_vs_Col[ir]=new TH1F(ss2b.str().c_str(),ss3b.str().c_str(),ncx,0.5,ncx+0.5);
        
        
        stringstream ss2c;
        ss2c << "Skew_vs_Col"<< ir+1  ;
        stringstream ss3c;
        ss3c << "Skewness of the Peak vs Column, for Row "<< ir+1 ;

        Skew_vs_Col[ir]=new TH1F(ss2c.str().c_str(),ss3c.str().c_str(),ncx,0.5,ncx+0.5);

    }
 
    int    NcolCnt[NROWS]={};
    for (int it=0; it<fL7TelescopeCount; it++){
    
        
    for (int ip=0; ip<NPIXELS; ip++){
        
        int row=ip%NROWS+1;
        int col=ip/NROWS+1+NCOLS*it;
        int itp = it*NPIXELS+ip;
        //int itp1 = xc[it]*NPIXELS+ip;
        
        Float_t Tmax1=hFineTmax->GetBinContent(itp+1);
        Float_t dTmax1=hFineTmax->GetBinError(itp+1);
        Float_t Qmax1=hFineQmax->GetBinContent(itp+1);
        Float_t dQmax1=hFineQmax->GetBinError(itp+1);
        Float_t sigmaT1=hFineSigmaT->GetBinContent(itp+1);
        Float_t dsigmaT1=hFineSigmaT->GetBinError(itp+1);
        Float_t skew1=hFineSkew->GetBinContent(itp+1);
        Float_t dskew1=hFineSkew->GetBinError(itp+1);

        // necessario un if per controllare cosa mettiamo e per contare i punti del fit
        if (Qmax1>300. && dQmax1/Qmax1<0.15 && sigmaT1 > 3. && dsigmaT1/sigmaT1<0.25 ){
            TimeMax_vs_Col[row-1]->SetBinContent(col,Tmax1);
            TimeMax_vs_Col[row-1]->SetBinError(col,dTmax1);
            NcolCnt[row-1]++;
            // 27/3/2017: sort them all with respect to LEFTMOST bay (index 'it' loops now on bays from the leftmost one)
            hFIT3Tmax->SetBinContent(itp+1,Tmax1);
            hFIT3Tmax->SetBinError(itp+1,dTmax1);
            hFIT3Qmax->SetBinContent(itp+1,Qmax1);
            hFIT3Qmax->SetBinError(itp+1,dQmax1);
            hFIT3SigmaT->SetBinContent(itp+1,sigmaT1);
            hFIT3SigmaT->SetBinError(itp+1,dsigmaT1);
            hFIT3Skew->SetBinContent(itp+1,skew1);
            hFIT3Skew->SetBinError(itp+1,dskew1);

        }
        
        NphMax_vs_Col[row-1]->SetBinContent(col,Qmax1);
        NphMax_vs_Col[row-1]->SetBinError(col,dQmax1);
        SigmaT_vs_Col[row-1]->SetBinContent(col,sigmaT1);
        SigmaT_vs_Col[row-1]->SetBinError(col,dsigmaT1);
        Skew_vs_Col[row-1]->SetBinContent(col,skew1);
        Skew_vs_Col[row-1]->SetBinError(col,dskew1);
 
      }
    }
    
    for (int it=0; it<fL7TelescopeCount; it++){
      for (int ip=0; ip<NPIXELS; ip++){
        
        int row=ip%NROWS+1;
        //int itp1 = xc[it]*NPIXELS+ip;
        int itp = it*NPIXELS+ip;
// 27/3/2017: sort them all with respect to LEFTMOST bay (index 'it' loops now on bays from the leftmost one)
        if (NcolCnt[row-1]<4) {
            hFIT3Tmax->SetBinContent(itp+1,0.);
            hFIT3Tmax->SetBinError(itp+1,0.);
         }
      }
    }
    
    TH1F *hrawAzi = new TH1F("hrawAzi","Raw Guess on Bolt Azimuth",100,-20.,80.);
 //   bool firstrow = true ;
    
    int NrowCnt = 0;
    
    for (int ir=0; ir<NROWS; ir++){
        if (NcolCnt[ir]<4) continue;
        TFitResultPtr r = TimeMax_vs_Col[ir]->Fit("pol2","S");
    // This is the fit to get the column with the first pixel in a row
        Double_t chi2   = r->Chi2(); // to retrieve the fit chi2
        Double_t a0   = r->Value(0); // retrieve the value for the parameter 0
        Double_t a1   = r->Value(1); // retrieve the value for the parameter 1
        Double_t a2   = r->Value(2); // retrieve the value for the parameter 2
        
        
        Double_t colphimin = -0.5*a1/a2+ 0.5*((ir+1)%2);
    // the term 0.5*((ir+1)%2 corrects for the staggering of the pixel center in each column
    // may be it is not necessary
        Double_t tmin  = a0-0.25*a1*a1/a2;
        
        cout << ThisGPSsec << " PREFIT row a0 a1 a2 col time " << ir+1 << " " << a0 << " " << a1 << " " << a2 << " " ;
        cout << colphimin << " " << tmin << endl;
        hPhiMinBoltFinder->SetBinContent(ir+1,colphimin);
        
        hChi2BoltFinder->SetBinContent(ir+1,chi2);
        
        if (chi2/NcolCnt[ir]<50.){
            hTMinBoltFinder->SetBinContent(ir+1,tmin);
            NrowCnt ++;
        }
        // unweighted filling
        hrawAzi->Fill(colphimin);

    }


    if (NrowCnt<4) {TVector3 rawbolt(1.,0.,0.); return rawbolt;}
    // this is the fit to get the row with the first pixel
    TFitResultPtr r2 = hTMinBoltFinder->Fit("pol2","S");
    
    Double_t chi2b   = r2->Chi2(); // to retrieve the fit chi2
    Double_t a0b   = r2->Value(0); // retrieve the value for the parameter 0
    Double_t a1b   = r2->Value(1); // retrieve the value for the parameter 1
    Double_t a2b   = r2->Value(2); // retrieve the value for the parameter 2
    Double_t rowzetmin = -0.5*a1b/a2b;
    Double_t tmin2  = a0b-0.25*a1b*a1b/a2b;
    
    cout << ThisGPSsec << " PREFIT2 rowmin tmin2 " << rowzetmin << " " <<  tmin2 << endl  ;

    // this fit averages all the rows to get the best azimuth: really necessary?
    TFitResultPtr r3 = hrawAzi->Fit("gaus","S");
    
    Double_t chi2c   = r3->Chi2(); // to retrieve the fit chi2
    Double_t a0c   = r3->Value(0); // retrieve the value for the parameter 0
    Double_t a1c   = r3->Value(1); // retrieve the value for the parameter 1
    Double_t a2c   = r3->Value(2); // retrieve the value for the parameter 2
    
    //cout << ThisGPSsec << " PREFIT3 Peak Mean=Ave? Sigma=RMS? " << a0c << " " <<  a1c ;
    //cout <<  "=" << hrawAzi->GetMean() <<  " ? " ;
    
    //cout << a2c << "="  << hrawAzi->GetRMS() << " ? " << endl  ;

    int colTmin = int(hrawAzi->GetMean());
    
    int tshift = 0;
    if (colTmin>NCOLS){ colTmin -= NCOLS; telTmin--; tshift=-1;}
    if (colTmin<1){ colTmin += NCOLS; telTmin++; tshift=1;}
    
    int rowTmin = int(rowzetmin);

    
    TCanvas *cBoltFinder = new TCanvas("cBoltFinder","Lightning Preliminary search", 1400, 1000);
    cBoltFinder -> Divide(2,1);
    cBoltFinder->cd(1);
    bool firstrow = true ;
    for (int ir=0; ir<NROWS; ir++){
        if (NcolCnt[ir]<4) continue;
        int rowcolor = ir+45 + 15*(ir%3);
        
        TimeMax_vs_Col[ir]->SetLineColor(rowcolor);
        TimeMax_vs_Col[ir]->SetMarkerStyle(21);
        TimeMax_vs_Col[ir]->SetMarkerColor(rowcolor);
        TimeMax_vs_Col[ir]->GetFunction("pol2")->SetLineColor(rowcolor);
        
        if (firstrow){ firstrow=false;  TimeMax_vs_Col[ir]->Draw(); }
        else TimeMax_vs_Col[ir]->Draw("same");

    }

    cBoltFinder->cd(2);

    hTMinBoltFinder->Draw();
    
    
    
    stringstream ss3;
    if (tshift==0)ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_BoltFinder_R" << rowTmin << "_C" << colTmin << "_T0.png";
    if (tshift==-1)ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_BoltFinder_R" << rowTmin << "_C" << colTmin << "_TR.png";
    if (tshift==1)ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_BoltFinder_R" << rowTmin << "_C" << colTmin << "_TL.png";
    
    cBoltFinder->SaveAs(ss3.str().c_str());
    
    delete cBoltFinder;
    

    
    TVector3 rawbolt=CalculateRawBoltLoc(colTmin,rowTmin,telTmin);
    

    return rawbolt;
}


TVector3 FdProcessElves::FdDoubleElvesBoltFinder(){
    /*
     * FdDoubleElvesBoltFinder: preliminary estimate of bolt location, assuming Hd=80 km and Hbolt=0km
     *  input: hFineTmax, split in rows
     */
    
    // in case of two bay event, the column number goes from 1 to 40
    //int xc[3]={0,0,0};
    
    int ncx = NCOLS*3 ;
    
    int telTmin = HitTel[0];
    
    //if (fL7TelescopeCount == 2 && HitTel[0] == HitTel[1]+1)xc[1]=1;
    //else if (fL7TelescopeCount == 2 && HitTel[1] == HitTel[0]+1){xc[0]=1; telTmin = HitTel[1];}

    // TODO1: what if the two Telescopes are not adjacent?
    // TODO2: what if fTelescopeCount == 3? (never seen sofar)
    
    //  TimeMax_vs_Col[row-1]  is going to be fitted
    // all other histos are used to clean the histos to be fit from spurious values
    // Cuts on NphMax, SigmaT, Skew
    
    
    for (int ir=0; ir<NROWS; ir++){
        stringstream ss2;
        ss2 << "TimeMax_vs_Col"<< ir+1  ;
        stringstream ss3;
        ss3 << "Peak Time vs Column, for Row "<< ir+1 ;
        
        TimeMax_vs_Col[ir]=new TH1F(ss2.str().c_str(),ss3.str().c_str(),ncx,0.5,ncx+0.5);
        
        stringstream ss2a;
        ss2a << "NphMax_vs_Col"<< ir+1  ;
        stringstream ss3a;
        ss3a << "NphMax vs Column, for Row "<< ir+1 ;
        
        NphMax_vs_Col[ir]=new TH1F(ss2a.str().c_str(),ss3a.str().c_str(),ncx,0.5,ncx+0.5);
        
        stringstream ss2b;
        ss2b << "SigmaT_vs_Col"<< ir+1  ;
        stringstream ss3b;
        ss3b << "Sigma Time Raising Edge vs Column, for Row "<< ir+1 ;
        
        
        SigmaT_vs_Col[ir]=new TH1F(ss2b.str().c_str(),ss3b.str().c_str(),ncx,0.5,ncx+0.5);
        
        
        stringstream ss2c;
        ss2c << "Skew_vs_Col"<< ir+1  ;
        stringstream ss3c;
        ss3c << "Skewness of the Peak vs Column, for Row "<< ir+1 ;
        
        Skew_vs_Col[ir]=new TH1F(ss2c.str().c_str(),ss3c.str().c_str(),ncx,0.5,ncx+0.5);
        
    }
    
    int    NcolCnt[NROWS]={};
    for (int it=0; it<fL7TelescopeCount; it++){
        
        
        for (int ip=0; ip<NPIXELS; ip++){
            
            int row=ip%NROWS+1;
            int col=ip/NROWS+1+NCOLS*it;
            int itp = it*NPIXELS+ip;
            //int itp1 = xc[it]*NPIXELS+ip;
            
            Float_t Tmax1=hFineTmax->GetBinContent(itp+1);
            Float_t dTmax1=hFineTmax->GetBinError(itp+1);
            Float_t Qmax1=hFineQmax->GetBinContent(itp+1);
            Float_t dQmax1=hFineQmax->GetBinError(itp+1);
            Float_t sigmaT1=hFineSigmaT->GetBinContent(itp+1);
            Float_t dsigmaT1=hFineSigmaT->GetBinError(itp+1);
            Float_t skew1=hFineSkew->GetBinContent(itp+1);
            Float_t dskew1=hFineSkew->GetBinError(itp+1);
            
            // necessario un if per controllare cosa mettiamo e per contare i punti del fit
            if (Qmax1>QMIN && dQmax1/Qmax1<0.15 && sigmaT1 > 3. && dsigmaT1/sigmaT1<0.25 ){
                TimeMax_vs_Col[row-1]->SetBinContent(col,Tmax1);
                TimeMax_vs_Col[row-1]->SetBinError(col,dTmax1);
                NcolCnt[row-1]++;
                hFIT3Tmax->SetBinContent(itp+1,Tmax1);
                hFIT3Tmax->SetBinError(itp+1,dTmax1);
                hFIT3Qmax->SetBinContent(itp+1,Qmax1);
                hFIT3Qmax->SetBinError(itp+1,dQmax1);
                hFIT3SigmaT->SetBinContent(itp+1,sigmaT1);
                hFIT3SigmaT->SetBinError(itp+1,dsigmaT1);
                hFIT3Skew->SetBinContent(itp+1,skew1);
                hFIT3Skew->SetBinError(itp+1,dskew1);
                
            }
            
            NphMax_vs_Col[row-1]->SetBinContent(col,Qmax1);
            NphMax_vs_Col[row-1]->SetBinError(col,dQmax1);
            SigmaT_vs_Col[row-1]->SetBinContent(col,sigmaT1);
            SigmaT_vs_Col[row-1]->SetBinError(col,dsigmaT1);
            Skew_vs_Col[row-1]->SetBinContent(col,skew1);
            Skew_vs_Col[row-1]->SetBinError(col,dskew1);
            
        }
    }
    
    for (int it=0; it<fL7TelescopeCount; it++){
        for (int ip=0; ip<NPIXELS; ip++){
            
            int row=ip%NROWS+1;
            int itp = it*NPIXELS+ip;
            
            if (NcolCnt[row-1]<4) {
                hFIT3Tmax->SetBinContent(itp,0.);
                hFIT3Tmax->SetBinError(itp,0.);
            }
        }
    }
    
    TH1F *hrawAzi = new TH1F("hrawAzi","Raw Guess on Bolt Azimuth",100,-20.,80.);
    //   bool firstrow = true ;
    
    int NrowCnt = 0;
    
    for (int ir=0; ir<NROWS; ir++){
        if (NcolCnt[ir]<4) continue;
        TFitResultPtr r = TimeMax_vs_Col[ir]->Fit("pol2","S");
        // This is the fit to get the column with the first pixel in a row
        Double_t chi2   = r->Chi2(); // to retrieve the fit chi2
        Double_t a0   = r->Value(0); // retrieve the value for the parameter 0
        Double_t a1   = r->Value(1); // retrieve the value for the parameter 1
        Double_t a2   = r->Value(2); // retrieve the value for the parameter 2
        
        
        Double_t colphimin = -0.5*a1/a2+ 0.5*((ir+1)%2);
        // the term 0.5*((ir+1)%2 corrects for the staggering of the pixel center in each column
        // may be it is not necessary
        Double_t tmin  = a0-0.25*a1*a1/a2;
        
        //cout << ThisGPSsec << " PREFIT row a0 a1 a2 col time " << ir+1 << " " << a0 << " " << a1 << " " << a2 << " " ;
        //cout << colphimin << " " << tmin << endl;
        hPhiMinBoltFinder->SetBinContent(ir+1,colphimin);
        
        hChi2BoltFinder->SetBinContent(ir+1,chi2);
        
        if (chi2/NcolCnt[ir]<50.){
            hTMinBoltFinder->SetBinContent(ir+1,tmin);
            NrowCnt ++;
        }
        // unweighted filling
        hrawAzi->Fill(colphimin);
        
    }
    
    
    if (NrowCnt<4) {TVector3 rawbolt(1.,0.,0.); return rawbolt;}
    // this is the fit to get the row with the first pixel
    TFitResultPtr r2 = hTMinBoltFinder->Fit("pol2","S");
    
    Double_t chi2b   = r2->Chi2(); // to retrieve the fit chi2
    Double_t a0b   = r2->Value(0); // retrieve the value for the parameter 0
    Double_t a1b   = r2->Value(1); // retrieve the value for the parameter 1
    Double_t a2b   = r2->Value(2); // retrieve the value for the parameter 2
    Double_t rowzetmin = -0.5*a1b/a2b;
    Double_t tmin2  = a0b-0.25*a1b*a1b/a2b;
    
    //cout << ThisGPSsec << " PREFIT2 rowmin tmin2 " << rowzetmin << " " <<  tmin2 << endl  ;
    
    // this fit averages all the rows to get the best azimuth: really necessary?
    TFitResultPtr r3 = hrawAzi->Fit("gaus","S");
    
    Double_t chi2c   = r3->Chi2(); // to retrieve the fit chi2
    Double_t a0c   = r3->Value(0); // retrieve the value for the parameter 0
    Double_t a1c   = r3->Value(1); // retrieve the value for the parameter 1
    Double_t a2c   = r3->Value(2); // retrieve the value for the parameter 2
    
    //cout << ThisGPSsec << " PREFIT3 Peak Mean=Ave? Sigma=RMS? " << a0c << " " <<  a1c ;
    //cout <<  "=" << hrawAzi->GetMean() <<  " ? " ;
    
    //cout << a2c << "="  << hrawAzi->GetRMS() << " ? " << endl  ;
    
    int colTmin = int(hrawAzi->GetMean());
    
    int tshift = 0;
    if (colTmin>NCOLS){ colTmin -= NCOLS; telTmin--; tshift=-1;}
    if (colTmin<1){ colTmin += NCOLS; telTmin++; tshift=1;}
    
    int rowTmin = int(rowzetmin);
    
    
    TCanvas *cBoltFinder = new TCanvas("cBoltFinder","Lightning Preliminary search", 1400, 1000);
    cBoltFinder -> Divide(2,1);
    cBoltFinder->cd(1);
    bool firstrow = true ;
    for (int ir=0; ir<NROWS; ir++){
        if (NcolCnt[ir]<4) continue;
        int rowcolor = ir+45 + 15*(ir%3);
        
        TimeMax_vs_Col[ir]->SetLineColor(rowcolor);
        TimeMax_vs_Col[ir]->SetMarkerStyle(21);
        TimeMax_vs_Col[ir]->SetMarkerColor(rowcolor);
        TimeMax_vs_Col[ir]->GetFunction("pol2")->SetLineColor(rowcolor);
        
        if (firstrow){ firstrow=false;  TimeMax_vs_Col[ir]->Draw(); }
        else TimeMax_vs_Col[ir]->Draw("same");
        
    }
    
    cBoltFinder->cd(2);
    
    hTMinBoltFinder->Draw();
    
    
    
    stringstream ss3;
    if (tshift==0)ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_BoltFinder_R" << rowTmin << "_C" << colTmin << "_T0.png";
    if (tshift==-1)ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_BoltFinder_R" << rowTmin << "_C" << colTmin << "_TR.png";
    if (tshift==1)ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_BoltFinder_R" << rowTmin << "_C" << colTmin << "_TL.png";
    
    cBoltFinder->SaveAs(ss3.str().c_str());
    
    delete cBoltFinder;
    
    
    
    TVector3 rawbolt=CalculateRawBoltLoc(colTmin,rowTmin,telTmin);
    
    
    return rawbolt;
}


TVector3 FdProcessElves::CalculateRawBoltLoc(int col, int row , int tel  ){

    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;

    const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
    const double ph = (11.66667 - row) * dPhi;
    
    const double mcosOO = -cos(oo);
    const double z = mcosOO * cos(ph);
    const double u = mcosOO * sin(ph);
    const double y = sin(oo);
    
    TVector3 XPixel(-u,-y,-z);
    XPixel.RotateY(73.5*degree);
    // tel is MIRROR NUMBER=1 to 6 , but can also be 0 or 7 if the ELVES center is outside view

    XPixel.RotateZ((backwall[site1]+15.+30.*(tel-1))*degree);
    
    //TVector3 PoloNord(0.,0.,1.);
    cout << ThisGPSsec << " ENTERING CalculateRawBoltLoc WITH col row tel " << col << " " << row << " " << tel << endl;
    
    cout << "SITE1LOCATION " << 90.-Site1Location.Theta()/degree << " " << Site1Location.Phi()/degree << endl;
    
    East1 = PoloNord.Cross(Site1Location);
    cout << ThisGPSsec << " EAST VECTOR  " << East1.Theta()/degree << " " << East1.Phi()/degree ;

    TVector3 AziView1 = East1;
    AziView1.Rotate(XPixel.Phi(),Site1Location);
    
    cout << "ROTATED " << XPixel.Phi()/degree << " degrees" << endl;
    TVector3 MaxCirc1 = Site1Location.Cross(AziView1);
    TVector3 rawbolt  = Site1Location;
    
    // this is the distance between the FD and the point of light emission,
    // at Hd altitude (default 80 km)
    
    double Rsite = Rearth;
    double EFD  = RootOf(1.,2.*Rsite*cos(XPixel.Theta()),Rsite*Rsite-pow(Rearth+Hd0,2),1);
    double cta = (EFD*EFD-pow(Rearth+Hd0,2)-Rsite*Rsite)/(2*Rsite*(Rearth+Hd0));
    if (cta>=1 || cta <=-1)cta = 1;
    double alpha = 2*acos(cta);

    rawbolt.Rotate(-alpha,MaxCirc1);

    return rawbolt;
}

int FdProcessElves::FdSingleElvesBoltFitter(TVector3 rawbolt, double t0  ){

    double boltLat = 90.-rawbolt.Theta()/degree;
    double boltLong = rawbolt.Phi()/degree;
    double hBolt = 0.;
    double tBolt = t0  ;
    
    cout << " FIT3DINPUT  t0, Lat, Long " << t0 << " " <<  boltLat  << " " << boltLong << endl;
    
    double Hd = Hd0;
   // corrected 00:05 Feb.6 : -15 degrees !
    double leftmostbaycenter = backwall[site1]+30*HitTel[0]-15.;
    
    if (fL7TelescopeCount == 2 && HitTel[1] == HitTel[0]+1)leftmostbaycenter +=30.;
    
    TF1 *fitT0bolt = new TF1("fitT0bolt",t0LinFitter ,0.,NPIXELS*fL7TelescopeCount*1. ,9);

// FIT NUMBER 1 : only tBolt,boltLat,boltLong ==========================================================================
    fitT0bolt -> SetParameters(tBolt,boltLat,boltLong,hBolt,Hd,hsite[site1],backwall[site1],latsite[site1],longsite[site1]);

    fitT0bolt -> FixParameter(3,0.);
    fitT0bolt -> FixParameter(4,Hd0);
    fitT0bolt -> FixParameter(5,hsite[site1]);
//    fitT0bolt -> FixParameter(6,backwall[site1]);
    fitT0bolt -> FixParameter(6,leftmostbaycenter);
    fitT0bolt -> FixParameter(7,latsite[site1]);
    fitT0bolt -> FixParameter(8,longsite[site1]);
    
    hFIT3Tmax->Fit("fitT0bolt","RB");

    // Refit if there are some points with very bad pulls
    if(CheckTimeFit5D_Output(fitT0bolt)>0)hFIT3Tmax->Fit("fitT0bolt","RB");
    for (int i=0; i<9; i++)BestFitParameter[i]=fitT0bolt->GetParameter(i);
   
    int fit5D_ndf[3];
    fit5D_ndf[0] =   PlotTimeFit5D_Output(fitT0bolt);
    //fitT0bolt -> ReleaseParameter(3);

    
    BoltLatitude = fitT0bolt->GetParameter(1);
    BoltLongitude = fitT0bolt->GetParameter(2);
    BoltAltitude = fitT0bolt->GetParameter(3);
    BoltGPStime = fitT0bolt->GetParameter(0);
    TLEAltitude = fitT0bolt->GetParameter(4);

    TMarker * mark1 = new TMarker(fitT0bolt->GetParameter(2), fitT0bolt->GetParameter(1),8);
    mark1->SetMarkerColor(sitecolor[site1]);
    mark1->SetMarkerSize(1.5);
    mark1->SetMarkerStyle(34);
    
    fSites2->Add(mark1);
    
    // FIT NUMBER 2 : tBolt,boltLat,boltLong  + Emission Layer Altitude  ===================================================================
    
    fitT0bolt -> ReleaseParameter(4);
    fitT0bolt -> SetParLimits(4,70.,110.);

    hFIT3Tmax->Fit("fitT0bolt","RB");

    for (int i=0; i<9; i++)BestFitParameter[i]=fitT0bolt->GetParameter(i);
    
    TMarker * mark1b = new TMarker(fitT0bolt->GetParameter(2), fitT0bolt->GetParameter(1),8);
    mark1b->SetMarkerColor(sitecolor[site1]);
    mark1b->SetMarkerSize(2.0);
    mark1b->SetMarkerStyle(34);
    
    fSites2->Add(mark1b);
    

    
        fit5D_ndf[1] = PlotTimeFit5D_Output(fitT0bolt);

    
    // FIT NUMBER 3 : tBolt,boltLat,boltLong  + Emission Layer Altitude + Bolt Altitude  ======================================================
    
    fitT0bolt -> ReleaseParameter(3);
    fitT0bolt -> SetParLimits(3,0.,20.);
    
    hFIT3Tmax->Fit("fitT0bolt","RB");
    for (int i=0; i<9; i++)BestFitParameter[i]=fitT0bolt->GetParameter(i);

    TMarker * mark1c = new TMarker(fitT0bolt->GetParameter(2), fitT0bolt->GetParameter(1),8);
    mark1c->SetMarkerColor(sitecolor[site1]);
    mark1c->SetMarkerSize(2.5);
    mark1c->SetMarkerStyle(34);
    
    fSites2->Add(mark1c);
    
  
    // Draw TWO circles at 100 and 200 km ArcDistance from the Bolt Vertical
    
    fB100km = new TGraph();
    fB200km = new TGraph();
    
    double Rsite=Rearth+hsite[site1];
    TVector3 BoltLocation(1.,1.,1.);
    BoltLocation.SetMag(Rsite+fitT0bolt->GetParameter(3));
    BoltLocation.SetPhi(fitT0bolt->GetParameter(2)*degree);
    BoltLocation.SetTheta((90.-fitT0bolt->GetParameter(1))*degree);

    TVector3 FDxB = Site1Location.Cross(BoltLocation);
    FDxB.SetMag(1.);
    TVector3 Bolt100km = BoltLocation;
    TVector3 Bolt200km = BoltLocation;
    Bolt100km.Rotate(100./Rearth,FDxB);
    Bolt200km.Rotate(200./Rearth,FDxB);
    
    
    for (int j=0; j<180; j++){
        Bolt100km.Rotate(2.*degree,BoltLocation);
        double Latpoint = (90. - Bolt100km.Theta()/degree);
        double Longpoint = Bolt100km.Phi()/degree;
        fB100km->SetPoint(j,Longpoint,Latpoint);
    }

    for (int j=0; j<180; j++){
        Bolt200km.Rotate(2.*degree,BoltLocation);
        double Latpoint = (90. - Bolt200km.Theta()/degree);
        double Longpoint = Bolt200km.Phi()/degree;
        fB200km->SetPoint(j,Longpoint,Latpoint);
    }

        fit5D_ndf[2] = PlotTimeFit5D_Output(fitT0bolt);
    return fit5D_ndf[0];
   
    
    /*BoltLatitude = boltLat;
    BoltLongitude = boltLong;
    BoltAltitude = hBolt;
    BoltGPStime = tBolt;
    TLEAltitude = Hd;  // or Hd0?
    
 
    BoltLatitude = fitT0bolt->GetParameter(1);
     BoltLongitude = fitT0bolt->GetParameter(2);
     BoltAltitude = fitT0bolt->GetParameter(3);
     BoltGPStime = fitT0bolt->GetParameter(0);
     TLEAltitude = fitT0bolt->GetParameter(4);  // or Hd0?
    
*/
    
   //cout << ThisGPSsec << " FIT3DOUTPUT2  t0, Lat, Long , Hb, Hd" << fitT0bolt->GetParameter(0) << " " << fitT0bolt->GetParameter(1) ;
    //cout << " " << fitT0bolt->GetParameter(2) << " " << fitT0bolt->GetParameter(3) <<" " << fitT0bolt->GetParameter(4) << endl;


}

int FdProcessElves::CheckTimeFit5D_Output(TF1* boltfit){
    
    int nbadpulls=0;
    for (int i=0; i<NPIXELS*fL7TelescopeCount; i++){
        
        int nnpix = i%NPIXELS;
        int bay=i/NPIXELS;
        int ir = nnpix%NROWS;
        double yt = hFIT3Tmax -> GetBinContent(i+1);
        double dyt = hFIT3Tmax -> GetBinError(i+1);
        if (yt>0.){
            double xt = boltfit->Eval(hFIT3Tmax->GetBinCenter(i+1));
            if (dyt>0.){
                double pull = ((yt-xt)/dyt);
                if (abs(pull)>10.){
                    nbadpulls++;
                    hFIT3Tmax -> SetBinContent(i+1,0.);
                    hFIT3Tmax -> SetBinError(i+1,0.);
            
                }
            }
            
        }
    }
    
    return nbadpulls;
}



int FdProcessElves::PlotTimeFit5D_Output(TF1* boltfit){

    Double_t siteLat= boltfit->GetParameter(7);
    Double_t siteLong= boltfit->GetParameter(8);

    TVector3 BoltLocation(1.,1.,1.);
    BoltLocation.SetMag(Rearth+boltfit->GetParameter(3));
    BoltLocation.SetPhi(boltfit->GetParameter(2)*degree);
    BoltLocation.SetTheta((90.-boltfit->GetParameter(1))*degree);

    double boltLat = boltfit->GetParameter(1);
    double boltLong = boltfit->GetParameter(2);
    
    double dlat  = boltLat-siteLat;
    double latm  = 0.5*(boltLat+siteLat);
    double dlong = boltLong-siteLong;
    double dN = Rearth*dlat*degree;
    double dE = Rearth*sin(latm*degree)*dlong*degree;
    
    
    double phibolt = atan2(dN,dE)/degree;
    double rawdist = sqrt(dN*dN+dE*dE);
    
    Double_t Rsite = Rearth+boltfit->GetParameter(5); //hsite[site1]

    double cph = cos(siteLong*degree);
    double sph = sin(siteLong*degree);
    double cth = cos((90.-siteLat)*degree);
    double sth = sin((90.-siteLat)*degree);
    
/*
    
    TVector3 BoltLocation(1.,1.,1.);
    BoltLocation.SetMag(Rearth+boltfit->GetParameter(3));
    BoltLocation.SetPhi(boltfit->GetParameter(2)*degree);
    BoltLocation.SetTheta((90.-boltfit->GetParameter(1))*degree);
    
    TVector3 FDxB = Site1Location.Cross(BoltLocation);
    TVector3 FDxBxFD = FDxB.Cross(Site1Location);
    //FDxBxFD.SetMag(1.);
    TVector3 AziBolt = East1.Cross(FDxBxFD);
    AziBolt.SetMag(1.);
    double sinPhi = (Site1Location*AziBolt)/Site1Location.Mag();
    double phibolt = asin(sinPhi)/degree;
   */
    
    cout << ThisGPSsec << " FIT5Dchisq/Ndof=" << int(boltfit->GetChisquare()+0.5) << "/" << boltfit->GetNDF() << endl;
    cout << ThisGPSsec << " FIT5DOUTPUT  t0, Lat, Long, Hbolt, Hd " << boltfit->GetParameter(0) << " " << boltfit->GetParameter(1) ;
    cout << " " << boltfit->GetParameter(2) << " " << boltfit->GetParameter(3)<< " " << boltfit->GetParameter(4) << " " << phibolt << endl;
    int fitndf = boltfit->GetNDF() ;
    int iHd = int(boltfit->GetParameter(4)+0.5);
    int iHb = int(boltfit->GetParameter(3)+0.5);
   
    TCanvas *cFIT3Plot = new TCanvas("cFIT3Plot","Fit Results", 2800, 2000);
    cFIT3Plot->Divide(2,1);
    
    
    stringstream ss;
    ss<< "GPS " << ThisGPSsec << "." << EventNanoTime << " EYE:" << site1+1 << " Fit #chi^{2}/Ndof=" ;
    ss << " " << int(boltfit->GetChisquare()+0.5) << "/" << boltfit->GetNDF() <<  "; Dist/c(#mu s); T(FD) (#mu s)" ;
    
    stringstream ss2;
    ss2 << "Tbolt=" << int(EventNanoTime/1000.+boltfit->GetParameter(0)) << " #Phi_{East}=" << phibolt;
    ss2 << " Lat=" << 0.01*int(100*boltfit->GetParameter(1)) << " Long=" << 0.01*int(100*boltfit->GetParameter(2));
    ss2 << " Hbolt=" << 0.1*int(10.*boltfit->GetParameter(3)+0.5) << " km  Hd=" << 0.1*int(10.*boltfit->GetParameter(4)+0.5) << " km ; Dist/c(#mu s); dT(FD) (#mu s)" ;
    

    stringstream ss3y;
    ss3y << "GPS " << ThisGPSsec << "." << EventNanoTime << " EYE:" << site1+1 << " Fit #chi^{2}/Ndof=" ;
    ss3y << " " << int(boltfit->GetChisquare()+0.5) << "/" << boltfit->GetNDF() <<  "; Dist_{BE}(km); dT(FD) (#mu s)" ;

    TGraphErrors *BestDTC_vs_DistB2E = new TGraphErrors();
    BestDTC_vs_DistB2E->SetTitle(ss3y.str().c_str());

    TGraphErrors *BestTC_vs_DistBolt = new TGraphErrors();
    BestTC_vs_DistBolt->SetTitle(ss.str().c_str());
    
    TGraphErrors *BestDTC_vs_DistBolt = new TGraphErrors();
    BestDTC_vs_DistBolt->SetTitle(ss2.str().c_str());
    
    int ntr[NROWS]={};int ntc[2*NCOLS]={}; int nq=0;
    for (int ir=0; ir<NROWS; ir++)InitializeRowMarkers(ir);
    for (int ic=0; ic<NCOLS*fL7TelescopeCount; ic++)InitializeColMarkers(ic);
    for (int i=0; i<NPIXELS*fL7TelescopeCount; i++){
        
        int nnpix = i%NPIXELS;
        int bay=i/NPIXELS;
        int ir = nnpix%NROWS;
        int ic = int(i/NROWS);
        dT_BFD[i] = PathFromBolt_to_FD(i,BestFitParameter)/SpeedOfLight;

        double yt = hFIT3Tmax -> GetBinContent(i+1);
        double ytc = yt-dT_BFD[i]-BestFitParameter[0];
        double dyt = hFIT3Tmax -> GetBinError(i+1);
        if (yt>0.){
            double xt = boltfit->Eval(hFIT3Tmax->GetBinCenter(i+1));
            
            double xtc = R_BE[i];
            
            double dumb=0.;
            if (dyt>0.){
             double pull = ((yt-xt)/dyt);
             if (pull>0)hFitPullsP->SetBinContent(i+1,pull);
             if (pull<0)hFitPullsN->SetBinContent(i+1,-pull);
            }
            BestTC_vs_DistBolt->SetPoint(nq,xt,yt);
            BestTC_vs_DistBolt->SetPointError(nq,dumb,dyt);
            BestDTC_vs_DistBolt->SetPoint(nq,xt,yt-xt);
            BestDTC_vs_DistBolt->SetPointError(nq,dumb,dyt);
            TC_vs_DistBoltRP[ir]->SetPoint(ntr[ir],xt,yt);
            DTC_vs_DistBoltRP[ir]->SetPoint(ntr[ir],xt,yt-xt);
            DTC_vs_DistBoltCP[ic]->SetPoint(ntc[ic],xt,yt-xt);
            
            BestDTC_vs_DistB2E->SetPoint(nq,xtc,ytc);
            BestDTC_vs_DistB2E->SetPointError(nq,dumb,dyt);
            DTC_vs_DistB2E_RP[ir]->SetPoint(ntr[ir],xtc,ytc);
            DTC_vs_DistB2E_CP[ic]->SetPoint(ntc[ic],xtc,ytc);
            
            
            nq++;ntr[ir]++;ntc[ic]++;
        }
    }
    
    
    
    
    /*
    cFIT3Plot-> cd(1);
    hFitPullsP -> SetFillColor(2);
    double maxPosPull = hFitPullsP->GetMaximum();
    double maxNegPull = hFitPullsN->GetMaximum();
    if (maxPosPull>maxNegPull){hFitPullsP->Draw();     hFitPullsN->Draw("same");}
    else {                     hFitPullsN->Draw();     hFitPullsP->Draw("same");}
    */
    /*
    cFIT3Plot-> cd(1);
    BestTC_vs_DistBolt->Draw("AWP");
    for (int ir=0; ir<NROWS; ir++){
        TC_vs_DistBoltRP[ir]->Draw("P");
    }
    */
    
    //cFIT3Plot-> cd(3);
    
    TGraph *dummyplot2 = new TGraph();
    
    dummyplot2->SetTitle(ss3y.str().c_str());
    dummyplot2->SetPoint(0,0,-35.);
    dummyplot2->SetPoint(1,500.,35.);

    
    cFIT3Plot-> cd(1);
    dummyplot2-> Draw("AWP");

    if (nq>0)BestDTC_vs_DistB2E->Draw("P");
    for (int ir=0; ir<NROWS; ir++){
        if (ntr[ir]>0)DTC_vs_DistB2E_RP[ir]->Draw("P");
    }
    for (int ic=0; ic<NCOLS; ic++){
        if (ntc[ic]>0)DTC_vs_DistB2E_CP[ic]->Draw("P");
    }
    
    
    //FdSingleElvesPull_CameraViewPlot(3);
    
    TGraph *dummyplot = new TGraph();
    
    dummyplot->SetTitle(ss2.str().c_str());
    dummyplot->SetPoint(0,0,-35.);
    dummyplot->SetPoint(1,300.+100*nPages,35.);
    
    
    
    cFIT3Plot-> cd(2);
    dummyplot-> Draw("AWP");
    if (nq>0)BestDTC_vs_DistBolt->Draw("P");
    
    
    
    
    TLegend *leg1r = new TLegend(0.75,0.15,0.9,0.9);
    leg1r->SetTextSize(0.025);
    
    for (int rp=0; rp<NROWS; rp++){
        if (ntr[rp]>0)DTC_vs_DistBoltRP[rp]->Draw("P");
        stringstream ssp;
        ssp<< "Row" << rp+1;
        stringstream sspn;
        sspn<< "TRow" << rp+1;
        leg1r->AddEntry(sspn.str().c_str(),ssp.str().c_str(),"p");
    }
    leg1r->Draw("P");
    
    TLegend *leg2c = new TLegend(0.6,0.15,0.75,0.9);
    leg2c->SetTextSize(0.025/fL7TelescopeCount);
    
    for (int cp=0; cp<NCOLS*fL7TelescopeCount; cp++){
        if (ntc[cp]>0)DTC_vs_DistBoltCP[cp]->Draw("P");
        stringstream ssp;
        ssp<< "Col" << cp+1;
        stringstream sspn;
        sspn<< "TCol" << cp+1;
        leg2c->AddEntry(sspn.str().c_str(),ssp.str().c_str(),"p");
    }
    leg2c->Draw("P");

    
    
    
    stringstream ss3;
    ss3<< "EVT" << ThisGPSsec << "." << EventNanoTime << "_HB_" << iHb << "_HD_" << iHd << ".png";
    
    cFIT3Plot->SaveAs(ss3.str().c_str());
    
    delete cFIT3Plot;
    
    
    if (rawdist>2000.)fitndf = 0;
    return fitndf;
}

void FdProcessElves::InitializeRowMarkers(int rp){
    stringstream ss;
    ss<< "TRow" << rp+1;
    stringstream ssd;
    ssd<< "TRow" << rp+1;

    int rowcolor = rp+45 + 15*(rp%3);
    
    TC_vs_DistBoltRP[rp] = new TGraph();
    TC_vs_DistBoltRP[rp]->SetMarkerSize(1.0);
//    TC_vs_DistBoltRP[rp]->SetMarkerColor(44*(rp%2+1)-2*rp+15);
    TC_vs_DistBoltRP[rp]->SetMarkerColor(rowcolor);
    TC_vs_DistBoltRP[rp]->SetMarkerStyle(20+rp%4);
    TC_vs_DistBoltRP[rp]->SetName(ss.str().c_str());

    DTC_vs_DistBoltRP[rp] = new TGraph();
    DTC_vs_DistBoltRP[rp]->SetMarkerSize(1.0);
 //   DTC_vs_DistBoltRP[rp]->SetMarkerColor(44*(rp%2+1)-2*rp+15);
    DTC_vs_DistBoltRP[rp]->SetMarkerColor(rowcolor);
    DTC_vs_DistBoltRP[rp]->SetMarkerStyle(20+rp%4);
    DTC_vs_DistBoltRP[rp]->SetName(ssd.str().c_str());

    DTC_vs_DistB2E_RP[rp] = new TGraph();
    DTC_vs_DistB2E_RP[rp]->SetMarkerSize(1.0);
    //   DTC_vs_DistBoltRP[rp]->SetMarkerColor(44*(rp%2+1)-2*rp+15);
    DTC_vs_DistB2E_RP[rp]->SetMarkerColor(rowcolor);
    DTC_vs_DistB2E_RP[rp]->SetMarkerStyle(20+rp%4);

    
}

void FdProcessElves::InitializeColMarkers(int cp){

    stringstream ssd;
    ssd<< "TCol" << cp+1;
    
    DTC_vs_DistBoltCP[cp] = new TGraph();
    float ms = 0.5+0.1*abs(cp-10*fL7TelescopeCount+0.5);
    
    DTC_vs_DistBoltCP[cp]->SetMarkerSize(ms);
    DTC_vs_DistBoltCP[cp]->SetMarkerColor(1);
    DTC_vs_DistBoltCP[cp]->SetMarkerStyle(32);
    if (cp<10*fL7TelescopeCount)DTC_vs_DistBoltCP[cp]->SetMarkerStyle(26);
    DTC_vs_DistBoltCP[cp]->SetName(ssd.str().c_str());

    DTC_vs_DistB2E_CP[cp] = new TGraph();

    DTC_vs_DistB2E_CP[cp]->SetMarkerSize(ms);
    DTC_vs_DistB2E_CP[cp]->SetMarkerColor(1);
    DTC_vs_DistB2E_CP[cp]->SetMarkerStyle(32);
    if (cp<10*fL7TelescopeCount)DTC_vs_DistB2E_CP[cp]->SetMarkerStyle(26);
    
}


void FdProcessElves::DrawSkewedDipole(double delta, double phi0){
    fBE1_100km = new TGraph();
    fBE1_200km = new TGraph();
    
    
    double RHD = Rearth+BestFitParameter[4];

    double Lat0 =BestFitParameter[1] + sin(phi0)*0.5*delta/RHD/degree;
    double Long0 = BestFitParameter[2] + cos(phi0)*0.5*delta/RHD/cos(Lat0*degree)/degree;
    
    TMarker * mark1c = new TMarker(Long0, Lat0,8);
    mark1c->SetMarkerColor(1);
    mark1c->SetMarkerSize(2.5);
    mark1c->SetMarkerStyle(33);
    
    fSites2->Add(mark1c);
    
    
    for (int j=0; j<180; j++){
        double AA = 100.; if(0.5*delta>AA)continue;
        double X = AA*cos(2*j*degree);
        double Y = AA*sqrt(1.-pow(0.5*delta/AA,2))*sqrt(1.-pow(X/AA,2));
        if (j>90)Y=-Y;
        double XP = X*cos(phi0)-Y*sin(phi0);
        double YP = Y*cos(phi0)+X*sin(phi0);
        double Longpoint = Long0 + XP/(RHD*cos(Lat0*degree))/degree;
        double Latpoint = Lat0 + YP/(RHD)/degree;
        
        cout  << " SKEWED_DIPOLE Long,Lat, Longpt, Latpt , RHD "  <<  Long0 << " " << Lat0 << " " << Longpoint << " " << Latpoint << " " << RHD << endl;
        fBE1_100km->SetPoint(j,Longpoint,Latpoint);
    }

for (int j=0; j<180; j++){
    double AA = 200.; if(0.5*delta>AA)continue;
    double X = AA*cos(2*j*degree);
    double Y = AA*sqrt(1.-pow(0.5*delta/AA,2))*sqrt(1.-pow(X/AA,2));
    if (j>90)Y=-Y;
    double XP = X*cos(phi0)-Y*sin(phi0);
    double YP = Y*cos(phi0)+X*sin(phi0);
    double Longpoint = Long0 + XP/(RHD*cos(Lat0*degree))/degree;
    double Latpoint = Lat0 + YP/(RHD)/degree;
    fBE1_200km->SetPoint(j,Longpoint,Latpoint);
}

                                   
    
}


Double_t FdProcessElves::t0LinFitter(Double_t *x, Double_t *par)
{
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    
    if (x[0]<=0 || x[0]>=2*NPIXELS) return 0;
    
    int npix = int(x[0]);
    
    
    
    Double_t boltLat = par[1];
    Double_t boltLong = par[2];
    Double_t hBolt = par[3];
    Double_t tBolt = par[0];
    Double_t Hd = par[4];
    Double_t Rsite = Rearth+par[5]; //hsite[site1]
    Double_t LeftBayCenter = par[6]; //backwall[site1]
    Double_t siteLat= par[7];
    Double_t siteLong= par[8];

    
    double cph = cos(siteLong*degree);
    double sph = sin(siteLong*degree);
    double cth = cos((90.-siteLat)*degree);
    double sth = sin((90.-siteLat)*degree);

    Double_t Xs1 = Rsite*sth*cph;
    Double_t Ys1 = Rsite*sth*sph;
    Double_t Zs1 = Rsite*cth;
    
    int nnpix = npix%NPIXELS;
    int bay=npix/NPIXELS;
    int row = nnpix%NROWS+1;
    int col = int(nnpix/NROWS)+1;
    const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
    const double ph = (11.66667 - row) * dPhi;
    
    const double mcosOO = -cos(oo);
    const double z = mcosOO * cos(ph);
    const double u = mcosOO * sin(ph);
    const double y = sin(oo);
    
    TVector3 XPixel(-u,-y,-z);
    XPixel.RotateY(73.5*degree);
    // tel is MIRROR NUMBER=1 to 6 , but can also be 0 or 7 if the ELVES center is outside view
    
//    Left side = Backwall+30*tel where tel is the left most telescope in the image
//    LeftBayCenter = Left side -15
    //XPixel.RotateZ((BWsite+15.+30.*(tel-1))*degree);
    XPixel.RotateZ((LeftBayCenter-30.*bay)*degree);
    TVector3 NorthPole(0.,0.,1.);
    TVector3 FDSiteLocation(Xs1,Ys1,Zs1);
    
    
    TVector3 East2 = NorthPole.Cross(FDSiteLocation);
    
    East2.SetMag(1.);

    
    
    TVector3 AziView1 = East2;
    AziView1.Rotate(XPixel.Phi(),FDSiteLocation);
    TVector3 MaxCirc1 = FDSiteLocation.Cross(AziView1);
    TVector3 PixelHdVert = FDSiteLocation;
    
    
    // this is the distance between the FD and the point of light emission,
    // at Hd altitude (default 90 km)
    double EFD  = RootOf(1.,2.*Rsite*cos(XPixel.Theta()),Rsite*Rsite-pow(Rearth+Hd,2),1);
    
    double cta = (EFD*EFD-pow(Rearth+Hd,2)-Rsite*Rsite)/(2*Rsite*(Rearth+Hd));
    if (cta>=1 || cta <=-1)cta = 1;
    double alpha_EFD = acos(cta);
    
    //R_EFD[npix] = (Rearth+Hd)*alpha_EFD;
    //dT_EFD[npix] = EFD/SpeedOfLight ;
    
    // The vertical at FD site is rotated on the max circle to the vertical at light emission
    // and Latitude+Longitude are calculated
    PixelHdVert.Rotate(-alpha_EFD,MaxCirc1);
    if (PixelHdVert.Angle(FDSiteLocation)>30.*degree) PixelHdVert = - PixelHdVert ;
    PixelHdVert.SetMag(Rearth+Hd);
    
    
    
    
    TVector3 BoltLocation(1.,1.,1.);
    BoltLocation.SetMag(Rearth+hBolt);
    BoltLocation.SetPhi(boltLong*degree);
    BoltLocation.SetTheta((90.-boltLat)*degree);
    
    TVector3 B2E = PixelHdVert-BoltLocation;
    double B2FDtime = (EFD+B2E.Mag())/SpeedOfLight+tBolt;
    
    return B2FDtime;
}




double FdProcessElves::PathFromBolt_to_FD(int npix, Double_t *par)
{
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    
    
    
    Double_t boltLat = par[1];
    Double_t boltLong = par[2];
    Double_t hBolt = par[3];
    Double_t tBolt = par[0];
    Double_t Hd = par[4];
    Double_t Rsite = Rearth+par[5]; //hsite[site1]
    Double_t LeftBayCenter = par[6]; //backwall[site1]
    Double_t siteLat= par[7];
    Double_t siteLong= par[8];
    
    
    int nnpix = npix%NPIXELS;
    int bay=npix/NPIXELS;
    int row = nnpix%NROWS+1;
    int col = int(nnpix/NROWS)+1;
    const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
    const double ph = (11.66667 - row) * dPhi;
    
    const double mcosOO = -cos(oo);
    const double z = mcosOO * cos(ph);
    const double u = mcosOO * sin(ph);
    const double y = sin(oo);
    
    TVector3 XPixel(-u,-y,-z);
    XPixel.RotateY(73.5*degree);
    // tel is MIRROR NUMBER=1 to 6 , but can also be 0 or 7 if the ELVES center is outside view
    
    //    Left side = Backwall+30*tel where tel is the left most telescope in the image
    //    LeftBayCenter = Left side -15
    //XPixel.RotateZ((BWsite+15.+30.*(tel-1))*degree);
    XPixel.RotateZ((LeftBayCenter-30.*bay)*degree);
    
    TVector3 AziView1 = East1;
    AziView1.Rotate(XPixel.Phi(),Site1Location);
    TVector3 MaxCirc1 = Site1Location.Cross(AziView1);
    TVector3 PixelHdVert = Site1Location;
    
    
    // this is the distance between the FD and the point of light emission,
    // at Hd altitude (default 90 km)
    double EFD  = RootOf(1.,2.*Rsite*cos(XPixel.Theta()),Rsite*Rsite-pow(Rearth+Hd,2),1);
    
    double cta = (EFD*EFD-pow(Rearth+Hd,2)-Rsite*Rsite)/(2*Rsite*(Rearth+Hd));
    if (cta>=1 || cta <=-1)cta = 1;
    double alpha_EFD = acos(cta);
    
    R_EFD[npix] = (Rearth+Hd)*(180.*degree-alpha_EFD); // This is a CURVED distance at Ionosphere Base
    dT_EFD[npix] = EFD/SpeedOfLight ;
    
    // The vertical at FD site is rotated on the max circle to the vertical at light emission
    // and Latitude+Longitude are calculated
    PixelHdVert.Rotate(-alpha_EFD,MaxCirc1);
    if (PixelHdVert.Angle(Site1Location)>30.*degree) PixelHdVert = - PixelHdVert ;
    PixelHdVert.SetMag(Rearth+Hd);
    
    
    
    
    TVector3 BoltLocation(1.,1.,1.);
    BoltLocation.SetMag(Rearth+hBolt);
    BoltLocation.SetPhi(boltLong*degree);
    BoltLocation.SetTheta((90.-boltLat)*degree);
    
    TVector3 B2E = PixelHdVert-BoltLocation;
    double alpha_BE = PixelHdVert.Angle(BoltLocation);
    //  double B2FDtime = (EFD+B2E.Mag())/SpeedOfLight+tBolt;
    double B2FD=EFD+B2E.Mag();
    R_BE[npix] = (Rearth+Hd)*alpha_BE;// This is a CURVED distance at Ionosphere Base
    Dist_BE[npix] = B2E.Mag(); // This is the LINEAR distance between Bolt and Emitting Pixel
    
    double theta_B2E = B2E.Angle(BoltLocation);
    SinTh_BE[npix] = sin(theta_B2E);
    
    TVector3 BoltAxisTilted = BoltLocation;
    TVector3 EastBolt = PoloNord.Cross(BoltLocation);
    TVector3 NorthBolt = BoltLocation.Cross(EastBolt);
    
    
    BoltAxisTilted.Rotate(10.*degree,NorthBolt);
    BoltAxisTilted.Rotate(-20.*degree,BoltLocation);
    
    double theta_B2Etilt = B2E.Angle(BoltAxisTilted);
    SinTh_BE[npix] = sin(theta_B2Etilt);
    
    
    
    
    
    
    return B2FD;
}



double FdProcessElves::CalculateEllipticDistance(int npix, double delta, double phi0)
{
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    
    
    
    /*
    Double_t boltLat = par[1];
    Double_t boltLong = par[2];
    Double_t hBolt = par[3];
    Double_t tBolt = par[0];
    Double_t Hd = par[4];
    Double_t Rsite = Rearth+par[5]; //hsite[site1]
    Double_t LeftBayCenter = par[6]; //backwall[site1]
    Double_t siteLat= par[7];
    Double_t siteLong= par[8];
    */
    double R_e = Rearth+BestFitParameter[4];
    TVector3 BoltVertical(1.,1.,1.);
    BoltVertical.SetMag(R_e);
    BoltVertical.SetPhi(BestFitParameter[2]*degree);
    BoltVertical.SetTheta((90.-BestFitParameter[1])*degree);
    
    //double Lat2rad=LatPixelVert[npix]*degree +    sin(phi0)*delta/(Rearth+BestFitParameter[4]);
    //double Long2rad=LongPixelVert[npix]*degree +  cos(phi0)*delta/(Rearth+BestFitParameter[4]) / cos(LatPixelVert[npix]*degree);
    
    double Lat2rad=BestFitParameter[1]*degree +    sin(phi0)*delta/(Rearth+BestFitParameter[4]);
    double Long2rad=BestFitParameter[2]*degree +  cos(phi0)*delta/(Rearth+BestFitParameter[4]) / cos(BestFitParameter[1]*degree);
    
    TVector3 PixelVertical(1.,1.,1.);
    PixelVertical.SetMag(R_e);
    PixelVertical.SetPhi(LongPixelVert[npix]*degree);
    PixelVertical.SetTheta((90.-LatPixelVert[npix])*degree);

    TVector3 F1Vertical(1.,1.,1.);
    F1Vertical.SetMag(R_e);
    F1Vertical.SetPhi(Long2rad);
    F1Vertical.SetTheta(90.*degree-Lat2rad);
    
    double ED = 0.5*R_e*(BoltVertical.Angle(PixelVertical)+F1Vertical.Angle(PixelVertical));
    double ED1 = 0.5*R_e*(BoltVertical.Angle(PixelVertical));
    double ED2 = 0.5*R_e*(F1Vertical.Angle(PixelVertical));
                          
    if (npix == 100) cout<< ThisGPSsec << " EllipticDistance " << npix << " " <<  delta << " " << phi0 << " " << Lat2rad/degree << " " << Long2rad/degree
       <<  " " << ED <<  " " << ED1 <<  " " << ED2 << endl;
    
    
    
    return ED;
}



    double FdProcessElves::AngDistFromFirstPixel(double elev)
{
/****    AngDistFromFirstPixel:
*     tes the preliminary distance (arc angle on earth surface)
*     from the lightning bolt.
*     Input: elevation of the first pixel, IN RADIANS!!!
*     Caveat: same altitude between source and observer
*     Default: Hd=Hd0
*****************************/
    
    double Rd=Rearth+Hd0;
    double A=pow(Rd,2)*(1./sin(elev));
    double B=-2*Rearth*Rd*cos(elev)/sin(elev);
    double C=pow(Rearth*sin(elev)/cos(elev),2)-pow(Rd,2);
    double cosalf=RootOf(A,B,C,1);
    double alf=2*acos(cosalf);

}




double FdProcessElves::CalculateGeomCorr( int npix, int sitex=0, int ltel05=5){
    
    
   
    TVector3 Site1Location(1.,1.,1.);
    Site1Location.SetMag(1.);
    Site1Location.SetPhi(longsite[sitex]*degree);
    Site1Location.SetTheta((90.-latsite[sitex])*degree);
    
    
    const double eta0 = 0.0261799387;
    const double dOmega = eta0;
    const double dPhi = sqrt(3.) * eta0/2;
    
    int nnpix = npix%NPIXELS;
    int row = nnpix%NROWS+1;
    int col = int(nnpix/NROWS)+1;
    int bay = npix/NPIXELS;
    
    
    const double oo = (col - ((row%2) ? 10. : 10.5)) * dOmega;
    const double ph = (11.66667 - row) * dPhi;
    
    const double mcosOO = -cos(oo);
    const double z = mcosOO * cos(ph);
    const double x = mcosOO * sin(ph);
    const double y = sin(oo);

    TVector3 XPixel(-x,-y,-z);
    
    XPixel.RotateY(73.5*degree);
    
    XPixel.RotateZ((backwall[sitex]+15.+30.*(ltel05-bay))*degree);
    
    East1 = PoloNord.Cross(Site1Location);
    
    East1.SetMag(1.);
    
    //cout << " SITE1LOCATION  LAT,LONG " << 90.-Site1Location.Theta()/degree << " " << Site1Location.Phi()/degree << endl;
    
    double Rsite = Rearth+hsite[sitex];
    double Hd = TLEAltitude;
    
    TVector3 AziView1 = East1;
    AziView1.Rotate(XPixel.Phi(),Site1Location);
    TVector3 MaxCirc1 = Site1Location.Cross(AziView1);
    TVector3 PixelHdVert = Site1Location;
    
    // this is the distance between the FD and the point of light emission,
    // at Hd altitude (default 90 km)
    double EFD  = RootOf(1.,2.*Rsite*cos(XPixel.Theta()),Rsite*Rsite-pow(Rearth+Hd,2),1);
    double PixelEFD = EFD;
    double cta = (EFD*EFD-pow(Rearth+Hd,2)-Rsite*Rsite)/(2*Rsite*(Rearth+Hd));
    if (cta>=1 || cta <=-1)cta = 1;
    double alpha_EFD = acos(cta);
    
    // The vertical at FD site is rotated on the max circle to the vertical at light emission
    // and Latitude+Longitude are calculated
    PixelHdVert.Rotate(-alpha_EFD,MaxCirc1);
    if (PixelHdVert.Angle(Site1Location)>30.*degree) PixelHdVert = - PixelHdVert ;
    PixelHdVert.SetMag(Rearth+Hd);
    
    double LatPixelHdV = 90.-PixelHdVert.Theta()/degree;
    double LongPixelHdV = PixelHdVert.Phi()/degree;
    
    LatPixelVert[npix] = LatPixelHdV;
    LongPixelVert[npix] = LongPixelHdV;
    if (row == 21 && (col == 0 || col ==19))
    cout << ThisGPSsec << " PIXELCHECK " << npix << " C" << col+1 << "R" << row+1 << " bay " << bay <<endl;
    
    
    
    //    PixelDirs[site1*NTELS+tel]->SetPoint(npix,XPixel.Phi()/degree,90.-XPixel.Theta()/degree);
    //    PixelDirsLL[site1*NTELS+tel]->SetPoint(npix,LongPixelHdV,LatPixelHdV);
    
    
    
    //PixelDirs2->SetPoint(npix,XPixel2.Phi()/degree,90.-XPixel2.Theta()/degree);
    double oo2[7], ph2[7];
    oo2[0] = oo+dOmega/2.; ph2[0] = ph-dPhi/3.;
    oo2[1] = oo+dOmega/2.; ph2[1] = ph+dPhi/3.;
    oo2[2] = oo;                     ph2[2] = ph+2.*dPhi/3.;
    oo2[3] = oo-dOmega/2.; ph2[3] = ph+dPhi/3.;
    oo2[4] = oo-dOmega/2.; ph2[4] = ph-dPhi/3.;
    oo2[5] = oo;                     ph2[5] = ph-2.*dPhi/3.;
    oo2[6] = oo2[0];               ph2[6] = ph2[0];
    PixelEdges[npix] = new TGraph();
    PixelEdgesLL[npix] = new TGraph();
    
    if (row == 21 && (col == 0 || col ==19)) cout << ThisGPSsec << "PIXRotation " << (backwall[sitex]+15.+30.*(ltel05-bay)) << " degrees " << endl;
    double HexArea = 0;
    for (int iV = 0; iV<7 ; iV++){
        
        double mcosOO2 = -cos(oo2[iV]);
        double z2 = mcosOO2 * cos(ph2[iV]);
        double x2 = mcosOO2 * sin(ph2[iV]);
        double y2 = sin(oo2[iV]);
        
        
        
        TVector3 *XPixEdge = new TVector3(-x2,-y2,-z2);
        XPixEdge->RotateY(73.5*degree);
        XPixEdge->RotateZ((backwall[sitex]+15.+30.*(ltel05-bay))*degree);


        PixelEdges[npix]->SetPoint(iV,XPixEdge->Phi()/degree,90.-XPixEdge->Theta()/degree);
        TVector3 AziView2 = East1;
        AziView2.Rotate(XPixEdge->Phi(),Site1Location);
        TVector3 MaxCirc2 = Site1Location.Cross(AziView2);
        TVector3 PixEdgeHdVert = Site1Location;
        // this is the distance between the FD and the point of light emission,
        // at Hd altitude (default 90 km)
        EFD = RootOf(1.,2.*Rsite*cos(XPixEdge->Theta()),Rsite*Rsite-pow(Rearth+Hd,2),1);
        cta = (EFD*EFD-pow(Rearth+Hd,2)-Rsite*Rsite)/(2*Rsite*(Rearth+Hd));
        if (cta>=1 || cta <=-1)cta = 1;
        alpha_EFD = acos(cta);
        
        PixEdgeHdVert.Rotate(-alpha_EFD,MaxCirc2);
        if (PixEdgeHdVert.Angle(Site1Location)>30.*degree) PixEdgeHdVert = - PixEdgeHdVert ;
        
        double LatPixEdgeHdV = 90.-PixEdgeHdVert.Theta()/degree;
        double LongPixEdgeHdV = PixEdgeHdVert.Phi()/degree;
        PixelEdgesLL[npix]->SetPoint(iV,LongPixEdgeHdV,LatPixEdgeHdV);
        
 //       if (row == 21 && (col == 0 || col ==19)) cout << "PIXEDGE"<< iV << " " << LongPixEdgeHdV << " " << LatPixEdgeHdV << endl;
        
        TVector3 RPixEdgeHdV = PixEdgeHdVert;
        RPixEdgeHdV.SetMag(Rearth+Hd);
        
        
        if (iV>0) {Double_t xLong, xLat;  PixelEdgesLL[npix]->GetPoint(iV-1,xLong,xLat);
            double xth = (90.-xLat)*degree;
            double xph = xLong*degree;
            
            TVector3 xPixEdge(sin(xth)*cos(xph),sin(xth)*sin(xph),cos(xth));
            double L0 = (Rearth+Hd)*xPixEdge.Angle(PixEdgeHdVert);
            double L1 = (Rearth+Hd)*PixelHdVert.Angle(PixEdgeHdVert);
            double L2 = (Rearth+Hd)*PixelHdVert.Angle(xPixEdge);
            double SP = (L0+L1+L2)/2.;
            double Area = sqrt(SP*(SP-L0)*(SP-L1)*(SP-L2));
            HexArea += Area;
        }
    }// end loop over Hexagon
    
    Area_Ionosphere[npix]=HexArea;
    double correction=(PixelEFD*PixelEFD)/Amirror/Area_Ionosphere[npix];
    return correction;
    
    
}



Double_t FdProcessElves::circumference(Double_t *x, Double_t *par){
    
    Float_t xx =x[0];
    
    Double_t a = -2*par[1];
    Double_t b = -2*par[2];
    Double_t c = TMath::Power(par[1],2) + TMath::Power(par[2],2) - TMath::Power(par[0],2);
    
    Double_t alpha = 1.;
    Double_t beta =  b;
    Double_t gamma = TMath::Power(xx,2) + a*xx + c;
    
    Double_t f = (- beta - TMath::Sqrt(TMath::Power(beta,2) - 4*alpha*gamma))/alpha;
    return f;
}

Double_t FdProcessElves::asymGaussian(Double_t *x, Double_t *par)
{
    Double_t y=0.,xr=0.;
    if (x[0]<=par[1]) {xr = (x[0]-par[1])/par[2]; y = par[0]*exp(-xr*xr);}
    if (x[0]>par[1]) {xr = (x[0]-par[1])/(par[2]*par[3]); y = par[0]*exp(-xr*xr);}
    return y;
}

Double_t FdProcessElves::P2exp(Double_t *x, Double_t *par)
{
    Double_t y=pow(par[0]*x[0],2)*exp(-2*x[0]/par[1]);
    return y;
}



Double_t FdProcessElves::MyFunc2D(Double_t *x, Double_t *par){
    
    Double_t xx =x[0]-par[2];
    Double_t yy =x[1]-par[3];
    Double_t zeta = TMath::ATan2(yy,xx);
    Double_t rr = par[0];
    Double_t ww = par[1];
    Double_t Norm = par[4];
    Double_t zeta0 = par[5];
    Double_t dzeta = par[6];
    Double_t radius=TMath::Sqrt(xx*xx+yy*yy);
    Double_t arg1=TMath::Power((radius-rr)/ww,2.)/2.;
    Double_t arg2=TMath::Power((zeta-zeta0)/dzeta,2.)/2.;
    
    
    Double_t f = Norm*exp(-arg1)*exp(-arg2);
    
    return f;
}




void FdProcessElves::Fit2(const fevt::Eye & eye){
    
    gROOT->SetStyle("Plain");
    
    // Color Palette creation
    int          palette[MaxColors];
    int          index;
    float        r, g, b;
    TColor       *color;
    
    for (int i=0 ; i<MaxColors ; i++) {
        index = 500 + i;
        palette[i] = index;
        
        b = 1.0 * (1.0 - float(MaxColors-i)/float(MaxColors));
        g = 1.0 * (1.0 - float(MaxColors-i)/float(MaxColors));
        r = 1.0 * (1.0 - float(MaxColors-i)/float(MaxColors));
        color = new TColor(index, r, g, b);
        
    }
    
    TColor::SetPalette(MaxColors,palette);
    gStyle->SetNumberContours(MaxColors);
    
    
    fPixels->Delete();
    
    // draw all  the pixels in the telescope onto a Canvas
    // the  ones with a pulse in different color
    
    //  int firstpixelid= fdet::Telescope::GetFirstPixelId();
    //  int lastpixelid = fdet::Telescope::GetLastPixelId();
    
    double minphi = 360.;
    double maxphi = 0;
    bool change_phi = false;
    
    Point eyeposition = det::Detector::GetInstance().GetFDetector().GetEye(eye).GetPosition();
    CoordinateSystemPtr eye2CS = fwk::LocalCoordinateSystem::Create(eyeposition);
    
    // eye2CS is centered in the Eye, with the X axis orientated towards East.
    // So, It's not parallel to the eye Back walll like the standar eye
    // coordinate system 'eyeCS'.
    
    fevt::Eye::ConstTelescopeIterator teliter;
    
    for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        const fevt::Telescope& tel = *teliter;
        
        // note the difference between det::Telescope and fevt::Telescope
        const fdet::Telescope & dettelescope = Detector::GetInstance().GetFDetector().GetTelescope(tel);
        int firstpixelid= dettelescope.GetFirstPixelId();
        int lastpixelid = dettelescope.GetLastPixelId();
        
        for (int i = firstpixelid; i<=lastpixelid; i++){
            const fdet::Pixel & detpixel= dettelescope.GetPixel(i);
            const Vector&  dir = detpixel.GetDirection();
            
            //double theta = dir.GetTheta(eye2CS)/ utl::degree;
            double phi =   dir.GetPhi(eye2CS) /  utl::degree ;  if (phi<0) phi+= 360.;
            
            
            if (phi > maxphi) maxphi = phi;
            if (phi < minphi) minphi = phi;
            
        }//for i
        
    }
    
    
    
    if ( maxphi - minphi > 180 ) {
        change_phi = true;
        maxphi = -99999;
        minphi = 999999;
    }
    
    
    if (change_phi == true ) {
        
        for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
            
            const fevt::Telescope& tel = *teliter;
            
            // note the difference between det::Telescope and fevt::Telescope
            const fdet::Telescope & dettelescope = Detector::GetInstance().GetFDetector().GetTelescope(tel);
            
            int firstpixelid= dettelescope.GetFirstPixelId();
            int lastpixelid = dettelescope.GetLastPixelId();
            
            for (int i = firstpixelid; i<=lastpixelid; i++){
                const fdet::Pixel & detpixel= dettelescope.GetPixel(i);
                const Vector&  dir = detpixel.GetDirection();
                
                //double theta = dir.GetTheta(eye2CS)/ utl::degree;
                double phi =   dir.GetPhi(eye2CS) /  utl::degree ;
                
                if (phi > 180) phi -= 360;
                
                if (phi > maxphi) maxphi = phi;
                if (phi < minphi) minphi = phi;
                
            }
        }
    }
    
    
    for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        const fevt::Telescope& tel = *teliter;
        
        // note the difference between det::Telescope and fevt::Telescope
        const fdet::Telescope & dettelescope = Detector::GetInstance().GetFDetector().GetTelescope(tel);
        int firstpixelid= dettelescope.GetFirstPixelId();
        int lastpixelid = dettelescope.GetLastPixelId();
        
        
        for (int i = firstpixelid; i<=lastpixelid; i++){
            const fdet::Pixel & detpixel= dettelescope.GetPixel(i);
            const Vector&  dir = detpixel.GetDirection();
            
            double theta = dir.GetTheta(eye2CS)/ utl::degree;
            double phi =   dir.GetPhi(eye2CS) /  utl::degree ;
            
            if (change_phi == false ) if (phi<0) phi+= 360.;
            
            if (change_phi == true ) if (phi > 180) phi -= 360;
            
            TMarker * mark = new TMarker(phi, 90-theta ,4); //empty dot
            mark->SetMarkerColor(15);
            fPixels->Add(mark);
            
        }//for i
    } // for teliter
    
    
    
    
    int FirstPulseStart=9999;
    double FirstPulseStartTheta=0., FirstPulseStartPhi=0.;
    
    fevt::EyeRecData::ConstPixelIterator pixiter;
    const fevt::EyeRecData& eyerec = eye.GetRecData();
    
    for (pixiter= eyerec.PulsedPixelsBegin(); pixiter!=eyerec.PulsedPixelsEnd(); ++pixiter){
        
        const fevt::Pixel & pixel = *pixiter;
        const fdet::Pixel & detpixel= Detector::GetInstance().GetFDetector().GetPixel(pixel);
        
        Vector dir = detpixel.GetDirection();
        dir.TransformTo(eye2CS);
        
        Vector sdp = eyerec.GetSDP();
        sdp.TransformTo(eye2CS);
        
        double theta = dir.GetTheta(eye2CS)/ utl::degree;
        double phi =   dir.GetPhi(eye2CS) /  utl::degree;
        
        if (change_phi == false ) if (phi<0) phi+= 360.;
        
        if (change_phi == true ) if (phi > 180) phi -= 360;
        
        const fevt::PixelRecData & pixelrecdata = pixel.GetRecData();
        
        
        if ( pixelrecdata.GetPulseStart() <= FirstPulseStart ){
            FirstPulseStartTheta = (90 - theta);
            FirstPulseStartPhi = phi;
            FirstPulseStart = pixelrecdata.GetPulseStart();
        }
        
    }
    
    cout<<"FirstPulseStart="<<FirstPulseStart<<"    FirstPulseStartTheta="<<FirstPulseStartTheta<<
				"    FirstPulseStartPhi="<<FirstPulseStartPhi<<endl;
    
    
    int size=0, nPixCount=0;
    
    for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        fevt::Telescope::ConstPixelIterator pixiter;
        const fevt::Telescope& tel = *teliter;
        const fdet::Camera &detCamera = det::Detector::GetInstance().GetFDetector().GetTelescope(*teliter).GetCamera();
        size = detCamera.GetFADCTraceLength();
        
        
        for (pixiter= tel.PixelsBegin(ComponentSelector::eHasData); pixiter!= tel.PixelsEnd(ComponentSelector::eHasData); ++pixiter){
            
            if (!pixiter->HasRecData()) continue;
            
            nPixCount++;
            
        } //pixels
        
    } // telescopes
    
    
    int eye_number = eye.GetId();
    stringstream s2;
    s2<<eye_number;
    s2>> fEye_number_string;
    
    TimeStamp	timeMy = eye.GetHeader().GetTimeStamp();
    
    stringstream ss;
    ss<<"EV_"<<timeMy.GetGPSSecond()<<"_"<<timeMy.GetGPSNanoSecond()<<"_SDP";
    
    TH2D * SDPPlot = new TH2D("SDP","SDP",300,minphi,maxphi,300,0.,40.);
    
    SDPPlot->SetStats(kFALSE);
    SDPPlot->GetXaxis()->SetTitle("#phi  [deg]");
    SDPPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]");
    SDPPlot->GetZaxis()->SetTitle("Amplitude");
    
    
    // retrieve  pixels with a reconstructed pulse
    
    //fevt::EyeRecData::ConstPixelIterator pixiter;
    //const fevt::EyeRecData& eyerec = eye.GetRecData();
    
    //searching the pixels minimum and maximum pulse
    double min_charge = 999999999.;
    double max_charge = 0.;
    double charge = 0.;
    double adding_charge = 0.;
    int nSlotInterval = 10;
    double pixel_slot[NPIXELS][(int)(size/nSlotInterval)];
    int counter = 0;
    
    
    
    for (teliter = eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        fevt::Telescope::ConstPixelIterator pixiter;
        const fevt::Telescope& tel = *teliter;
        
        counter=0;
        
        for (pixiter= tel.PixelsBegin(ComponentSelector::eHasData); pixiter!= tel.PixelsEnd(ComponentSelector::eHasData); ++pixiter){
            
            charge=0.;
            adding_charge=0.;
            
            const fevt::Pixel & pixel = *pixiter;
            
            if (!pixiter->HasRecData()) continue;
            
            const fevt::PixelRecData & pixelrecdata = pixel.GetRecData();
            const TraceD & trace = pixelrecdata.GetPhotonTrace();
            
            double y[size];
            std::copy(trace.Begin(), trace.End(),y);
            
            for ( int i=0; i<size; i++ ){
                
                if( y[i] < -1000. ) adding_charge = 0.;
                else adding_charge = y[i];
                
                charge = charge + adding_charge;
                
                if( i%nSlotInterval != nSlotInterval-1 ) continue;
                
                else{
                    pixel_slot[counter][((i-nSlotInterval+1)/nSlotInterval)] = charge;
                    //cout<<"counter="<<counter<<
                    //			"   ((i-nSlotInterval+1)/nSlotInterval)="<<((i-nSlotInterval+1)/nSlotInterval)<<
                    //			"   i="<<i<<
                    //			"   nSlotInterval="<<nSlotInterval<<endl;
                    if (charge < min_charge) min_charge = charge;
                    if (charge > max_charge) max_charge = charge;
                    charge=0.;
                }
            }// for y[i]
            counter++;
            
            
        } //pixels
        
        
    } // telescopes
    
    cout<<"nPixels="<<nPixCount<<"  (size/nSlotInterval)="<<(int)(size/nSlotInterval)<<endl;
    cout<<"max_charge="<<max_charge<<endl<<"min_charge="<<min_charge<<endl;
    
    int counter2=0;
    double max_charge2=0., min_charge2=99999999.;
    vector<double> threshold_charges;
    //TH1F *histo = new TH1F();
    
    fOut->cd();
    fOut->mkdir("Fit2");
    
    for ( int slottime = 0; slottime < (int)(size/nSlotInterval); slottime++){
        
        //cout<<(slottime+1)*nSlotInterval<<endl;
        
        counter2=0;
        counter=0;
        vector<double> charges;
        max_charge2=0.;
        min_charge2=99999999.;
        
        //stringstream ss1;
        //ss1<<"fitgaus_"<<(nSlotInterval*(1+slottime));
        //TCanvas *c1 = new TCanvas(ss1.str().c_str(),ss1.str().c_str());
        
        for (teliter = eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
            
            fevt::Telescope::ConstPixelIterator pixiter;
            const fevt::Telescope& tel = *teliter;
            
            for (pixiter= tel.PixelsBegin(ComponentSelector::eHasData); pixiter!= tel.PixelsEnd(ComponentSelector::eHasData); ++pixiter){
                
                const fevt::Pixel & evtpixel = *pixiter;
                
                if (!pixiter->HasRecData()) continue;
                
                const fdet::Pixel & detpixel= Detector::GetInstance().GetFDetector().GetPixel(evtpixel);
                
                Vector dir = detpixel.GetDirection();
                dir.TransformTo(eye2CS);
                
                Vector sdp = eyerec.GetSDP();
                sdp.TransformTo(eye2CS);
                
                double theta = dir.GetTheta(eye2CS)/ utl::degree;
                double phi =   dir.GetPhi(eye2CS) /  utl::degree;
                
                if (change_phi == false ) if (phi<0) phi+= 360.;
                if (change_phi == true ) if (phi > 180) phi -= 360;
                
                
                if ( (90 - theta) <= FirstPulseStartTheta ){
                    charges.push_back(pixel_slot[counter][slottime]);
                    //cout<<pixel_slot[counter][slottime]<<"  ";
                    if( pixel_slot[counter][slottime] <= min_charge2 ) min_charge2 = pixel_slot[counter][slottime];
                    if( pixel_slot[counter][slottime] >= max_charge2 ) max_charge2 = pixel_slot[counter][slottime];
                    counter2++;
                }
                counter++;
                
            } // for pixiter
        }// for teliter
        
        
        
        TH1F *histo = new TH1F("charge","charge",counter2,min_charge2,max_charge2);
        
        for (unsigned int i=0; i<charges.size(); i++) histo->Fill(charges[i]);
        
        TF1 *f1 = new TF1("f1", "landau", min_charge2,max_charge2);
        histo->Fit("f1","Q");
        
        //cout<<f1->GetParameter(1)<<"   "<<f1->GetParameter(2)<<endl;
        
        threshold_charges.push_back(f1->GetParameter(1) + 7*(f1->GetParameter(2)));
        
        //c1->cd();
        //histo->Draw();
        
        //c1->Write();
        
        delete histo;
        
    } // slotime
    
    cout<<"minphi="<<minphi<<"   maxphi="<<maxphi<<endl;
    
    SDPPlot->Fill(210.,30.,0);
    
    
    fOut->cd();
    fOut->cd("Fit2");
    
    
    double maxphi2=0., minphi2=0., thetamaxphi2=0., thetaminphi2=0., charge_sum=0.;
    
    for ( int slottime = 0; slottime < (int)(size/nSlotInterval); slottime++){
        
        vector<double> thetaFittedPixels, phiFittedPixels, chargeFittedPixels, weightFittedPixels; //charges_sum;
        
        charge_sum=0.;
        counter=0;
        
        TCanvas *c2 = new TCanvas("name","title",2);
        //c2->SetCanvasSize(800,800);
        //c2->SetWindowSize(15+c2->GetWw(),30+c2->GetWh());
        
        stringstream ss2;
        ss2.str("");
        if( (nSlotInterval*(1+slottime)) >= 0 &&(nSlotInterval*(1+slottime)) < 10)
            ss2<<ss.str()<<"_AllPixels_SlotTime_000"<<(nSlotInterval*(1+slottime));
        else if( (nSlotInterval*(1+slottime)) >= 10 &&(nSlotInterval*(1+slottime)) < 100)
            ss2<<ss.str()<<"_AllPixels_SlotTime_00"<<(nSlotInterval*(1+slottime));
        else if( (nSlotInterval*(1+slottime)) > 99 && (nSlotInterval*(1+slottime)) < 1000)
            ss2<<ss.str()<<"_AllPixels_SlotTime_0"<<(nSlotInterval*(1+slottime));
        else
            ss2<<ss.str()<<"_AllPixels_SlotTime_"<<(nSlotInterval*(1+slottime));
        
        c2->SetName(ss2.str().c_str());
        c2->SetTitle(ss2.str().c_str());
        c2->SetFrameFillColor(1);
        
        stringstream ss3;
        ss3<<(nSlotInterval*(1+slottime));
        
        SDPPlot->SetTitle(ss3.str().c_str());
        
        c2->cd();
        SDPPlot->GetZaxis()->SetRangeUser(min_charge,max_charge);
        SDPPlot->Draw("ZCOL");
        gPad->Update();
        TPaletteAxis * pal = (TPaletteAxis*)SDPPlot->GetListOfFunctions()->FindObject("palette");
        
        pal->SetLabelColor(1);
        pal->SetLabelFont(62);
        pal->SetLabelOffset(0.005);
        pal->SetLabelSize(0.02);
        pal->SetTitleOffset(0.8);
        pal->SetTitleSize(0.04);
        pal->SetFillColor(1);
        pal->SetFillStyle(1001);
        //cout<<"x1="<<pal->GetX1NDC()<<"  x2="<<pal->GetX2NDC()<<"  y1="<<pal->GetY1NDC()<<"  y2="<<pal->GetY2NDC()<<endl;
        pal->SetX2NDC(0.928);
        
        TIter next(fPixels);
        
        while( TMarker * mark = (TMarker *) next() ) {
            c2->cd();
            mark->Draw("SAME");
        }
        
        fOut2.open("FittingPixels.txt",ios::out);
        
        for (teliter = eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
            
            fevt::Telescope::ConstPixelIterator pixiter;
            const fevt::Telescope& tel = *teliter;
            
            for (pixiter= tel.PixelsBegin(ComponentSelector::eHasData); pixiter!= tel.PixelsEnd(ComponentSelector::eHasData); ++pixiter){
                
                const fevt::Pixel & evtpixel = *pixiter;
                
                if (!pixiter->HasRecData()) continue;
                
                const fdet::Pixel & detpixel= Detector::GetInstance().GetFDetector().GetPixel(evtpixel);
                
                Vector dir = detpixel.GetDirection();
                dir.TransformTo(eye2CS);
                
                Vector sdp = eyerec.GetSDP();
                sdp.TransformTo(eye2CS);
                
                double theta = dir.GetTheta(eye2CS)/ utl::degree;
                double phi =   dir.GetPhi(eye2CS) /  utl::degree;
                
                if (change_phi == false ) if (phi<0) phi+= 360.;
                if (change_phi == true ) if (phi > 180) phi -= 360;
                
                TMarker * mark2 = new TMarker(phi, 90-theta,8);
                
                for(int j = 0 ; j < MaxColors ; j++){
                    
                    if (j!= MaxColors) {
                        if (pixel_slot[counter][slottime] < (min_charge+(max_charge-min_charge)*(j+1)/MaxColors) &&
                            pixel_slot[counter][slottime] >= (min_charge+(max_charge-min_charge)*j/MaxColors) ) {
                            mark2->SetMarkerColor(palette[j]);
                            mark2->SetMarkerSize(1.7+(float)j/(float)MaxColors);
                            c2->cd();
                            mark2->Draw("SAME");
                        }
                        else continue;
                    }
                    else {
                        if( pixel_slot[counter][slottime] <= (min_charge+(max_charge-min_charge)*(j+1)/MaxColors) &&
                           pixel_slot[counter][slottime] >= (min_charge+(max_charge-min_charge)*j/MaxColors) ) {
                            mark2->SetMarkerColor(palette[j]);
                            mark2->SetMarkerSize(1.7+(float)j/(float)MaxColors);
                            c2->cd();
                            mark2->Draw("SAME");
                        }
                        else continue;
                    }
                    
                } // for j
                
                if ( (90 - theta) <= FirstPulseStartTheta && pixel_slot[counter][slottime] >= threshold_charges[slottime]){
                    thetaFittedPixels.push_back(90-theta);
                    phiFittedPixels.push_back(phi);
                    chargeFittedPixels.push_back(TMath::Sqrt(pixel_slot[counter][slottime]));
                    charge_sum += pixel_slot[counter][slottime];
                    fOut2<<90-theta<<"   "<<phi<<"   "<<pixel_slot[counter][slottime]<<endl;
                }
                counter++;
                
            } // for pixiter
        }// for teliter
        
        
        fOut2.close();
        
        
        if ( ((slottime+1)*nSlotInterval) >= 400/*FirstPulseStart*/ ){
            
            //charges_sum.push_back(charge_sum);
            
            for (unsigned int i=0; i<chargeFittedPixels.size(); i++){
                weightFittedPixels.push_back(TMath::Sqrt(charge_sum/chargeFittedPixels[i]));
            }
            
            minphi2=9999.;
            maxphi2=0.;
            thetaminphi2=0.;
            thetamaxphi2=0.;
            
            for (unsigned int i=0; i<phiFittedPixels.size(); i++){
                if (phiFittedPixels[i] <= minphi2) {minphi2 = phiFittedPixels[i]; thetaminphi2 = thetaFittedPixels[i];}
                if (phiFittedPixels[i] >= maxphi2) {maxphi2 = phiFittedPixels[i]; thetamaxphi2 = thetaFittedPixels[i];}
                
            }
            
            
            Double_t x1=FirstPulseStartPhi, y1=FirstPulseStartTheta, x2=maxphi2, y2=thetamaxphi2;
            Double_t r= TMath::Sqrt(TMath::Power(x1-x2,2)+TMath::Power(y1-y2,2));
            
            cout<<(slottime+1)*nSlotInterval<<endl;
            //cout<<"X1=("<<x1<<" ; "<<y1<<")"<<"    X2=("<<x2<<" ; "<<y2<<")"<<endl;
            //cout<<"r="<<r<<endl;
            
            /*
             TCanvas *c3 = new TCanvas("Chi2","Chi2");
             
             Double_t chi2 = 0., min_chi2=99999999., max_chi2=0.;
             
             for (unsigned int j=0; j<thetaFittedPixels.size(); j++){
             chi2 = TMath::Power(thetaFittedPixels[j]- y1 - TMath::Sqrt(TMath::Power(r,2)-TMath::Power(phiFittedPixels[j]-x1,2)),2);
             if(chi2 <= min_chi2) min_chi2 = chi2;
             if(chi2 >= max_chi2) max_chi2 = chi2;
             }
             
             TH1F *histoChi2 = new TH1F("Chi2","Chi2",thetaFittedPixels.size(),min_chi2,max_chi2);
             
             for (unsigned int j=0; j<thetaFittedPixels.size(); j++){
             chi2 = TMath::Power(thetaFittedPixels[j]- y1 - TMath::Sqrt(TMath::Power(r,2)-TMath::Power(phiFittedPixels[j]-x1,2)),2);
             histoChi2->Fill(chi2);
             }
             
             c3->cd();
             histoChi2->Draw();
             c3->Write();
             */
            
            TGraphErrors *graph = new
            TGraphErrors(thetaFittedPixels.size(),&phiFittedPixels.front(),&thetaFittedPixels.front(),&weightFittedPixels.front(),&weightFittedPixels.front());
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(0.75);
            graph->SetMarkerColor(2);
            graph->SetLineColor(2);
            
            
            TF1 *fa = new TF1("fa",FdProcessElves::circumference,minphi2,maxphi2,3);
            fa->SetLineColor(3);
            fa->SetParameters(r,FirstPulseStartPhi,FirstPulseStartTheta);
            //fa->SetParLimits(0,0.,2*r);
            //fa->SetParLimits(1,300.,315.);
            //fa->SetParLimits(2,10.,25.);
            
            fa->FixParameter(1,FirstPulseStartPhi);
            fa->FixParameter(2,FirstPulseStartTheta);
            
            graph->Fit("fa","RQ");
            
            fa->SetParameter(0,fa->GetParameter(0));
            fa->ReleaseParameter(1);
            fa->ReleaseParameter(2);
            fa->FixParameter(0,fa->GetParameter(0));
            
            graph->Fit("fa","RQ");
            
            fa->SetParameter(1,fa->GetParameter(1));
            fa->SetParameter(2,fa->GetParameter(2));
            fa->ReleaseParameter(0);
            
            graph->Fit("fa","R");
            
            c2->cd();
            graph->Draw("PX");
            
            TGraphErrors *graph2 = new TGraphErrors(1);
            
            graph2->SetLineColor(3);
            graph2->SetMarkerColor(3);
            graph2->SetMarkerStyle(8);
            graph2->SetPoint(0,fa->GetParameter(1),fa->GetParameter(2));
            graph2->SetPointError(0,fa->GetParError(1),fa->GetParError(2));
            
            c2->cd();
            graph2->Draw("P");
            
            /*
             
             TMinuit minuit(3);
             minuit.SetFCN(FdSDPFinder::MinuitFitFunc);
             
             minuit.SetPrintLevel(-1);
             
             double vstrt[3];
             double vstp[3];
             double bmin[3];
             double bmax[3];
             
             vstrt[0] = 0.; vstrt[1] = FirstPulseStartPhi; vstrt[2] = FirstPulseStartTheta;
             vstp[0] = 0.1; vstp[1] = 0.01; vstp[2] = 0.01;
             bmin[0] = 0.; bmin[1] = 0.; bmin[2] = 0.;
             bmax[0] = 0.; bmax[1] = 0.; bmax[2] = 0.;
             
             double arglist[10];
             int ierflg = 0;
             
             arglist[0]=1;
             minuit.mnexcm("SET ERR", arglist,1,ierflg);
             
             minuit.mnparm(0,"radius",vstrt[0],vstp[0],bmin[0],bmax[0],ierflg);
             minuit.mnparm(1,"phi0",vstrt[1],vstp[1],bmin[1],bmax[1],ierflg);
             minuit.mnparm(2,"theta0",vstrt[2],vstp[2],bmin[2],bmax[2],ierflg);
             
             arglist[0] = 1;
             minuit.mnexcm("SET STR",arglist,1,ierflg);
             
             minuit.FixParameter(1);
             minuit.FixParameter(2);
             
             //arglist[0] = TMath::Power(10,-9);
             //minuit.mnexcm("SET EPS",arglist,1,ierflg);
             
             arglist[0] = 500;
             minuit.mnexcm("MIGRAD",arglist,1,ierflg);
             */
            
        }
        
        c2->SetCanvasSize(800,800);
        c2->SetWindowSize(15+c2->GetWw(),30+c2->GetWh());
        c2->Write();
        
        //ss3<<"_test";
        //TCanvas *c3 = new TCanvas(ss3.str().c_str(),ss3.str().c_str(),2);
        //fOut->cd();
        //fOut->cd("test");
        //OutTest->cd();
        //c3->Write();
        //delete c3;
        
        bool gif = false;
        if(gif){
            char str1[62];
            Int_t sys_err;
            c2->Print("sdp1.eps");
            sys_err=system("pstopnm -ppm -xborder 0 -yborder 0 -portrait sdp1.eps");
            sys_err=system("ppmtogif sdp1.eps001.ppm > sdp1.gif");
            stringstream ss5;
            ss5<<"/home/epsi/mussa/exotica/FR4/gif/Fit/"
            <<timeMy.GetGPSSecond()<<"/FixedRadius/"<<ss2.str();
            string tit = ss5.str();
            tit += ".gif";
            strcpy (str1,"mv sdp1.gif ");
            strcat (str1,tit.c_str());
            sys_err = system(str1);
        }
        
    } // slotime
    
    
}


void   FdProcessElves::MinuitFitFunc(int  &npar, double *gin, double &f, double *par, int iflag){
    
    ifstream infile;
    
    double phi, theta, charge;
    vector<double> phis, thetas, charges;
    
    infile.open("FittingPixels.txt", ios::in);
    
    while (!infile.eof()) {
        
        infile >>theta>>phi>>charge;
        
        if (infile.eof()) // condizione di end of file (eof) per non ripetere due volte l'ultima riga di dati
            break;
        
        phis.push_back(phi);
        thetas.push_back(theta);
        charges.push_back(charge);
        //cout<<"phi="<<phi<<"    theta="<<theta<<"    charge="<<charge<<endl;
    }
    
    double chisq =0;
    
    double& radius = par[0];
    double& centerPhi = par[1];
    double& centerTheta = par[2];
    
    // iterate over pixels with a pulse
    unsigned int ndof =0;
    double total_charge=0;
    
    for (unsigned int i=0; i<charges.size(); i++){
        
        double tmp = thetas[i] - (centerTheta - TMath::Sqrt(radius*radius - (phis[i] - centerPhi)*(phis[i] - centerPhi)));
        
        chisq += tmp*tmp  * charges[i];
        ndof++;
        total_charge += charges[i];
        
    }// for charges
    
    // normalize chisq
    chisq /= total_charge;
    chisq *= ndof;
    //
    f= chisq;
    
    infile.close();
    
}




void FdProcessElves::Fit3D(const fevt::Eye & eye){
    
    gROOT->SetStyle("Plain");
    
    // Color Palette creation
    int          palette[MaxColors];
    int          index;
    float        r, g, b;
    TColor       *color;
    
    for (int i=0 ; i<MaxColors ; i++) {
        index = 500 + i;
        palette[i] = index;
        
        b = 1.0 * (1.0 - float(MaxColors-i)/float(MaxColors));
        g = 1.0 * (1.0 - float(MaxColors-i)/float(MaxColors));
        r = 1.0 * (1.0 - float(MaxColors-i)/float(MaxColors));
        color = new TColor(index, r, g, b);
        
    }
    
    TColor::SetPalette(MaxColors,palette);
    gStyle->SetNumberContours(MaxColors);
    
    
    fPixels->Delete();
    
    // draw all  the pixels in the telescope onto a Canvas
    // the  ones with a pulse in different color
    
    //  int firstpixelid= fdet::Telescope::GetFirstPixelId();
    //  int lastpixelid = fdet::Telescope::GetLastPixelId();
    
    double minphi = 360.;
    double maxphi = 0;
    bool change_phi = false;
    
    Point eyeposition = det::Detector::GetInstance().GetFDetector().GetEye(eye).GetPosition();
    CoordinateSystemPtr eye2CS = fwk::LocalCoordinateSystem::Create(eyeposition);
    
    // eye2CS is centered in the Eye, with the X axis orientated towards East.
    // So, It's not parallel to the eye Back walll like the standar eye
    // coordinate system 'eyeCS'.
    
    fevt::Eye::ConstTelescopeIterator teliter;
    
    for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        const fevt::Telescope& tel = *teliter;
        
        // note the difference between det::Telescope and fevt::Telescope
        const fdet::Telescope & dettelescope = Detector::GetInstance().GetFDetector().GetTelescope(tel);
        int firstpixelid= dettelescope.GetFirstPixelId();
        int lastpixelid = dettelescope.GetLastPixelId();
        
        for (int i = firstpixelid; i<=lastpixelid; i++){
            const fdet::Pixel & detpixel= dettelescope.GetPixel(i);
            const Vector&  dir = detpixel.GetDirection();
            
            //double theta = dir.GetTheta(eye2CS)/ utl::degree;
            double phi =   dir.GetPhi(eye2CS) /  utl::degree ;  if (phi<0) phi+= 360.;
            
            
            if (phi > maxphi) maxphi = phi;
            if (phi < minphi) minphi = phi;
            
        }//for i
        
    }
    
    
    
    if ( maxphi - minphi > 180 ) {
        change_phi = true;
        maxphi = -99999;
        minphi = 999999;
    }
    
    
    if (change_phi == true ) {
        
        for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
            
            const fevt::Telescope& tel = *teliter;
            
            // note the difference between det::Telescope and fevt::Telescope
            const fdet::Telescope & dettelescope = Detector::GetInstance().GetFDetector().GetTelescope(tel);
            int firstpixelid= dettelescope.GetFirstPixelId();
            int lastpixelid = dettelescope.GetLastPixelId();
            
            
            for (int i = firstpixelid; i<=lastpixelid; i++){
                const fdet::Pixel & detpixel= dettelescope.GetPixel(i);
                const Vector&  dir = detpixel.GetDirection();
                
                //double theta = dir.GetTheta(eye2CS)/ utl::degree;
                double phi =   dir.GetPhi(eye2CS) /  utl::degree ;
                
                if (phi > 180) phi -= 360;
                
                if (phi > maxphi) maxphi = phi;
                if (phi < minphi) minphi = phi;
                
            }
        }
    }
    
    
    for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        const fevt::Telescope& tel = *teliter;
        
        // note the difference between det::Telescope and fevt::Telescope
        const fdet::Telescope & dettelescope = Detector::GetInstance().GetFDetector().GetTelescope(tel);
        int firstpixelid= dettelescope.GetFirstPixelId();
        int lastpixelid = dettelescope.GetLastPixelId();
        
        for (int i = firstpixelid; i<=lastpixelid; i++){
            const fdet::Pixel & detpixel= dettelescope.GetPixel(i);
            const Vector&  dir = detpixel.GetDirection();
            
            double theta = dir.GetTheta(eye2CS)/ utl::degree;
            double phi =   dir.GetPhi(eye2CS) /  utl::degree ;
            
            if (change_phi == false ) if (phi<0) phi+= 360.;
            
            if (change_phi == true ) if (phi > 180) phi -= 360;
            
            TMarker * mark = new TMarker(phi, 90-theta ,4); //empty dot
            mark->SetMarkerColor(15);
            fPixels->Add(mark);
            
        }//for i
    } // for teliter
    
    
    
    
    int FirstPulseStart=9999;
    double FirstPulseStartTheta=0., FirstPulseStartPhi=0.;
    
    fevt::EyeRecData::ConstPixelIterator pixiter;
    const fevt::EyeRecData& eyerec = eye.GetRecData();
    
    for (pixiter= eyerec.PulsedPixelsBegin(); pixiter!=eyerec.PulsedPixelsEnd(); ++pixiter){
        
        const fevt::Pixel & pixel = *pixiter;
        const fdet::Pixel & detpixel= Detector::GetInstance().GetFDetector().GetPixel(pixel);
        
        Vector dir = detpixel.GetDirection();
        dir.TransformTo(eye2CS);
        
        Vector sdp = eyerec.GetSDP();
        sdp.TransformTo(eye2CS);
        
        double theta = dir.GetTheta(eye2CS)/ utl::degree;
        double phi =   dir.GetPhi(eye2CS) /  utl::degree;
        
        if (change_phi == false ) if (phi<0) phi+= 360.;
        if (change_phi == true ) if (phi > 180) phi -= 360;
        
        const fevt::PixelRecData & pixelrecdata = pixel.GetRecData();
        
        
        if ( pixelrecdata.GetPulseStart() <= FirstPulseStart ){
            FirstPulseStartTheta = (90 - theta);
            FirstPulseStartPhi = phi;
            FirstPulseStart = pixelrecdata.GetPulseStart();
        }
        
    }
    
    cout<<"FirstPulseStart="<<FirstPulseStart<<"    FirstPulseStartTheta="<<FirstPulseStartTheta<<
				"    FirstPulseStartPhi="<<FirstPulseStartPhi<<endl;
    
    
    int size=0, nPixCount=0;
    
    for (teliter =  eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        fevt::Telescope::ConstPixelIterator pixiter;
        const fevt::Telescope& tel = *teliter;
        const fdet::Camera &detCamera = det::Detector::GetInstance().GetFDetector().GetTelescope(*teliter).GetCamera();
        size = detCamera.GetFADCTraceLength();
        
        
        for (pixiter= tel.PixelsBegin(ComponentSelector::eHasData); pixiter!= tel.PixelsEnd(ComponentSelector::eHasData); ++pixiter){
            
            if (!pixiter->HasRecData()) continue;
            
            nPixCount++;
            
        } //pixels
        
    } // telescopes
    
    
    int eye_number = eye.GetId();
    stringstream s2; 
    s2<<eye_number;
    s2>> fEye_number_string;
    
    TimeStamp	timeMy = eye.GetHeader().GetTimeStamp();
    
    stringstream ss;
    ss<<"EV_"<<timeMy.GetGPSSecond()<<"_"<<timeMy.GetGPSNanoSecond()<<"_SDP";
    
    TH2D * SDPPlot = new TH2D("SDP","SDP",300,minphi,maxphi,300,0.,32.); 
    
    SDPPlot->SetStats(kFALSE);
    SDPPlot->GetXaxis()->SetTitle("#phi  [deg]"); 
    SDPPlot->GetYaxis()->SetTitle("pixel elevation angle   [deg]"); 
    SDPPlot->GetZaxis()->SetTitle("Amplitude");
    
    
    // retrieve  pixels with a reconstructed pulse
    
    //fevt::EyeRecData::ConstPixelIterator pixiter;
    //const fevt::EyeRecData& eyerec = eye.GetRecData();
    
    //searching the pixels minimum and maximum pulse
    double min_charge = 999999999.;
    double max_charge = 0.;
    double charge = 0.;
    double adding_charge = 0.;
    int nSlotInterval = 10;
    double pixel_slot[NPIXELS][(int)(size/nSlotInterval)];
    double Pixel_X[NPIXELS];
    double Pixel_Y[NPIXELS];
    int counter = 0;
    
    
    TH2F *hQ2D = new TH2F("QcamT","Q in camera at time slot 18",720,-360,360,66,0,33);
    
    for (teliter = eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
        
        fevt::Telescope::ConstPixelIterator pixiter;
        const fevt::Telescope& tel = *teliter;
        
        counter=0;
        
        for (pixiter= tel.PixelsBegin(ComponentSelector::eHasData); pixiter!= tel.PixelsEnd(ComponentSelector::eHasData); ++pixiter){
            
            charge=0.;
            adding_charge=0.;
            
            const fevt::Pixel & pixel = *pixiter;        
            
            if (!pixiter->HasRecData()) continue;
            
            const fevt::PixelRecData & pixelrecdata = pixel.GetRecData();
            const TraceD & trace = pixelrecdata.GetPhotonTrace();	
            
            double y[size];
            std::copy(trace.Begin(), trace.End(),y);
            
            for ( int i=0; i<size; i++ ){
                
                if( y[i] < -1000. ) adding_charge = 0.;
                else adding_charge = y[i];
                
                charge = charge + adding_charge;
                
                if( i%nSlotInterval != nSlotInterval-1 ) continue;
                
                else{
                    pixel_slot[counter][((i-nSlotInterval+1)/nSlotInterval)] = charge;
                    //cout<<"i="<<i<<endl;
                    if (charge < min_charge) min_charge = charge;
                    if (charge > max_charge) max_charge = charge;
                    charge=0.;
                }
            }// for y[i]
            counter++;
            
            
        } //pixels
        
        
    } // telescopes
    
    
    cout<<"max_charge="<<max_charge<<endl<<"min_charge="<<min_charge<<endl;
    
    double mean_charge=0., sigma2_charge=0.;
    int counter2=0;
    //vector<double> means_charge, sigma2s_charge;
    double max_charge2=0., min_charge2=99999999.;
    vector<double> threshold_charges;
    
    fOut->cd();
    fOut->mkdir("Fit");
    //fOut->cd("Fit");
    
    for ( int slottime = 0; slottime < (int)(size/nSlotInterval); slottime++){
        
        //mean_charge=0.;
        //sigma2_charge=0.;
        counter2=0;			
        counter=0;
        vector<double> charges;
        max_charge2=0.;
        min_charge2=99999999.;
        
        //stringstream ss1;
        //ss1<<"fitgaus_"<<(nSlotInterval*(1+slottime));		
        //TCanvas *c1 = new TCanvas(ss1.str().c_str(),ss1.str().c_str());		
        
        for (teliter = eye.TelescopesBegin(ComponentSelector::eHasData); teliter != eye.TelescopesEnd(ComponentSelector::eHasData); ++teliter  ){
            
            fevt::Telescope::ConstPixelIterator pixiter;
            const fevt::Telescope& tel = *teliter;
            
            for (pixiter= tel.PixelsBegin(ComponentSelector::eHasData); pixiter!= tel.PixelsEnd(ComponentSelector::eHasData); ++pixiter){
                
                const fevt::Pixel & evtpixel = *pixiter;
                
                if (!pixiter->HasRecData()) continue;		
                
                const fdet::Pixel & detpixel= Detector::GetInstance().GetFDetector().GetPixel(evtpixel);
                
                Vector dir = detpixel.GetDirection();
                dir.TransformTo(eye2CS);
                
                Vector sdp = eyerec.GetSDP();
                sdp.TransformTo(eye2CS);
                
                double theta = dir.GetTheta(eye2CS)/ utl::degree; 
                double phi =   dir.GetPhi(eye2CS) /  utl::degree; 
                
                if (change_phi == false ) if (phi<0) phi+= 360.; 
                if (change_phi == true ) if (phi > 180) phi -= 360;    
                
                Pixel_X[counter] = phi;
                Pixel_Y[counter] = theta;
                
                if ( (90 - theta) <= FirstPulseStartTheta ){
                    charges.push_back(pixel_slot[counter][slottime]);
                    if( pixel_slot[counter][slottime] <= min_charge2 ) min_charge2 = pixel_slot[counter][slottime];
                    if( pixel_slot[counter][slottime] >= max_charge2 ) max_charge2 = pixel_slot[counter][slottime];
                    counter2++;
                    //mean_charge = mean_charge*(counter2-1)/counter2 + pixel_slot[counter][slottime]/counter2;
                    //sigma2_charge = sigma2_charge*(counter2-1)/counter2 + TMath::Power(pixel_slot[counter][slottime] - mean_charge,2)/counter2;
                }			
                if (slottime >=0 and slottime <nslots/3) {
                    cout << "PIXEL_X_Y_Q_at_slottime " << slottime << " === " <<  Pixel_X[counter];
                    cout             << " " << 90.-Pixel_Y[counter];
                    cout             << " " << pixel_slot[counter][slottime];
                    
                    cout << endl;
                    
                    hQ2D->Fill(Pixel_X[counter],90.-Pixel_Y[counter],pixel_slot[counter][slottime]);
                    
                }
                
                counter++;			
                
            } // for pixiter
        }// for teliter
        
        cout<<(slottime+1)*nSlotInterval<<endl;
        
        TH1F *histo = new TH1F("charge","charge",counter2,min_charge2,max_charge2);		
        
        for (unsigned int i=0; i<charges.size(); i++) histo->Fill(charges[i]);
        
        TF1 *f1 = new TF1("f1", "landau", min_charge2,max_charge2);
        histo->Fit("f1");
        
        //cout<<f1->GetParameter(1)<<"   "<<f1->GetParameter(2)<<endl;
        
        threshold_charges.push_back(f1->GetParameter(1) + 7*(f1->GetParameter(2)));
        
        //c1->cd();
        //histo->Draw();
        
        //c1->Write();
        /*
         if ( ((slottime+1)*nSlotInterval) >= FirstPulseStart ){
         means_charge.push_back(mean_charge);
         sigma2s_charge.push_back(sigma2_charge);
         }
         */
        
        
    } // slotime	
    
    fOut->cd();		
    fOut->mkdir("Slot18");
    fOut->cd("Slot18");
    
    TCanvas* fitCanvas = new TCanvas("fitTest","FitTest");
    fitCanvas->cd();
    hQ2D->Draw("colz");
    fitCanvas->Write();
    
}

double FdProcessElves::RootOf(double a, double b, double c, int i) {
    
    double d = b*b-4*a*c;
    double r = 9.999e99;
    if (d >= 0){ 
        r = (-b-sqrt(d))/(2*a) ;
        if (i == 1) r = (-b+sqrt(d))/(2*a);
    }
    return r;
}


// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "make -C .. FdSDPFinderOG/FdSDPFinder.o  -k"
// End:
