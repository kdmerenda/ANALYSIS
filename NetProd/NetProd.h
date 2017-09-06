#ifndef _NetProd_NetProd_h_
#define _NetProd_NetProd_h_

/**
 * \file
 * \author Kevin-Druis Merenda
 * \date 10 Aug 2017
 */

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"

/*
 * Avoid using using namespace declarations in your headers,
 * doing so makes all symbols from each namespace visible
 * to the client which includes your header.
 */

namespace NetProd {

  class NetProd : 
      public boost::noncopyable,
      public fwk::VModule {
  public:
    NetProd();
    ~NetProd();
    VModule::ResultFlag Init();
    VModule::ResultFlag Run(evt::Event& e);
    VModule::ResultFlag Finish();
    VModule::ResultFlag GlueTrace(int, TH1F*, int, int);
    VModule::ResultFlag DumpTraces();
    std::string exec(const char*);
  private:

    
    int fNPages;
    const static int fNPixels = 440;
    const static int fNTels = 6;
    const static int fNRows = 22;
    const static int fNColumns = 20;
    TH1F* fhRawPixel[fNTels][fNPixels];
    TH1F* fhRawPixelMerge[fNTels*fNPixels];
    TString hnameBase[fNTels][fNPixels];
    TH1F* hdeltaTMedian;
    TH1F* hdeltaTMean;
    TH1F* hdeltaTMeanCenter;
    int fBinning;
    int fNBinTrace;
    int fMaxBinTrace;
    Int_t fNBinPage;
    int fMinBinTrace;
    int fEventCounter;
    int fEventTracker;
    int fDataEventCounter;
    TString* fOutFileName;
    double fEventTime;
    double fFirstEventTime;
    double fFirstEventTracker;
    bool fData;
    TString fFirstEventId;
    REGISTER_MODULE("NetProd",NetProd);
  };
}

#endif // _NetProd_NetProd_h_
