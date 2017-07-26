#ifndef _DoubleELVESAnalysis_DoubleELVESAnalysis_h_
#define _DoubleELVESAnalysis_DoubleELVESAnalysis_h_

/**
 * \file
 * \author Kevin-Druis Merenda
 * \date 17 May 2017
 */

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include "TH1F.h"
#include <vector>
/*
 * Avoid using using namespace declarations in your headers,
 * doing so makes all symbols from each namespace visible
 * to the client which includes your header.
 */

namespace DoubleELVESAnalysis {

  class DoubleELVESAnalysis : 
      public boost::noncopyable,
      public fwk::VModule {
  public:
    DoubleELVESAnalysis();
    ~DoubleELVESAnalysis();
    VModule::ResultFlag Init();
    VModule::ResultFlag Run(evt::Event& e);
    VModule::ResultFlag Finish();
    VModule::ResultFlag GlueTrace( int,TH1F*, int, int);
    double CalcMHWScore(std::vector<double>);

  private:
    // Declare down here your data following this conventions:
    // fFieldName  : members.
    // fgFieldName : static data members.
    
    const static int fNPixels = 440;
    const static int fNTels = 6;
    const static int fNRows = 22;
    const static int fNColumns = 20;
    TH1F* fhRawPixel[fNTels][fNPixels];//added a dimension for sim vs data

    TH1F* hdeltaTMedian;
    TH1F* hdeltaTMean;
    TH1F* hdeltaTMeanCenter;
    int fNBinTrace;
    int fMaxBinTrace;
    Int_t fNBinPage;
    int fMinBinTrace;
    int fEventCounter;
    int fDataEventCounter;


    // Declare down here your private functions like this:

    //
    // /// Brief description.
    // /*! Detailed description that can span several 
    //     lines.
    // */
    // type foo();

    // This goes at the end.
    REGISTER_MODULE("DoubleELVESAnalysis",DoubleELVESAnalysis);
  };
}

#endif // _DoubleELVESAnalysis_DoubleELVESAnalysis_h_
