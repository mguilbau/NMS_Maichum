#ifndef MULTICUMULANTS_TOYMCFLOWDISTGENERATOR_H
#define MULTICUMULANTS_TOYMCFLOWDISTGENERATOR_H

#ifndef __CINT__
	#include "vendor/loguru/loguru.hpp"
#endif

#include "ToyMCDistGenerator.h"

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include <string>

/**
 * Class to define Sets and subsets
 * 
 */
namespace toymc
{

class ToyMCFlowDistGenerator : public ToyMCDistGenerator
{
public:
	virtual const char* name() const { return "ToyMCFlowDistGenerator"; }
	virtual const char* classname() const { return "ToyMCFlowDistGenerator"; }
	
	/**
	 * Constructor - 
	 */
	ToyMCFlowDistGenerator( )
	{
                LOG_SCOPE_FUNCTION( INFO );
                formula_ = "x/TMath::Power([0],2)*TMath::Exp(-1.*(x*x + [1]*[1])/2./TMath::Power([0],2))*TMath::BesselI0(x*[1]/TMath::Power([0],2))";
	}

	ToyMCFlowDistGenerator( int nbins, double min, double max, std::string title )
	{
                LOG_SCOPE_FUNCTION( INFO );
                formula_ = "x/TMath::Power([0],2)*TMath::Exp(-1.*(x*x + [1]*[1])/2./TMath::Power([0],2))*TMath::BesselI0(x*[1]/TMath::Power([0],2))";
	}

        ~ToyMCFlowDistGenerator()
        {
                LOG_SCOPE_FUNCTION( 9 );
        }

private:
        std::string formula_;

};

}

#endif
