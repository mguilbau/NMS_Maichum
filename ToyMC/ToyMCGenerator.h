#ifndef MULTICUMULANTS_TOYMCGENERATOR_H
#define MULTICUMULANTS_TOYMCGENERATOR_H

#include <string>

#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>

#include "ToyMC/ToyMCEvent.h"

/**
 * Class to define Sets and subsets
 * 
 */
namespace toymc
{

class ToyMCGenerator
{
public:
	virtual const char* name() const { return "ToyMCGenerator"; }
	virtual const char* classname() const { return "ToyMCGenerator"; }
	
	/**
	 * Constructor - 
	 */
	ToyMCGenerator()
        {
           this->_nevts = 100;
           
           //To be move in a config class
           this->_multSpectrum = TF1("multSpectrum","exp(15.4-x*0.00243)",0,500);
           this->_ptSpectrum   = TF1("ptSpectrum","x*exp(-1*(sqrt(pow(0.13957,2)+pow(x,2)))/(0.44))",0.,10.);
           this->_etaSpectrum  = TF1("etaSpectrum","gaus",-10.,10.);
           this->_phiSpectrum  = TF1("multSpectrum","exp(0)",0.,2.*TMath::Pi());

           this->_events = new TTree("eventTree", "eventTree");
        }

	ToyMCGenerator(int nevts)
        {
           this->_nevts = nevts;

           //To be move in a config class
           this->_multSpectrum = TF1("multSpectrum","exp(15.4-x*0.00243)",0,500);
           this->_ptSpectrum   = TF1("ptSpectrum","x*exp(-1*(sqrt(pow(0.13957,2)+pow(x,2)))/(0.44))",0.,10.);
           this->_etaSpectrum  = TF1("etaSpectrum","gaus",-10.,10.);
           this->_phiSpectrum  = TF1("multSpectrum","exp(0)",0.,2.*TMath::Pi());

           this->_events = new TTree("eventTree", "eventTree");
        }

	~ToyMCGenerator()
        {
           delete _events;
        }

        void generate()
        {

           ToyMCEvent* mcEvt = new ToyMCEvent;
           this->_events->Branch("Event", "ToyMCEvent", &mcEvt, 32000, 3);
           for(int ievt = 0; ievt<this->_nevts; ++ievt)
           {
             int mult = static_cast<int>(this->_multSpectrum.GetRandom());
             mcEvt->set(mult);

             for(int ipart = 0; ipart < mult; ++ipart)
             {
                double pt  = this->_ptSpectrum.GetRandom();
                double eta = this->_etaSpectrum.GetRandom();
                double phi = this->_phiSpectrum.GetRandom();
                
                ToyMCParticle mcPart;
                mcPart.set(pt, eta, phi);
                mcEvt->add(mcPart);
             }
             //Fill TTree
             this->_events->Fill(); 
           }

           //Write TTree in TFile
           TFile* fout = new TFile("ToyMCTTree.root", "RECREATE");
           this->_events->Write();
           fout->Close();
           delete fout;

           delete mcEvt;
        }
	
	/**
	 * builds the string representation of the object
	 *
	 * @return String representation of this
	 */
	std::string toString(){
		std::string s = "";
		s += this->classname(); 
		return s; 
	}

protected:
        int _nevts;
        TTree* _events;
        TF1 _multSpectrum;
        TF1 _ptSpectrum;
        TF1 _etaSpectrum;
        TF1 _phiSpectrum;
};
}
#endif
