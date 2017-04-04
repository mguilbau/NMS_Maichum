#ifndef MULTICUMULANTS_TOYMCGENERATOR_H
#define MULTICUMULANTS_TOYMCGENERATOR_H

#include <string>


// ROOT
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TObject.h>
#include <TRandom3.h>

#include "ToyMC/ToyMCEvent.h"
#include "ToyMC/ToyMCParticle.h"

#include "ToyMC/BranchWriter.h"
#include "ToyMC/TClonesArrayWriter.h"



/**
 * Class to define Sets and subsets
 * 
 */
namespace toymc
{

class ToyMCGenerator : public TObject
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

		commonSetup();

		makeTree();
	}

	ToyMCGenerator(int nevts)
	{
		this->_nevts = nevts;
		commonSetup();

		makeTree();
	}

	~ToyMCGenerator()
	{
		this->_tree->Write();
		this->_output->Close();
	}

	void commonSetup(){
		//To be move in a config class
		this->_multSpectrum = TF1("multSpectrum","exp(15.4-x*0.00243)",0,500);
		this->_ptSpectrum   = TF1("ptSpectrum","x*exp(-1*(sqrt(pow(0.13957,2)+pow(x,2)))/(0.44))",0.,10.);
		this->_etaSpectrum  = TF1("etaSpectrum","gaus",-10.,10.);
		this->_phiSpectrum  = TF1("multSpectrum","exp(0)",0.,2.*TMath::Pi());

		
	}

	void makeTree(){

		// TODO: config output filename
		this->_output = new TFile( "ToyMC.root", "RECREATE" );
		// if you open the file first then the TTree can do buffered outputs

		this->_tree = new TTree("ToyMCTree", "Toy MC Tree");
		this->_eventWriter.createBranch( this->_tree, "Event" );
		this->_plcsWriter.createBranch( this->_tree, "Particles" );
	}

	void generate( int nevts = -1)
	{
		if ( nevts > 0 )
			this->_nevts = nevts;

		for(int ievt = 0; ievt<this->_nevts; ++ievt)
		{
			int mult = static_cast<int>(this->_multSpectrum.GetRandom());
			// mcEvt->set(mult);

			this->_event.set( gRandom->Uniform( -5, 5 ), gRandom->Uniform( 0, 40 ) );
			this->_eventWriter.set( this->_event );

			this->_plcsWriter.reset();
			for(int ipart = 0; ipart < mult; ++ipart)
			{
				// make sure we start fresh
				// reusing the smae object is a huge performance boost inside loops
				this->_plc.reset();

				double pt  = this->_ptSpectrum.GetRandom();
				double eta = this->_etaSpectrum.GetRandom();
				double phi = this->_phiSpectrum.GetRandom();
				
				this->_plc.set(pt, eta, phi);
				this->_plcsWriter.add( this->_plc );
			}
			//Fill TTree
			this->_tree->Fill(); 
		}
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
	TFile *_output;
	TTree* _tree;
	ToyMCEvent _event;
	ToyMCParticle _plc;

	#ifndef __CINT__
	BranchWriter<ToyMCEvent> _eventWriter;
	TClonesArrayWriter<ToyMCParticle> _plcsWriter;
	#endif

	TF1 _multSpectrum;
	TF1 _ptSpectrum;
	TF1 _etaSpectrum;
	TF1 _phiSpectrum;

	ClassDef( ToyMCGenerator, 1 )
};
}	// toymc NAMESPACE
#endif
