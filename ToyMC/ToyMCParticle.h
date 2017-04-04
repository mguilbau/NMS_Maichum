#ifndef MULTICUMULANTS_TOYMCPARTICLE_H
#define MULTICUMULANTS_TOYMCPARTICLE_H

#include <string>

/**
 * Class to define Sets and subsets
 * 
 */
namespace toymc
{

class ToyMCParticle
{
public:
	virtual const char* name() const { return "ToyMCParticle"; }
	virtual const char* classname() const { return "ToyMCParticle"; }
	
	/**
	 * Constructor - 
	 */
	ToyMCParticle()
        {
           set(0.,0.,0.);
        }

	~ToyMCParticle(){}

        void set(double pt, double eta, double phi)
        {
          this->_phi = phi;
          this->_pt  = pt;
          this->_eta = eta;
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
        double _phi;
        double _pt;
        double _eta;
};
}
#endif

