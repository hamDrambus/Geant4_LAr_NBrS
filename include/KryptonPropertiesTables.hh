#ifndef KryptonPropertiesTables_h
#define KryptonPropertiesTables_h

#include <string>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "utilities/FunctionTable.hh"
#include "utilities/GlobalUtilities.hh"
#include "utilities/IntegrationInterval.hh"
#include "utilities/PolynomialFit.hh"
#include "GlobalParameters.hh"
#include "MediumPropertiesTables.hh"

class KryptonPropertiesTables : public MediumPropertiesTables {
public:
	KryptonPropertiesTables() { classname = "KryptonPropertiesTables"; }
	virtual void InitializeFields(void) override;
};

#endif // KryptonPropertiesTables_h
