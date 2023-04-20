#ifndef XenonPropertiesTables_h
#define XenonPropertiesTables_h

#include <string>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "utilities/FunctionTable.hh"
#include "utilities/GlobalUtilities.hh"
#include "utilities/IntegrationInterval.hh"
#include "utilities/PolynomialFit.hh"
#include "GlobalParameters.hh"
#include "MediumPropertiesTables.hh"

class XenonPropertiesTables : public MediumPropertiesTables {
public:
	XenonPropertiesTables() { classname = "XenonPropertiesTables"; }
  virtual void InitializeFields(void) override;
};

#endif // XenonPropertiesTables_h
