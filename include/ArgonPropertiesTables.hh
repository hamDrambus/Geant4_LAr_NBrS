#ifndef ArgonPropertiesTables_h
#define ArgonPropertiesTables_h

#include <string>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "utilities/FunctionTable.hh"
#include "utilities/GlobalUtilities.hh"
#include "utilities/IntegrationInterval.hh"
#include "utilities/PolynomialFit.hh"
#include "GlobalParameters.hh"
#include "MediumPropertiesTables.hh"

class ArgonPropertiesTables : public MediumPropertiesTables {
public:
  ArgonPropertiesTables() { classname = "ArgonPropertiesTables"; }
  virtual void InitializeFields(void) override;
};

#endif // ArgonPropertiesTables_h
