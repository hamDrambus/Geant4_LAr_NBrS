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
  virtual void Initialize(void) override;
protected:
  // Integration intervals are selected manually after looking at integrated functions behavior.
  // In case of another materials, input XSs or domains these may need to be changed.
  virtual IntegrationRange GetIntervalElectronDistributions(double field) const override;
  virtual IntegrationRange GetIntervalFDistributions(double field) const override;
};

#endif // ArgonPropertiesTables_h
