#include <chrono>
#include <thread>

#include "KryptonPropertiesTables.hh"

void KryptonPropertiesTables::InitializeFields()
{
	const double kVcm = kilovolt / cm;
	double sfc = gPars::medium_props.field_step_factor; // Step factor to test results dependencies on step size.
	interval_field_NBrS = IntegrationInterval(1000 * kVcm, 1300.0 * kVcm, sfc * 20 * kVcm);
	interval_field_NBrS += IntegrationInterval(200 * kVcm, 1000 * kVcm, sfc * 10 * kVcm);
	interval_field_NBrS += IntegrationInterval(80 * kVcm, 200 * kVcm, sfc * 5 * kVcm);
	interval_field_NBrS += IntegrationInterval(10.0 * kVcm, 80.0 * kVcm, sfc * 1 * kVcm);
	interval_field_NBrS += IntegrationInterval(3.0 * kVcm, 12.0 * kVcm, sfc * 0.2 * kVcm);

	interval_field_drift = IntegrationInterval(0.01 * kVcm, 1.0 * kVcm, sfc * 0.02 * kVcm);
	interval_field_drift += IntegrationInterval(1 * kVcm, 1.5 * kVcm, sfc * 0.05 * kVcm);
	interval_field_drift += IntegrationInterval(1.5 * kVcm, 3.0 * kVcm, sfc * 0.1 * kVcm);
	interval_field_drift += interval_field_NBrS;

	low_high_threshold_field = interval_field_NBrS.min();
}
