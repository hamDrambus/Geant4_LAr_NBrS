/*  Utility class for non uniform integration ranges.
 *  This is useful when integrated functions have regions of interest
 *  where they have large contribution to integral and thus small integration
 *  step dx is required but computations would be very slow if small dx is
 *  unnecessarily applied to whole domain.
 */

#ifndef INTEGRATION_INTERVAL_H
#define INTEGRATION_INTERVAL_H

#include <iostream>
#include <vector>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctgmath>
#include <float.h>
#include <boost/lexical_cast.hpp>

class IntegrationRange;
class IntegrationInterval
{
protected:
	long double step_;
	long double left_;
	long double right_;
public:
	IntegrationInterval(long double left, long double right, long double step);
	IntegrationRange operator += (const IntegrationRange &r);
	IntegrationRange operator += (const IntegrationInterval &r);
	IntegrationInterval& operator*= (const long double& r);

	friend class IntegrationRange;
	friend IntegrationRange operator+ (const IntegrationRange &l, const IntegrationInterval& r);
	friend IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationRange& r);
	friend IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationInterval& r);
};

class IntegrationRange
{
protected:
	std::vector<IntegrationInterval> arr_;
public:
	IntegrationRange(void);
	IntegrationRange(const IntegrationInterval &inter);

	long int NumOfIndices(void) const;
	long double Value (long int index) const;
	void Trim (long double left, long double right);
	void Print(std::ostream & str);
	long double max(void) const;
	long double min(void) const;

	IntegrationRange& operator += (const IntegrationRange &r);
	IntegrationRange& operator += (const IntegrationInterval &r);

	friend IntegrationRange operator+ (const IntegrationRange &l, const IntegrationRange& r);
	friend IntegrationRange operator+ (const IntegrationRange &l, const IntegrationInterval& r);
	friend IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationRange& r);
	friend IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationInterval& r);
	friend IntegrationRange operator* (const IntegrationRange &l, const double & r);
	friend IntegrationRange operator/ (const IntegrationRange &l, const double & r);
};

IntegrationRange operator+ (const IntegrationRange &l, const IntegrationRange& r);
IntegrationRange operator+ (const IntegrationRange &l, const IntegrationInterval& r);
IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationRange& r);
IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationInterval& r);
IntegrationRange operator* (const IntegrationRange &l, const long double & r);
IntegrationRange operator/ (const IntegrationRange &l, const long double & r);

#endif // INTEGRATION_INTERVAL_H
