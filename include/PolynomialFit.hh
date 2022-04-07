#ifndef POLYNOMIAL_FIT_H
#define POLYNOMIAL_FIT_H

#include <iostream>
#include <fstream>
#include <string>
#include <cfloat>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>

#include "GlobalUtilities.hh"

//TVectorD parameters are [0]+[1]*x+[2]*x^2+...
class PolynomialFit {
protected:
	std::size_t _order;
public:
	PolynomialFit(std::size_t order);
	~PolynomialFit();
	void setOrder(std::size_t n) {
		_order = n;
	}
	std::size_t getOrder(void) const {
		return _order;
	}

	//in_x0 - relative to what point carry out fit. Automatic value is set if boost::none is passed.
	std::vector<double> operator ()(const std::vector<std::pair<double, double>> &vals_in, boost::optional<double> &in_x0) const;
	//Fit only part of a vector. offset+N_points-1 must in the range of the vector
	std::vector<double> operator ()(const std::vector<std::pair<double, double>> &vals_in, int offset, int N_points, boost::optional<double> &in_x0) const;
};

//Wraps PolynomialFit: stores raw data, N points used in every fit and last region (cache_n_from, cache_n_to) in which fit/interpolation took place.
//The latter is required for optimization, because polynomial coefficients are not updated unless necessary (x moved from previous region)
class DataVector {
	//TODO: add thread-local cache for getX_indices and getY_indices
protected:
	std::vector<std::pair<double, double>> xys;
	PolynomialFit fitter;
	std::size_t N_used;

	bool use_left, use_right;
	boost::optional<double> left_value, right_value;
public:
	DataVector(std::size_t fit_order = 1, std::size_t N_used = 2);
	DataVector(std::vector<double> &xx, std::vector<double> &yy, std::size_t fit_order, std::size_t N_used);
	DataVector(std::string data_file);
	virtual ~DataVector();

	void initialize(std::vector<double> &xx, std::vector<double> &yy, std::size_t fit_order, std::size_t N_used);

	void setOrder(std::size_t ord) {
		fitter.setOrder(ord);
	}
	std::size_t getOrder(void) const {
		return fitter.getOrder();
	}
	void setNused(std::size_t N) {
		N_used = std::max((std::size_t)1, N);
	}
	std::size_t getNused(void) const {
		return N_used;
	}
	//precedence goes to use_left-/right-most methods.
	void use_leftmost(bool use) {
		use_left = use;
	}
	void use_rightmost(bool use) {
		use_right = use;
	}
	void set_leftmost(double val) {
		left_value = val;
	}
	void unset_leftmost(void) {
		left_value = boost::none;
	}
	void set_rightmost(double val) {
		right_value = val;
	}
	void unset_rightmost(void) {
		right_value = boost::none;
	}
	void set_out_value(double val) {
		set_leftmost(val);
		set_rightmost(val);
	}
	void unset_out_value(void) {
		unset_leftmost();
		unset_rightmost();
	}
	void erase(std::size_t n) {
		xys.erase(xys.begin() + n);
	}
	void erase() {
		xys.clear();
	}
	void clear(void) {
		xys.clear();
	}
	std::size_t size(void) const {
		return xys.size();
	}
	double getX(std::size_t n) const {
		return xys[n].first;
	}
	double getY(std::size_t n) const {
		return xys[n].second;
	}
	void scaleX(double factor) {
	  for (auto i : xys)
	    i.first *= factor;
	}
	void scaleY(double factor) {
	  for (auto i : xys)
      i.second *= factor;
	}
	void scaleXY(double factorX, double factorY) {
	  scaleX(factorX);
	  scaleY(factorY);
	}
	std::vector<double> get_Xs(double scaling = 1.0) const;
	std::vector<double> get_Ys(double scaling = 1.0) const;
	std::pair<double, double> getXY(std::size_t n) const {
		return xys[n];
	}
	std::pair<double, double>& operator[](std::size_t ind) {
		return xys[ind];
	}
	void resize(std::size_t sz) {
		xys.resize(sz);
	}
	bool isValid(void) const {
		return (N_used > fitter.getOrder()) && N_used <= size();
	}

	double operator()(double X_point, boost::optional<double> x0 = boost::none) const; //x0 = point is recommended to use. At least x0 must be close to point, or there will be large errors otherwise
	void insert(double x, double y);
	void push_back (double x, double y);

	//save/load full state except cache from file
	void read(std::ifstream& str);
	void write(std::string fname, std::string comment = "") const;
	void write(std::ofstream& str, std::string comment = "") const;
	//Warning! These functions have defined behaviour only when X/Y values are sorted in the ascending order.
	//DataVector's X values are supposed to be always sorted, but there is no guarantee about Ys;
	boost::optional<std::pair<std::size_t, std::size_t>> getX_indices(double X_point) const; //[n_min, n_max] are used, not [n_min, n_max).N_used == n_max - n_min + 1 >= order + 1
	boost::optional<std::pair<std::size_t, std::size_t>> getY_indices(double Y_point) const; //[n_min, n_max] are used, not [n_min, n_max).N_used == n_max - n_min + 1 >= order + 1
protected:
	double polynomial_value(double x, double x0, const std::vector<double>& coefs) const;
public:
	friend class FunctionTable;
};

struct pdf_data {
	double x;
	double pdf;//===y
	double cdf;//===Int(-inf, x){y(x)dx}
};

//TODO: create common parent for PDF_routine and DataVector
//does not support indefinite domain.
class PDF_routine { //Probability Density Function routine - load not normalized distribution from file,
	// construct cumulative DF and use it to generate values.
	//PDF's order is fixed to 1 and Nused to 2 (linear interpolation).
protected:
	std::vector<pdf_data> vals;
	bool cdf_ready;
	void normalize_cdf(void);
public:
	PDF_routine();
	PDF_routine(std::vector<double> &pdf_xx, std::vector<double> &pdf_yy);
	~PDF_routine();
	void read(std::ifstream& str);
	bool read(std::string fname);
	void write(std::string fname, std::string comment = "") const;
	void write(std::ofstream& str, std::string comment = "") const;
	double operator()(double Rand) const; //generates random number according to pdf. Has no internal random engine
	void pdf_to_cdf(void);
	//void initialize(std::vector<double> &xx, std::vector<double> &yy);

	void erase(std::size_t n, bool recalculate_cdf = true) {
		vals.erase(vals.begin() + n);
		cdf_ready = false;
		if (recalculate_cdf)
			pdf_to_cdf();
	}
	void erase() {
		cdf_ready = false;
		vals.clear();
	}
	void clear(void) {
		cdf_ready = false;
		vals.clear();
	}
	std::size_t size(void) const {
		return vals.size();
	}
	double getX(std::size_t n) const {
		return vals[n].x;
	}
	double getPDF(std::size_t n) const {
		return vals[n].pdf;
	}
	double getCDF(std::size_t n) const {
		return vals[n].cdf;
	}
	pdf_data get(std::size_t n) const {
		return vals[n];
	}
	pdf_data& operator[](std::size_t ind) {
		cdf_ready = false;
		return vals[ind];
	}
	void resize(std::size_t sz) {
		cdf_ready = false;
		vals.resize(sz);
	}
	bool isValid(void) const {
		return (1 <= size()) && cdf_ready;
	}

	void insert(double x, double y); //Preserves sorting, recalculates pdf, quite expensive when insering not to the end
	void push_back(double x, double y, bool recalculate_cdf = true); //Faster version, not checking whether sorting is preserved.
	//void set_back(double x, double y, std::size_t ind);
protected:
	//Warning! These functions have defined behaviour only when X/Y values are sorted in the ascending order.
	boost::optional<std::pair<std::size_t, std::size_t>> getX_indices(double X_point) const; //[n_min, n_max] are used, not [n_min, n_max).N_used == n_max - n_min + 1 >= order + 1
	boost::optional<std::pair<std::size_t, std::size_t>> getY_indices(double Y_point) const; //[n_min, n_max] are used, not [n_min, n_max).N_used == n_max - n_min + 1 >= order + 1
};

#endif
