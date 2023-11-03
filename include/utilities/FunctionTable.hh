#ifndef FUNCTION_TABLE_H
#define FUNCTION_TABLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>

#include "PolynomialFit.hh"

//TODO: add profiling function along x and y (from DataTableExtractor' main.cc)

//find_E_indexes_by_value works only for monotone rising function!
//F(E1,Ey)>F(E2,Ey) if E1>E2
class FunctionTable {
protected:
	std::deque<DataVector> ys_; //DataVector is supposed to be sorted by its xs (y for FunctionTable). Search of Y (getY()) by value works only when DataVector's ys are monotonically increasing.
	std::deque<double> xs_; //data is supposed to be sorted by xs;
	//std::pair<long int, long int> find_E_indexes (double E, std::size_t Ey_index) const;
	//std::pair<long int, long int> find_E_indexes_by_value (double val, std::size_t Ey_index) const;
public:
	FunctionTable(void);
	~FunctionTable();

	boost::optional<std::pair<std::size_t, std::size_t>> getX_indices(double X) const;
	void resize(std::size_t n) {
		xs_.resize(n);
		ys_.resize(n);
	}
	void clear(void) {
		xs_.clear();
		ys_.clear();
	}
	std::size_t size(void) const {
		return xs_.size();
	}
	DataVector& operator[](std::size_t X_index) {
		return ys_[X_index];
	}
	const DataVector& getY_data(std::size_t X_index) const {
		return ys_[X_index];
	}
	double getX(std::size_t X_index) const {
	  return xs_[X_index];
	}

	void push (double X, double Y, double val);
	void push (double X, const DataVector& vals);
	void scale (double x_scale, double y_scale, double z_scale) {
    for (auto &x : xs_)
      x *= x_scale;
    for (auto &yz : ys_)
      yz.scaleXY(y_scale, z_scale);
	}

	double operator()(double X, double Y) const;
	double getY(double X, double val) const;

	void read(std::ifstream& str);
	void write(std::string fname) const;
	void write (std::ofstream& str) const;
	bool is_empty(void) const {
		return xs_.empty();
	}


	boost::optional<std::pair<double, double>> getRangeAtX(double X) const;
	boost::optional<std::size_t> getSizeAtX(double X) const;
};

#endif // FUNCTION_TABLE_H
