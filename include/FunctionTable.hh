#ifndef FUNCTION_TABLE_H
#define FUNCTION_TABLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include "PolynomialFit.hh"

struct ftable_point {
	double X;
	double Y;
	double VAL;
};

//find_E_indexes_by_value works only for monotone rising function!
//F(E1,Ey)>F(E2,Ey) if E1>E2
class FunctionTable {
protected:
	std::deque<DataVector> ys_; //DataVector is supposed to be sorted by its xs (y for FunctionTable). Search of Y (getY()) by value works only when DataVector's ys are monotonically increasing.
	std::deque<double> xs_; //data is supposed to be sorted by xs;
	boost::optional<std::pair<std::size_t, std::size_t>> getX_indices(double X) const;
	//std::pair<long int, long int> find_E_indexes (double E, std::size_t Ey_index) const;
	//std::pair<long int, long int> find_E_indexes_by_value (double val, std::size_t Ey_index) const;
public:
	FunctionTable(void);
	~FunctionTable();

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

	void push (double X, double Y, double val);

	double operator()(double X, double Y) const;
	double getY(double X, double val) const;

	void read(std::ifstream& str);
	void write(std::string fname) const;
	void write (std::ofstream& str) const;
	bool is_empty(void) const {
		return xs_.empty();
	}
};

#endif
