#include <utilities/FunctionTable.hh>

FunctionTable::FunctionTable(void)
{}

FunctionTable::~FunctionTable()
{}

boost::optional<std::pair<double, double>> FunctionTable::getRangeAtX(double X) const
{
	boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(X);
	if (!inds) {
		return boost::none;
	}
	const DataVector *vec1 = &getY_data(inds->first);
	const DataVector *vec2 = &getY_data(inds->second);
	if (!vec1->size() && !vec2->size())
		return boost::none;
	double y_min, y_max;
	if (vec1->size() && vec2->size()) {
		y_min = std::min(vec1->getX(0), vec2->getX(0));
		y_max = std::max(vec1->getX(vec1->size()-1), vec2->getX(vec2->size()-1));
		return std::pair<double, double>(y_min, y_max);
	}
	if (!vec1->size())
		vec1 = vec2;
	y_min = vec1->getX(0);
	y_max = vec1->getX(vec1->size()-1);
	return std::pair<double, double>(y_min, y_max);
}

boost::optional<std::size_t> FunctionTable::getSizeAtX(double X) const
{
	boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(X);
	if (!inds) {
		return boost::none;
	}
	const DataVector *vec1 = &getY_data(inds->first);
	const DataVector *vec2 = &getY_data(inds->second);
	if (!vec1->size() && !vec2->size())
		return boost::none;
	if (vec1->size() && vec2->size()) {
		return std::max(vec1->size(), vec2->size());
	}
	if (!vec1->size())
		vec1 = vec2;
	return vec1->size();
}

boost::optional<std::pair<std::size_t, std::size_t>> FunctionTable::getX_indices(double X) const
{
	return find_bounds(xs_, X);
}

double FunctionTable::operator ()(double X, double Y) const
{
	boost::optional<std::pair<std::size_t, std::size_t>> X_inds = getX_indices(X);
	if (!X_inds)
		return 0;
	if (X_inds->second >= xs_.size())
		X_inds->second = X_inds->first;
	if (X_inds->first == X_inds->second) {
		boost::optional<std::pair<std::size_t, std::size_t>> Y_inds = ys_[X_inds->first].getX_indices(Y);
		if (!Y_inds)
			return 0;
		if (Y_inds->first == Y_inds->second) {
			return ys_[X_inds->first].getY(Y_inds->first);
		} else {
			return ys_[X_inds->first](Y);
		}
	} else {
		boost::optional<std::pair<std::size_t, std::size_t>> Y_inds_left = ys_[X_inds->first].getX_indices(Y);
		boost::optional<std::pair<std::size_t, std::size_t>> Y_inds_right = ys_[X_inds->second].getX_indices(Y);
		if (Y_inds_left)
			if (Y_inds_left->second >= ys_[X_inds->first].size())
				Y_inds_left->second = Y_inds_left->first;
		if (Y_inds_right)
			if (Y_inds_right->second >= ys_[X_inds->second].size())
				Y_inds_right->second = Y_inds_right->first;

		if ((!Y_inds_right)&&(!Y_inds_left))
			return 0;
		if (!Y_inds_right) {
			if (Y_inds_left->first == Y_inds_left->second) {
				return ys_[X_inds->first].getY(Y_inds_left->first);
			} else {
				double z0 = ys_[X_inds->first].getY(Y_inds_left->first), z1 = ys_[X_inds->first].getY(Y_inds_left->second);
				double y0 = ys_[X_inds->first].getX(Y_inds_left->first), y1 = ys_[X_inds->first].getX(Y_inds_left->second);
				return z0 + (z1 - z0) * (Y - y0) / (y1 - y0);
			}
		}
		if (!Y_inds_left) {
			if (Y_inds_right->first == Y_inds_right->second) {
				return ys_[X_inds->first].getY(Y_inds_right->first);
			} else {
				double z0 = ys_[X_inds->first].getY(Y_inds_right->first), z1 = ys_[X_inds->first].getY(Y_inds_right->second);
				double y0 = ys_[X_inds->first].getX(Y_inds_right->first), y1 = ys_[X_inds->first].getX(Y_inds_right->second);
				return z0 + (z1 - z0) * (Y - y0) / (y1 - y0);
			}
		}
		if ((Y_inds_left->first == Y_inds_left->second)&&(Y_inds_right->first == Y_inds_right->second)) {
			double z0 = ys_[X_inds->first].getY(Y_inds_left->first), z1 = ys_[X_inds->second].getY(Y_inds_right->first);
			double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
			return z0 + (z1 - z0) * (X-x0) / (x1-x0);
		}
		if (Y_inds_left->first == Y_inds_left->second) {//three point linear interpolation
			double z0_0 = ys_[X_inds->first].getY(Y_inds_left->first), z1_0 = ys_[X_inds->second].getY(Y_inds_right->first), z1_1 = ys_[X_inds->second].getY(Y_inds_right->second);
			double y0_0 = ys_[X_inds->first].getX(Y_inds_left->first), y1_0 = ys_[X_inds->second].getX(Y_inds_right->first), y1_1 = ys_[X_inds->second].getX(Y_inds_right->second);
			double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
			//first interpolate along X axis, then Y axis
			double y0 = y0_0 + (X-x0)*(y1_0-y0_0)/(x1-x0);
			double y1 = y0_0 + (X-x0)*(y1_1-y0_0)/(x1-x0);
			double z0 = z0_0 + (X-x0)*(z1_0-z0_0)/(x1-x0);
			double z1 = z0_0 + (X-x0)*(z1_1-z0_0)/(x1-x0);
			return z0 + (z1-z0)*(Y-y0)/(y1-y0);
		}
		if (Y_inds_right->first == Y_inds_right->second) {//three point linear interpolation
			double z0_0 = ys_[X_inds->first].getY(Y_inds_left->first), z1_0 = ys_[X_inds->second].getY(Y_inds_right->first), z0_1 = ys_[X_inds->first].getY(Y_inds_left->second);
			double y0_0 = ys_[X_inds->first].getX(Y_inds_left->first), y1_0 = ys_[X_inds->second].getX(Y_inds_right->first), y0_1 = ys_[X_inds->first].getX(Y_inds_left->second);
			double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
			//first interpolate along X axis, then Y axis
			double y0 = y0_0 + (X-x0)*(y1_0-y0_0)/(x1-x0);
			double y1 = y0_1 + (X-x0)*(y1_0-y0_1)/(x1-x0);
			double z0 = z0_0 + (X-x0)*(z1_0-z0_0)/(x1-x0);
			double z1 = z0_1 + (X-x0)*(z1_0-z0_1)/(x1-x0);
			return z0 + (z1-z0)*(Y-y0)/(y1-y0);
		}
		//4 point linear interpolation
		double z0_0 = ys_[X_inds->first].getY(Y_inds_left->first), z0_1 = ys_[X_inds->first].getY(Y_inds_left->second),
				z1_0 = ys_[X_inds->second].getY(Y_inds_right->first), z1_1 = ys_[X_inds->second].getY(Y_inds_right->second);
		double y0_0 = ys_[X_inds->first].getX(Y_inds_left->first), y0_1 = ys_[X_inds->first].getX(Y_inds_left->second),
				y1_0 = ys_[X_inds->second].getX(Y_inds_right->first), y1_1 = ys_[X_inds->second].getX(Y_inds_right->second);
		double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
		//first interpolate along X axis, then Y axis
		double y0 = y0_0 + (X-x0)*(y1_0-y0_0)/(x1-x0);
		double y1 = y0_1 + (X-x0)*(y1_1-y0_1)/(x1-x0);
		double z0 = z0_0 + (X-x0)*(z1_0-z0_0)/(x1-x0);
		double z1 = z0_1 + (X-x0)*(z1_1-z0_1)/(x1-x0);
		return z0 + (z1-z0)*(Y-y0)/(y1-y0);
	}
	return 0;
}

double FunctionTable::getY(double X, double val) const
{
	boost::optional<std::pair<std::size_t, std::size_t>> X_inds = getX_indices(X);
	if (!X_inds)
		return -1;
	if (X_inds->second >= xs_.size())
		X_inds->second = X_inds->first;
	if (X_inds->first == X_inds->second) {
		boost::optional<std::pair<std::size_t, std::size_t>> Y_inds = ys_[X_inds->first].getY_indices(val);
		if (!Y_inds)
			return -1;
		if (Y_inds->first == Y_inds->second) {
			return ys_[X_inds->first].getX(Y_inds->first);
		} else {
			double z0 = ys_[X_inds->first].getY(Y_inds->first), z1 = ys_[X_inds->first].getY(Y_inds->second);
			double y0 = ys_[X_inds->first].getX(Y_inds->first), y1 = ys_[X_inds->first].getX(Y_inds->second);
			return y0 + (y1 - y0) * (val - z0) / (z1 - z0);
		}
	} else {
		boost::optional<std::pair<std::size_t, std::size_t>> Y_inds_left = ys_[X_inds->first].getY_indices(val);
		boost::optional<std::pair<std::size_t, std::size_t>> Y_inds_right = ys_[X_inds->second].getY_indices(val);
		if (Y_inds_left)
			if (Y_inds_left->second >= ys_[X_inds->first].size())
				Y_inds_left->second = Y_inds_left->first;
		if (Y_inds_right)
			if (Y_inds_right->second >= ys_[X_inds->second].size())
				Y_inds_right->second = Y_inds_right->first;

		if ((!Y_inds_right) && (!Y_inds_left))
			return 0;
		if (!Y_inds_right) {
			if (Y_inds_left->first == Y_inds_left->second) {
				return ys_[X_inds->first].getX(Y_inds_left->first);
			} else {
				double z0 = ys_[X_inds->first].getY(Y_inds_left->first), z1 = ys_[X_inds->first].getY(Y_inds_left->second);
				double y0 = ys_[X_inds->first].getX(Y_inds_left->first), y1 = ys_[X_inds->first].getX(Y_inds_left->second);
				return y0 + (y1 - y0) * (val - z0) / (z1 - z0);
			}
		}
		if (!Y_inds_left) {
			if (Y_inds_right->first == Y_inds_right->second) {
				return ys_[X_inds->first].getX(Y_inds_right->first);
			} else {
				double z0 = ys_[X_inds->first].getY(Y_inds_right->first), z1 = ys_[X_inds->first].getY(Y_inds_right->second);
				double y0 = ys_[X_inds->first].getX(Y_inds_right->first), y1 = ys_[X_inds->first].getX(Y_inds_right->second);
				return y0 + (y1 - y0) * (val - z0) / (z1 - z0);
			}
		}
		if ((Y_inds_left->first == Y_inds_left->second) && (Y_inds_right->first == Y_inds_right->second)) {
			double y0 = ys_[X_inds->first].getX(Y_inds_left->first), y1 = ys_[X_inds->second].getX(Y_inds_right->first);
			double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
			return y0 + (y1 - y0) * (X - x0) / (x1 - x0);
		}
		if (Y_inds_left->first == Y_inds_left->second) {//three point linear interpolation
			double z0_0 = ys_[X_inds->first].getY(Y_inds_left->first), z1_0 = ys_[X_inds->second].getY(Y_inds_right->first), z1_1 = ys_[X_inds->second].getY(Y_inds_right->second);
			double y0_0 = ys_[X_inds->first].getX(Y_inds_left->first), y1_0 = ys_[X_inds->second].getX(Y_inds_right->first), y1_1 = ys_[X_inds->second].getX(Y_inds_right->second);
			double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
			//first interpolate along X axis, then Y axis
			double y0 = y0_0 + (X - x0)*(y1_0 - y0_0) / (x1 - x0);
			double y1 = y0_0 + (X - x0)*(y1_1 - y0_0) / (x1 - x0);
			double z0 = z0_0 + (X - x0)*(z1_0 - z0_0) / (x1 - x0);
			double z1 = z0_0 + (X - x0)*(z1_1 - z0_0) / (x1 - x0);
			return y0 + (y1 - y0)*(val - z0) / (z1 - z0);
		}
		if (Y_inds_right->first == Y_inds_right->second) {//three point linear interpolation
			double z0_0 = ys_[X_inds->first].getY(Y_inds_left->first), z1_0 = ys_[X_inds->second].getY(Y_inds_right->first), z0_1 = ys_[X_inds->first].getY(Y_inds_left->second);
			double y0_0 = ys_[X_inds->first].getX(Y_inds_left->first), y1_0 = ys_[X_inds->second].getX(Y_inds_right->first), y0_1 = ys_[X_inds->first].getX(Y_inds_left->second);
			double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
			//first interpolate along X axis, then Y axis
			double y0 = y0_0 + (X - x0)*(y1_0 - y0_0) / (x1 - x0);
			double y1 = y0_1 + (X - x0)*(y1_0 - y0_1) / (x1 - x0);
			double z0 = z0_0 + (X - x0)*(z1_0 - z0_0) / (x1 - x0);
			double z1 = z0_1 + (X - x0)*(z1_0 - z0_1) / (x1 - x0);
			return y0 + (y1 - y0)*(val - z0) / (z1 - z0);
		}
		//4 point linear interpolation
		double z0_0 = ys_[X_inds->first].getY(Y_inds_left->first), z0_1 = ys_[X_inds->first].getY(Y_inds_left->second),
			z1_0 = ys_[X_inds->second].getY(Y_inds_right->first), z1_1 = ys_[X_inds->second].getY(Y_inds_right->second);
		double y0_0 = ys_[X_inds->first].getX(Y_inds_left->first), y0_1 = ys_[X_inds->first].getX(Y_inds_left->second),
			y1_0 = ys_[X_inds->second].getX(Y_inds_right->first), y1_1 = ys_[X_inds->second].getX(Y_inds_right->second);
		double x0 = xs_[X_inds->first], x1 = xs_[X_inds->second];
		//first interpolate along X axis, then Y axis
		double y0 = y0_0 + (X - x0)*(y1_0 - y0_0) / (x1 - x0);
		double y1 = y0_1 + (X - x0)*(y1_1 - y0_1) / (x1 - x0);
		double z0 = z0_0 + (X - x0)*(z1_0 - z0_0) / (x1 - x0);
		double z1 = z0_1 + (X - x0)*(z1_1 - z0_1) / (x1 - x0);
		return y0 + (y1 - y0)*(val - z0) / (z1 - z0);
	}
	return -1;
}

void FunctionTable::push(double X, const DataVector& vals)
{
  std::size_t X_index = 0;
  if (xs_.empty()) {
    xs_.push_back(X);
    ys_.push_back(vals);
  } else {
    if (X > xs_.back()) {
      xs_.push_back(X);
      ys_.push_back(vals);
      return;
    }
    if (X < xs_.front()) {
      X_index = 0;
      xs_.insert(xs_.begin(), X);
      ys_.insert(ys_.begin(), vals);
      return;
    }
    boost::optional<std::pair<std::size_t, std::size_t>> x_inds = getX_indices(X);
    if (x_inds->first == x_inds->second) {
      X_index = x_inds->first; // Merging present data with new vals
      for (std::size_t i = 0, i_end_ = vals.size(); i!=i_end_; ++i) {
        ys_[X_index].insert(vals.getX(i), vals.getY(i));
      }
    }
    xs_.insert(xs_.begin() + x_inds->second, X);
    ys_.insert(ys_.begin() + x_inds->second, vals);
  }
}

void FunctionTable::push(double X, double Y, double val)
{
	std::size_t X_index = 0;
	if (xs_.empty()) {
		xs_.push_back(X);
		ys_.push_back(DataVector(1, 2));
	} else {
		if (X > xs_.back()) {
			xs_.push_back(X);
			ys_.push_back(DataVector(1, 2));
			X_index = xs_.size()-1;
			goto second;
		}
		if (X < xs_.front()) {
			X_index = 0;
			xs_.insert(xs_.begin(), X);
			ys_.insert(ys_.begin(), DataVector(1, 2));
			goto second;
		}
		boost::optional<std::pair<std::size_t, std::size_t>> x_inds = getX_indices(X);
		if (x_inds->first == x_inds->second) {
			X_index = x_inds->first;
			goto second;
		}
		xs_.insert(xs_.begin() + x_inds->second, X);
		ys_.insert(ys_.begin() + x_inds->second, DataVector(1, 2));
		X_index = x_inds->second;
	}
second:
	ys_[X_index].insert(Y, val);//not a push_back, insert/replaces value preserving order
}

void FunctionTable::read(std::ifstream& str)
{
	clear();
	double val, val2;
	size_t size, size_E;
	std::size_t counter = 0, counter_E;
	bool valid = true;
	while (!str.eof()) {
		str.read((char*)&size,sizeof(std::size_t));
		ys_.resize(size, DataVector(1, 2));
		xs_.resize(size);
		counter = 0;
		while (!str.eof() && counter != size) {
			str.read((char*)&val, sizeof(double));
			xs_[counter] = val;
			if (str.eof()) {
				valid = false;
				break;
			}
			str.read((char*)&size_E, sizeof(std::size_t));
			ys_[counter].resize(size_E);
			counter_E = 0;
			while (!str.eof() && counter_E != size_E) {
				str.read((char*)&val,sizeof(double));
				if (str.eof()) {
					valid = false;
					break;
				}
				str.read((char*)&val2,sizeof(double));
				ys_[counter][counter_E] = std::pair<double, double>(val, val2);
				++counter_E;
			}
			++counter;
		}
		break;
	}
	if ((counter_E!=size_E)||(counter!=size)||!valid) {
		std::cout<<"FunctionTable: Error while reading file"<<std::endl;
		clear();
	}
}

void FunctionTable::write(std::string fname) const
{
	std::ofstream str;
	str.open(fname, std::ios_base::trunc | std::ios_base::binary);
	write(str);
	str.close();
}

void FunctionTable::write (std::ofstream& str) const
{
	std::size_t real_size= xs_.size();
	str.write((char*)&real_size, sizeof(std::size_t));
	for (std::size_t Ey_ind = 0, Ey_ind_end_= real_size; Ey_ind != Ey_ind_end_; ++Ey_ind) {
		str.write((char*)&xs_[Ey_ind], sizeof(double));
		std::size_t size = ys_[Ey_ind].size();
		str.write((char*)&size, sizeof(std::size_t));
		for (std::size_t E = 0, E_end_= size; E != E_end_; ++E) {
			std::pair<double, double> yz = ys_[Ey_ind].getXY(E);
			str.write((char*)&yz.first, sizeof(double));
			str.write((char*)&yz.second, sizeof(double));
		}
	}
}
