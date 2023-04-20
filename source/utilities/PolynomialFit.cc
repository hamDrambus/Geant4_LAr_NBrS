#include <utilities/PolynomialFit.hh>

namespace uBLAS = boost::numeric::ublas;

PolynomialFit::PolynomialFit(std::size_t order)
{
	setOrder(order);
}

PolynomialFit::~PolynomialFit()
{}

std::vector<double> PolynomialFit::operator ()(const std::vector<std::pair<double, double>> &vals_in, boost::optional<double> &in_x0) const {
	return (*this)(vals_in, 0, vals_in.size(), in_x0);
}

std::vector<double> PolynomialFit::operator ()(const std::vector<std::pair<double, double>> &vals,
	int offset, int N_points, boost::optional<double> &in_x0) const//only for a part of a vector
{
	std::vector<double> out;
	if ((vals.size()-offset) < N_points) {
		std::cout<<"PolynomialFit::operator(): Error: N points is out of range:"<<std::endl;
		std::cout<<"\tx.size()="<< vals.size()<<" offset="<<offset<<" N_points="<<N_points<<std::endl;
		return out;
	}
	if (offset < 0) {
		std::cout<<"PolynomialFit::operator(): Error: offset is out of range:"<<std::endl;
		std::cout<<"\tx.size()="<< vals.size()<<" offset="<<offset<<" N_points="<<N_points<<std::endl;
		return out;
	}
	if (N_points < (_order + 1)) {
		std::cout<<"PolynomialFit::operator(): Error: no enough N points for fit:"<<std::endl;
		std::cout<<"\torder="<<_order<<" N_points="<<N_points<<std::endl;
		return out;
	}
	in_x0 = (in_x0 ? *in_x0 : vals[offset].first); //It is bad to set x0 to some fixed value (e.g. 0) because
	//interpolating too far from it will result in unstable results due to limited precision.
	//Ideally x0 should be set to the point at which we interpolate the data.
	if (1 == _order) {
		out.resize(2);
		out[1] = (vals[offset + 1].second - vals[offset].second) / (vals[offset + 1].first - vals[offset].first);
		out[0] = vals[offset].second + (*in_x0 - vals[offset].first)*out[1];
		//^value at in_x0 point
	} else {
		uBLAS::matrix<double> mat(N_points, _order + 1);
		for (int col = 0, col_end_ = mat.size2(); col < col_end_; ++col)
			for (int row = 0, row_end_ = mat.size1(); row < row_end_; ++row)
				mat(row, col) = pow(vals[offset + row].first - *in_x0, col);
		uBLAS::vector<double> Y(N_points);
		for (int row = 0, row_end_ = Y.size(); row < row_end_; ++row)
			Y[row] = vals[offset + row].second;
		//Solve the equation mat^T*mat*X = mat^T*Y for X via LU decomposition (mat is generally not diagonal)
		Y = uBLAS::prod(uBLAS::trans(mat), Y);
		mat = uBLAS::prod(uBLAS::trans(mat), mat);
		int res = uBLAS::lu_factorize(mat);
		if (res != 0)
			return out;
		uBLAS::inplace_solve(mat, Y, uBLAS::unit_lower_tag());
		uBLAS::inplace_solve(mat, Y, uBLAS::upper_tag());
		out.resize(Y.size());
		std::copy(Y.begin(), Y.end(), out.begin());
	}
	if (out.size() != (_order+1)) {
		out.resize(0);
		return out;
	}
	return out;
}

//=========================================================

DataVector::DataVector(std::size_t fit_order, std::size_t N_used_) :
		fitter(fit_order), use_left(false), use_right(false)
{
	setNused(N_used_);
}

DataVector::DataVector(std::vector<double> &xx, std::vector<double> &yy, std::size_t fit_order, std::size_t N_used_) : DataVector(fit_order, N_used_)
{
	initialize(xx, yy, fit_order, N_used_);
}

DataVector::DataVector(std::string fname) : DataVector(1, 2)
{
	std::ifstream str;
	str.open(fname);
	if (!str.is_open() || !str.good()) {
		std::cout << "DataVector::DataVector:: Warning: Could not open file \"" << fname<<"\"" << std::endl;
	}
	read(str);
	if (!isValid()) {
		std::cout << "DataVector::DataVector:: Warning: file \"" << fname << "\" contains invalid data" << std::endl;
	}
}

DataVector::~DataVector() {}

void DataVector::initialize(std::vector<double> &xx, std::vector<double> &yy, std::size_t fit_order, std::size_t N_used_)
{
	setNused(N_used_);
	fitter.setOrder(fit_order);
	if (xx.size()!=yy.size()) {
		std::cout<<"DataVector::initialize(): Error: x and y data size mismatch!"<<std::endl;
		return;
	}
	std::size_t i_end_ = xx.size();
	xys.resize(i_end_);
	for (std::size_t i = 0; i != i_end_; ++i)
		xys[i] = std::pair<double, double>(xx[i], yy[i]);
	std::sort(xys.begin(), xys.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b)->bool {
		return a.first < b.first;
	});
}

void DataVector::insert(double x, double y) //does not disrupt order
{
	std::size_t sz = xys.size();
	if (0 == sz) {
		xys.push_back(std::pair<double, double>(x, y));
		return;
	}
	if (x > xys.back().first) {
		xys.push_back(std::pair<double, double>(x, y));
		return;
	}
	if (x < xys.front().first) {
		xys.insert(xys.begin(), std::pair<double, double>(x, y));
		return;
	}
	boost::optional<std::pair<std::size_t, std::size_t>> inds = getX_indices(x);
	if (inds->first == inds->second) //do not insert points with equal x, replace only
		xys[inds->first].second = y;
	else
		xys.insert(xys.begin() + inds->second, std::pair<double, double>(x, y));
}

void DataVector::push_back (double x, double y)//faster version not checking that the new array is ordered.
{
	xys.push_back(std::pair<double, double> (x, y));
}

void DataVector::push_back (std::pair<double, double> xy)
{
	xys.push_back(xy);
}

double DataVector::integrate(const IntegrationRange &range)
{
  double out = 0.0;
  if (range.NumOfIndices()<=0)
    return out;
  double X_prev = range.Value(0);
  double Y_prev = 0;
  for (long int i = 0, i_end_ = range.NumOfIndices(); i != i_end_; ++i) {
    double X = range.Value(i);
    double Y = operator ()(X);
    if ((Y == DBL_MAX || Y == -DBL_MAX)&&(Y_prev == DBL_MAX || Y_prev == -DBL_MAX)) {
      out = DBL_MAX;
      return out;
    }
    if (Y == DBL_MAX || Y == -DBL_MAX) {
      out += Y_prev*(X - X_prev);
      X_prev = X;
      Y_prev = Y;
      continue;
    }
    if (Y_prev == DBL_MAX || Y_prev == -DBL_MAX) {
      out += Y*(X - X_prev);
      X_prev = X;
      Y_prev = Y;
      continue;
    }
    out += 0.5*(Y + Y_prev)*(X - X_prev);
    X_prev = X;
    Y_prev = Y;
  }
  return out;
}

void DataVector::integrate(void)
{
  if (xys.empty())
    return;
  double X_prev = xys.front().first;
  double Y_prev = 0;
  for (long int i = 0, i_end_ = xys.size(); i != i_end_; ++i) {
    double X = xys[i].first;
    double Y = xys[i].second;
    if (i != 0)
      xys[i].second = xys[i-1].second + 0.5*(Y + Y_prev)*(X - X_prev);//Integral
    else
      xys[i].second = 0;
    X_prev = X;
    Y_prev = Y;
  }
}

std::vector<double> DataVector::get_Xs(double scaling) const
{
	std::vector<double> out(xys.size());
	for (std::size_t i = 0, i_end_ = xys.size(); i!=i_end_; ++i) {
		out[i] = xys[i].first * scaling;
	}
	return out;
}

std::vector<double> DataVector::get_Ys(double scaling) const
{
	std::vector<double> out(xys.size());
	for (std::size_t i = 0, i_end_ = xys.size(); i!=i_end_; ++i) {
		out[i] = xys[i].second * scaling;
	}
	return out;
}

void DataVector::read(std::ifstream& str)
{
	clear();
	std::string line, word;
	int line_n = 0;
	while (!str.eof() && str.is_open()) {
		std::getline(str, line);
		++line_n;
		bool invalid_header_is_data = true;
		if (1 == line_n) {
			//parse "//Order	N_used	use_left use_right is_set_left is_set_right left_value right_value"
			int word_n = 0;
			try {
				double dval;
				std::size_t ival;
				bool is_set_right;
				bool is_set_left;
				if (line.size() < 2)
					throw std::runtime_error("Header line has wrong format (too small)");
				if ((line[0] != '/') || (line[1] != '/'))
					throw std::runtime_error("Header line has wrong format (does not start with \"//\")");
				invalid_header_is_data = false;
				line.erase(line.begin(), line.begin() + 2);
				word = strtoken(line, "\t \r");
				++word_n;
				ival = boost::lexical_cast<std::size_t>(word);
				setOrder(ival);

				word = strtoken(line, "\t \r");
				++word_n;
				ival = boost::lexical_cast<std::size_t>(word);
				setNused(ival);

				word = strtoken(line, "\t \r");
				++word_n;
				ival = std::stoi(word);
				use_leftmost(ival);

				word = strtoken(line, "\t \r");
				++word_n;
				ival = std::stoi(word);
				use_rightmost(ival);

				word = strtoken(line, "\t \r");
				++word_n;
				ival = std::stoi(word);
				is_set_left = ival;

				word = strtoken(line, "\t \r");
				++word_n;
				ival = std::stoi(word);
				is_set_right = ival;

				word = strtoken(line, "\t \r");
				++word_n;
				dval = boost::lexical_cast<double>(word);
				if (is_set_left)
					left_value = dval;
				else
					left_value = boost::none;
				word = strtoken(line, "\t \r");
				++word_n;
				dval = boost::lexical_cast<double>(word);
				if (is_set_right)
					right_value = dval;
				else
					right_value = boost::none;
				continue;
			} catch (boost::bad_lexical_cast &e) {
				use_rightmost(false);
				use_leftmost(false);
				setOrder(1);
				setNused(2);
				set_out_value(0);
				if (!invalid_header_is_data) {
          std::cout << "DataVector::read: Error on line " << line_n << ". Can't convert word #"<<word_n<<" \""<<word<<"\" to numerical value" << std::endl;
          std::cout << e.what() << std::endl;
          std::cout << "DataVector::read: bad header, using default"<<std::endl;
					continue;
				}
			} catch (std::exception &e) {
				use_rightmost(false);
				use_leftmost(false);
				setOrder(1);
				setNused(2);
				set_out_value(0);
				if (!invalid_header_is_data) {
          std::cout << "DataVector::read: Unforeseen exception on line " << line_n << " word #"<<word_n<<":" << std::endl;
          std::cout << e.what() << std::endl;
          std::cout << "DataVector::read: bad header, using default"<<std::endl;
					continue;
				}
			}
		}
		if (line.size() >= 2) //Ignore simple c style comment
			if ((line[0] == '/') && (line[1] == '/'))
				continue;
		try {
			word = strtoken(line, "\t \r");
			double x = boost::lexical_cast<double>(word);
			word = strtoken(line, "\t \r");
			double val = boost::lexical_cast<double>(word);
			insert(x, val);
		} catch (boost::bad_lexical_cast &e) {
			continue;
		} catch (std::exception &e) {
			std::cerr << "DataVector::read: Unforeseen exception on line " << line_n << std::endl;
			std::cerr << e.what() << std::endl;
			return;
		}
	}
}

void DataVector::write(std::string fname, std::string comment) const
{
	std::ofstream str;
	str.open(fname, std::ios_base::trunc);
	write(str, comment);
	str.close();
}

void DataVector::write(std::ofstream& str, std::string comment) const
{
	//"//Order	N_used	use_left use_right is_set_left is_set_right left_value right_value"
	str << "//" << boost::lexical_cast<std::string>(getOrder())
		<< "\t" << boost::lexical_cast<std::string>(N_used)
		<< "\t" << (use_left ? 1 : 0) << "\t" << (use_right ? 1 : 0)
		<< "\t" << (left_value ? 1 : 0) << "\t" << (right_value ? 1 : 0)
		<< "\t" << boost::lexical_cast<std::string>(left_value ? *left_value : 0)
		<< "\t" << boost::lexical_cast<std::string>(right_value ? *right_value : 0 )<< std::endl;
	if (!comment.empty())
		str << "//" << comment << std::endl;
	for (std::size_t i = 0, i_end_ = xys.size(); i != i_end_; ++i) {
		str << boost::lexical_cast<std::string>(xys[i].first) << "\t"
			<< boost::lexical_cast<std::string>(xys[i].second) << std::endl;
	}
}

double DataVector::operator()(double X_point, boost::optional<double> x0) const
{
	if (!isValid())
		return DBL_MAX;
	if (X_point < xys.front().first) {
		if (use_left)
			return xys.front().second; //TODO: maybe add scaling
		if (left_value)
			return *left_value;
	}
	if (X_point > xys.back().first) {
		if (use_right)
			return xys.back().second;
		if (right_value)
			return *right_value;
	}
	boost::optional<std::pair<std::size_t, std::size_t>> indices = getX_indices(X_point);
	if (boost::none == indices)
		return DBL_MAX;
	//expand indices to [n_min, n_max] are used, not [n_min, n_max). N_used == n_max - n_min + 1 >= order + 1
	std::size_t span = (N_used - 1) / 2;  //asymmetrical interpolation range in the case of odd order.
	if (indices->first < span) { //first is too low
		indices = std::pair<std::size_t, std::size_t>(0, N_used - 1);
	} else {
		indices->first = indices->first - span;
		indices->second = indices->first + (N_used - 1);
		std::size_t sz = xys.size();
		if (indices->second >= sz) {
			indices = std::pair<std::size_t, std::size_t>(sz - N_used, sz - 1);
		}
	}
	std::vector<double> coefs = fitter(xys, indices->first, indices->second - indices->first + 1, x0); //i_max-i_min+1==N_used
	if (0 != coefs.size())
		return polynomial_value(X_point, *x0, coefs);
	return DBL_MAX;
}

boost::optional<std::pair<std::size_t, std::size_t>> DataVector::getX_indices(double x) const
{
	auto picker = [] (const std::vector<std::pair<double, double>> & xys, std::size_t index) -> double {
		return xys[index].first;
	};
	return find_bounds(xys, x, picker);
}
//I chose to copy getX_indices code here instead of using parameter or lambda value picker function in the fear that it will reduce the performance. I did not test that it would.
boost::optional<std::pair<std::size_t, std::size_t>> DataVector::getY_indices(double y) const
{
	auto picker = [] (const std::vector<std::pair<double, double>> & xys, std::size_t index) -> double {
		return xys[index].second;
	};
	return find_bounds(xys, y, picker);
}

double DataVector::polynomial_value(double x, double x0, const std::vector<double>& coefs) const
{
	std::size_t order = coefs.size();
	double out_ = 0;
	for (std::size_t o_O = 0; o_O < order; ++o_O)
		out_ += coefs[o_O]*pow(x-x0, o_O);
	return out_;
}

//=========================================================

PDF_routine::PDF_routine(): cdf_ready(false)
{}

PDF_routine::PDF_routine(std::vector<double> &pdf_xx, std::vector<double> &pdf_yy) : cdf_ready(false)
{
	if (pdf_xx.size() != pdf_yy.size()) {
		std::cout << "PDF_routine::PDF_routine(): Error: x and y data size mismatch!" << std::endl;
		return;
	}
	std::size_t i_end_ = pdf_xx.size();
	vals.resize(i_end_);
	pdf_data temp;
	for (std::size_t i = 0; i != i_end_; ++i) {
		temp.x = pdf_xx[i];
		temp.pdf = pdf_yy[i];
		temp.cdf = 0;
		vals[i] = temp;
	}
	std::sort(vals.begin(), vals.end(), [](const pdf_data &a, const pdf_data &b)->bool {
		return a.x < b.x;
	});
	pdf_to_cdf();
}

void PDF_routine::normalize_cdf(void)
{
	std::size_t i_end_ = vals.size();
	if (0 == i_end_)
		return;
	if (1 == i_end_) {
		vals[0].cdf = 1;
		return;
	}
	double denom = 1.0 / vals[i_end_ - 1].cdf;
	for (std::size_t i = 0; i != i_end_; ++i)
		vals[i].cdf *= denom;//multiplication is cheaper than division
}

void PDF_routine::pdf_to_cdf(void)
{
	double X_prev = vals.front().x;
	double Y_prev = 0;
	for (long int i = 0, i_end_ = vals.size(); i != i_end_; ++i) {
		double X = vals[i].x;
		double Y = vals[i].pdf;
		if (i != 0)
			vals[i].cdf = vals[i - 1].cdf + 0.5*(Y + Y_prev)*(X - X_prev);//Integral
		else
			vals[i].cdf = 0;
		X_prev = X;
		Y_prev = Y;
	}
	normalize_cdf();
	cdf_ready = true;
}

PDF_routine::~PDF_routine()
{}

void PDF_routine::insert(double x, double y) //Preserves sorting, recalculates pdf, quite expensive when insering not to the end
{
	pdf_data temp;
	for (auto i = vals.begin(), i_end_ = vals.end(); i != i_end_; ++i) {
		if (x == i->x) {
			i->pdf = y; //do not insert points with equal x, replace only
			pdf_to_cdf();
			return;
		}
		if (i->x > x) {
			temp.x = x;
			temp.pdf = y;
			temp.cdf = 0;
			vals.insert(i, temp);
			pdf_to_cdf();
			return;
		}
	}
	temp.x = x;
	temp.pdf = y;
	std::size_t sz = vals.size();
	if (0 == sz) {
		temp.cdf = 0;
		vals.insert(vals.end(), temp);
	} else {
		--sz;
		temp.cdf = vals[sz].cdf + 0.5*(y + vals[sz].pdf)*(x - vals[sz].x); //Integral
		vals.insert(vals.end(), temp);
	}
	normalize_cdf();
}

void PDF_routine::push_back(double x, double y, bool recalculate_cdf) //Faster version, not checking whether sorting is preserved.
{
	pdf_data temp;
	temp.x = x;
	temp.pdf = y;
	std::size_t sz = vals.size();
	if (sz == 0) {
		temp.cdf = 1;
		vals.insert(vals.end(), temp);
	} else {
		if (vals[sz - 1].x == x) { //check only equality with the last value
			--sz;
			vals[sz].pdf = y;
			if (recalculate_cdf) {
				if (sz > 0) {
					vals[sz].cdf = vals[sz - 1].cdf + 0.5*(vals[sz - 1].pdf + vals[sz].pdf)*(vals[sz].x - vals[sz - 1].x); //Integral
					normalize_cdf();
					return;
				} else {
					normalize_cdf();
					return;
				}
			}
			cdf_ready = false;
			return;
		}
		if (recalculate_cdf) {
			--sz;
			temp.cdf = vals[sz].cdf + 0.5*(y + vals[sz].pdf)*(x - vals[sz].x); //Integral
			vals.insert(vals.end(), temp);
			normalize_cdf();
			cdf_ready = true;
		} else {
			cdf_ready = false;
			temp.cdf = 0;
			vals.insert(vals.end(), temp);
		}
	}
}

void PDF_routine::write(std::string fname, std::string comment) const
{
	std::ofstream str;
	str.open(fname, std::ios_base::trunc);
	write(str, comment);
	str.close();
}

void PDF_routine::write(std::ofstream& str, std::string comment) const
{
	if (!comment.empty())
		str << "//" << comment << std::endl;
	for (std::size_t i = 0, i_end_ = vals.size(); i != i_end_; ++i) {
		str << boost::lexical_cast<std::string>(vals[i].x) << "\t"
			<< boost::lexical_cast<std::string>(vals[i].pdf) << std::endl;
	}
}

bool PDF_routine::read(std::string fname)
{
	std::ifstream str;
	str.open(fname);
	if (!str.is_open() || !str.good()) {
		std::cout << "PDF_routine::read:: Warning: Could not open file \"" << fname<<"\"" << std::endl;
		return false;
	}
	read(str);
	if (!isValid()) {
		std::cout << "PDF_routine::read:: Warning: file \"" << fname << "\" contains invalid data" << std::endl;
		return false;
	}
	return true;
}

void PDF_routine::read(std::ifstream& str)
{
	clear();
	std::vector<double> ixx, iyy;//it's better to push whole read array instead of element-wise filling
	std::string line, word;
	int line_n = 0;
	while (!str.eof() && str.is_open()) {
		std::getline(str, line);
		++line_n;
		if (line.size() >= 2) //Ignore simple c style comment
			if ((line[0] == '/') && (line[1] == '/'))
				continue;
		try {
			double val, x;
			word = strtoken(line, "\t \r");
			x = boost::lexical_cast<double>(word);
			word = strtoken(line, "\t \r");
			val = boost::lexical_cast<double>(word);
			insert(x, val);
		}  catch (boost::bad_lexical_cast &e) {
			continue;
		} catch (std::exception &e) {
			std::cerr << "PDF_routine::read: Unforeseen exception on line " << line_n << std::endl;
			std::cerr << e.what() << std::endl;
			pdf_to_cdf();
			return;
		}
	}
	pdf_to_cdf();
}

//has no internal random engine
double PDF_routine::operator()(double Rand) const
{
	for (std::size_t i = 0, i_end_ = vals.size(); i != i_end_; ++i) {
		if (Rand <= vals[i].cdf) {
			if (0 == i)
				return vals[i].x;
			double x0 = vals[i - 1].x, x1 = vals[i].x;
			double y0 = vals[i - 1].cdf, y1 = vals[i].cdf;
			return x0 + (x1 - x0)*(Rand - y0) / (y1 - y0);
		}
	}
	return (vals.size() ? vals.back().x : 0);
}
