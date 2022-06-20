#include <utilities/IntegrationInterval.hh>

IntegrationInterval::IntegrationInterval(long double left, long double right, long double step):
		left_(std::min(left, right)), right_(std::max(left, right)), step_(std::fabs(step))
{}

IntegrationRange IntegrationInterval::operator += (const IntegrationRange &r)
{
  IntegrationRange out(*this);
	out = out + r;
	return out;
}

IntegrationRange IntegrationInterval::operator += (const IntegrationInterval &r)
{
  IntegrationRange out(*this);
	out = out + r;
	return out;
}

IntegrationInterval& IntegrationInterval::operator*= (const long double& r)
{
  left_ *= r;
  right_ *= r;
  step_ *= r;
  return *this;
}

void IntegrationInterval::Rescale(const long double&& new_left, const long double&& new_right)
{
  long double l = std::min(new_left, new_right);
  long double r = std::max(new_left, new_right);
  long double factor = (r - l)/(right_ - left_);
  left_ = l;
  right_ = r;
  step_ *= factor;
}


IntegrationRange::IntegrationRange(void) {}
IntegrationRange::IntegrationRange(const IntegrationInterval &inter)
{
  if (inter.left_!=inter.right_)
    arr_.push_back(inter);
}

long int IntegrationRange::NumOfIndices(void) const
{
	long int N = 0;
	for (int j=0, end_ = arr_.size(); j!=end_; ++j) {
		N+=std::max(std::lround((arr_[j].right_-arr_[j].left_)/arr_[j].step_), (long int)1);
		if (j==(end_-1)) {
			N+=1;
		} else {
			if (arr_[j+1].left_!=arr_[j].right_)
				N+=1; //intervals are not adjacent
		}
	}
	return N;
}

long double IntegrationRange::Value (long int index) const
{
	long int Nsum = 0, Noff = 0, N = 0;
	int j=0, end_ = arr_.size();
	bool closed = false;
	for (; j!=end_; ++j) {
		N = std::max(std::lround((arr_[j].right_-arr_[j].left_)/arr_[j].step_), (long int)1);
		if (j==(end_-1)) {
			N+=1;
			closed= true;
		} else {
			closed = false;
			if (arr_[j+1].left_!=arr_[j].right_) {
				N+=1; //intervals are not adjacent
				closed = true;
			}
		}
		Nsum+=N;
		if (index<Nsum) {
			Noff = Nsum - N;
			break;
		}
	}
	if (j==end_)
		return DBL_MAX;
	return arr_[j].left_+(index-Noff)*(arr_[j].right_-arr_[j].left_)/(N- (closed ?  1 : 0));
}

void IntegrationRange::Trim (long double left, long double right)
{
	if (right<left) {
		long double temp = left;
		left = right;
		right = temp;
	}
	std::vector<IntegrationInterval> new_arr;
	new_arr.reserve(arr_.size());
  for (int j=0, end_ = arr_.size(); j!=end_; ++j) {
    if ((left<=arr_[j].left_)&&(right>=arr_[j].right_)) {
      new_arr.push_back(arr_[j]);
      continue;
    }
    if (arr_[j].right_<left || arr_[j].left_>right)
      continue;
    if (arr_[j].right_ == left) {
      if (j != (end_ - 1))
        if (arr_[j+1].left_ <= left)
          continue; // Edge point is present in the next interval, not duplicating
      //new_arr.push_back(IntegrationInterval(left, left, 1.0)); //leave one point intersection
      continue;
    }
    if (arr_[j].left_ == right) {
      if (j != (0))
        if (arr_[j-1].right_ >= right)
          continue; // Edge point is present in the next interval, not duplicating
      //new_arr.push_back(IntegrationInterval(right, right, 1.0)); //leave one point intersection
      continue;
    }
    new_arr.push_back(IntegrationInterval(std::max(arr_[j].left_, left), std::min(arr_[j].right_, right), arr_[j].step_));
  }
  new_arr.erase(std::remove_if(new_arr.begin(), new_arr.end(),
      [](IntegrationInterval e){return e.left_ == e.right_;}),
      new_arr.end());
  arr_ = new_arr;
}

long double IntegrationRange::max(void) const
{
  long double out = -DBL_MAX;
  for (int j=0, end_ = arr_.size(); j!=end_; ++j)
    out = std::max(out, arr_[j].right_);
  return out;
}

long double IntegrationRange::min(void) const
{
  long double out = DBL_MAX;
  for (int j=0, end_ = arr_.size(); j!=end_; ++j)
    out = std::min(out, arr_[j].left_);
  return out;
}

void IntegrationRange::Print(std::ostream & str)
{
	for (int j=0, end_ = arr_.size(); j!=end_; ++j) {
		str << "["<<boost::lexical_cast<std::string>(arr_[j].left_) << "; "
			<< boost::lexical_cast<std::string>(arr_[j].step_) << "; "
			<< boost::lexical_cast<std::string>(arr_[j].right_) << "]" << std::endl;
	}
}

IntegrationRange& IntegrationRange::operator += (const IntegrationRange &r)
{
	*this = *this + r;
	return *this;
}

IntegrationRange& IntegrationRange::operator += (const IntegrationInterval &r)
{
	*this = *this + r;
	return *this;
}

IntegrationRange operator+ (const IntegrationRange &l, const IntegrationRange& r)
{
  IntegrationRange out = l;
	for (int j=0, end_=r.arr_.size(); j!=end_; ++j) {
		out += r.arr_[j];
	}
	return out;
}

IntegrationRange operator+ (const IntegrationRange &l, const IntegrationInterval& r)
{
	std::vector<IntegrationInterval> new_l;
	bool is_inside_r = false;
	if (l.arr_.empty()) {
	  IntegrationRange out(r);
		return out;
	}
	for (int j=0, end_=l.arr_.size(); j!=end_; ++j) {
		if (is_inside_r) {
			if (r.right_<l.arr_[j].left_) {
				new_l.push_back(IntegrationInterval(l.arr_[j-1].right_, r.right_, r.step_));
				new_l.insert(new_l.end(),l.arr_.begin()+j,l.arr_.end());
				break;
			}
			if (r.right_<l.arr_[j].right_) {
				new_l.push_back(IntegrationInterval(l.arr_[j-1].right_, l.arr_[j].left_, r.step_));
				new_l.push_back(IntegrationInterval(l.arr_[j].left_, r.right_, std::min(r.step_, l.arr_[j].step_)));
				new_l.push_back(IntegrationInterval(r.right_, l.arr_[j].right_, l.arr_[j].step_));
				new_l.insert(new_l.end(),l.arr_.begin()+j+1,l.arr_.end());
				break;
			}
			//!(r.right_<l.arr_[j].right_)
			new_l.push_back(IntegrationInterval(l.arr_[j-1].right_, l.arr_[j].left_, r.step_));
			new_l.push_back(IntegrationInterval(l.arr_[j].left_, l.arr_[j].right_, std::min(r.step_, l.arr_[j].step_)));
		} else {
			if (r.right_<l.arr_[j].left_) {
				new_l.push_back(IntegrationInterval(r.left_, r.right_, r.step_));
				new_l.insert(new_l.end(),l.arr_.begin()+j,l.arr_.end());
				break;
			}
			if ((r.right_<l.arr_[j].right_)&&(r.left_<l.arr_[j].left_)) {
				new_l.push_back(IntegrationInterval(r.left_, l.arr_[j].left_, r.step_));
				new_l.push_back(IntegrationInterval(l.arr_[j].left_, r.right_, std::min(r.step_, l.arr_[j].step_)));
				new_l.push_back(IntegrationInterval(r.right_, l.arr_[j].right_, l.arr_[j].step_));
				new_l.insert(new_l.end(),l.arr_.begin()+j+1,l.arr_.end());
				break;
			}
			if ((r.right_<l.arr_[j].right_)&&!(r.left_<l.arr_[j].left_)) {
				new_l.push_back(IntegrationInterval(l.arr_[j].left_, r.left_, l.arr_[j].step_));
				new_l.push_back(IntegrationInterval(r.left_, r.right_, std::min(r.step_, l.arr_[j].step_)));
				new_l.push_back(IntegrationInterval(r.right_, l.arr_[j].right_, l.arr_[j].step_));
				new_l.insert(new_l.end(),l.arr_.begin()+j+1,l.arr_.end());
				break;
			}
			if (!(r.left_<l.arr_[j].right_)) {
				new_l.push_back(IntegrationInterval(l.arr_[j].left_, l.arr_[j].right_, l.arr_[j].step_));
			} else { //!(r.right_<l.arr_[j].right_) && (r.left_<l.arr_[j].right_)
				is_inside_r = true;
				if (r.left_<l.arr_[j].left_) {
					new_l.push_back(IntegrationInterval(r.left_, l.arr_[j].left_, r.step_));
					new_l.push_back(IntegrationInterval(l.arr_[j].left_, l.arr_[j].right_, std::min(r.step_, l.arr_[j].step_)));
				} else {
					new_l.push_back(IntegrationInterval(l.arr_[j].left_, r.left_, l.arr_[j].step_));
					new_l.push_back(IntegrationInterval(r.left_, l.arr_[j].right_, std::min(r.step_, l.arr_[j].step_)));
				}
			}
		}
		if ((end_-1) == j) { //special conditions for edges
			if (is_inside_r) {
				new_l.push_back(IntegrationInterval(l.arr_[j].right_, r.right_, r.step_));
			} else { //(r.left_>=l.arr_[j].right_)
				new_l.push_back(IntegrationInterval(r.left_, r.right_, r.step_));
			}
		}
	}
	IntegrationRange out;
	for (int j =0, end_ = new_l.size(); j!=end_; ++j) {
		if (out.arr_.empty()) {
			if (new_l[j].right_>new_l[j].left_)
				out.arr_.push_back(new_l[j]);
			continue;
		}
		if ((out.arr_.back().step_==new_l[j].step_)&&(out.arr_.back().right_==new_l[j].left_))
			out.arr_.back().right_ = new_l[j].right_; //merging
		else {
			if (new_l[j].right_>new_l[j].left_)
				out.arr_.push_back(new_l[j]);
		}
	}
	out.arr_.erase(std::remove_if(out.arr_.begin(), out.arr_.end(),
      [](IntegrationInterval e){return e.left_ == e.right_;}),
	    out.arr_.end());
	return out;
}

IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationRange& r)
{
	return r + l;
}

IntegrationRange operator+ (const IntegrationInterval &l, const IntegrationInterval& r)
{
  IntegrationRange left (l);
	return left + r;
}

IntegrationRange operator* (const IntegrationRange &l, const double & r)
{
  IntegrationRange out = l;
  for (int j=0, end_=out.arr_.size(); j!=end_; ++j) {
    out.arr_[j] *= r;
  }
  return out;
}

IntegrationRange operator/ (const IntegrationRange &l, const double & r)
{
  return l * (1.0 / r);
}
