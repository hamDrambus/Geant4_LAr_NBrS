#ifndef FIND_BOUNDS_CONTAINER_HPP_
#define FIND_BOUNDS_CONTAINER_HPP_

#include <boost/optional.hpp>

//Warning! These functions have defined behavior only when container values are sorted in the ascending order.
template <class Value, class Container, class Picker>
boost::optional<std::pair<std::size_t, std::size_t>> find_bounds(const Container & cont, Value x, Picker picker) {
	boost::optional<std::pair<std::size_t, std::size_t>> out;
	std::size_t sz = cont.size();
	if (0 == sz)
		return out;
	if (x <= picker(cont, 0)) {
		out = std::pair<std::size_t, std::size_t>(0, 0);
		return out;
	}
	if (x >= picker(cont, sz - 1)) {
		out = std::pair<std::size_t, std::size_t>(sz - 1, sz - 1);
		return out;
	}
	//find first x which is not less that X_point. That is index bounding X_point: xs[first] <= X_point < xs[first + 1]
	//See std::lower_bound and std::upper_bound
	std::size_t count = sz;
	std::size_t first = 0;
	//std::lower_bound(xys.begin(), xys.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b)->bool{
	//	return a.first<b.first;
	//});
	while (count > 0) {
		std::size_t step = count / 2;
		std::size_t ind = first + step;
		if (!(x < picker(cont, ind))) {
			first = ++ind;
			count -= step + 1;
		} else
			count = step;
	}
	//first is such, that x>=xs[first-1] and x<xs[first]
	//first always != 0 here
	--first;
	if (x == picker(cont, first)) {
		out = std::pair<std::size_t, std::size_t>(first, first);
		return out;
	}
	out = std::pair<std::size_t, std::size_t>(first, first + 1);
	return out;
}

template <class Value, class Container>
boost::optional<std::pair<std::size_t, std::size_t>> find_bounds(const Container & cont, Value x)
{
	auto value_picker = [](const Container & cont, std::size_t index) -> Value {
		return cont[index];
	};
	return find_bounds(cont, x, value_picker);
}

#endif // FIND_BOUNDS_CONTAINER_HPP_
