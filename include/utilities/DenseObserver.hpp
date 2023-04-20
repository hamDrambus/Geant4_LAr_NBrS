#ifndef DENSE_OBSERVER_HPP_
#define DENSE_OBSERVER_HPP_

#include <vector>
#include <boost/type_traits/is_same.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>

#include "FindBoundsContainer.hpp"

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

// Special observer for dense steppers in boost::odeint.
// Its purpose is to call usual observer function not only at the step ends, but
// also at intermediate points so that solution at the step can be interpolated linearly
// This class is created in custom integrate functions (because errors are private members of steppers)
// TODO: conceptually this is not a very good solution, it is better to make special stepper based on dense_stepper instead.
// But it that case, it should be implemented as in the boost:odeint: with tag dispatching and factories (generation functions).
template <class Stepper, class Observer>
class dense_observer {
public:
	typedef typename Stepper::value_type value_type;
	typedef typename Stepper::state_type state_type;
	typedef typename Stepper::time_type time_type;
	typedef typename Stepper::algebra_type algebra_type;
	typedef typename Stepper::operations_type operations_type;
	typedef typename odeint::unwrap_reference< Stepper >::type& stepper_type_tmp;
	typedef Stepper stepper_type;
	typedef Observer observer_type;
	dense_observer(value_type abs_err = static_cast< value_type >( 1.0e-6 ),
			value_type rel_err = static_cast< value_type >( 1.0e-6 ), observer_type obs = null_observer())
			: m_abs_err(abs_err < static_cast< value_type >(0) ? -abs_err : abs_err),
				m_rel_err(rel_err < static_cast< value_type >(0) ? -rel_err : rel_err),
				m_observer(obs), m_algebra(), is_first(true) {}
	state_type interpolate_linear(const std::vector<time_type> & ts, const std::vector<state_type> & xs, time_type t) const {
		boost::optional<std::pair<std::size_t, std::size_t> > ind = find_bounds(ts, t);
		if (ind->first == ind->second)
			return xs[ind->first];
		const time_type t0 = ts[ind->first];
		const time_type t1 = ts[ind->second];
		const time_type dt = ( t1 - t0 );
		const value_type theta = ( t - t0 ) / dt;
		const value_type mtheta_1 = static_cast< value_type >( 1 ) - theta;
		state_type out;
		m_algebra.for_each3( out , xs[ind->first] , xs[ind->second] ,
				typename operations_type::template scale_sum2< value_type , time_type >( mtheta_1 , theta ) );
		return out;
	}

	void operator() (state_type x, time_type t) {
		m_observer(x, t);
	}

	void operator() (stepper_type stepper, dense_output_stepper_tag) {
		if (!is_first) { // first call of the observer is before the first step is made
			state_type x = stepper.current_state();
			time_type t = stepper.current_time();
			if (t == t_prev) {
				x_prev = t;
				m_observer(x_prev, t_prev);
				return;
			}
			// Determining points suitable for linear interpolation of the solution.
			std::size_t Nsteps = decltype(stepper)::stepper_type::order_value;
			bool first_pass = false;
			bool increasing = false;
			std::vector<state_type> points_x, points_prev_x;
			std::vector<time_type> points_t, points_prev_t;

			points_prev_t.push_back(t_prev);
			points_prev_x.push_back(x_prev);
			points_prev_t.push_back(t);
			points_prev_x.push_back(x);
			while (Nsteps > 1) {
				points_t.reserve(Nsteps + 1);
				points_x.reserve(Nsteps + 1);
				points_t.push_back(t_prev);
				points_x.push_back(x_prev);
				value_type max_abs_error = static_cast< value_type >(0);
				value_type max_rel_error = static_cast< value_type >(0);
				time_type dt = (t - t_prev) / Nsteps;
				if (dt <= std::numeric_limits<double>::epsilon())
					break;
				for (std::size_t i = 1; i != Nsteps; ++i) {
					time_type t_inter = (i < Nsteps / 2) ?  t_prev + i * dt : t - (Nsteps - i) * dt;
					state_type x_inter;
					stepper.calc_state(t_inter, x_inter);
					state_type x_lin = interpolate_linear(points_prev_t, points_prev_x, t_inter);
					points_t.push_back(t_inter);
					points_x.push_back(x_inter);
					state_type dx;
					m_algebra.for_each3( dx , x_inter , x_lin ,
								typename operations_type::template scale_sum2< value_type , time_type >(
										static_cast< value_type >( 1 ) , static_cast< value_type >( -1 ) ) );
					value_type abs_err = m_algebra.norm_inf( dx );
					value_type rel_err = abs_err / m_algebra.norm_inf( x_inter );
					max_abs_error = std::max(max_abs_error, abs_err);
					max_rel_error = std::max(max_rel_error, rel_err);
				}
				points_t.push_back(t);
				points_x.push_back(x);
				bool condition = (m_abs_err == static_cast< value_type >(0) ?
						false : max_abs_error / m_abs_err < static_cast< value_type >(1) + std::numeric_limits<value_type>::epsilon()) ||
						(m_rel_err == static_cast< value_type >(0) ?
						false : max_rel_error / m_rel_err < static_cast< value_type >(1) + std::numeric_limits<value_type>::epsilon());
				if (!condition) {
					if (!first_pass)
						increasing = true;
					else {
						if (!increasing) { // Nsteps was decreasing, but has to increase this step: optimum found at the previous iteration.
							points_t = points_prev_t;
							points_x = points_prev_x;
							break;
						}
					}
					Nsteps *= 2;
				} else {
					if (!first_pass) {
						if (increasing)	// Nsteps was increasing, but has to decrease this step: optimum found at this iteration.
							break;
					} else {
						increasing = false;
					}
					Nsteps = (Nsteps % 2 == 0 ? Nsteps / 2 : (Nsteps + 1) / 2);
				}
				points_prev_t = points_t;
				points_prev_x = points_x;
				points_t.clear();
				points_x.clear();
				first_pass = false;
			}
			if (points_t.size() > 2) {
				for (std::size_t i = 1, i_end_ = points_t.size() - 1; i != i_end_; ++i) {
					m_observer(points_x[i], points_t[i]);
				}
			}
		} else {
			is_first = m_abs_err < std::numeric_limits<value_type>::epsilon() && m_rel_err < std::numeric_limits<value_type>::epsilon();
		}
		x_prev = stepper.current_state();
		t_prev = stepper.current_time();
		m_observer(x_prev, t_prev);
	}

protected:
	value_type m_abs_err;
	value_type m_rel_err;
	observer_type m_observer;
	algebra_type m_algebra;
	bool is_first;
	state_type x_prev;
	time_type t_prev;
};

} // namespace detail
} // odeint
} // numeric
} // boost

#endif // DENSE_OBSERVER_HPP_
