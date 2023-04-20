#ifndef INTEGRATION_CONDITIONAL_H
#define INTEGRATION_CONDITIONAL_H

// Extension of boost/numeric/odeint/integrate/integrate_adaptive.hpp
// for adaptive integration of ODEs. Instead of fixed time t1, integration is
// performed until some condition is satisfied (e.g. integral (state) >= max_value, time >= t1, etc.)

#include <boost/type_traits/is_same.hpp>

#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>
#include <boost/numeric/odeint.hpp>
#include "DenseObserver.hpp"

namespace boost {
namespace numeric {
namespace odeint {

namespace detail {

/* TODO: Stopper concept also requires Trimming concept (or method?) so that
 * the stopping condition is satisfied exactly at the end point.
 * For this saving the previous state (step) and the current one is required.
 */


// forward declaration
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_const(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt ,
        Observer observer , stepper_tag );

/*
 * integrate adaptive for controlled stepper
 */
template< class Stepper , class System , class State , class Time , class Stopper , class Observer >
size_t integrate_adaptive_conditional(
        Stepper stepper , System system , State &start_state ,
        Time &start_time , Stopper stopper , Time &dt ,
        Observer observer , controlled_stepper_tag
)
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;
    typename odeint::unwrap_reference< Stopper >::type &stop = stopper;

    failed_step_checker fail_checker;  // to throw a runtime_error if step size adjustment fails
    size_t count = 0;
    while( !stop(start_state , start_time , dt ) )
    {
        obs( start_state , start_time );
        //if( less_with_sign( end_time , static_cast<Time>(start_time + dt) , dt ) )
        //{
        //    dt = end_time - start_time;
        //}

        controlled_step_result res;
        do
        {
            res = st.try_step( system , start_state , start_time , dt );
            fail_checker();  // check number of failed steps
        }
        while( res == fail );
        fail_checker.reset();  // if we reach here, the step was successful -> reset fail checker

        ++count;
    }
    obs( start_state , start_time );
    return count;
}

/*
 * integrate adaptive for dense output steppers.
 * the whole stepper is passed to the observer in case intermediate values should be recorded.
 *
 * step size control is used if the stepper supports it
 */
template< class Stepper , class System , class State , class Time , class Stopper , class Observer >
size_t integrate_adaptive_conditional(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Stopper stopper , Time dt ,
        Observer observer , dense_output_stepper_tag )
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;
    typename odeint::unwrap_reference< Stopper >::type &stop = stopper;

    size_t count = 0;
    st.initialize( start_state , start_time , dt );

    while( !stop(st.current_state() , st.current_time() , st.current_time_step() ) )
    {
    	//while( less_eq_with_sign( static_cast<Time>(st.current_time() + st.current_time_step()) ,
			//			 end_time ,
			//			 st.current_time_step() ) )
			//{   //make sure we don't go beyond the end_time
					obs(st.current_state(), st.current_time());
					st.do_step( system );
					++count;
			//}
			// calculate time step to arrive exactly at end time
			//st.initialize( st.current_state() , st.current_time() , static_cast<Time>(end_time - st.current_time()) );
    }
    obs(st.current_state(), st.current_time());
    // overwrite start_state with the final point
    boost::numeric::odeint::copy( st.current_state() , start_state );
    return count;
}

/*
 * integrate adaptive for dense output steppers.
 * the whole stepper is passed to the observer in case intermediate values should be recorded.
 *
 * step size control is used if the stepper supports it
 */
template< class Stepper , class System , class State , class Time , class Stopper , class Observer >
size_t integrate_adaptive_conditional_dense(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Stopper stopper , Time dt ,
        Observer observer , dense_output_stepper_tag )
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;
    typename odeint::unwrap_reference< Stopper >::type &stop = stopper;

    size_t count = 0;
    st.initialize( start_state , start_time , dt );

    while( !stop(st.current_state() , st.current_time() , st.current_time_step() ) )
    {
    	//while( less_eq_with_sign( static_cast<Time>(st.current_time() + st.current_time_step()) ,
			//			 end_time ,
			//			 st.current_time_step() ) )
			//{   //make sure we don't go beyond the end_time
					obs(st, dense_output_stepper_tag());
					st.do_step( system );
					++count;
			//}
			// calculate time step to arrive exactly at end time
			//st.initialize( st.current_state() , st.current_time() , static_cast<Time>(end_time - st.current_time()) );
    }
    obs(st, dense_output_stepper_tag());
    // overwrite start_state with the final point
    boost::numeric::odeint::copy( st.current_state() , start_state );
    return count;
}

/*
 * integrate adaptive for dense output steppers.
 * the whole stepper is passed to the observer in case intermediate values should be recorded.
 *
 * step size control is used if the stepper supports it
 */
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive_dense(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt ,
        Observer observer , dense_output_stepper_tag )
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;

    size_t count = 0;
    st.initialize( start_state , start_time , dt );

    obs( st, dense_output_stepper_tag() );
    while( less_with_sign( st.current_time() , end_time , st.current_time_step() ) )
		{
				while( less_eq_with_sign( static_cast<Time>(st.current_time() + st.current_time_step()) ,
							 end_time ,
							 st.current_time_step() ) )
				{   //make sure we don't go beyond the end_time
						st.do_step( system );
						obs( st, dense_output_stepper_tag() );
						++count;
				}
				// calculate time step to arrive exactly at end time
				st.initialize( st.current_state() , st.current_time() , static_cast<Time>(end_time - st.current_time()) );
		}
    // overwrite start_state with the final point
    boost::numeric::odeint::copy( st.current_state() , start_state );
    return count;
}

} // namespace detail

/*
 * the two overloads are needed in order to solve the forwarding problem
 */
template< class Stepper , class System , class State , class Time , class Stopper , class Observer >
size_t integrate_adaptive_conditional(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Stopper stopper , Time dt ,
        Observer observer )
{
    typedef typename odeint::unwrap_reference< Stepper >::type::stepper_category stepper_category;
    return detail::integrate_adaptive_conditional(
            stepper , system , start_state ,
            start_time , stopper , dt ,
            observer , stepper_category() );
}

/**
 * \brief Second version to solve the forwarding problem,
 * can be called with Boost.Range as start_state.
 */
template< class Stepper , class System , class State , class Time , class Stopper , class Observer >
size_t integrate_adaptive_conditional(
        Stepper stepper , System system , const State &start_state ,
        Time start_time , Stopper stopper , Time dt ,
        Observer observer )
{
    typedef typename odeint::unwrap_reference< Stepper >::type::stepper_category stepper_category;
    return detail::integrate_adaptive(
            stepper , system , start_state ,
            start_time , stopper , dt ,
            observer , stepper_category() );
}

/*
 * the two overloads are needed in order to solve the forwarding problem
 */
template< class Stepper , class System , class State , class Time , class Stopper , class Observer >
size_t integrate_adaptive_conditional_dense(
        Stepper stepper , typename Stepper::value_type abs_tol, typename Stepper::value_type rel_tol,
				System  system , State &start_state , Time start_time , Stopper stopper , Time dt ,
        Observer observer )
{
	typedef typename result_of::make_dense_output< Stepper >::type new_stepper_type;
	typedef detail::dense_observer<new_stepper_type, Observer> observer_wrapper;
	new_stepper_type step = make_dense_output(abs_tol, rel_tol, stepper);
	observer_wrapper obs_wrapper(abs_tol, rel_tol, observer);
	typedef typename odeint::unwrap_reference< new_stepper_type >::type::stepper_category stepper_category;
	return detail::integrate_adaptive_conditional_dense(
					step , system , start_state ,
					start_time , stopper , dt ,
					obs_wrapper, stepper_category() );
}

/**
 * \brief Second version to solve the forwarding problem,
 * can be called with Boost.Range as start_state.
 */
template< class Stepper , class System , class State , class Time , class Stopper , class Observer >
size_t integrate_adaptive_conditional_dense(
        Stepper stepper , typename Stepper::value_type abs_tol, typename Stepper::value_type rel_tol,
				System  system , const State &start_state , Time start_time , Stopper stopper , Time dt ,
        Observer observer )
{
	typedef typename result_of::make_dense_output< Stepper >::type new_stepper_type;
	typedef detail::dense_observer<new_stepper_type, Observer> observer_wrapper;
	new_stepper_type step = make_dense_output(abs_tol, rel_tol, stepper);
	observer_wrapper obs_wrapper(abs_tol, rel_tol, observer);
	typedef typename odeint::unwrap_reference< new_stepper_type >::type::stepper_category stepper_category;
	return detail::integrate_adaptive_conditional_dense(
					step , system , start_state ,
					start_time , stopper , dt ,
					obs_wrapper, stepper_category() );
}

/*
 * the two overloads are needed in order to solve the forwarding problem
 */
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive_dense(
				Stepper stepper , typename Stepper::value_type abs_tol, typename Stepper::value_type rel_tol,
				System  system , State &start_state , Time start_time , Time end_time , Time dt ,
				Observer observer )
{
		typedef typename result_of::make_dense_output< Stepper >::type new_stepper_type;
		typedef detail::dense_observer<new_stepper_type, Observer> observer_wrapper;
		new_stepper_type step = make_dense_output(abs_tol, rel_tol, stepper);
		observer_wrapper obs_wrapper(abs_tol, rel_tol, observer);
    typedef typename odeint::unwrap_reference< new_stepper_type >::type::stepper_category stepper_category;
    return detail::integrate_adaptive_dense(
    				step , system , start_state ,
            start_time , end_time , dt ,
						obs_wrapper , stepper_category() );
}

/**
 * \brief Second version to solve the forwarding problem,
 * can be called with Boost.Range as start_state.
 */
template< class Stepper , class System , class State , class Time , class Observer >
size_t integrate_adaptive_dense(
				Stepper stepper , typename Stepper::value_type abs_tol, typename Stepper::value_type rel_tol,
				System  system , const State &start_state , Time start_time , Time end_time , Time dt ,
				Observer observer )
{
		typedef typename result_of::make_dense_output< Stepper >::type new_stepper_type;
		typedef detail::dense_observer<new_stepper_type, Observer> observer_wrapper;
		new_stepper_type step = make_dense_output(abs_tol, rel_tol, stepper);
		observer_wrapper obs_wrapper(abs_tol, rel_tol, observer);
    typedef typename odeint::unwrap_reference< new_stepper_type >::type::stepper_category stepper_category;
    return detail::integrate_adaptive_dense(
    				step , system , start_state ,
            start_time , end_time , dt ,
						obs_wrapper , stepper_category() );
}

/**
 * \brief integrate_adaptive without an observer.
 */
template< class Stepper , class System , class State , class Time , class Stopper >
size_t integrate_adaptive_conditional(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Stopper stopper , Time dt )
{
    return integrate_adaptive_conditional( stepper , system , start_state , start_time , stopper , dt , null_observer() );
}

/**
 * \brief Second version to solve the forwarding problem,
 * can be called with Boost.Range as start_state.
 */
template< class Stepper , class System , class State , class Time , class Stopper >
size_t integrate_adaptive_conditional(
        Stepper stepper , System system , const State &start_state ,
        Time start_time , Stopper stopper , Time dt )
{
    return integrate_adaptive_conditional( stepper , system , start_state , start_time , stopper , dt , null_observer() );
}

} // namespace odeint
} // namespace numeric
} // namespace boost

#endif // INTEGRATION_CONDITIONAL_H
