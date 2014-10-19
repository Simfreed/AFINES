#include "motor.h"
#include "motor_ensemble.h"
#include "actin_ensemble.h"
#define BOOST_TEST_MODULE motor_ensemble_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    actin_ensemble ae = actin_ensemble();
    motor m = motor(1, 1, 3.1415926535897/2, 1, &ae, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "blue");
    double tol = 0.001;

    BOOST_CHECK_CLOSE( m.get_hx()[0],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hx()[1],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[0], 0.5, tol);                // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[1], 1.5, tol);               
    BOOST_CHECK_EQUAL( m.get_color(), "blue");
    BOOST_CHECK_EQUAL( m.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m.get_states()[1], 0);
    
    actin * a = new actin(0, 0, 0, 1, 0, 0, 0, 0, 0);
    ae.add_monomer(a, 0);

    motor m2 = motor(1, 1, 0, 1, &ae, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "blue");
    BOOST_CHECK_EQUAL( m2.get_states()[0], 1);
    BOOST_CHECK_EQUAL( m2.get_states()[1], 0);
    
    motor m3 = motor(1, 1, 0, 1, &ae, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "blue");
    BOOST_CHECK_EQUAL( m3.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m3.get_states()[1], 1);
    
} 

BOOST_AUTO_TEST_CASE( friction_test )
{
}

BOOST_AUTO_TEST_CASE( force_test)
{
}

BOOST_AUTO_TEST_CASE( direction_test)
{
}

BOOST_AUTO_TEST_CASE( start_end_test)
{
}

BOOST_AUTO_TEST_CASE( get_intpoint_test)
{
}
// EOF
