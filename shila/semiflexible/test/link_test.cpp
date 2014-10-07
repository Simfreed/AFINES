#include "Link.h"
#include "actin_ensemble.h"
#define BOOST_TEST_MODULE link_test
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    actin * a1 = new actin(-1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin( 1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a1, 0);
    ae.add_monomer(a2, 0);
    Link l( 1, 100, 100, &ae, 0, 1, "green");

    double tol = 0.001;

    BOOST_CHECK_EQUAL( l.get_kb(), 100);                  // 1 //
    BOOST_CHECK_EQUAL( l.get_hx()[0], -0.5);                   // 1 //
    BOOST_CHECK_EQUAL( l.get_hx()[1],  0.5);                   // 1 //
    BOOST_CHECK_EQUAL( l.get_hy()[0], 2);                // 1 //
    BOOST_CHECK_EQUAL( l.get_hy()[1], 2);               
    BOOST_CHECK_EQUAL( l.get_color(), "green");

    BOOST_CHECK_CLOSE( l.get_posx(), 0, tol);
    BOOST_CHECK_CLOSE( l.get_posy(), 2, tol);
    
} 

BOOST_AUTO_TEST_CASE( actin_update )
{
    actin * a1 = new actin(-1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin( 1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a1, 0);
    ae.add_monomer(a2, 0);
    Link l( 1, 100, 100, &ae, 0, 1, "green");
    l.actin_update();

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
