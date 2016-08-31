#include "globals.h"

#define BOOST_TEST_MODULE globals_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( event_test )
{

} 

BOOST_AUTO_TEST_CASE( rij_periodic_test )
{
    //double dx=0.5, dy=5, xl=10, yl=7;
//array<double, 2> = rij_periodic(dx, dy, xl, yl);

}


BOOST_AUTO_TEST_CASE( coord2quad_test)
{
    double X=50, x=24.6, n=100;
    BOOST_CHECK_EQUAL(coord2quad_floor(X,n,x),99);
    BOOST_CHECK_EQUAL(coord2quad_ceil(X,n,x),0);


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

BOOST_AUTO_TEST_CASE( set_get_xcm )
{
}
// EOF
