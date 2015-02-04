#include "Link.h"
#include "filament.h"
#define BOOST_TEST_MODULE link_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    actin * a1 = new actin(-1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin( 1, 2, 0, 1, 0, 0, 0, 0, 0);
    double link_length = 0;
    double stretching_stiffness = 100;
    double bending_stiffness = 1000;
    filament f;
    f.add_rod(a1, link_length, stretching_stiffness, bending_stiffness);
  
    f.add_rod(a2, link_length, stretching_stiffness, bending_stiffness);

    Link l( 1, 100, 100, &f, 0, 1);

    double tol = 0.001;

    BOOST_CHECK_EQUAL( l.get_kb(), 100);                  // 1 //
    BOOST_CHECK_EQUAL( l.get_hx()[0], -0.5);                   // 1 //
    BOOST_CHECK_EQUAL( l.get_hx()[1],  0.5);                   // 1 //
    BOOST_CHECK_EQUAL( l.get_hy()[0], 2);                // 1 //
    BOOST_CHECK_EQUAL( l.get_hy()[1], 2);               

    BOOST_CHECK_CLOSE( l.get_xcm(), 0, tol);
    BOOST_CHECK_CLOSE( l.get_ycm(), 2, tol);
    
    delete a1;
    delete a2;
} 

BOOST_AUTO_TEST_CASE( filament_update )
{
    
    actin * a1 = new actin(-1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin( 1, 2, 0, 1, 0, 0, 0, 0, 0);
    double link_length = 0;
    double stretching_stiffness = 100;
    double bending_stiffness = 1000;
    filament f;
    f.add_rod(a1, link_length, stretching_stiffness, bending_stiffness);
    f.add_rod(a2, link_length, stretching_stiffness, bending_stiffness);
    
    Link l( 1, 100, 100, &f, 0, 1);
    l.filament_update();
    
    delete a1;
    delete a2;

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
