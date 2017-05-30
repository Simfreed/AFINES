#include "actin.h"
//#include "globals.h"
#define BOOST_TEST_MODULE actin_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    actin * a = new actin(0.25, 0.5, 0.75, 1.25);
    BOOST_CHECK_EQUAL( a->get_xcm(), 0.25);                  // 1 //
    BOOST_CHECK_EQUAL( a->get_ycm(), 0.5);                   // 1 //
    BOOST_CHECK_EQUAL( a->get_length(), 0.75);                // 1 //
    BOOST_CHECK_EQUAL( a->get_viscosity(), 1.25);               
    BOOST_CHECK_EQUAL( a->get_velocity()[0], 0);               
    BOOST_CHECK_EQUAL( a->get_velocity()[1], 0);               
    BOOST_CHECK_EQUAL( a->get_force()[0], 0);               
    BOOST_CHECK_EQUAL( a->get_force()[1], 0);               
    BOOST_CHECK_CLOSE( a->get_friction(), 17.6715, 0.001);               
    delete a;
} 

BOOST_AUTO_TEST_CASE( force_test)
{
    actin a(0.25, 0.5, 0.75, 1.25);
    array<double, 2> f1 = {3.5,-6.7};
    array<double, 2> f2 = {-2,1};
    array<double, 2> ftot = {f1[0]+f2[0],f1[1]+f2[1]};

    a.update_force(f1[0], f1[1]);
    BOOST_CHECK_EQUAL( a.get_force()[0] , f1[0]); 
    BOOST_CHECK_EQUAL( a.get_force()[1] , f1[1]); 
    a.update_force(f2[0], f2[1]);
    BOOST_CHECK_EQUAL( a.get_force()[0] , ftot[0]); 
    BOOST_CHECK_EQUAL( a.get_force()[1] , ftot[1]); 
    
    a.reset_force();
    BOOST_CHECK_EQUAL( a.get_force()[0], 0);               
    BOOST_CHECK_EQUAL( a.get_force()[1], 0);               
}

BOOST_AUTO_TEST_CASE( velocity_test)
{
    actin a(0.25, 0.5, 0.75, 1.25);
    array<double, 2> v1 = {3.5,-6.7};
    array<double, 2> v2 = {-2,1};
    array<double, 2> vtot = {v1[0]+v2[0],v1[1]+v2[1]};
    double vtotsq = dot(vtot[0], vtot[1], vtot[0], vtot[1]);

    a.update_velocity(v1[0], v1[1]);
    BOOST_CHECK_EQUAL( a.get_velocity()[0] , v1[0]); 
    BOOST_CHECK_EQUAL( a.get_velocity()[1] , v1[1]); 
    a.update_velocity(v2[0], v2[1]);
    BOOST_CHECK_EQUAL( a.get_velocity()[0] , vtot[0]); 
    BOOST_CHECK_EQUAL( a.get_velocity()[1] , vtot[1]); 
    BOOST_CHECK_EQUAL( a.get_vsquared() , vtotsq); 
    
    a.reset_velocity();
    BOOST_CHECK_EQUAL( a.get_velocity()[0], 0);               
    BOOST_CHECK_EQUAL( a.get_velocity()[1], 0);               
}

BOOST_AUTO_TEST_CASE( pos_test)
{
    actin a(0.25, 0.5, 0.75, 1);
    a.set_xcm(-0.45);
    a.set_ycm(5.67);
    
    BOOST_CHECK_EQUAL( a.get_xcm(), -0.45);
    BOOST_CHECK_EQUAL( a.get_ycm() , 5.67);
}
// EOF
