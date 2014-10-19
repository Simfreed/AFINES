#define BOOST_TEST_MODULE actin_test

#include "actin.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    actin * a = new actin(0.25, 0.5, 0.75, 1, 2, 2, 4, 4, 0.2); 
    BOOST_CHECK_EQUAL( a->get_xcm(), 0.25);                  // 1 //
    BOOST_CHECK_EQUAL( a->get_ycm(), 0.5);                   // 1 //
    BOOST_CHECK_EQUAL( a->get_angle(), 0.75);                // 1 //
    BOOST_CHECK_EQUAL( a->get_length(), 1);               
    delete a;
} 

BOOST_AUTO_TEST_CASE( friction_test )
{
    actin a(0.25, 0.5, 0.75, 1, 2, 2, 4, 4, 0.2);
    double * fric = a.get_friction();
    double tol = 0.001; //%
    BOOST_CHECK_CLOSE( fric[0] , 0.340655, tol);
    BOOST_CHECK_CLOSE( fric[1] , 0.681311, tol);
    BOOST_CHECK_CLOSE( fric[2] , 0.085164, tol);
    delete[] fric;
}

BOOST_AUTO_TEST_CASE( force_test)
{
    actin a(0.25, 0.5, 0.75, 1, 2, 2, 4, 4, 0.2);
    BOOST_CHECK_EQUAL( a.get_forces()[0] , 0); 
    BOOST_CHECK_EQUAL( a.get_forces()[1] , 0); 
    BOOST_CHECK_EQUAL( a.get_forces()[2] , 0); 

    a.update_force(11,22,33);
    BOOST_CHECK_EQUAL( a.get_forces()[0] , 11); 
    BOOST_CHECK_EQUAL( a.get_forces()[1] , 22); 
    BOOST_CHECK_EQUAL( a.get_forces()[2] , 33); 

}

BOOST_AUTO_TEST_CASE( direction_test)
{
    actin a(0.25, 0.5, 0.75, 1, 2, 2, 4, 4, 0.2);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( a.get_direction()[0] , 0.731689 , tol); 
    BOOST_CHECK_CLOSE( a.get_direction()[1] , 0.681639 , tol); 
}

BOOST_AUTO_TEST_CASE( start_end_test)
{
    actin a(0, 0, 0, 1, 0, 0, 0, 0, 0);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( a.get_start()[0], -0.500000, tol); 
    BOOST_CHECK_CLOSE( a.get_start()[1],  0.000000, tol); 
    BOOST_CHECK_CLOSE( a.get_end()[0],  0.500000, tol); 
    BOOST_CHECK_CLOSE( a.get_end()[1],  0.000000, tol); 
    
}

BOOST_AUTO_TEST_CASE( get_intpoint_test)
{
    actin a(0, 0, 0, 1, 0, 0, 0, 0, 0);
    double x = 0.1, y = 1;
    double tol = 0.001;
    double * intpoint = a.get_intpoint(x,y);
    BOOST_CHECK_CLOSE( intpoint[0], 0.1, tol);
    BOOST_CHECK_CLOSE( intpoint[1], 0, tol);
    delete[] intpoint;
}

BOOST_AUTO_TEST_CASE( set_get_xcm )
{
    actin a(1, 0, 0, 1, 0, 0, 0, 0, 0);
}
// EOF
