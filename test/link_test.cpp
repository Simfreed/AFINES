#include "Link.h"
#include "filament.h"
#define BOOST_TEST_MODULE link_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nactin              = 2;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double tol = 0.001;
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    actin a; 
   
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    
    Link * l = f->get_link(0);
    BOOST_CHECK_EQUAL( l->get_kl(), 100);                  
    BOOST_CHECK_EQUAL( l->get_hx()[0], 0);             
    BOOST_CHECK_EQUAL( l->get_hx()[1],  1);            
    BOOST_CHECK_EQUAL( l->get_hy()[0], 0);              
    BOOST_CHECK_EQUAL( l->get_hy()[1], 0);               

    BOOST_CHECK_CLOSE( l->get_xcm(), 0.5, tol);
    BOOST_CHECK_CLOSE( l->get_ycm(), 0, tol);
    
    delete f;
} 

BOOST_AUTO_TEST_CASE( get_distance_test)
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nactin              = 2;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double tol = 0.001;
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    actin a; 
   
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    
    Link * l = f->get_link(0);
    
    BOOST_CHECK_CLOSE(l->get_distance(0.75, 2), 2, tol);
    BOOST_CHECK_CLOSE(l->get_distance(4, -4), 5, tol);
    BOOST_CHECK_CLOSE(l->get_distance(-5, 12), 13, tol);
    
    delete f;
}

BOOST_AUTO_TEST_CASE( get_quadrants_test )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nactin              = 2;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    //TEST 1
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
   
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    
    Link * l = f->get_link(0);
    l->quad_update();
    vector<array<int,2> > quads = l->get_quadrants();
    vector<array<int,2> > expected_quads;
    expected_quads.push_back({0,0});
    expected_quads.push_back({1,0});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nExpected Quadrants : don't equal Link Quadrants : \n");
    cout<<"\nLink Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;

    //TEST 2
    link_length = 1.5;
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    
    l = f->get_link(0);
    l->quad_update();
    quads = l->get_quadrants();
    expected_quads.push_back({0,0});
    expected_quads.push_back({1,0});
    expected_quads.push_back({2,0});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nExpected Quadrants : don't equal Link Quadrants : \n");
    cout<<"\nLink Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;
    
    //TEST 3
    startphi = pi/4;
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    
    l = f->get_link(0);
    l->quad_update();
    quads = l->get_quadrants();
    expected_quads.push_back({0,0});
    expected_quads.push_back({0,1});
    expected_quads.push_back({0,2});
    expected_quads.push_back({1,0});
    expected_quads.push_back({1,1});
    expected_quads.push_back({1,2});
    expected_quads.push_back({2,0});
    expected_quads.push_back({2,1});
    expected_quads.push_back({2,2});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nExpected Quadrants : don't equal Link Quadrants : \n");
    cout<<"\nLink Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;

    //TEST 4
    startphi =0; startx = 1; link_length = 1; nactin = 2;
   
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    
    l = f->get_link(0);
    l->quad_update();
    
    quads = l->get_quadrants();
    expected_quads.push_back({2,0});
    expected_quads.push_back({3,0});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nExpected Quadrants : don't equal Link Quadrants : \n");
    cout<<"\nLink Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;
    
    //TEST 5
    startphi = pi/2; startx = -0.25; starty = 0.2; link_length = 1; nactin = 2;
   
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    
    l = f->get_link(0);
    l->quad_update();
    
    quads = l->get_quadrants();
    expected_quads.push_back({-1,0});
    expected_quads.push_back({-1,1});
    expected_quads.push_back({-1,2});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nExpected Quadrants : don't equal Link Quadrants : \n");
    cout<<"\nLink Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;
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
