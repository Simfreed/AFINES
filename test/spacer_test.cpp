#include "spacer.h"
#include "filament_ensemble.h"
#include "globals.h"
#define BOOST_TEST_MODULE spacer_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{

    double tol = 0.001;
    double actin_density = 0;
    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2}, state = {0,0}, findex = {0,0}, lindex = {0, 0};
    vector<array<double, 3> > pos_sets;
    int nactin = 2;

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double actin_rad = 0.5, link_len = 1, spacer_len = 1;
    double v0 = 0;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 0, kend = 0, fstall = 3.85, fbreak = 3.85, ebind = 0.04;
    double frac_force = 0;
    
    
    filament_ensemble * f = new filament_ensemble(actin_density, fov, nq, dt, temp, 
            actin_rad, vis, nactin, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed);
    spacer  m = spacer(array<double, 3>{1, 1, 3.1416/2}, spacer_len, f, state, 
            findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, ebind, vis, bc);

    BOOST_CHECK_CLOSE( m.get_hx()[0],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hx()[1],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[0], 0.5, tol);                // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[1], 1.5, tol);               
    BOOST_CHECK_EQUAL( m.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m.get_states()[1], 0);
    
    delete f;
    
    actin_density = nactin/(fov[0]*fov[1]);
    pos_sets.push_back({-0.5,0,0});


    f = new filament_ensemble(actin_density, fov, nq, dt, temp, actin_rad, vis, nactin, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed);
    state = {1,0};

    spacer m2 = spacer(array<double, 3>{1, 1, 0}, spacer_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, ebind, vis, bc);
    BOOST_CHECK_EQUAL( m2.get_states()[0], 1);
    BOOST_CHECK_EQUAL( m2.get_states()[1], 0);
    
    state = {0,1};
    spacer m3 = spacer(array<double, 3>{1, 1, 0}, spacer_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, ebind, vis, bc);
    BOOST_CHECK_EQUAL( m3.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m3.get_states()[1], 1);
   
    delete f;
} 


BOOST_AUTO_TEST_CASE( update_bending_one_attached)
{
    //Filament ENSEMBLE
    array<double, 2> fov = {4,4};
    vector<vector<double> > actin_sets;
    array<int, 2> nq = {8,8}, state = {1,0}, findex = {0,-1}, lindex = {0, 0};

    double dt = 0.0001, temp = 0, vis = 0.001;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double stretching = 1, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, fbreak = 3.85, ebind = 0.04;
    double frac_force = 100000;
//    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={0, 0, actin_rad,0}, pos2={1, 0,actin_rad,0};
    
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    
    /*****************************
     * TEST 1: motor ang = pi /4 *
     * **************************/

    filament_ensemble * f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc);
    
    double mx = 0.5+sqrt(2)/4, my = sqrt(2)/4, mang = pi/4, mlen = 1, mstiff=1;
    spacer m = spacer(array<double, 3>{mx, my, mang}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, ebind, vis, bc);
    m.set_bending(0.04,pi/2); 
    
    double ang0, ang1;
    int tf = 1000;
    ang0 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.actin_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        ang1 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);

        BOOST_CHECK_MESSAGE(abs(angBC(pi/2-ang0, 2*pi)) >= abs(angBC(pi/2-ang1, 2*pi)), "ang0 = "<<ang0<<" is closer to pi/2 than ang1 = "<<ang1<<" at time "<<t);
        ang0 = ang1;
    }

    // delete f;
    
    /*****************************
     * TEST 2: motor ang = -pi /4 *
     * **************************/

    f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc);
    
    my = -sqrt(2)/4, mang = -pi/4.0;
    m = spacer(array<double, 3>{mx, my, mang}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, ebind, vis, bc);
    m.set_bending(0.04,pi/2); 
    tf =1000; 
    ang0 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.actin_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        ang1 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);

        BOOST_CHECK_MESSAGE(abs(angBC(-pi/2-ang0, 2*pi)) >= abs(angBC(-pi/2-ang1, 2*pi)), "TEST 2: ang0 = "<<ang0<<" is closer to -pi/2 than ang1 = "<<ang1<<" at time "<<t);
        ang0 = ang1;
    }

    delete f;
    
    /*****************************
     * TEST 3: motor ang = -3pi /4 *
     * **************************/

    f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc);
    
    mx = 0.5-sqrt(2)/4, mang = -3*pi/4.0;
    m = spacer(array<double, 3>{mx, my, mang}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, ebind, vis, bc);
    m.set_bending(0.04,pi/2); 
    tf =1000; 
    ang0 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.actin_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        ang1 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);

        BOOST_CHECK_MESSAGE(abs(angBC(-pi/2-ang0, 2*pi)) >= abs(angBC(-pi/2-ang1, 2*pi)), "TEST 3: ang0 = "<<ang0<<" is closer to -pi/2 than ang1 = "<<ang1<<" at time "<<t);
        ang0 = ang1;
    }

    delete f;
    
    /*****************************
     * TEST 4: motor ang = 3pi /4 *
     * **************************/

    f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc);
    
    my = sqrt(2)/4, mang = 3*pi/4.0;
    m = spacer(array<double, 3>{mx, my, mang}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, ebind, vis, bc);
    m.set_bending(0.04,pi/2); 
    tf =1000; 
    ang0 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.actin_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        ang1 = angBC(m.get_angle() - f->get_angle(0,0), 2*pi);

        BOOST_CHECK_MESSAGE(abs(angBC(pi/2-ang0, 2*pi)) >= abs(angBC(pi/2-ang1, 2*pi)), "TEST 4: ang0 = "<<ang0<<" is closer to pi/2 than ang1 = "<<ang1<<" at time "<<t);
        ang0 = ang1;
    }

    delete f;
}   

/* Functions to test : 
        void attach(int hd);

        void brownian(double t, double gamma);

        void step_twoheads();

        void actin_update();

        void update_shape();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd);

        void update_pos_a_end(int hd, double pos);

*/
// EOF
