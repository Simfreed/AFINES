#include "motor.h"
#include "filament_ensemble.h"
#include "globals.h"
#define BOOST_TEST_MODULE motor_test
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
    
    double actin_rad = 0.5, link_len = 1, motor_len = 1;
    double v0 = 0;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 0, kend = 0;
    double frac_force = 0;
    
    
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_density, fov, nq, dt, temp, 
            actin_rad, vis, nactin, link_len, pos_sets, stretching, bending, frac_force, bc, seed);
    motor<ATfilament_ensemble>  m = motor<ATfilament_ensemble>({1, 1, 3.1416/2}, motor_len, f, state, 
            findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);

    BOOST_CHECK_CLOSE( m.get_hx()[0],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hx()[1],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[0], 0.5, tol);                // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[1], 1.5, tol);               
    BOOST_CHECK_EQUAL( m.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m.get_states()[1], 0);
    
    actin_density = nactin/(fov[0]*fov[1]);
    pos_sets.push_back({-0.5,0,0});

    delete f;

    f = new ATfilament_ensemble(actin_density, fov, nq, dt, temp, actin_rad, vis, nactin, link_len, pos_sets, stretching, bending, frac_force, bc, seed);
    state = {1,0};

    motor<ATfilament_ensemble> m2 = motor<ATfilament_ensemble>({1, 1, 0}, motor_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL( m2.get_states()[0], 1);
    BOOST_CHECK_EQUAL( m2.get_states()[1], 0);
    
    state = {0,1};
    motor<ATfilament_ensemble> m3 = motor<ATfilament_ensemble>({1, 1, 0}, motor_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL( m3.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m3.get_states()[1], 1);
   
    delete f;
} 


BOOST_AUTO_TEST_CASE( step_onehead )
{

    //Filament ENSEMBLE
    double tol = 0.001;
    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2}, state = {1,0}, findex = {0,-1}, lindex = {0, -1};
    vector<array<double, 3> > pos_sets;
    int nactin = 2;
    double actin_density = nactin/(fov[0]*fov[1]);

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 2, kend = 0;
    double frac_force = 0;
    
    pos_sets.push_back({-0.4,0,0});
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_density, fov, nq, dt, temp, actin_rad, vis, nactin, link_len, pos_sets, stretching, bending, frac_force, bc, seed);

    //MOTOR
    double mx = 0.5, my = 0.5, mang = pi/2, mlen = 1;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);

    //IF detach check that:
    //(a) new state of head is correct (was 1, now 0)
    //(b) position of head relative to end of rod is correct (0)
    //(c) position of head is correct
    //
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.1, tol);
    m.step_onehead(0); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[0], 0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], my - 0.5*mlen*sin(mang), tol);

    // IF no detach
    // (a) new position is correct
    koff = -1;
    kend = -1;
    m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.1, tol);
    m.step_onehead(0); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.35, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.25, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);

    delete f;
}

BOOST_AUTO_TEST_CASE( move_end_detach )
{
    //Filament ENSEMBLE
    double tol = 0.001;
    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2}, state = {1,0}, findex = {0,-1}, lindex = {0, -1};
    vector<array<double, 3> > pos_sets;
    int nactin = 4;
    double actin_density = nactin/(fov[0]*fov[1]);

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 2, kend = 2;
    double frac_force = 0;
    
    pos_sets.push_back({0,0,0});
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_density, fov, nq, dt, temp, actin_rad, vis, nactin, link_len, pos_sets, stretching, bending, frac_force, bc, seed);

    //MOTOR
    double mx = 0.125, my = 0.5, mang = pi/2, mlen = 1;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);

    // IF position is greater than filament length
    double pos = 1.125;
    //  If last rod in filament
    //      try to detach
    //      if detach:
    //          (a) new state of head is correct (was 1, now 0)
    //          (b) position of head relative to end of rod is correct (0)
    //          (c) position of head is correct

    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[0], 0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], my - 0.5*mlen*sin(mang), tol);
    
    //      else (don't detach): 
    //          stay
    kend = -1;
    m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.125, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    //  Else (not last rod in filament)
    mx = 1.125;
    lindex = {1, -1};
    m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    
    //      (a) motor moves to next rod in filament
    //      (b) motor moves appropriate distance
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.125, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.875, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    // Else ( new position isn't more than filament length
    //  (a) motor moves to appropriate new position
    pos = 0.375;
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], pos, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.625, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    delete f;
}

BOOST_AUTO_TEST_CASE( attach )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2}, state = {0,0}, findex = {-1,-1}, lindex = {-1, -1};
    vector<array<double, 3> > pos_sets;
    int nactin = 4;
    int nfil = 1;//3;
    double actin_density = nfil*nactin/(fov[0]*fov[1]);

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0;
    double frac_force = 0;
    
    //pos_sets.push_back({1,1,pi/2});
    pos_sets.push_back({0,0,0});
    //pos_sets.push_back({-2,3,pi});
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_density, fov, nq, dt, temp, actin_rad, vis, nactin, link_len, pos_sets, stretching, bending, frac_force, bc, seed);
    f->quad_update(); 
    //MOTOR
    double mx = 2.35, my = 0.5, mang = pi/2, mlen = 1;
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {f_index, l_index} = {1, 2}, because the distance between head 0 and that link is 0
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    m.attach(1); //shouldn't attach because the distance is greater than max binding distance = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    
    kon = -1;
    m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {f_index, l_index} = {0, 3}, because the distance between head 0 and that link is 0
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    m.attach(1); //shouldn't attach because the distance is greater than max binding distance = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
}   

//moving across a boundary

BOOST_AUTO_TEST_CASE( step_onehead_periodic )
{

    //Filament ENSEMBLE
    double tol = 0.001;
    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2}, state = {1,0}, findex = {0,-1}, lindex = {1, -1};
    vector<vector<double> > actin_sets;

    double dt = 1, temp = 0, vis = 0;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 2, kend = 0;
    double frac_force = 0;
    
    vector<double> pos1={23.6,0,actin_rad,0},pos2={24.6,0,actin_rad,0},pos3={-24.4,0,actin_rad,0};
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    actin_sets.push_back(pos3);
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, bending, frac_force, bc);
//    cout<<f->get_filament(0)->to_string();
    //MOTOR
    double mx = -24.875, my = 0.5, mang = pi/2, mlen = 1;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    cout<<m.to_string();
    //IF detach check that:
    //(a) new state of head is correct (was 1, now 0)
    //(b) position of head relative to end of rod is correct (0)
    //(c) position of head is correct
    //
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.475, tol);
    m.step_onehead(0); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[0], 0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], my - 0.5*mlen*sin(mang), tol);

    // IF no detach
    // (a) new position is correct
    koff = -1;
    kend = -1;
    m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.475, tol);
    m.step_onehead(0); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.725, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 24.875, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);

    delete f;
}
//Check for attaching across a boundary

BOOST_AUTO_TEST_CASE( attach_periodic )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2}, state = {0,0}, findex = {-1,-1}, lindex = {-1, -1};
    vector<vector<double> > actin_sets;

    double dt = 1, temp = 0, vis = 0;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0;
    double frac_force = 0;
    
    //pos_sets.push_back({0,0,0});
    vector<double> pos1={0,24.875,actin_rad,0},pos2={1,24.875,actin_rad,0},pos3={2,24.875,actin_rad,0}, pos4={3,24.875,actin_rad,0};
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    actin_sets.push_back(pos3);
    actin_sets.push_back(pos4);
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, bending, frac_force, bc);
    f->quad_update(); 
    
    //MOTOR
    double mx = 2.35, my = -24.625, mang = pi/2, mlen = 1;
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {f_index, l_index} = {1, 2}, because the distance between head 0 and that link is 0
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    m.attach(1); //shouldn't attach because the distance is greater than max binding distance = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    
    kon = -1;
    m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {f_index, l_index} = {0, 3}, because the distance between head 0 and that link is 0
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    m.attach(1); //shouldn't attach because the distance is greater than max binding distance = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
}   

//stretch is as expected

BOOST_AUTO_TEST_CASE( attach_twoheads_periodic )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2}, state = {0,0}, findex = {-1,-1}, lindex = {-1, -1};
    vector<vector<double> > actin_sets;

    double dt = 1, temp = 0, vis = 0;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0;
    double frac_force = 0;
    
    //pos_sets.push_back({0,0,0});
    vector<double> pos1={0, 24.875, actin_rad,0}, pos2={1,24.875,actin_rad,0}, pos3={2,24.875,actin_rad,0}, pos4={3,24.875,actin_rad,0};
    vector<double> pos5={1, -24, actin_rad,1},    pos6={2,-24,actin_rad,1},    pos7={3,-24,actin_rad,1},    pos8={4,-24,actin_rad,1};
    
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    actin_sets.push_back(pos3);
    actin_sets.push_back(pos4);
    actin_sets.push_back(pos5);
    actin_sets.push_back(pos6);
    actin_sets.push_back(pos7);
    actin_sets.push_back(pos8);
    
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, bending, frac_force, bc);
    f->quad_update(); 
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    //MOTOR
    double mx = 2.35, my = -24.625, mang = pi/2, mlen = 1;
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    //cout<<m.to_string();
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {f_index, l_index} = {0, 2}, because the distance between head 0 and that link is 0
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    //cout<<m.to_string();
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    m.attach(1); //should attach to {f_index, l_index} = {1, 1}, because kon * dt > 1
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 1);
    
    m.update_force();
    //m.step_twoheads();
    m.actin_update();
    
    array<double, 2> forcea, forceb, /*bottom 2*/ forcec, forced /*top 2*/;
    forcea = f->get_filament(0)->get_actin(2)->get_force();
    forceb = f->get_filament(0)->get_actin(3)->get_force();
    forcec = f->get_filament(1)->get_actin(1)->get_force();
    forced = f->get_filament(1)->get_actin(2)->get_force();
    
    double tol = 0.001;
    BOOST_CHECK_CLOSE( forcea[0]+1,1, tol); 
    BOOST_CHECK_CLOSE( forceb[0]+1,1, tol); 
    BOOST_CHECK_CLOSE( forcec[0]+1,1, tol); 
    BOOST_CHECK_CLOSE( forced[0]+1,1, tol); 
    
    BOOST_CHECK_CLOSE( forcea[1],mstiff*0.65*0.125, tol); 
    BOOST_CHECK_CLOSE( forceb[1],mstiff*0.35*0.125, tol); 
    BOOST_CHECK_CLOSE( forcec[1],mstiff*-0.65*0.125, tol); 
    BOOST_CHECK_CLOSE( forced[1],mstiff*-0.35*0.125, tol); 
    
    m.step_onehead(0);
    m.step_onehead(1);
    
    //koff = kend = 0 so shouldn't detach
    //filaments and motor are perpendicular so velocity = v0
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.9, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 2.1, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 24.875, tol);
    
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 0.9, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[1], 2.1, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[1], -24, tol);

    
}   

BOOST_AUTO_TEST_CASE( step_twoheads )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {50,50};
    vector<vector<double> > actin_sets;
    array<int, 2> nq = {100,100}, state = {1,1}, findex = {0,1}, lindex = {2, 2};

    double dt = 1, temp = 0, vis = 0;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={3, 0, actin_rad,0}, pos2={2, 0,actin_rad,0}, pos3={1, 0,actin_rad,0}, pos4={0, 0,actin_rad,0};
    vector<double> pos5={0, 3, actin_rad,1}, pos6={0, 2,actin_rad,1}, pos7={0, 1,actin_rad,1}, pos8={0, 0,actin_rad,1};
    
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    actin_sets.push_back(pos3);
    actin_sets.push_back(pos4);
    actin_sets.push_back(pos5);
    actin_sets.push_back(pos6);
    actin_sets.push_back(pos7);
    actin_sets.push_back(pos8);
    
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, bending, frac_force, bc);
    f->quad_update(); 
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    //MOTOR
    double mx = 0.25, my = sqrt(3)/4, mang = (pi/2.+pi/6.), mlen = 1;
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 2*mx, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 2*mx, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 2*my, tol);
    BOOST_CHECK_SMALL(m.get_hx()[1], zero);
    BOOST_CHECK_CLOSE(m.get_hy()[1], 2*my, tol);
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 2);
    
    /*##################################*/
    m.update_force();
    m.step_onehead(0);
    m.step_onehead(1);
    /*##################################*/

    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 2*mx+v0*dt, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 2*mx+v0*dt, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 2*my+v0*dt - 1, tol);
    BOOST_CHECK_SMALL(m.get_hx()[1], zero);
    BOOST_CHECK_CLOSE(m.get_hy()[1], 2*my+v0*dt, tol);
    
    m.actin_update();
//    array<double, 2> forcea, forceb, /*bottom 2*/ forcec, forced /*top 2*/;
    
    m.step_onehead(0);
    m.step_onehead(1);
    
}   

BOOST_AUTO_TEST_CASE( force_attached )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {50,50};
    vector<vector<double> > actin_sets;
    array<int, 2> nq = {100,100}, state = {1,1}, findex = {0,1}, lindex = {2, 2};

    double dt = 1, temp = 0, vis = 0;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={3, 0, actin_rad,0}, pos2={2, 0,actin_rad,0}, pos3={1, 0,actin_rad,0}, pos4={0, 0,actin_rad,0};
    vector<double> pos5={0, 3, actin_rad,1}, pos6={0, 2,actin_rad,1}, pos7={0, 1,actin_rad,1}, pos8={0, 0,actin_rad,1};
    
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    actin_sets.push_back(pos3);
    actin_sets.push_back(pos4);
    actin_sets.push_back(pos5);
    actin_sets.push_back(pos6);
    actin_sets.push_back(pos7);
    actin_sets.push_back(pos8);
    
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, bending, frac_force, bc);
    f->quad_update(); 
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    //MOTOR
    double mx = 0.25, my = sqrt(3)/4, mang = (pi/2.+pi/6.), mlen = 1;
    motor<ATfilament_ensemble> m = motor<ATfilament_ensemble>({mx, my, mang}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 2*mx, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 2*mx, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 2*my, tol);
    BOOST_CHECK_SMALL(m.get_hx()[1], zero);
    BOOST_CHECK_CLOSE(m.get_hy()[1], 2*my, tol);
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 2);
    
    /*##################################*/
    m.update_force();
    m.step_onehead(0);
    m.step_onehead(1);
    /*##################################*/

    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 2*mx+v0*dt, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 2*mx+v0*dt, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 2*my+v0*dt - 1, tol);
    BOOST_CHECK_SMALL(m.get_hx()[1], zero);
    BOOST_CHECK_CLOSE(m.get_hy()[1], 2*my+v0*dt, tol);
    
    m.actin_update();
//    array<double, 2> forcea, forceb, /*bottom 2*/ forcec, forced /*top 2*/;
    
    m.step_onehead(0);
    m.step_onehead(1);
    
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

        void move_end_detach(int hd, double pos);

*/
// EOF
