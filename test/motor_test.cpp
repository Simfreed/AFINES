#include "motor.h"
#include "filament_ensemble.h"
#include "globals.h"
#define BOOST_TEST_MODULE motor_test
#include <boost/test/unit_test.hpp>

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
    
    
    ATfilament_ensemble * f = new ATfilament_ensemble(actin_density, fov, nq, dt, temp, actin_rad, vis, nactin, link_len, pos_sets, stretching, bending, frac_force, bc, seed);
    motor<ATfilament_ensemble>  m = motor<ATfilament_ensemble>({1, 1, 3.1416/2}, motor_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, kon, koff, kend, actin_rad, vis, bc);

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
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
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

        inline void reflect(double t, double gamma, double x1, double x2, double y1, double y2);
        
        inline void periodic(double t, double gamma, double x1, double x2, double y1, double y2);

*/
// EOF
