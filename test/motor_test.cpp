#include "motor.h"
#include "filament_ensemble.h"
#define BOOST_TEST_MODULE motor_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{

    double tol = 0.001;
    double bead_density = 1;
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{0,0}}, findex = {{0,0}}, lindex = {{-1, -1}};
    vector<array<double, 3> > pos_sets;
    int nbead = 3;

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double bead_rad = 0.5, spring_len = 1, motor_len = 1;
    double v0 = 0;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    
    filament_ensemble * f = new filament_ensemble(bead_density, fov, nq, dt, temp, 
            bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);
    cout<<"\nDEBUG: nfilaments = "<<f->get_nfilaments();
    motor  m = motor(array<double, 3>{{1, 1, 3.1416/2}}, motor_len, f, state, 
            findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    cout<<"\nDEBUG: created motors";
    BOOST_CHECK_CLOSE( m.get_hx()[0],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hx()[1],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[0], 0.5, tol);                // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[1], 1.5, tol);               
    BOOST_CHECK_EQUAL( m.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m.get_states()[1], 0);
    
    //delete f;

    //bead_density = nbead/(fov[0]*fov[1]);
    //pos_sets.push_back({{-0.5,0,0}});

    //f = new filament_ensemble(bead_density, fov, nq, dt, temp, bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);
    state = {{1,0}};
    lindex = {{0,-1}};
    f->get_filament(0)->get_spring(0);
    cout<<"\nDEBUG: got to line 50";

    motor m2 = motor(array<double, 3>{{1, 1, 0}}, motor_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    BOOST_CHECK_EQUAL( m2.get_states()[0], 1);
    BOOST_CHECK_EQUAL( m2.get_states()[1], 0);
    
    state = {{0,1}};
    lindex = {{-1,0}};
    motor m3 = motor(array<double, 3>{{1, 1, 0}}, motor_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    BOOST_CHECK_EQUAL( m3.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m3.get_states()[1], 1);
   
    delete f;
} 


BOOST_AUTO_TEST_CASE( step_onehead )
{

    //Filament ENSEMBLE
    double tol = 0.001;
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{1,0}}, findex = {{0,-1}}, lindex = {{0, -1}};
    vector<array<double, 3> > pos_sets;
    int nbead = 2;
    double bead_density = nbead/(fov[0]*fov[1]);

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 2, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    pos_sets.push_back({{-0.4,0,0}});
    filament_ensemble * f = new filament_ensemble(bead_density, fov, nq, dt, temp, bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);

    //MOTOR
    double mx = 0.5, my = 0.5, mang = pi/2, mlen = 1;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);

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
/*
    // IF no detach
    // (a) new position is correct
    koff = -1;
    kend = -1;
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
*/
    delete f;
}

BOOST_AUTO_TEST_CASE( update_pos_a_end )
{
    //Filament ENSEMBLE
    double tol = 0.001, zero = 1e-10;   
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{1,0}}, findex = {{0,-1}}, lindex = {{0, -1}};
    vector<array<double, 3> > pos_sets;
    int nbead = 4;
    double bead_density = nbead/(fov[0]*fov[1]);

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 2, kend = 2, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    pos_sets.push_back({{0,0,0}});
    filament_ensemble * f = new filament_ensemble(bead_density, fov, nq, dt, temp, bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);

    //MOTOR
    double mx = 0.125, my = 0.5, mang = pi/2, mlen = 1;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);

    // IF position is greater than filament length
    double pos = 1.125;
    //  If last rod in filament and pos is greater than length of spring, go to barbed end

    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.update_pos_a_end(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[0], spring_len);
    
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], my - 0.5*mlen*sin(mang), tol);
    
    m.update_position_attached(0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    //  Else (not last rod in filament)
    kend = -1;
    mx = 1.125;
    lindex = {{1, -1}};
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
    //      (a) motor moves to next rod in filament
    //      (b) motor moves appropriate distance
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.update_pos_a_end(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.125, tol);
    
    m.update_position_attached(0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.875, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    // Else ( new position isn't more than filament length
    //  (a) motor moves to appropriate new position
    pos = 0.375;
    m.update_pos_a_end(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], pos, tol);
    
    m.update_position_attached(0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.625, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
   
    //Test 2 : Motor on filament facing left, so, moves right
    pos_sets.clear();
    pos_sets.push_back({{0,0,pi}});
    f = new filament_ensemble(bead_density, fov, nq, dt, temp, bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);

    //MOTOR
    mx = -0.125;
    findex = {{0,-1}}, lindex = {{0, -1}}; 
    kon = 0, koff = 2, kend = 2; 

    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);

    // IF position is greater than filament length
    pos = 1.125;
    
    //  If last rod in filament and pos is greater than length, move to barbed end
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.update_pos_a_end(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[0], 1);
    
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], my - 0.5*mlen*sin(mang), tol);
    
    m.update_position_attached(0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    //  Else (not last rod in filament)
    mx = -1.125;
    lindex = {{1, -1}};
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
    //      (a) motor moves to next rod in filament
    //      (b) motor moves appropriate distance
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.update_pos_a_end(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.125, tol);
    
    m.update_position_attached(0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], -0.875, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    // Else ( new position isn't more than filament length
    //  (a) motor moves to appropriate new position
    pos = 0.375;
    m.update_pos_a_end(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], pos, tol);
    
    m.update_position_attached(0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], -0.625, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    //TEST 3, OTHER HEAD
    mx = 0.125, my =-0.5;
    state = {{0,1}},  findex = {{-1,0}}, lindex = {{-1, 0}}; 
    kon = 0, koff = 2, kend = 2; 
    pos_sets.clear();
    pos_sets.push_back({{0,0,0}});
    f = new filament_ensemble(bead_density, fov, nq, dt, temp, bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);

    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);

    // IF position is greater than filament length
    pos = 1.125;
    //  If last rod in filament, move to barbed end

    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 0.875, tol);
    m.update_pos_a_end(1, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 0);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[1], spring_len);
    
    BOOST_CHECK_CLOSE(m.get_hx()[1], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[1], my + 0.5*mlen*sin(mang), tol);
    
    m.update_position_attached(1);
    BOOST_CHECK_CLOSE(m.get_hx()[1], 0, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[1], 0, tol);
    
    //  Else (not last rod in filament)
    kend = -1;
    mx = 1.125;
    lindex = {{-1, 1}};
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
    //      (a) motor moves to next rod in filament
    //      (b) motor moves appropriate distance
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 0.875, tol);
    m.update_pos_a_end(1, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], 0.125, tol);
    
    m.update_position_attached(1);
    BOOST_CHECK_CLOSE(m.get_hx()[1], 0.875, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[1], 0, tol);
    
    // Else ( new position isn't more than filament length
    //  (a) motor moves to appropriate new position
    pos = 0.375;
    m.update_pos_a_end(1, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[1], pos, tol);
    
    m.update_position_attached(1);
    BOOST_CHECK_CLOSE(m.get_hx()[1], 0.625, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[1], 0, tol);
    delete f;
}

BOOST_AUTO_TEST_CASE( attach )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<array<double, 3> > pos_sets;
    int nbead = 5;
    int nfil = 1;//3;
    double bead_density = nfil*nbead/(fov[0]*fov[1]);

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    double seed = -1;
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    //pos_sets.push_back({{1,1,pi/2}});
    pos_sets.push_back({{0,0,0}});
    //pos_sets.push_back({{-2,3,pi}});
    filament_ensemble * f = new filament_ensemble(bead_density, fov, nq, dt, temp, bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);
    f->quad_update_serial(); 
    //MOTOR
    double mx = 2.35, my = 0.5, mang = pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
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
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {{f_index, l_index}} = {{0, 3}}, because the distance between head 0 and that spring is 0
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
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{1,0}}, findex = {{0,-1}}, lindex = {{1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 2, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    vector<double> pos1={{23.6,0,bead_rad,0}},pos2={{24.6,0,bead_rad,0}},pos3={{-24.4,0,bead_rad,0}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
//    cout<<f->get_filament(0)->to_string();
    //MOTOR
    double mx = -24.875, my = 0.5, mang = pi/2, mlen = 1;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{100,100}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    //pos_sets.push_back({{0,0,0}});
    vector<double> pos1={{0,24.875,bead_rad,0}},pos2={{1,24.875,bead_rad,0}},pos3={{2,24.875,bead_rad,0}}, pos4={{3,24.875,bead_rad,0}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    //MOTOR
    double mx = 2.35, my = -24.625, mang = pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
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
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {{f_index, l_index}} = {{0, 3}}, because the distance between head 0 and that spring is 0
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
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.189;
    double frac_force = 0;
    
    //pos_sets.push_back({{0,0,0}});
    vector<double> pos1={{0, 24.875, bead_rad,0}}, pos2={{1,24.875,bead_rad,0}}, pos3={{2,24.875,bead_rad,0}}, pos4={{3,24.875,bead_rad,0}};
    vector<double> pos5={{1, -24, bead_rad,1}},    pos6={{2,-24,bead_rad,1}},    pos7={{3,-24,bead_rad,1}},    pos8={{4,-24,bead_rad,1}};
    
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    bead_sets.push_back(pos5);
    bead_sets.push_back(pos6);
    bead_sets.push_back(pos7);
    bead_sets.push_back(pos8);
    
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    //MOTOR
    double mx = 2.35, my = -24.625, mang = pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    //cout<<m.to_string();
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    m.attach(0); //should attach to {{f_index, l_index}} = {{0, 2}}, because the distance between head 0 and that spring is 0
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_states()[1], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], -1);
    //cout<<m.to_string();
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    m.attach(1); //should attach to {{f_index, l_index}} = {{1, 1}}, because kon * dt > 1
    BOOST_CHECK_EQUAL(m.get_states()[1], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[1], 1);
    BOOST_CHECK_EQUAL(m.get_l_index()[1], 1);
    
    m.update_angle();
    m.update_force();
    //m.step_twoheads();
    m.filament_update();
    
    array<double, 2> forcea, forceb, /*bottom 2*/ forcec, forced /*top 2*/;
    forcea = f->get_filament(0)->get_bead(2)->get_force();
    forceb = f->get_filament(0)->get_bead(3)->get_force();
    forcec = f->get_filament(1)->get_bead(1)->get_force();
    forced = f->get_filament(1)->get_bead(2)->get_force();
    
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
    array<double, 2> fov = {{50,50}};
    vector<vector<double> > bead_sets;
    array<int, 2> nq = {{100,100}}, state = {{1,1}}, findex = {{0,1}}, lindex = {{2, 2}};

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={{3, 0, bead_rad,0}}, pos2={{2, 0,bead_rad,0}}, pos3={{1, 0,bead_rad,0}}, pos4={{0, 0,bead_rad,0}};
    vector<double> pos5={{0, 3, bead_rad,1}}, pos6={{0, 2,bead_rad,1}}, pos7={{0, 1,bead_rad,1}}, pos8={{0, 0,bead_rad,1}};
    
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    bead_sets.push_back(pos5);
    bead_sets.push_back(pos6);
    bead_sets.push_back(pos7);
    bead_sets.push_back(pos8);
    
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    //MOTOR
    double mx = 0.25, my = sqrt(3)/4, mang = (pi/2.+pi/6.), mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
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
    m.update_angle();
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
    
    m.filament_update();
//    array<double, 2> forcea, forceb, /*bottom 2*/ forcec, forced /*top 2*/;
    
    m.step_onehead(0);
    m.step_onehead(1);
    
}   

BOOST_AUTO_TEST_CASE( force_attached )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    vector<vector<double> > bead_sets;
    array<int, 2> nq = {{100,100}}, state = {{1,1}}, findex = {{0,1}}, lindex = {{2, 2}};

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={{3, 0, bead_rad,0}}, pos2={{2, 0,bead_rad,0}}, pos3={{1, 0,bead_rad,0}}, pos4={{0, 0,bead_rad,0}};
    vector<double> pos5={{0, 3, bead_rad,1}}, pos6={{0, 2,bead_rad,1}}, pos7={{0, 1,bead_rad,1}}, pos8={{0, 0,bead_rad,1}};
    
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    bead_sets.push_back(pos5);
    bead_sets.push_back(pos6);
    bead_sets.push_back(pos7);
    bead_sets.push_back(pos8);
    
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    //cout<<f->get_filament(0)->to_string(); 
    //cout<<f->get_filament(1)->to_string(); 
    //MOTOR
    double mx = 0.25, my = sqrt(3)/4, mang = (pi/2.+pi/6.), mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
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
    m.update_angle();
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
    
    m.filament_update();
//    array<double, 2> forcea, forceb, /*bottom 2*/ forcec, forced /*top 2*/;
    
    m.step_onehead(0);
    m.step_onehead(1);
    
}   

BOOST_AUTO_TEST_CASE( dead_head )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    vector<vector<double> > bead_sets;
    array<int, 2> nq = {{100,100}}, state = {{1,-1}}, findex = {{0,-1}}, lindex = {{2, -1}};

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={{3, 0, bead_rad,0}}, pos2={{2, 0,bead_rad,0}}, pos3={{1, 0,bead_rad,0}}, pos4={{0, 0,bead_rad,0}};
    
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    //MOTOR
    double mx = 0.45, my = 0.5, mang = pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    //m.kill_head(1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], mx, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/

    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
   
    double h0x1 = mx+v0*dt;
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], h0x1, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], h0x1, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx = -mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*cos(atan2(mlen, v0*dt));
    double mfy =  mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*sin(atan2(mlen, v0*dt));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy, tol);

//    f->update_positions();
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], mfx*h0x1, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[0], mfx*(1-h0x1), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], mfy*h0x1, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[1], mfy*(1-h0x1), tol);
    
// PURPOSELY NOT UPDATING POSITIONS OF ACTINS, SO THAT I DON'T HAVE TO FIGURE OUT 
// TOO MANY THINGS FOR THE SLOWED DOWN MYOSIN
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/
    
    double fmax = 3.85;
    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    
    double h0x2 = mx + v0*dt + v0*(1+mfx/fmax)*dt;
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], h0x2 , tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], h0x2, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx2 = -mstiff*(sqrt((h0x2-mx)*(h0x2-mx)+ mlen*mlen)-mlen)*cos(atan2(mlen, h0x2-mx));
    double mfy2 =  mstiff*(sqrt((h0x2-mx)*(h0x2-mx)+ mlen*mlen)-mlen)*sin(atan2(mlen, h0x2-mx));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx2, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy2, tol);

    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], mfx*h0x1 + mfx2*h0x2, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[0], mfx*(1-h0x1) + mfx2*(1-h0x2), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], mfy*h0x1 + mfy2*h0x2, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[1], mfy*(1-h0x1) + mfy2*(1-h0x2), tol);

    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/
    
    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    
    double h0x3 = mx + v0*dt + v0*(1+mfx/fmax)*dt + v0*(1+mfx2/fmax)*dt;
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], h0x3 - 1, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], h0x3, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx3 = -mstiff*(sqrt((h0x3-mx)*(h0x3-mx)+ mlen*mlen)-mlen)*cos(atan2(mlen, h0x3-mx));
    double mfy3 =  mstiff*(sqrt((h0x3-mx)*(h0x3-mx)+ mlen*mlen)-mlen)*sin(atan2(mlen, h0x3-mx));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx3, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy3, tol);

    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], mfx*h0x1 + mfx2*h0x2 + mfx3*(2-h0x3), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[0], mfx3*(h0x3-1), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], mfy*h0x1 + mfy2*h0x2 + mfy3*(2-h0x3), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[1], mfy3*(h0x3-1), tol);

}   

BOOST_AUTO_TEST_CASE( dead_head_upside_down )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    vector<vector<double> > bead_sets;
    array<int, 2> nq = {{100,100}}, state = {{1,-1}}, findex = {{0,-1}}, lindex = {{2, -1}};

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={{3, 0, bead_rad,0}}, pos2={{2, 0,bead_rad,0}}, pos3={{1, 0,bead_rad,0}}, pos4={{0, 0,bead_rad,0}};
    
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    //MOTOR
    double mx = 0.25, my = -0.5, mang = -pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    //m.kill_head(1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], mx, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/

    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], mx+v0*dt, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx+v0*dt, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx = -mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*cos(atan2(mlen, v0*dt));
    double mfy = -mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*sin(atan2(mlen, v0*dt));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy, tol);

//    f->update_positions();
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], mfx/2, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[0], mfx/2, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], mfy/2, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[1], mfy/2, tol);

// PURPOSELY NOT UPDATING POSITIONS OF ACTINS, SO THAT I DON'T HAVE TO FIGURE OUT 
// TOO MANY THINGS FOR THE SLOWED DOWN MYOSIN
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/
    
    double fmax = 3.85;
    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    
    double h0x = mx + v0*dt + v0*(1+mfx/fmax)*dt;
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], h0x , tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], h0x, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx2 = -mstiff*(sqrt((h0x-mx)*(h0x-mx)+ mlen*mlen)-mlen)*cos(atan2(mlen, h0x-mx));
    double mfy2 = -mstiff*(sqrt((h0x-mx)*(h0x-mx)+ mlen*mlen)-mlen)*sin(atan2(mlen, h0x-mx));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx2, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy2, tol);

    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], mfx/2 + mfx2*h0x, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[0], mfx/2 + mfx2*(1-h0x), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], mfy/2 + mfy2*h0x, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 3)[1], mfy/2 + mfy2*(1-h0x), tol);

}   

BOOST_AUTO_TEST_CASE( dead_head_bwd )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    vector<vector<double> > bead_sets;
    array<int, 2> nq = {{100,100}}, state = {{1,-1}}, findex = {{0,-1}}, lindex = {{1, -1}};

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={{0, 0, bead_rad,0}}, pos2={{1, 0,bead_rad,0}}, pos3={{2, 0,bead_rad,0}}, pos4={{3, 0,bead_rad,0}};
    
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    //MOTOR
    double mx = 1.5, my = 0.5, mang = pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    //m.kill_head(1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.5, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/

    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.75, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx-v0*dt, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx =  mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*cos(atan2(mlen, v0*dt));
    double mfy =  mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*sin(atan2(mlen, v0*dt));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy, tol);

//    f->update_positions();
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[0], 0.75*mfx, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], 0.25*mfx, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[1], 0.75*mfy, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], 0.25*mfy, tol);

// PURPOSELY NOT UPDATING POSITIONS OF ACTINS, SO THAT I DON'T HAVE TO FIGURE OUT 
// TOO MANY THINGS FOR THE SLOWED DOWN MYOSIN
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/
    
    double fmax = 3.85;
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    
    double h0x2 = mx - v0*dt - v0*(1-mfx/fmax)*dt;
    double h0pos_a_end = 2-h0x2; //1 - (h0x-1)
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], h0pos_a_end , tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], h0x2, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx2 = mstiff*(sqrt((h0x2-mx)*(h0x2-mx)+ mlen*mlen)-mlen)*cos(atan2(mlen, mx - h0x2));
    double mfy2 = mstiff*(sqrt((h0x2-mx)*(h0x2-mx)+ mlen*mlen)-mlen)*sin(atan2(mlen, mx - h0x2));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx2, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy2, tol);

    BOOST_CHECK_CLOSE(f->get_force(0, 1)[0], 0.75*mfx + mfx2*h0pos_a_end, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], 0.25*mfx + mfx2*(1-h0pos_a_end), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[1], 0.75*mfy + mfy2*h0pos_a_end, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], 0.25*mfy + mfy2*(1-h0pos_a_end), tol);

    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 0);
    
    double h0x3 = mx - v0*dt - v0*(1-mfx/fmax)*dt - v0*(1-mfx2/fmax)*dt;
    //double h0pos_a_end = 2-h0x; //1 - (h0x-1)
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 1-h0x3 , tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], h0x3, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx3 = mstiff*(sqrt((h0x3-mx)*(h0x3-mx)+ mlen*mlen)-mlen)*cos(atan2(mlen, mx - h0x3));
    double mfy3 = mstiff*(sqrt((h0x3-mx)*(h0x3-mx)+ mlen*mlen)-mlen)*sin(atan2(mlen, mx - h0x3));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx3, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy3, tol);

    BOOST_CHECK_CLOSE(f->get_force(0, 1)[0], 0.75*mfx + mfx2*h0pos_a_end + mfx3*h0x3, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 0)[0], (1-h0x3)*mfx3, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[1], 0.75*mfy + mfy2*h0pos_a_end + mfy3*h0x3, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 0)[1], (1-h0x3)*mfy3, tol);
}   

BOOST_AUTO_TEST_CASE( dead_head_bwd_upside_down )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    vector<vector<double> > bead_sets;
    array<int, 2> nq = {{100,100}}, state = {{1,-1}}, findex = {{0,-1}}, lindex = {{1, -1}};

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={{0, 0, bead_rad,0}}, pos2={{1, 0,bead_rad,0}}, pos3={{2, 0,bead_rad,0}}, pos4={{3, 0,bead_rad,0}};
    
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    //MOTOR
    double mx = 1.5, my = -0.5, mang = 3*pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    //m.kill_head(1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.5, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/

    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.75, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx-v0*dt, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx =  mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*cos(atan2(mlen, v0*dt));
    double mfy =  -mstiff*(sqrt(v0*v0*dt*dt+ mlen*mlen)-mlen)*sin(atan2(mlen, v0*dt));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy, tol);

//    f->update_positions();
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[0], 0.75*mfx, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], 0.25*mfx, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[1], 0.75*mfy, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], 0.25*mfy, tol);

// PURPOSELY NOT UPDATING POSITIONS OF ACTINS, SO THAT I DON'T HAVE TO FIGURE OUT 
// TOO MANY THINGS FOR THE SLOWED DOWN MYOSIN
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/
    
    double fmax = 3.85;
    //force is near 0, so each head just move v0*dt = 0.25
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 1);
    
    double h0x = mx - v0*dt - v0*(1-mfx/fmax)*dt;
    double h0pos_a_end = 2-h0x; //1 - (h0x-1)
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], h0pos_a_end , tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], h0x, tol);
    BOOST_CHECK_SMALL(m.get_hy()[0], zero);

    double mfx2 = mstiff*(sqrt((h0x-mx)*(h0x-mx)+ mlen*mlen)-mlen)*cos(atan2(mlen, mx - h0x));
    double mfy2 = -mstiff*(sqrt((h0x-mx)*(h0x-mx)+ mlen*mlen)-mlen)*sin(atan2(mlen, mx - h0x));

    BOOST_CHECK_CLOSE(m.get_force()[0], mfx2, tol);
    BOOST_CHECK_CLOSE(m.get_force()[1], mfy2, tol);

    BOOST_CHECK_CLOSE(f->get_force(0, 1)[0], 0.75*mfx + mfx2*h0pos_a_end, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[0], 0.25*mfx + mfx2*(1-h0pos_a_end), tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 1)[1], 0.75*mfy + mfy2*h0pos_a_end, tol);
    BOOST_CHECK_CLOSE(f->get_force(0, 2)[1], 0.25*mfy + mfy2*(1-h0pos_a_end), tol);

}   

//Check for attaching at different parts of a filament

BOOST_AUTO_TEST_CASE( attach_difft_spots )
{
    //Filament ENSEMBLE
    
    array<double, 2> fov = {{33,2}};
    array<int, 2> nq = {{33,2}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0;
    
    double mstiff = 1, stretching = 0, bending = 0; //spring constants
    double kon = 0.5, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    //pos_sets.push_back({{0,0,0}});
//    vector<double> pos1={{0,0,bead_rad,0}},pos2={{1,0,bead_rad,0}},pos3={{2,0,bead_rad,0}}, pos4={{3,0,bead_rad,0}};
    int nbead = 16;
    vector<double> pos;
    for (int i =0; i<nbead; i++){
        pos = {{double(i), 0, bead_rad, 0}};
        bead_sets.push_back(pos);
    }
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    //attachment
    double mx = 0, my = 0.075, mang = pi/2, mlen = 0.15, posx, posy;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 15> num_attached; 
    int nevents = 1000;
    int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i < nbead-1; i++){
        mx = i + 0.35;
        num_attached[i]=0;
        posx=0;
        posy=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_attached[i] += m.get_states()[0];
            if (m.get_states()[0] == 1){
                posx += m.get_hx()[0];
                posy += m.get_hy()[0];
            }
        }
        cout<<"\nmx = "<<mx<<";\tnum_attached = "<<num_attached[i]<<";\t(<x>, <y>) = ( "<<posx/num_attached[i]<<" , "<<posy/num_attached[i]<<" )";
        BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    }
 
    my=0.075; mlen=0.15;
    kon=0.05;
    nevents = 10000;
    threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i < nbead-1; i++){
        mx = i + 0.5;
        num_attached[i]=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_attached[i] += m.get_states()[0];
        }
        BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    }
    
    kon=0.9;
    nevents = 10000;
    threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i < nbead-1; i++){
        mx = i + 0.5;
        num_attached[i]=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_attached[i] += m.get_states()[0];
        }
        BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    }
    
    //detachment
    array<int,15> num_detached;
    koff = 0.1;
    kend = 0.1;
    nevents = 10000;
    threesig = int(3*sqrt(koff*(1-koff)*double(nevents)));
    state = {{1,0}};
    findex = {{0,-1}};
    for (int i = 0; i < nbead-1; i++){
        mx = i + 0.5;
        lindex = {{i, -1}};
        num_detached[i]=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.step_onehead(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_detached[i] += pr(m.get_states()[0]);
        }
        BOOST_CHECK_SMALL(num_detached[i] - int(nevents*koff), threesig); 
    }
    
    koff = 0.5;
    kend = 0.5;
    nevents = 20000;
    threesig = int(3*sqrt(koff*(1-koff)*double(nevents)));
    state = {{1,0}};
    findex = {{0,-1}};
    for (int i = 0; i < nbead-1; i++){
        mx = i + 0.5;
        lindex = {{i, -1}};
        num_detached[i]=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.step_onehead(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_detached[i] += pr(m.get_states()[0]);
        }
        cout<<"\nmx = "<<mx<<";\tnum_attached = "<<num_detached[i];//<<";\t(<x>, <y>) = ( "<<posx/num_detached[i]<<" , "<<posy/num_detached[i]<<" )";
        BOOST_CHECK_SMALL(num_detached[i] - int(nevents*koff), threesig); 
    }
    
    koff = 0.85;
    kend = 0.85;
    nevents = 10000;
    threesig = int(3*sqrt(koff*(1-koff)*double(nevents)));
    state = {{1,0}};
    findex = {{0,-1}};
    for (int i = 0; i < nbead-1; i++){
        mx = i + 0.5;
        lindex = {{i, -1}};
        num_detached[i]=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.step_onehead(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_detached[i] += pr(m.get_states()[0]);
        }
        BOOST_CHECK_SMALL(num_detached[i] - int(nevents*koff), threesig); 
    }
    
}   

BOOST_AUTO_TEST_CASE( attach_difft_spots_rig )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{25,2}};
    array<int, 2> nq = {{25,2}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 5, spring_len = 10;
    double v0 = 0;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 0.5, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    set_seed(12); 
    //pos_sets.push_back({{0,0,0}});
    vector<double> pos1={{-5,0,bead_rad,0}},pos2={{5,0,bead_rad,0}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    cout<<f->get_filament(0)->to_string();
    
    //attachment
    double mx = 0, my = 0.075, mang = pi/2, mlen = 0.15, posx, posy;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 11> num_attached; 
    int nevents = 10000;
    int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i <= int(spring_len); i++){
        mx = i - 5;
        num_attached[i]=0;
        posx=0;
        posy=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.attach(0);
            num_attached[i] += m.get_states()[0];
            if (m.get_states()[0] == 1){
                posx += m.get_hx()[0];
                posy += m.get_hy()[0];
            }
        }
        cout<<"\nmx = "<<mx<<";\tnum_attached = "<<num_attached[i]<<";\t(<x>, <y>) = ( "<<posx/num_attached[i]<<" , "<<posy/num_attached[i]<<" )";
        BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    }
    
    kon=0.05;
    nevents = 10000;
    threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i <= int(spring_len); i++){
        mx = i - 5;
        num_attached[i]=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_attached[i] += m.get_states()[0];
        }
        BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    }
    
    kon=0.9;
    nevents = 10000;
    threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i <= int(spring_len); i++){
        mx = i - 5;
        num_attached[i]=0;
        for (int j = 0; j < nevents; j++){
            m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_attached[i] += m.get_states()[0];
        }
        BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    }

}

BOOST_AUTO_TEST_CASE( attach_two_options )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{25,2}};
    array<int, 2> nq = {{1,1}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 5, spring_len = 10;
    double v0 = 0;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 0.5, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    set_seed(12); 
    //pos_sets.push_back({{0,0,0}});
    vector<double> pos1={{-5,0,bead_rad,0}},pos2={{5,0,bead_rad,0}};
    vector<double> pos3={{5,0,bead_rad,1}},pos4={{-5,0,bead_rad,1}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
//    cout<<f->get_filament(0)->to_string();
    
    //attachment
    double mx = -4.5, my = 0, mang = 0, mlen = 0.15;//, posx, posy
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 2> num_attached; 
    int nevents = 10000;
    //int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    num_attached[0]=0;
    num_attached[1]=0;
    for (int j = 0; j < nevents; j++){
        //cout<<"\niter "<<j<<endl;
        m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
        m.attach(0);
        if (m.get_states()[0] == 1){
            num_attached[m.get_f_index()[0]]+=1;
            //posx += m.get_hx()[0];
            //posy += m.get_hy()[0];
        }
    }
    cout<<"\nnum_attached[0] = "<<num_attached[0]<<";\tnum_attached[1] = "<<num_attached[1]<<endl;
        //BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    delete f;
    
    //Reverse Filament Order
    bead_sets.clear();
    pos1={{5,0,bead_rad,0}},pos2={{-5,0,bead_rad,0}};
    pos3={{-5,0,bead_rad,1}},pos4={{5,0,bead_rad,1}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    
    f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    num_attached[0]=0;
    num_attached[1]=0;
    for (int j = 0; j < nevents; j++){
        //cout<<"\niter "<<j<<endl;
        m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
        m.attach(0);
        if (m.get_states()[0] == 1){
            num_attached[m.get_f_index()[0]]+=1;
            //posx += m.get_hx()[0];
            //posy += m.get_hy()[0];
        }
    }
    cout<<"\nnum_attached[0] = "<<num_attached[0]<<";\tnum_attached[1] = "<<num_attached[1]<<endl;
   delete f; 
}

BOOST_AUTO_TEST_CASE( brownian_attach)
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{25,2}};
    array<int, 2> nq = {{1,1}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 5, spring_len = 10;
    double v0 = 0;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 0.5, koff = 0, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    set_seed(12); 
    //pos_sets.push_back({{0,0,0}});
    vector<double> pos1={{-5,0,bead_rad,0}},pos2={{5,0,bead_rad,0}};
    vector<double> pos3={{5,0,bead_rad,1}},pos4={{-5,0,bead_rad,1}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    
    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
//    cout<<f->get_filament(0)->to_string();
    
    //attachment
    double mx = -4.5, my = 0, mang = 0, mlen = 0.15;//, posx, posy
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 2> num_attached; 
    int nevents = 10000;
    //int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    num_attached[0]=0;
    num_attached[1]=0;
    for (int j = 0; j < nevents; j++){
        //cout<<"\niter "<<j<<endl;
        m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
        m.attach(0);
        if (m.get_states()[0] == 1){
            num_attached[m.get_f_index()[0]]+=1;
            //posx += m.get_hx()[0];
            //posy += m.get_hy()[0];
        }
    }
    cout<<"\nnum_attached[0] = "<<num_attached[0]<<";\tnum_attached[1] = "<<num_attached[1]<<endl;
        //BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    delete f;
    
    //Reverse Filament Order
    bead_sets.clear();
    pos1={{5,0,bead_rad,0}},pos2={{-5,0,bead_rad,0}};
    pos3={{-5,0,bead_rad,1}},pos4={{5,0,bead_rad,1}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    
    f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    num_attached[0]=0;
    num_attached[1]=0;
    for (int j = 0; j < nevents; j++){
        //cout<<"\niter "<<j<<endl;
        m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
        m.attach(0);
        if (m.get_states()[0] == 1){
            num_attached[m.get_f_index()[0]]+=1;
            //posx += m.get_hx()[0];
            //posy += m.get_hy()[0];
        }
    }
    cout<<"\nnum_attached[0] = "<<num_attached[0]<<";\tnum_attached[1] = "<<num_attached[1]<<endl;
   delete f; 
}

BOOST_AUTO_TEST_CASE( motor_slow_down_stall )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    vector<vector<double> > bead_sets;
    array<int, 2> nq = {{100,100}}, state = {{1,-1}}, findex = {{0,-1}}, lindex = {{14, -1}};

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 1, rcut = 0.063;
    double frac_force = 0;
    double zero = 1e-10;   

    vector<double> pos;
    for (int i = 15; i >=0; i--)
    {
        pos = {{double(i),0,bead_rad,0}};
        bead_sets.push_back(pos);
    }

    filament_ensemble * f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    
    //MOTOR
    double mx = 0.5, my = 0.5, mang = pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
    /*##################################*/
    m.step_onehead(0);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/

    double x = m.get_hx()[0];
    double vel = (x - mx)/dt;
    while(vel > zero ){
        //cout<<"\nDEBUG: hx = "<<x<<"\tmotor speed  = "<<vel;
        m.step_onehead(0);
        m.update_angle();
        m.update_force();
        m.filament_update();
        BOOST_CHECK((m.get_hx()[0]-x)/dt < vel);
        vel = (m.get_hx()[0] - x)/dt;
        x = m.get_hx()[0];
    }
    
    //OTHER MOTOR HEAD
    mang = 3*pi/2;
    state = {{-1,1}}, findex = {{-1,0}}, lindex = {{-1, 14}};
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
    /*##################################*/
    m.step_onehead(1);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/

    x = m.get_hx()[1];
    vel = (x - mx)/dt;
    while(fabs(vel) > zero ){
        //cout<<"\nDEBUG: hx = "<<x<<"\tmotor speed  = "<<vel;
        m.step_onehead(1);
        m.update_angle();
        m.update_force();
        m.filament_update();
        BOOST_CHECK((m.get_hx()[1]-x)/dt < vel);
        vel = (m.get_hx()[1] - x)/dt;
        x = m.get_hx()[1];
    }
    
    //Filament facing the other direction
    bead_sets.clear();
    for (int i = 0; i <= 15; i++)
    {
        pos = {{double(i),0,bead_rad,0}};
        bead_sets.push_back(pos);
    }

    f = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    f->quad_update_serial(); 
    mx = 14.75;
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
    /*##################################*/
    m.step_onehead(1);
    m.update_angle();
    m.update_force();
    m.filament_update();
    /*##################################*/

    x = m.get_hx()[1];
    vel = (x - mx)/dt;
    while(fabs(vel) > zero ){
        //cout<<"\nDEBUG: hx = "<<x<<"\tmotor speed  = "<<vel;
        m.step_onehead(1);
        m.update_angle();
        m.update_force();
        m.filament_update();
        BOOST_CHECK(fabs((m.get_hx()[1]-x)/dt) < fabs(vel));
        vel = (m.get_hx()[1] - x)/dt;
        x = m.get_hx()[1];
    }
}   

BOOST_AUTO_TEST_CASE( attach_detach )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<array<double, 3> > pos_sets;
    int nbead = 5;
    int nfil = 1;//3;
    double bead_density = nfil*nbead/(fov[0]*fov[1]);

    double dt = 1, temp = 0.004, vis = 0, tol=1e-6;
    string bc = "PERIODIC";
    double seed = -1;
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 10000, koff = 10000, kend = 0, fstall = 3.85, rcut = 0.189;
    double frac_force = 0;
    
    pos_sets.push_back({{0,0,0}});
    filament_ensemble * f = new filament_ensemble(bead_density, fov, nq, dt, temp, bead_rad, vis, nbead, spring_len, pos_sets, stretching, 1, bending, frac_force, bc, seed, 0, 0);
    f->quad_update_serial(); 
    
    //MOTOR
    double mx = 2.35, my = 0.6, mang = pi/2, mlen = 1;
    motor m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);

    BOOST_CHECK_EQUAL(m.get_hx()[0], 2.35);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0.1, tol);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    
    //Attach
    m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because kon is saturated
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], 0);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], 2);
    BOOST_CHECK_EQUAL(m.get_hx()[0], 2.35);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);

    //Detach
    m.step_onehead(0);
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_f_index()[0], -1);
    BOOST_CHECK_EQUAL(m.get_l_index()[0], -1);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 2.35, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0.1, tol);
    
    //MOTOR: attach, walk to end, detach perpendicularly
    v0 = 0.5;
    kon = 10000, koff = 0, kend = 10000;
    mx = 2.5, my = 0.6, mang = pi/2, mlen = 1;
    m = motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);

    //Attach
    m.attach(0); //should attach to 

    //walk
    for(int i=0; i<6; i++){
        BOOST_CHECK_EQUAL(m.get_states()[0], 1);
        m.step_onehead(0);
        m.update_angle();
        m.update_force();
    }
    
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0.1, tol);

    delete f;
}   

/* Functions to test : 
        void attach(int hd);

        void brownian(double t, double gamma);

        void step_twoheads();

        void filament_update();

        void update_shape();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd);

        void update_pos_a_end(int hd, double pos);

*/
// EOF
