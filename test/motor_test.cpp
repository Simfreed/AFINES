#include "motor.h"
#include "actin_ensemble.h"
#include "link_ensemble.h"
#include "globals.h"
#define BOOST_TEST_MODULE motor_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    actin_ensemble ae = actin_ensemble();
    motor m = motor(1, 1, 3.1415926535897/2, 1, &ae, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "blue");
    double tol = 0.001;

    BOOST_CHECK_CLOSE( m.get_hx()[0],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hx()[1],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[0], 0.5, tol);                // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[1], 1.5, tol);               
    BOOST_CHECK_EQUAL( m.get_color(), "blue");
    BOOST_CHECK_EQUAL( m.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m.get_states()[1], 0);
    
    actin * a = new actin(0, 0, 0, 1, 0, 0, 0, 0, 0);
    ae.add_monomer(a, 0);

    motor m2 = motor(1, 1, 0, 1, &ae, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "blue");
    BOOST_CHECK_EQUAL( m2.get_states()[0], 1);
    BOOST_CHECK_EQUAL( m2.get_states()[1], 0);
    
    motor m3 = motor(1, 1, 0, 1, &ae, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "blue");
    BOOST_CHECK_EQUAL( m3.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m3.get_states()[1], 1);
    
} 


BOOST_AUTO_TEST_CASE( step_onehead )
{

    //ACTIN ENSEMBLE
    double startx = 0.5, starty = 0, phi0 = 0;
    int pol_index = 0, nmon = 3;
    double fovx = 20, fovy = 20;
    std::string col = "";

    actin_ensemble ae = actin_ensemble();
    ae.set_straight_filaments(true);
    ae.set_ld(1);
    ae.set_visc(0);
    ae.set_fov(fovx, fovy);
    ae.set_nq(2*fovx, 2*fovy);
    ae.add_polymer(startx, starty, phi0, pol_index, nmon);
    
    double lnk_len = 0, kl = 0, kb = 0;
    link_ensemble lks;
    
    ae.connect_polymers( &lks, lnk_len, kl, kb, col );

    //MOTOR
    double mx = 0.5, my = 0.5, mang = pi/2, mlen = 1;
    int state0 = 1, state1 = 0, aindex0 = 0, aindex1 = -1;
    double temp = 0, v0 = 0.25, stiffness = 0;
    double ron = 0, roff = 2, rend = 0, delta_t = 1;
    double actin_len = 0, vis = 0;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor m = motor(mx, my, mang, mlen, &ae, state0, state1, aindex0, aindex1, fovx, fovy, delta_t, temp,
            v0, stiffness, ron, roff, rend, actin_len, vis, col);

    double tol = 0.001;
    
    //IF detach that:
    //(a) new state of head is correct (was 1, now 0)
    //(b) position of head relative to end of rod is correct (0)
    //(c) position of head is correct
    //
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.5, tol);
    m.step_onehead(0); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], -1);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[0], 0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], my - 0.5*mlen*sin(mang), tol);

    // IF no detach
    // (a) new position is correct
    roff = -1;
    rend = -1;
    m = motor(mx, my, mang, mlen, &ae, state0, state1, aindex0, aindex1, fovx, fovy, delta_t, temp,
            v0, stiffness, ron, roff, rend, actin_len, vis, col);
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.5, tol);
    m.step_onehead(0); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.75, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.25, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);

}

BOOST_AUTO_TEST_CASE( move_end_detach )
{
    //ACTIN ENSEMBLE
    double startx = 0.5, starty = 0, phi0 = 0;
    int pol_index = 0, nmon = 3;
    double fovx = 20, fovy = 20;
    std::string col = "";

    actin_ensemble ae = actin_ensemble();
    ae.set_straight_filaments(true);
    ae.set_ld(1);
    ae.set_visc(0);
    ae.set_fov(fovx, fovy);
    ae.set_nq(2*fovx, 2*fovy);
    ae.add_polymer(startx, starty, phi0, pol_index, nmon);
    
    double lnk_len = 0, kl = 0, kb = 0;
    link_ensemble lks;
    
    ae.connect_polymers( &lks, lnk_len, kl, kb, col );

    //MOTOR
    double mx = 0.125, my = 0.5, mang = pi/2, mlen = 1;
    int state0 = 1, state1 = 0, aindex0 = 0, aindex1 = -1;
    double temp = 0, v0 = 0.25, stiffness = 0;
    double ron = 0, roff = 2, rend = 2, delta_t = 1;
    double actin_len = 0, vis = 0;
    
    /* Note:
     * To control "event(rate, timestep)" function:
     *  rate*timestep > 1 ==> "event" returns 1
     *  rate*timestep < 0 ==> "event" returns 0
     */
    motor m = motor(mx, my, mang, mlen, &ae, state0, state1, aindex0, aindex1, fovx, fovy, delta_t, temp,
            v0, stiffness, ron, roff, rend, actin_len, vis, col);

    double tol = 0.001;
    
    // IF position is greater than filament length
    double pos = 1.125;
    //  If last rod in filament
    //      try to detach
    //      if detach:
    //          (a) new state of head is correct (was 1, now 0)
    //          (b) position of head relative to end of rod is correct (0)
    //          (c) position of head is correct

    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 0);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], -1);
    BOOST_CHECK_EQUAL(m.get_pos_a_end()[0], 0);
    BOOST_CHECK_CLOSE(m.get_hx()[0], mx - 0.5*mlen*cos(mang), tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], my - 0.5*mlen*sin(mang), tol);
    
    //      else (don't detach): 
    //          stay
    rend = -1;
    m = motor(mx, my, mang, mlen, &ae, state0, state1, aindex0, aindex1, fovx, fovy, delta_t, temp,
            v0, stiffness, ron, roff, rend, actin_len, vis, col);
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.125, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    //  Else (not last rod in filament)
    mx = 1.125;
    aindex0 = 1;
    m = motor(mx, my, mang, mlen, &ae, state0, state1, aindex0, aindex1, fovx, fovy, delta_t, temp,
            v0, stiffness, ron, roff, rend, actin_len, vis, col);
    
    //      (a) motor moves to next rod in filament
    //      (b) motor moves appropriate distance
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 1);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.875, tol);
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], 0.125, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.875, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);
    
    // Else ( new position isn't more than filament length
    pos = 0.375;
    m.move_end_detach(0, pos); 
    BOOST_CHECK_EQUAL(m.get_states()[0], 1);
    BOOST_CHECK_EQUAL(m.get_aindex()[0], 0);
    BOOST_CHECK_CLOSE(m.get_pos_a_end()[0], pos, tol);
    BOOST_CHECK_CLOSE(m.get_hx()[0], 0.625, tol);
    BOOST_CHECK_CLOSE(m.get_hy()[0], 0, tol);

    //  (a) motor moves to appropriate new position
}
// EOF
