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
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{0,0}}, findex = {{0,0}}, lindex = {{0, 0}};
    vector<array<double, 3> > pos_sets;
    int nactin = 2;

    double dt = 1, temp = 0, vis = 0;
    string bc = "REFLECTIVE";
    double seed = -1;
    
    double actin_rad = 0.5, link_len = 1, spacer_len = 1;
    double v0 = 0;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 0, kend = 0, fstall = 3.85, fbreak = 3.85;
    double frac_force = 0;
    
    
    filament_ensemble * f = new filament_ensemble(actin_density, fov, nq, dt, temp, 
            actin_rad, vis, nactin, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed,0,0);
    spacer  m = spacer(array<double, 3>{{1, 1, 3.1416/2}}, spacer_len, f, state, 
            findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, vis, bc);

    BOOST_CHECK_CLOSE( m.get_hx()[0],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hx()[1],  1, tol);                   // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[0], 0.5, tol);                // 1 //
    BOOST_CHECK_CLOSE( m.get_hy()[1], 1.5, tol);               
    BOOST_CHECK_EQUAL( m.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m.get_states()[1], 0);
    
    delete f;
    
    actin_density = nactin/(fov[0]*fov[1]);
    pos_sets.push_back({{-0.5,0,0}});


    f = new filament_ensemble(actin_density, fov, nq, dt, temp, actin_rad, vis, nactin, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed,0,0);
    state = {{1,0}};

    spacer m2 = spacer(array<double, 3>{{1, 1, 0}}, spacer_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, vis, bc);
    BOOST_CHECK_EQUAL( m2.get_states()[0], 1);
    BOOST_CHECK_EQUAL( m2.get_states()[1], 0);
    
    state = {{0,1}};
    spacer m3 = spacer(array<double, 3>{{1, 1, 0}}, spacer_len, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, vis, bc);
    BOOST_CHECK_EQUAL( m3.get_states()[0], 0);
    BOOST_CHECK_EQUAL( m3.get_states()[1], 1);
   
    delete f;
} 


BOOST_AUTO_TEST_CASE( update_bending_one_attached)
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{4,4}};
    vector<vector<double> > actin_sets;
    array<int, 2> nq = {{8,8}}, state = {{1,0}}, findex = {{0,-1}}, lindex = {{0, 0}};

    double dt = 0.0001, temp = 0, vis = 0.001;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, link_len = 1;
    double v0 = 0.25;
    
    double stretching = 1, bending = 0; //spring constants
    double kon = 2, koff = 0, kend = 0, fstall = 3.85, fbreak = 3.85;
    double frac_force = 100000;
//    double tol = 0.001, zero = 1e-10;   

    vector<double> pos1={{0, 0, actin_rad,0}}, pos2={{1, 0,actin_rad,0}};
    
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    
    /*****************************
     * TEST 1: spacer ang = pi /4 *
     * **************************/

    filament_ensemble * f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc,0,0);
    
    double mx = 0.5+sqrt(2)/4, my = sqrt(2)/4, mang = pi/4, mlen = 1, mstiff=1;
    spacer m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, vis, bc);
    m.set_bending(0.04,pi/2); 
    
    double eng0, eng1;
    int tf = 1000;
    eng0 = m.get_bending_energy()[0];
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.filament_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        eng1 = m.get_bending_energy()[0];

        BOOST_CHECK_MESSAGE(eng1 <= eng0, "eng0 = "<<eng0<<" is smaller than eng1 = "<<eng1<<" at time "<<t);
        eng0 = eng1;
    }

    // delete f;
    
    /*****************************
     * TEST 2: spacer ang = -pi /4 *
     * **************************/

    f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc,0,0);
    
    my = -sqrt(2)/4, mang = -pi/4.0;
    m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, vis, bc);
    m.set_bending(0.04,pi/2); 
    tf =1000; 
    eng0 = m.get_bending_energy()[0];
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.filament_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        eng1 = m.get_bending_energy()[0];

        BOOST_CHECK_MESSAGE(eng1 <= eng0, "eng0 = "<<eng0<<" is smaller than eng1 = "<<eng1<<" at time "<<t);
        eng0 = eng1;
    }

    delete f;
    
    /*****************************
     * TEST 3: spacer ang = -3pi /4 *
     * **************************/

    f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc,0,0);
    
    mx = 0.5-sqrt(2)/4, mang = -3*pi/4.0;
    m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, vis, bc);
    m.set_bending(0.04,pi/2); 
    tf =1000; 
    eng0 = m.get_bending_energy()[0];
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.filament_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        eng1 = m.get_bending_energy()[0];

        BOOST_CHECK_MESSAGE(eng1 <= eng0, "eng0 = "<<eng0<<" is smaller than eng1 = "<<eng1<<" at time "<<t);
        eng0 = eng1;
    }

    delete f;
    
    /*****************************
     * TEST 4: spacer ang = 3pi /4 *
     * **************************/

    f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, link_len, stretching, 1, bending, frac_force, bc,0,0);
    
    my = sqrt(2)/4, mang = 3*pi/4.0;
    m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, 
            fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, fbreak, vis, bc);
    m.set_bending(0.04,pi/2); 
    tf =1000; 
    eng0 = m.get_bending_energy()[0];
    for (int t = 0; t < tf; t++){

        /*##################################*/
        m.update_angle();
        m.update_force();
        m.step_onehead(0);
        m.brownian_relax(1);
        m.filament_update();
        f->update();
        f->clear_broken();
        /*##################################*/
    
        eng1 = m.get_bending_energy()[0];

        BOOST_CHECK_MESSAGE(eng1 <= eng0, "eng0 = "<<eng0<<" is smaller than eng1 = "<<eng1<<" at time "<<t);
        eng0 = eng1;
    }

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
    spacer m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
    m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
    spacer m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 15> num_attached; 
    int nevents = 1000;
    int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i < nbead-1; i++){
        mx = i + 0.35;
        num_attached[i]=0;
        posx=0;
        posy=0;
        for (int j = 0; j < nevents; j++){
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.step_onehead(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_detached[i] += pr(m.get_states()[0]);
        }
        BOOST_CHECK_SMALL(num_detached[i] - int(nevents*koff), threesig); 
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
    spacer m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 2> num_attached; 
    int nevents = 10000;
    //int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    num_attached[0]=0;
    num_attached[1]=0;
    for (int j = 0; j < nevents; j++){
        //cout<<"\niter "<<j<<endl;
        m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
        m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
    spacer m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 2> num_attached; 
    int nevents = 10000;
    //int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    num_attached[0]=0;
    num_attached[1]=0;
    for (int j = 0; j < nevents; j++){
        //cout<<"\niter "<<j<<endl;
        m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
        m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
    spacer m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    array<int, 11> num_attached; 
    int nevents = 10000;
    int threesig = int(3*sqrt(kon*(1-kon)*double(nevents)));
    for (int i = 0; i <= int(spring_len); i++){
        mx = i - 5;
        num_attached[i]=0;
        posx=0;
        posy=0;
        for (int j = 0; j < nevents; j++){
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
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
            m = spacer(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
            m.attach(0); //should attach to {{f_index, l_index}} = {{1, 2}}, because the distance between head 0 and that spring is 0
            num_attached[i] += m.get_states()[0];
        }
        BOOST_CHECK_SMALL(num_attached[i] - int(nevents*kon), threesig); 
    }

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
