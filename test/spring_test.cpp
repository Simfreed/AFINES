#include "filament.h"
#include "filament_ensemble.h"
#define BOOST_TEST_MODULE spring_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nbead              = 2;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  spring_length         = 1;
    double  bead_length        = spring_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    bead a; 
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    
    spring * l = f->get_spring(0);
    BOOST_CHECK_EQUAL( l->get_kl(), 100);                  
    BOOST_CHECK_EQUAL( l->get_hx()[0], 0);             
    BOOST_CHECK_EQUAL( l->get_hx()[1],  1);            
    BOOST_CHECK_EQUAL( l->get_hy()[0], 0);              
    BOOST_CHECK_EQUAL( l->get_hy()[1], 0);               

    delete f;
} 

BOOST_AUTO_TEST_CASE( get_distance_test)
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nbead              = 2;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  spring_length         = 1;
    double  bead_length        = spring_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double tol = 0.001;
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    bead a; 
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    
    spring * l = f->get_spring(0);
    
    array<double, 2> p1 = {{0.75,2}}, p2 = {{4, -4}}, p3 = {{-5, 12}};
    double d1 = 2, d2 = 5, d3 = 13;

    l->calc_intpoint(bc, 0, p1[0], p1[1]);
    BOOST_CHECK_CLOSE(l->get_distance_sq(bc, 0, p1[0], p1[1]), d1*d1, tol);

    l->calc_intpoint(bc, 0, p2[0], p2[1]);
    BOOST_CHECK_CLOSE(l->get_distance_sq(bc, 0, p2[0], p2[1]), d2*d2, tol);
    
    l->calc_intpoint(bc, 0, p3[0], p3[1]);
    BOOST_CHECK_CLOSE(l->get_distance_sq(bc, 0, p3[0], p3[1]), d3*d3, tol);
    
    delete f;
}

BOOST_AUTO_TEST_CASE( get_quadrants_test )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nbead              = 2;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  spring_length         = 1;
    double  bead_length        = spring_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "PERIODIC";
    
    //TEST 1
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    
    spring * l = f->get_spring(0);
    l->quad_update(bc, 0);
    vector<array<int,2> > quads = l->get_quadrants();
    vector<array<int,2> > expected_quads;
    expected_quads.push_back({{xgrid/2 , ygrid/2}});
    expected_quads.push_back({{xgrid/2 + 1, ygrid/2}});
    expected_quads.push_back({{xgrid/2 + 2, ygrid/2}});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 1 : Expected Quadrants : don't equal spring Quadrants : \n");
    cout<<"\nspring Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;

    //TEST 2
    spring_length = 1.5;
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    
    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    expected_quads.push_back({{0+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{1+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{2+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{3+xgrid/2,0+ygrid/2}});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 2 : Expected Quadrants : don't equal spring Quadrants : \n");
    cout<<"\nspring Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;
    
    //TEST 3
    startphi = pi/4;
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    
    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    expected_quads.push_back({{0+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{0+xgrid/2,1+ygrid/2}});
    expected_quads.push_back({{0+xgrid/2,2+ygrid/2}});
    expected_quads.push_back({{0+xgrid/2,3+ygrid/2}});
    expected_quads.push_back({{1+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{1+xgrid/2,1+ygrid/2}});
    expected_quads.push_back({{1+xgrid/2,2+ygrid/2}});
    expected_quads.push_back({{1+xgrid/2,3+ygrid/2}});
    expected_quads.push_back({{2+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{2+xgrid/2,1+ygrid/2}});
    expected_quads.push_back({{2+xgrid/2,2+ygrid/2}});
    expected_quads.push_back({{2+xgrid/2,3+ygrid/2}});
    expected_quads.push_back({{3+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{3+xgrid/2,1+ygrid/2}});
    expected_quads.push_back({{3+xgrid/2,2+ygrid/2}});
    expected_quads.push_back({{3+xgrid/2,3+ygrid/2}});
    
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 3: Expected Quadrants : don't equal spring Quadrants : \n");
    cout<<"\nspring Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;

    //TEST 4
    startphi =0; startx = 1; spring_length = 1; nbead = 2;
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    
    l = f->get_spring(0);
    l->quad_update(bc, 0);
    
    quads = l->get_quadrants();
    expected_quads.push_back({{2+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{3+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{4+xgrid/2,0+ygrid/2}});
    
    
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 4: Expected Quadrants : don't equal spring Quadrants : \n");
    cout<<"\nspring Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    delete f;
    
    //TEST 5
    startphi = pi/2; startx = -0.26; starty = 0.2; spring_length = 1; nbead = 2;
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    
    l = f->get_spring(0);
    l->quad_update(bc, 0);
    
    quads = l->get_quadrants();
    expected_quads.push_back({{-1+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{-1+xgrid/2,1+ygrid/2}});
    expected_quads.push_back({{-1+xgrid/2,2+ygrid/2}});
    expected_quads.push_back({{-1+xgrid/2,3+ygrid/2}});
    expected_quads.push_back({{ 0+xgrid/2,0+ygrid/2}});
    expected_quads.push_back({{ 0+xgrid/2,1+ygrid/2}});
    expected_quads.push_back({{ 0+xgrid/2,2+ygrid/2}});
    expected_quads.push_back({{ 0+xgrid/2,3+ygrid/2}});


    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 5: Expected Quadrants : don't equal spring Quadrants : \n");
    cout<<"\nspring Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    expected_quads.clear();
    
    
    //TEST 6: Different angles, but same set of quads
    spring_length = 3*sqrt(2) - 0.01; 
    nbead = 2;
    set<array<int, 2> > e_quads;
    for(int x = xgrid/2; x <= xgrid/2+6; x++)
        for(int y = ygrid/2; y <= ygrid/2+6; y++)
            e_quads.insert({{{x,y}}});
    
    //A
    startphi = pi/4; startx = 0; starty = 0; 
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    
    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    set<array<int, 2> > a_quads(quads.begin(), quads.end());
    
    BOOST_CHECK_MESSAGE(a_quads == e_quads, "\nTEST 6A: Expected Quadrants : don't equal spring Quadrants : \n");
    if (a_quads != e_quads){
        cout<<"\nspring Quadrants:"; 
        for_each(a_quads.begin(), a_quads.end(), intarray_printer);
        cout<<"\nExpected Quadrants:"; 
        for_each(e_quads.begin(), e_quads.end(), intarray_printer);
    }
    a_quads.clear(); 
    
    //B
    startphi = 7*pi/4; startx = 0; starty = 2.999;
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    copy(quads.begin(), quads.end(), inserter(a_quads, a_quads.end()));
    
    BOOST_CHECK_MESSAGE(a_quads == e_quads, "\nTEST 6B: Expected Quadrants : don't equal spring Quadrants : \n");
    
    if (a_quads != e_quads){
        cout<<"\nspring Quadrants:"; 
        for_each(a_quads.begin(), a_quads.end(), intarray_printer);
        cout<<"\nExpected Quadrants:"; 
        for_each(e_quads.begin(), e_quads.end(), intarray_printer);
    }
    a_quads.clear(); 
    
    //C
    startphi = -3*pi/4; startx = 2.999; starty = 2.999;
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    copy(quads.begin(), quads.end(), inserter(a_quads, a_quads.end()));
    
    BOOST_CHECK_MESSAGE(a_quads == e_quads, "\nTEST 6C: Expected Quadrants : don't equal spring Quadrants : \n");
    
    if (a_quads != e_quads){
        cout<<"\nspring Quadrants:"; 
        for_each(a_quads.begin(), a_quads.end(), intarray_printer);
        cout<<"\nExpected Quadrants:"; 
        for_each(e_quads.begin(), e_quads.end(), intarray_printer);
    }
    a_quads.clear(); 
    
    //D
    startphi = -5*pi/4; startx = 2.999; starty = 0;
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    copy(quads.begin(), quads.end(), inserter(a_quads, a_quads.end()));
    
    BOOST_CHECK_MESSAGE(a_quads == e_quads, "\nTEST 6D: Expected Quadrants : don't equal spring Quadrants : \n");
    
    if (a_quads != e_quads){
        cout<<"\nspring Quadrants:"; 
        for_each(a_quads.begin(), a_quads.end(), intarray_printer);
        cout<<"\nExpected Quadrants:"; 
        for_each(e_quads.begin(), e_quads.end(), intarray_printer);
    }
    a_quads.clear(); 
    delete f;
    
    //TEST 7: periodic boundaries
    //A : right
    startx = 24.501, starty = 0, startphi = 0, spring_length=1, bead_length=0.5;
    
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    
    expected_quads.push_back({{xgrid - 1 , ygrid/2}});
    expected_quads.push_back({{0         , ygrid/2}});
    expected_quads.push_back({{1         , ygrid/2}});
    expected_quads.push_back({{2         , ygrid/2}});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 7A : Expected Quadrants : don't equal spring Quadrants : \n");
    quads.clear();
    delete f;
    
    //B : left
    startx = -24.499, starty = 0, startphi = pi;
    
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 7B : Expected Quadrants : don't equal spring Quadrants : \n");
    quads.clear();
    expected_quads.clear();
    delete f;
    
    //B : top
    startx = 0, starty = 24.501, startphi = pi/2.0;
    
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    
    expected_quads.push_back({{xgrid/2 , 99}});
    expected_quads.push_back({{xgrid/2 , 0}});
    expected_quads.push_back({{xgrid/2 , 1}});
    expected_quads.push_back({{xgrid/2 , 2}});
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 7C : Expected Quadrants : don't equal spring Quadrants : \n");
    cout<<"\nspring Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    delete f;
    
    startx = 0, starty = -24.499, startphi = 3*pi/2.0;
    
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    l = f->get_spring(0);
    l->quad_update(bc, 0);
    quads = l->get_quadrants();
    //BOOST_CHECK_MESSAGE(quads == expected_quads, quads_error_message("7D",expected_quads,quads));
    BOOST_CHECK_MESSAGE(quads == expected_quads, "\nTEST 7D : Expected Quadrants : don't equal spring Quadrants : \n");
    cout<<"\nspring Quadrants:"; 
    for_each(quads.begin(), quads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for_each(expected_quads.begin(), expected_quads.end(), intarray_printer);
    quads.clear();
    delete f;
}

BOOST_AUTO_TEST_CASE( add_mot_test )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > actin_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double actin_rad = 0.5, spring_len = 1;
    double v0 = 0.25;
    
    double mstiff = 0, stretching = 0, bending = 0; //spring constants
    double kon = 0, koff = 2, kend = 0, fstall = 3.85, rcut = 0.063;
    double frac_force = 0;
    
    vector<double> pos1={{23.6,0,actin_rad,0}},pos2={{24.6,0,actin_rad,0}},pos3={{-24.4,0,actin_rad,0}};
    actin_sets.push_back(pos1);
    actin_sets.push_back(pos2);
    actin_sets.push_back(pos3);
    filament_ensemble * f = new filament_ensemble(actin_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc, 0, 0);
    
    spring * l0 = f->get_filament(0)->get_spring(0);
    spring * l1 = f->get_filament(0)->get_spring(1);
    
    double mx = -24.875, my = 0.5, mang = pi/2, mlen = 1;
    
    motor * m = new motor(array<double, 3>{{mx, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    
//    map<motor *, int> l0->get_mots() = l0->get_mots();
//    map<motor *, int> l1->get_mots() = l1->get_mots();

    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor(s) on spring 0 that haven't been added");
    BOOST_CHECK_MESSAGE(l1->get_mots().size() == 0, "\n"+std::to_string(l1->get_mots().size())+" motor(s) on spring 1 that haven't been added");
    m->set_f_index(0,0);
    m->set_l_index(0,0);
    BOOST_CHECK_MESSAGE(l0->get_mots()[m] == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 1, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l1->get_mots().size() == 0, "\nmotor on spring 1 without being added");
    m->set_l_index(0,1);
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l1->get_mots()[m] == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l1->get_mots().size() == 1, "\nmotor not on spring 1 after added");
    m->set_l_index(0,-1);
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l1->get_mots().size() == 0, "\nmotor not added to spring correctly");

    //increment l_index test
    m->set_f_index(0,0);
    m->set_l_index(0,0);
    BOOST_CHECK_MESSAGE(l0->get_mots()[m] == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 1, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l1->get_mots().size() == 0, "\nmotor on spring 1 without being added");
    m->set_l_index(0,1);
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor spring index not incremented");
    BOOST_CHECK_MESSAGE(l1->get_mots()[m] == 0, "\nmotor spring index not incremented");
    BOOST_CHECK_MESSAGE(l1->get_mots().size() == 1, "\nmotor spring index not incremented");
    m->set_l_index(0,-1);
  
    //add / remove motor test
    
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor(s) on spring 0 that haven't been added");
    l0->add_mot(m, 0);
    BOOST_CHECK_MESSAGE(l0->get_mots()[m] == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 1, "\nmotor not added to spring correctly");
    l0->remove_mot(m);
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor not added to spring correctly");
    
    //add / remove motor test
    motor * m2 = new motor(array<double, 3>{{mx+1, my, mang}}, mlen, f, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor(s) on spring 0 that haven't been added");
    l0->add_mot(m, 0);
    l0->add_mot(m2, 0);
    BOOST_CHECK_MESSAGE(l0->get_mots()[m] == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l0->get_mots().at(m2) == 0, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 2, "\nmotor not added to spring correctly");
    l0->remove_mot(m);
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 1, "\nmotor not added to spring correctly");
    BOOST_CHECK_MESSAGE(l0->get_mots().at(m2) == 0, "\nmotor not added to spring correctly");
    l0->remove_mot(m2);
    BOOST_CHECK_MESSAGE(l0->get_mots().size() == 0, "\nmotor not added to spring correctly");

    //delete l0;
    //delete l1;
    delete m;
    delete f;


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
