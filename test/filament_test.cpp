#define BOOST_TEST_MODULE filament_test

#include "globals.h"
#include "Link.h"
#include "actin.h"
#include "filament.h"
#include <boost/test/unit_test.hpp>

/*
       
        ~filament();
    
        void set_shear(double g);

        void update_positions(double t);
        
        
        vector<filament *> update_stretching();

        void update_shear(); 
        
        actin * get_actin(int i);

        vector<vector<vector<int> > > get_quadrants();
        
        string write_rods();
        
        string write_links();
        
        string to_string();
        
        bool operator==(const filament& that);    
    
        vector<actin *> get_actins(unsigned int first, unsigned int last);
        
        void update_forces(int index, double f1, double f2, double f3);
*/

BOOST_AUTO_TEST_CASE( constructors_test )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nrod                = 10;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    Link l; 
    actin a; 
   
    f = new filament({startx, starty, startphi}, nrod, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    //Expect to have actins  at (0, 0), (1, 0), (2, 0), ..., (9, 0)
    //Expect to have links at (0.5, 0), (1.5, 0), (2.5, 0), ...., (8.5, 0)
    
    for (int i = 0; i < nrod - 1; i++){
        a = actin( i, 0, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_actin(i)), "\n" + f->get_actin(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid});
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    a = actin( (nrod - 1), 0, actin_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_actin((nrod - 1))), "\n" + f->get_actin((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

    delete f;
    
    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament({startx, starty, startphi}, nrod, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    //Expect to have actin monomers at (-1, 0), (-1, -1), (-1, -2), ..., (-1, -9)
    for (int i = 0; i < nrod-1; i++){
        a = actin( -1 , -i, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_actin(i)), "\n" + f->get_actin(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid} );
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    a = actin( -1, -(nrod - 1), actin_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_actin((nrod - 1))), "\n" + f->get_actin((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

    delete f;

} 

BOOST_AUTO_TEST_CASE( fracture_constructor )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nrod                = 10;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    double  shear               = 0;
    string  bc                  = "REFLECTIVE";
    
    filament * f; 
    Link l; 
    actin a;
    vector<actin *> rodvec;

    f = new filament(rodvec, {xrange, yrange}, {xgrid, ygrid}, link_length, stretching_stiffness, bending_stiffness, dt, temp, fracture_force, shear, bc);
    delete f; 
    
    for (int i = 0; i < nrod; i++){
        rodvec.push_back(new actin( i, 0, actin_length, viscosity ));
    }

    f = new filament(rodvec, {xrange, yrange}, {xgrid, ygrid}, link_length, stretching_stiffness, bending_stiffness, dt, temp, fracture_force, shear, bc);
    //Expect to have rods  at (0, 0), (2, 0), (4, 0), ..., (18, 0)
    //Expect to have links at (-1, 0), (1, 0), (3, 0), ...., (19, 0)
    
    for (int i = 0; i < nrod-1; i++){
        a = actin( i, 0, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_actin(i)), "\n" + f->get_actin(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid} );
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    
    a = actin( (nrod - 1), 0, actin_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_actin((nrod - 1))), "\n" + f->get_actin((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );
    
    delete f;
    for(unsigned int i=0; i<rodvec.size(); i++)
        delete rodvec[i];
    rodvec.clear();
} 

BOOST_AUTO_TEST_CASE( get_quadrants_test )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nactin                = 5;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    Link l; 
    actin a; 
   
    f = new filament({startx, starty, startphi}, nactin, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    vector<vector<array<int,2> > > quads, expected_quads;
    quads = f->get_quadrants(0);

    for(int i = 0; i < nactin - 1; i++){
        expected_quads.push_back(vector<array<int,2> >() );
        expected_quads[i].push_back({2*i, 0});
        expected_quads[i].push_back({2*i+1, 0});
    }

    cout<<"\nFilament Quadrants:"; 
    for (unsigned int i=0; i < quads.size(); i++) 
        for_each(quads[i].begin(), quads[i].end(), intarray_printer);
    cout<<"\nExpected Quadrants:"; 
    for (unsigned int i=0; i < expected_quads.size(); i++) 
        for_each(expected_quads[i].begin(), expected_quads[i].end(), intarray_printer);

    BOOST_CHECK_MESSAGE(expected_quads == quads, "\nExpected quadrants not equal to filament quadrants");
    delete f;  
}
BOOST_AUTO_TEST_CASE( get_actins_test )
{
/*vector<actin *> get_actins(unsigned int first, unsigned int last);*/      
    
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nrod                = 10;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double startx = 0, starty = 0, startphi = 0;
    filament * f; 
    Link l; 
    actin a; 

    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament({startx, starty, startphi}, nrod, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    vector<actin *> rods0to0  = f->get_actins(0,0);
    vector<actin *> rods0to3  = f->get_actins(0, 3);
    vector<actin *> rods5to7  = f->get_actins(5, 7);
    vector<actin *> rods9to13 = f->get_actins(9, 13);
    vector<actin *> rodsAll   = f->get_actins(0, nrod);
    
    BOOST_CHECK_EQUAL(rods0to0.size(), 0);
    BOOST_CHECK_EQUAL(rods0to3.size(), 3);
    BOOST_CHECK_EQUAL(rods5to7.size(), 2);
    BOOST_CHECK_EQUAL(rods9to13.size(), 1);
    BOOST_CHECK_EQUAL(rodsAll.size(), 10);


    //check that everythings the same with the original filament
    //Expect to have actin monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
    for (int i = 0; i < nrod-1; i++){
        a = actin( -1 , -i, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_actin(i)), "\n" + f->get_actin(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid} );
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    a = actin( -1, -(nrod - 1), actin_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_actin((nrod - 1))), "\n" + f->get_actin((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

    //check the rods arrays
    for (int i = 0; i < 3; i++){
        a = actin( -1 , -i, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods0to3[i]), "\n" + rods0to3[i]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 5; i < 7; i++){
        a = actin( -1 , -i, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods5to7[i-5]), "\n" + rods5to7[i-5]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 9; i < nrod; i++){
        a = actin( -1 , -i, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods9to13[i-9]), "\n" + rods9to13[i-9]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 0; i < nrod; i++){
        a = actin( -1 , -i, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rodsAll[i]), "\n" + rodsAll[i]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
     
    delete f; 
    
    for (unsigned int i=0; i<rods0to0.size(); i++)
        delete rods0to0[i];
    for (unsigned int i=0; i<rods0to3.size(); i++)
        delete rods0to3[i];
    for (unsigned int i=0; i<rods5to7.size(); i++)
        delete rods5to7[i];
    for (unsigned int i=0; i<rods9to13.size(); i++)
        delete rods9to13[i];
    for (unsigned int i=0; i<rodsAll.size(); i++)
        delete rodsAll[i];
    rods0to0.clear();
    rods0to3.clear(); //this was an array of pointers, and those pointers were deleted with f
    rods5to7.clear();
    rods9to13.clear();
    rodsAll.clear();
}

BOOST_AUTO_TEST_CASE( fracture_test)
{
    // vector<filament *> fracture(int node);
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nrod                = 10;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    double startx = 0, starty = 0, startphi = 0;
    
    filament *f, *f1, *f2, *f2a, *f2b, *fcheck; 
    Link l; 
    actin a; 

    f      = new filament({startx, starty, startphi}, nrod, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    fcheck = new filament({startx, starty, startphi}, nrod, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    vector<filament *> newfils = f->fracture(5);
    
    BOOST_CHECK_MESSAGE( *f==*fcheck, "\nFilament not the same after it was fractured.\n"<<fcheck->to_string()<<"\ndoes not equal\n"<<f->to_string());
    BOOST_CHECK_EQUAL(newfils.size(), 2);

    f1 = new filament({startx, starty, startphi}, 6, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    f2 = new filament({startx + 6, starty, startphi}, 4, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    //check that the new filaments are the ones we expect
 
    BOOST_CHECK_MESSAGE(*(newfils[0]) == *f1, "\n" + newfils[0]->to_string() + "\ndoes not equal\n" + f1->to_string());
    BOOST_CHECK_MESSAGE(*(newfils[1]) == *f2, "\n" + newfils[1]->to_string() + "\ndoes not equal\n" + f2->to_string());

    f2a = new filament({startx + 6, starty, startphi}, 3, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    f2b = new filament({startx + 9, starty, startphi}, 1, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    vector<filament *> newfils2 = newfils[1]->fracture(2);
    BOOST_CHECK_MESSAGE(*(newfils2[0]) == *f2a, "\n" + newfils2[0]->to_string() + "\ndoes not equal\n" + f2a->to_string());
    BOOST_CHECK_MESSAGE(*(newfils2[1]) == *f2b, "\n" + newfils2[1]->to_string() + "\ndoes not equal\n" + f2b->to_string());
 
    // Try updating all new filaments
    f->update_positions(1);
    newfils[0]->update_positions(1);
    newfils[1]->update_positions(1);

    delete f;
    delete f1;
    delete f2;
    delete fcheck;
    // newfils.clear();
    for (unsigned int i = 0; i<newfils.size(); i++){
        delete newfils[i];
    }
    newfils.clear();
    
    delete f2a; 
    delete f2b;
    for (unsigned int i = 0; i<newfils2.size(); i++){
        delete newfils2[i];
    }
    newfils2.clear();
    
}

BOOST_AUTO_TEST_CASE( stretching_test_two_beads )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double  viscosity               = 0.5;
    double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 50;
    double  fovy                    = 50;
    int     nx                      = 100;
    int     ny                      = 100;
    double  tfinal                  = 10;
   
    cout<<"\n----------STRETCHING TEST, TWO BEADS----------\n"; 
    cout<<"\nINITIALLY STRETCHED SPRING\n";
    filament * f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 

    actin * act1 = new actin(0, 0, actin_length, viscosity);
    actin * act2 = new actin(2, 0, actin_length, viscosity); 
    
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);

    cout<<"\nfilament configuration (t = 0)\n"<<f->to_string();
    
    double e0 = numeric_limits<double>::max();
    double e1;
    vector<filament *> newfils;
    for (double t = 0; t < tfinal; t+=dt){
        //cout<<"\nStretching energy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        //f->update_bending();
        f->update_positions(t);
        e1 = f->get_stretching_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\ntime: " + to_string(t) + "e0: "+ to_string(e0) + "e1: " + to_string(e1));
        e0 = e1;
        //cfg= "\nfilament configuration (t = "+to_string(double(t)*dt)+")\n"+f->to_string();
    }
    
    cout<<"\nfilament configuration (t = "<<tfinal<<")\n"<<f->to_string();
    delete f;
    delete act1;
    delete act2;

    cout<<"\nINITIALLY SQUISHED SPRING\n";
    
    f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 

    act1 = new actin(0, 0, actin_length, viscosity);
    act2 = new actin(0.5, 0, actin_length, viscosity); 
    
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);

    cout<<"\nfilament configuration (t = 0)\n"<<f->to_string();
    
 
    e0 = numeric_limits<double>::max();
    for (double t = 0; t < tfinal; t+=dt){
        //cout<<"\nStretching energy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        //f->update_bending();
        f->update_positions(t);
        e1 = f->get_stretching_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\ntime: " + to_string(t) + "e0: "+ to_string(e0) + "e1: " + to_string(e1));
        e0 = e1;
    //    cfg= "\nfilament configuration (t = "+to_string(double(t)*dt)+")\n"+f->to_string();
    }
    
    cout<<"\nfilament configuration (t = "<<tfinal<<")\n"<<f->to_string();
    
    delete f;
    delete act1;
    delete act2;
    
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
}

BOOST_AUTO_TEST_CASE( stretching_test_three_beads )
{
    //Check if the stretching energy decreases after some small number of integration steps
    double  viscosity               = 0.5;
    double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 50;
    double  fovy                    = 50;
    int     nx                      = 100;
    int     ny                      = 100;
    int     nsteps                  = 10000;
   
    
    cout<<"\n----------STRETCHING TEST, THREE BEADS----------\n";
    filament * f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 
    actin * act1 = new actin(0, 0, actin_length, viscosity);
    actin * act2 = new actin(0, 1.5, actin_length, viscosity);
    actin * act3 = new actin(0, 2, actin_length, viscosity);
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->add_actin(act3, link_length, stretching_stiffness);
    
    double e0 = numeric_limits<double>::max(), e1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending(t);
        newfils = f->update_stretching(t);
        f->update_positions(t);
        e1 = f->get_stretching_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nEnergy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
    delete act1;
    delete act2;
    delete act3;
}    

BOOST_AUTO_TEST_CASE( bending_test_three_beads )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    //Check if the stretching energy decreases after some small number of integration steps
    double  viscosity               = 0.5;
    double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 50;
    double  fovy                    = 50;
    int     nx                      = 100;
    int     ny                      = 100;
    int     nsteps                  = 10000;
   
    cout<<"\n----------BENDING TEST, THREE BEADS----------\n";
    
    filament * f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 
    actin * act1 = new actin(0, 0, actin_length, viscosity);
    actin * act2 = new actin(1, 0, actin_length, viscosity);
    actin * act3 = new actin(1, 1, actin_length, viscosity);
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->add_actin(act3, link_length, stretching_stiffness);

    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    double e0 = f->get_bending_energy(), e1;
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions(t);
        e1 = f->get_bending_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nBending energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
    delete act1;
    delete act2;
    delete act3;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
} 

BOOST_AUTO_TEST_CASE( total_energy_test_three_beads )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    //Check if the stretching energy decreases after some small number of integration steps
    double  viscosity               = 0.5;
    double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 50;
    double  fovy                    = 50;
    int     nx                      = 100;
    int     ny                      = 100;
    int     nsteps                  = 10000;
   
    cout<<"\n----------BENDING AND STRETCHING TEST, THREE BEADS----------\n";
    
    filament * f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 
    actin * act1 = new actin(0, 0, actin_length, viscosity);
    actin * act2 = new actin(0, 1.5, actin_length, viscosity);
    actin * act3 = new actin(0.5, 1.5, actin_length, viscosity);
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->add_actin(act3, link_length, stretching_stiffness);

    vector<filament *> newfils = f->update_stretching(0);
    f->update_bending(0);
    f->update_positions(0);
    
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    double e0 = f->get_total_energy(), e1;
    for (int t = 1; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions(t);
        e1 = f->get_total_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nPotential energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
    delete act1;
    delete act2;
    delete act3;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
} 

BOOST_AUTO_TEST_CASE(energy_test_12_rods)
{

    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 0.5;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 1;
    double viscosity = 0.5;
    double stretching_stiffness = 0.4;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    int nrod = 12;
    filament * f; 
    cout<<"\n----------BENDING AND STRETCHING TEST, TWELVE BEADS----------\n";
    double startx = 0, starty = 0, startphi = 0;
    
    f = new filament({startx, starty, startphi}, nrod, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, false, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    
    vector<filament *> newfils = f->update_stretching(0);
    f->update_bending(0);
    f->update_positions(0);
    
    double e0 = f->get_total_energy(), e1;
    int nsteps = 1000;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 1; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nTotal energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending(t);
        newfils = f->update_stretching(t);
        f->update_positions(t);
        e1 = f->get_total_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nTotal Energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    delete f; 
   
}

/*BOOST_AUTO_TEST_CASE( shearing_test_two_beads )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double  viscosity               = 0.5;
    double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 50;
    double  fovy                    = 50;
    int     nx                      = 100;
    int     ny                      = 100;
    
    double  shear_freq              = 10; //Hz
    double  shear_dt                = 1.0/shear_freq;
    double  tfinal                  = 10;
    double  strain_rate             = 0.15/(tfinal*shear_freq);

    cout<<"\n----------SHEARING TEST, TWO BEADS----------\n"; 
    cout<<"\nSHEARING SPRING\n";
    filament * f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 

    actin * act1 = new actin(0, 0, actin_length, viscosity);
    actin * act2 = new actin(0, 1, actin_length, viscosity); 
    
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);

    cout<<"\nfilament configuration (t = 0)\n"<<f->to_string();
    
    vector<filament *> newfils = f->update_stretching(0);
    f->update_positions(0);
    double e0 = f->get_stretching_energy();
    double e1;
    
    for (double t = dt; t < tfinal; t+=dt){
//        cout<<"\nStretching energy at time "<<t<<" : "<<e0<<"\n";
        if (fmod(t, shear_dt) < dt){
            f->update_delrx(strain_rate*t);
            f->update_shear(t);
        }
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_positions(t);
        e1 = f->get_stretching_energy();
        //BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\ntime: " + to_string(t) + "\te0: "+ to_string(e0) + "\te1: " + to_string(e1));
        e0 = e1;
        //cfg= "\nfilament configuration (t = "+to_string(double(t)*dt)+")\n"+f->to_string();
    }
    delete f; 
    delete act1;
    delete act2;
    
}*/

BOOST_AUTO_TEST_CASE( reflective_bnd_cnd_test )
{
    // vector<filament *> fracture(int node);
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    double  dt                  = 1;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 1/(4*pi*actin_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
   
    filament * f = new filament({xrange,yrange}, {xgrid,ygrid}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    actin * act1 = new actin(24, 0, actin_length, viscosity);
    f->add_actin(act1, link_length, stretching_stiffness);

    f->update_forces(0,3,0);
    f->update_positions(0);
    array<double, 2> pos = f->get_bead_position(0);

    BOOST_CHECK_EQUAL(pos[0], 21);
    BOOST_CHECK_EQUAL(pos[1], 0);

    f->update_forces(0,0,10);
    f->update_positions(0);
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], 21);
    BOOST_CHECK_EQUAL(pos[1], 10);
    
    f->update_forces(0,0,16);
    f->update_positions(0);
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], 21);
    BOOST_CHECK_EQUAL(pos[1], -6);

    delete f;
    delete act1;

}

BOOST_AUTO_TEST_CASE( periodic_bnd_cnd_pos_test )
{
    // vector<filament *> fracture(int node);
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    double  dt                  = 1;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 1/(4*pi*actin_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "PERIODIC";
   
    filament * f = new filament({xrange,yrange}, {xgrid,ygrid}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    actin * act1 = new actin(24, 0, actin_length, viscosity);
    f->add_actin(act1, link_length, stretching_stiffness);

    f->update_forces(0,3,0);
    f->update_positions(0);
    array<double, 2> pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], -23);
    BOOST_CHECK_EQUAL(pos[1], 0);

    f->update_forces(0,0,10);
    f->update_positions(0);
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], -23);
    BOOST_CHECK_EQUAL(pos[1], 10);
    
    f->update_forces(0,0,16);
    f->update_positions(0);
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], -23);
    BOOST_CHECK_EQUAL(pos[1], -24);

    delete f;
    delete act1;

}

BOOST_AUTO_TEST_CASE( reflective_bnd_cnd_stretch_test )
{
    // vector<filament *> fracture(int node);
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    double  dt                  = 1;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 1/(4*pi*actin_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
  
    double tol                  = 0.001;

    filament * f = new filament({xrange,yrange}, {xgrid,ygrid}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    actin * act1 = new actin(24, 1, actin_length, viscosity);
    actin * act2 = new actin(-23, 1, actin_length, viscosity);
    
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(47-link_length)*(47-link_length));
    f->update_positions(0);
    array<double, 2> pos = f->get_bead_position(0);
    BOOST_CHECK_CLOSE(pos[0], -22,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);
    
    pos = f->get_bead_position(1);
    BOOST_CHECK_CLOSE(pos[0], 23,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);

    delete f;
    delete act1;
    delete act2;

}

BOOST_AUTO_TEST_CASE( periodic_bnd_cnd_stretch_test )
{
    // vector<filament *> fracture(int node);
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    double  dt                  = 1;
    double  temp                = 0;
    double  link_length         = 1;
    double  actin_length        = link_length/2;
    double  viscosity           = 1/(4*pi*actin_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "PERIODIC";
    double  tol                 = 0.001;

    filament * f = new filament({xrange,yrange}, {xgrid,ygrid}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    actin * act1 = new actin(24, 1, actin_length, viscosity);
    actin * act2 = new actin(-23, 1, actin_length, viscosity);
    
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(3-link_length)*(3-link_length));
    f->update_positions(0);
    array<double, 2> pos = f->get_bead_position(0);
    BOOST_CHECK_CLOSE(pos[0], -24,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);
    
    pos = f->get_bead_position(1);
    BOOST_CHECK_CLOSE(pos[0], -25,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);

    f = new filament({xrange,yrange}, {xgrid,ygrid}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    act1 = new actin(1, 24, actin_length, viscosity);
    act2 = new actin(1, -23, actin_length, viscosity);
    
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(3-link_length)*(3-link_length));
    f->update_positions(0);
    pos = f->get_bead_position(0);
    BOOST_CHECK_CLOSE(pos[1], -24,tol);
    BOOST_CHECK_CLOSE(pos[0], 1,tol);
    
    pos = f->get_bead_position(1);
    BOOST_CHECK_CLOSE(pos[1], -25,tol);
    BOOST_CHECK_CLOSE(pos[0], 1,tol);

    delete f;
    delete act1;
    delete act2;

}

BOOST_AUTO_TEST_CASE( bending_test_three_beads_periodic )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    //Check if the stretching energy decreases after some small number of integration steps
    double  viscosity               = 0.5;
    double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 50;
    double  fovy                    = 50;
    int     nx                      = 100;
    int     ny                      = 100;
    int     nsteps                  = 10000;
   
    cout<<"\n----------BENDING TEST, THREE BEADS----------\n";
    
    filament * f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "PERIODIC"); 
    actin * act1 = new actin(24.5, 0, actin_length, viscosity);
    actin * act2 = new actin(-24.5, 0, actin_length, viscosity);
    actin * act3 = new actin(-24.5, 1, actin_length, viscosity);
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->add_actin(act3, link_length, stretching_stiffness);

    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    double e0 = f->get_bending_energy(), e1;
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions(t);
        e1 = f->get_bending_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nBending energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
    delete act1;
    delete act2;
    delete act3;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
} 

BOOST_AUTO_TEST_CASE( total_energy_test_three_beads_periodic )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    //Check if the stretching energy decreases after some small number of integration steps
    double  viscosity               = 0.5;
    double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 50;
    double  fovy                    = 50;
    int     nx                      = 100;
    int     ny                      = 100;
    int     nsteps                  = 10000;
   
    cout<<"\n----------BENDING AND STRETCHING TEST, THREE BEADS----------\n";
    
    filament * f = new filament({fovx,fovy}, {nx,ny}, dt, temp, 0, frac, bending_stiffness, "PERIODIC"); 
    actin * act1 = new actin(0, 24, actin_length, viscosity);
    actin * act2 = new actin(-1, -24.5, actin_length, viscosity);
    actin * act3 = new actin(1.5, -24.0, actin_length, viscosity);
    f->add_actin(act1, link_length, stretching_stiffness);
    f->add_actin(act2, link_length, stretching_stiffness);
    f->add_actin(act3, link_length, stretching_stiffness);

    vector<filament *> newfils = f->update_stretching(0);
    f->update_bending(0);
    f->update_positions(0);
    
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    double e0 = f->get_total_energy(), e1;
    for (int t = 1; t < nsteps; t++){
    //    cout<<"\nEnergy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions(t);
        e1 = f->get_total_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nPotential energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
    delete act1;
    delete act2;
    delete act3;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
} 
// EOF
