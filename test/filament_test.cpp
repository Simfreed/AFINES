#define BOOST_TEST_MODULE filament_test

#include "filament.h"
#include "filament_ensemble.h"
#include "motor.h"
#include <boost/test/unit_test.hpp>

/*
       
        ~filament();
    
        void set_shear(double g);

        void update_positions(double t);
        
        
        vector<filament *> update_stretching();

        void update_shear(); 
        
        bead * get_bead(int i);

        vector<vector<vector<int> > > get_quadrants();
        
        string write_rods();
        
        string write_springs();
        
        string to_string();
        
        bool operator==(const filament& that);    
    
        vector<bead *> get_beads(unsigned int first, unsigned int last);
        
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
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    spring l; 
    bead a; 
   
    f = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    //Expect to have beads  at (0, 0), (1, 0), (2, 0), ..., (9, 0)
    //Expect to have springs at (0.5, 0), (1.5, 0), (2.5, 0), ...., (8.5, 0)
    
    for (int i = 0; i < nrod - 1; i++){
        a = bead( i, 0, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_bead(i)), "\n" + f->get_bead(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = spring( spring_len, stretching_stiffness, 1, f, {{i, i+1}}, {{xrange, yrange}}, {{xgrid, ygrid}});
        BOOST_CHECK_MESSAGE( l == *(f->get_spring(i)), "\nspring : " + f->get_spring(i)->to_string() + "\ndoes not equal\nspring : " + l.to_string() );
    }
    a = bead( (nrod - 1), 0, bead_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_bead((nrod - 1))), "\n" + f->get_bead((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

    delete f;
    
    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    //Expect to have bead monomers at (-1, 0), (-1, -1), (-1, -2), ..., (-1, -9)
    for (int i = 0; i < nrod-1; i++){
        a = bead( -1 , -i, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_bead(i)), "\n" + f->get_bead(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = spring( spring_len, stretching_stiffness, 1, f, {{i, i+1}}, {{xrange, yrange}}, {{xgrid, ygrid}} );
        BOOST_CHECK_MESSAGE( l == *(f->get_spring(i)), "\nspring : " + f->get_spring(i)->to_string() + "\ndoes not equal\nspring : " + l.to_string() );
    }
    a = bead( -1, -(nrod - 1), bead_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_bead((nrod - 1))), "\n" + f->get_bead((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

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
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    double  shear               = 0;
    string  bc                  = "REFLECTIVE";
    
    filament * f; 
    spring l; 
    bead a;
    vector<bead *> rodvec;

    f = new filament(rodvec, {{xrange, yrange}}, {{xgrid, ygrid}}, spring_len, stretching_stiffness, 1, bending_stiffness, dt, temp, fracture_force, shear, bc);
    delete f; 
    
    for (int i = 0; i < nrod; i++){
        rodvec.push_back(new bead( i, 0, bead_length, viscosity ));
    }

    f = new filament(rodvec, {{xrange, yrange}}, {{xgrid, ygrid}}, spring_len, stretching_stiffness, 1, bending_stiffness, dt, temp, fracture_force, shear, bc);
    //Expect to have rods  at (0, 0), (2, 0), (4, 0), ..., (18, 0)
    //Expect to have springs at (-1, 0), (1, 0), (3, 0), ...., (19, 0)
    
    for (int i = 0; i < nrod-1; i++){
        a = bead( i, 0, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_bead(i)), "\n" + f->get_bead(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = spring( spring_len, stretching_stiffness, 1, f, {{i, i+1}}, {{xrange, yrange}}, {{xgrid, ygrid}} );
        BOOST_CHECK_MESSAGE( l == *(f->get_spring(i)), "\nspring : " + f->get_spring(i)->to_string() + "\ndoes not equal\nspring : " + l.to_string() );
    }
    
    a = bead( (nrod - 1), 0, bead_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_bead((nrod - 1))), "\n" + f->get_bead((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );
    
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
    
    int     nbead              = 5;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double  startx = 0.001, starty = 0, startphi = 0;
    
    filament * f; 
    spring l; 
    bead a; 
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    vector<vector<array<int,2> > > quads, expected_quads;
    quads = f->get_quadrants();

    for(int i = 0; i < nbead - 1; i++){
        expected_quads.push_back(vector<array<int,2> >() );
        expected_quads[i].push_back({{2*i+xgrid/2, ygrid/2}});
        expected_quads[i].push_back({{2*i+1+xgrid/2, ygrid/2}});
        expected_quads[i].push_back({{2*i+2+xgrid/2, ygrid/2}});
        expected_quads[i].push_back({{2*i+3+xgrid/2, ygrid/2}});
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
BOOST_AUTO_TEST_CASE( get_beads_test )
{
/*vector<bead *> get_beads(unsigned int first, unsigned int last);*/      
    
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nrod                = 10;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double startx = 0, starty = 0, startphi = 0;
    filament * f; 
    spring l; 
    bead a; 

    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    vector<bead *> rods0to0  = f->get_beads(0,0);
    vector<bead *> rods0to3  = f->get_beads(0, 3);
    vector<bead *> rods5to7  = f->get_beads(5, 7);
    vector<bead *> rods9to13 = f->get_beads(9, 13);
    vector<bead *> rodsAll   = f->get_beads(0, nrod);
    
    BOOST_CHECK_EQUAL(int(rods0to0.size()), 0);
    BOOST_CHECK_EQUAL(int(rods0to3.size()), 3);
    BOOST_CHECK_EQUAL(int(rods5to7.size()), 2);
    BOOST_CHECK_EQUAL(int(rods9to13.size()), 1);
    BOOST_CHECK_EQUAL(int(rodsAll.size()), 10);


    //check that everythings the same with the original filament
    //Expect to have bead monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
    for (int i = 0; i < nrod-1; i++){
        a = bead( -1 , -i, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_bead(i)), "\n" + f->get_bead(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = spring( spring_len, stretching_stiffness, 1, f, {{i, i+1}}, {{xrange, yrange}}, {{xgrid, ygrid}} );
        BOOST_CHECK_MESSAGE( l == *(f->get_spring(i)), "\nspring : " + f->get_spring(i)->to_string() + "\ndoes not equal\nspring : " + l.to_string() );
    }
    a = bead( -1, -(nrod - 1), bead_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_bead((nrod - 1))), "\n" + f->get_bead((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

    //check the rods arrays
    for (int i = 0; i < 3; i++){
        a = bead( -1 , -i, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods0to3[i]), "\n" + rods0to3[i]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 5; i < 7; i++){
        a = bead( -1 , -i, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods5to7[i-5]), "\n" + rods5to7[i-5]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 9; i < nrod; i++){
        a = bead( -1 , -i, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods9to13[i-9]), "\n" + rods9to13[i-9]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 0; i < nrod; i++){
        a = bead( -1 , -i, bead_length, viscosity );
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
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    double startx = 0, starty = 0, startphi = 0;
    
    filament *f, *f1, *f2, *f2a, *f2b, *fcheck; 
    spring l; 
    bead a; 

    f      = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    fcheck = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    vector<filament *> newfils = f->fracture(5);
    
    BOOST_CHECK_MESSAGE( *f==*fcheck, "\nFilament not the same after it was fractured.\n"<<fcheck->to_string()<<"\ndoes not equal\n"<<f->to_string());
    BOOST_CHECK_EQUAL(int(newfils.size()), 2);

    f1 = new filament({{startx, starty, startphi}}, 6, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    f2 = new filament({{startx + 6, starty, startphi}}, 4, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    //check that the new filaments are the ones we expect
 
    BOOST_CHECK_MESSAGE(*(newfils[0]) == *f1, "\n" + newfils[0]->to_string() + "\ndoes not equal\n" + f1->to_string());
    BOOST_CHECK_MESSAGE(*(newfils[1]) == *f2, "\n" + newfils[1]->to_string() + "\ndoes not equal\n" + f2->to_string());

    f2a = new filament({{startx + 6, starty, startphi}}, 3, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    f2b = new filament({{startx + 9, starty, startphi}}, 1, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    vector<filament *> newfils2 = newfils[1]->fracture(2);
    BOOST_CHECK_MESSAGE(*(newfils2[0]) == *f2a, "\n" + newfils2[0]->to_string() + "\ndoes not equal\n" + f2a->to_string());
    BOOST_CHECK_MESSAGE(*(newfils2[1]) == *f2b, "\n" + newfils2[1]->to_string() + "\ndoes not equal\n" + f2b->to_string());
 
    // Try updating all new filaments
    f->update_positions();
    newfils[0]->update_positions();
    newfils[1]->update_positions();

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
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 

    bead * act1 = new bead(0, 0, bead_length, viscosity);
    bead * act2 = new bead(2, 0, bead_length, viscosity); 
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);

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
        f->update_positions();
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
    
    f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 

    act1 = new bead(0, 0, bead_length, viscosity);
    act2 = new bead(0.5, 0, bead_length, viscosity); 
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);

    cout<<"\nfilament configuration (t = 0)\n"<<f->to_string();
    
 
    e0 = numeric_limits<double>::max();
    for (double t = 0; t < tfinal; t+=dt){
        //cout<<"\nStretching energy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        //f->update_bending();
        f->update_positions();
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
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 
    bead * act1 = new bead(0, 0, bead_length, viscosity);
    bead * act2 = new bead(0, 1.5, bead_length, viscosity);
    bead * act3 = new bead(0, 2, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->add_bead(act3, spring_len, stretching_stiffness,1);
    
    double e0 = numeric_limits<double>::max(), e1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending(t);
        newfils = f->update_stretching(t);
        f->update_positions();
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
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 
    bead * act1 = new bead(0, 0, bead_length, viscosity);
    bead * act2 = new bead(1, 0, bead_length, viscosity);
    bead * act3 = new bead(1, 1, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->add_bead(act3, spring_len, stretching_stiffness,1);

    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    f->init_ubend(); 
    double e0 = f->get_bending_energy(), e1;
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions();
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
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 
    bead * act1 = new bead(0, 0, bead_length, viscosity);
    bead * act2 = new bead(0, 1.5, bead_length, viscosity);
    bead * act3 = new bead(0.5, 1.5, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->add_bead(act3, spring_len, stretching_stiffness,1);

    vector<filament *> newfils = f->update_stretching(0);
    f->update_bending(0);
    f->update_positions();
    
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    double e0 = f->get_potential_energy(), e1;
    for (int t = 1; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions();
        e1 = f->get_potential_energy();
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
    double bead_length = 0.5;
    double dt = 1e-3;
    double temp = 0;
    double spring_len = 1;
    double viscosity = 0.5;
    double stretching_stiffness = 0.4;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    int nrod = 12;
    filament * f; 
    cout<<"\n----------BENDING AND STRETCHING TEST, TWELVE BEADS----------\n";
    double startx = 0, starty = 0, startphi = 0;
    
    f = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, false, bead_length, spring_len, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    
    vector<filament *> newfils = f->update_stretching(0);
    f->update_bending(0);
    f->update_positions();
    
    double e0 = f->get_potential_energy(), e1;
    int nsteps = 1000;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 1; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nTotal energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending(t);
        newfils = f->update_stretching(t);
        f->update_positions();
        e1 = f->get_potential_energy();
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
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "REFLECTIVE"); 

    bead * act1 = new bead(0, 0, bead_length, viscosity);
    bead * act2 = new bead(0, 1, bead_length, viscosity); 
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);

    cout<<"\nfilament configuration (t = 0)\n"<<f->to_string();
    
    vector<filament *> newfils = f->update_stretching(0);
    f->update_positions();
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
        f->update_positions();
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
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 1/(6*pi*bead_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
   
    filament * f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    bead * act1 = new bead(24, 0, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);

    f->update_forces(0,3,0);
    f->update_positions();
    array<double, 2> pos = f->get_bead_position(0);

    BOOST_CHECK_EQUAL(pos[0], 21);
    BOOST_CHECK_EQUAL(pos[1], 0);

    f->update_forces(0,0,10);
    f->update_positions();
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], 21);
    BOOST_CHECK_EQUAL(pos[1], 10);
    
    f->update_forces(0,0,16);
    f->update_positions();
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
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 1/(6*pi*bead_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "PERIODIC";
   
    filament * f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    bead * act1 = new bead(24, 0, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);

    f->update_forces(0,3,0);
    f->update_positions();
    array<double, 2> pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], -23);
    BOOST_CHECK_EQUAL(pos[1], 0);

    f->update_forces(0,0,10);
    f->update_positions();
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], -23);
    BOOST_CHECK_EQUAL(pos[1], 10);
    
    f->update_forces(0,0,16);
    f->update_positions();
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], -23);
    BOOST_CHECK_EQUAL(pos[1], -24);

    delete f;
    delete act1;

}

BOOST_AUTO_TEST_CASE( lees_edwards_bnd_cnd_pos_test )
{
    // vector<filament *> fracture(int node);
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    double  dt                  = 1;
    double  temp                = 0;
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 1/(6*pi*bead_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "LEES-EDWARDS";
    double  strain              = 0.2;

    filament * f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 

    bead * act1 = new bead(0, -23, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    
    f->update_forces(0,0,-3);
    f->update_positions();
    array<double, 2> pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], 0);
    BOOST_CHECK_EQUAL(pos[1], 24);
    
    f->update_delrx(strain*xrange/2);
    
    f->update_forces(0,0,3);
    f->update_positions();
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0], -5);
    BOOST_CHECK_EQUAL(pos[1], -23);
    
    f->update_forces(0,-21,0);
    f->update_positions();
    pos = f->get_bead_position(0);
    BOOST_CHECK_EQUAL(pos[0],  24);
    BOOST_CHECK_EQUAL(pos[1], -23);
    
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
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 1/(6*pi*bead_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
  
    double tol                  = 0.001;

    filament * f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    bead * act1 = new bead(24, 1, bead_length, viscosity);
    bead * act2 = new bead(-23, 1, bead_length, viscosity);
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(47-spring_len)*(47-spring_len));
    f->update_positions();
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
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 1/(6*pi*bead_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "PERIODIC";
    double  tol                 = 0.001;

    filament * f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    bead * act1 = new bead(24, 1, bead_length, viscosity);
    bead * act2 = new bead(-23, 1, bead_length, viscosity);
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(3-spring_len)*(3-spring_len));
    f->update_positions();
    array<double, 2> pos = f->get_bead_position(0);
    BOOST_CHECK_CLOSE(pos[0], -24,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);
    
    pos = f->get_bead_position(1);
    BOOST_CHECK_CLOSE(pos[0], -25,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);

    f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    act1 = new bead(1, 24, bead_length, viscosity);
    act2 = new bead(1, -23, bead_length, viscosity);
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(3-spring_len)*(3-spring_len));
    f->update_positions();
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

BOOST_AUTO_TEST_CASE( lees_edwards_bnd_cnd_stretch_test )
{
    // vector<filament *> fracture(int node);
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    double  dt                  = 1;
    double  temp                = 0;
    double  spring_len         = 1;
    double  bead_length        = spring_len/2;
    double  viscosity           = 1/(6*pi*bead_length);
    double  stretching_stiffness= 1;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "LEES-EDWARDS";
    double  tol                 = 0.001;
    double  strain              = 0.16; //this makes the initial distance sqrt[3^2 + 4^2] = 5

    filament * f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    f->update_delrx(strain*xrange/2);

    bead * act1 = new bead(1, 24, bead_length, viscosity);
    bead * act2 = new bead(1,-23, bead_length, viscosity);
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->update_stretching(0);
    BOOST_CHECK_CLOSE(f->get_stretching_energy(), 0.5*stretching_stiffness*(5-spring_len)*(5-spring_len),tol);
    
    f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    f->update_delrx(strain*xrange/2);

    act1 = new bead(1, 24, bead_length, viscosity);
    act2 = new bead(-3,-23, bead_length, viscosity);
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(3-spring_len)*(3-spring_len));
    f->update_positions();
    array<double, 2> pos = f->get_bead_position(0);
    BOOST_CHECK_CLOSE(pos[0], -3,tol);
    BOOST_CHECK_CLOSE(pos[1], -24,tol);
    
    pos = f->get_bead_position(1);
    BOOST_CHECK_CLOSE(pos[0], -3, tol);
    BOOST_CHECK_CLOSE(pos[1], -25,tol);
    
    f = new filament({{xrange,yrange}}, {{xgrid,ygrid}}, dt, temp, 0, fracture_force, bending_stiffness, bc); 
    f->update_delrx(strain*xrange/2);
    
    act1 = new bead(24, 1, bead_length, viscosity);
    act2 = new bead(-23, 1, bead_length, viscosity);
    
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->update_stretching(0);
    BOOST_CHECK_EQUAL(f->get_stretching_energy(), 0.5*stretching_stiffness*(3-spring_len)*(3-spring_len));
    f->update_positions();
    pos = f->get_bead_position(0);
    BOOST_CHECK_CLOSE(pos[0], -24,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);
    
    pos = f->get_bead_position(1);
    BOOST_CHECK_CLOSE(pos[0], -25,tol);
    BOOST_CHECK_CLOSE(pos[1], 1,tol);
    
    delete f;
    delete act1;
    delete act2;

}
BOOST_AUTO_TEST_CASE( bending_test_three_beads_periodic )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    //Check if the stretching energy decreases after some small number of integration steps
    double  viscosity               = 0.5;
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "PERIODIC"); 
    bead * act1 = new bead(24.5, 0, bead_length, viscosity);
    bead * act2 = new bead(-24.5, 0, bead_length, viscosity);
    bead * act3 = new bead(-24.5, 1, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->add_bead(act3, spring_len, stretching_stiffness,1);

    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    f->init_ubend(); 
    double e0 = f->get_bending_energy(), e1;
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions();
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
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "PERIODIC"); 
    bead * act1 = new bead(0, 24, bead_length, viscosity);
    bead * act2 = new bead(-1, -24.5, bead_length, viscosity);
    bead * act3 = new bead(1.5, -24.0, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->add_bead(act3, spring_len, stretching_stiffness,1);

    vector<filament *> newfils = f->update_stretching(0);
    f->update_bending(0);
    f->update_positions();
    
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    double e0 = f->get_potential_energy(), e1;
    for (int t = 1; t < nsteps; t++){
    //    cout<<"\nEnergy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions();
        e1 = f->get_potential_energy();
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

BOOST_AUTO_TEST_CASE( total_energy_test_three_beads_lees_edwards )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    //Check if the stretching energy decreases after some small number of integration steps
    double  viscosity               = 0.5;
    double  bead_length            = 0.5;
    double  spring_len             = 1;
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
    double  strain                  = 0.5;
    
    cout<<"\n----------BENDING AND STRETCHING TEST, THREE BEADS, LE, "<<strain<<" Strain----------\n";
    
    filament * f = new filament({{fovx, fovy}}, {{nx, ny}}, dt, temp, 0, frac, bending_stiffness, "LEES-EDWARDS"); 
    f->update_delrx(strain*fovx/2);
    f->update_shear(0);

    bead * act1 = new bead(0, 24, bead_length, viscosity);
    bead * act2 = new bead(-1, -24.5, bead_length, viscosity);
    bead * act3 = new bead(1.5, -24.0, bead_length, viscosity);
    f->add_bead(act1, spring_len, stretching_stiffness,1);
    f->add_bead(act2, spring_len, stretching_stiffness,1);
    f->add_bead(act3, spring_len, stretching_stiffness,1);

    vector<filament *> newfils = f->update_stretching(0);
    f->update_bending(0);
    f->update_positions();
    
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    
    double e0 = f->get_potential_energy(), e1;
    for (int t = 1; t < nsteps; t++){
    //    cout<<"\nEnergy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching(t);
        if (newfils.size()>0)
            break;
        f->update_bending(t);
        f->update_positions();
        e1 = f->get_potential_energy();
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

BOOST_AUTO_TEST_CASE( grow_test_l0)
{
    // vector<filament *> fracture(int node);
    double  xrange = 50, yrange = 50;
    int     xgrid = (int)(2*xrange), ygrid = (int)(2*yrange), nrod = 10;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  spring_length         = 1;
    double  bead_length        = spring_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    double startx = 0, starty = 0, startphi = 0;
    double l0max = 2*spring_length, lgrow = 0.5;

    filament * f;
    spring l; 
    bead a; 

    f      = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    
    f->set_l0_max(l0max);
    f->grow(lgrow);

    BOOST_CHECK_MESSAGE( f->get_spring(0)->get_l0()==spring_length+lgrow, "\nbarbed end didn't grow");
    BOOST_CHECK_MESSAGE( f->get_nsprings() == nrod-1, "\ndifferent # of springs than expected");
    BOOST_CHECK_MESSAGE( f->get_nbeads() == nrod, "\ndifferent # of beads than expected");

    delete f;
    
}

BOOST_AUTO_TEST_CASE( grow_test_add_bead)
{
    // vector<filament *> fracture(int node);
    double  xrange = 50, yrange = 50;
    int     xgrid = (int)(2*xrange), ygrid = (int)(2*yrange), nbead = 3;
    double  dt                  = 1e-3;
    double  temp                = 0;
    double  spring_length         = 1;
    double  bead_length        = spring_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "PERIODIC";
    double startx = 0, starty = 0, startphi = 0;
    double l0max = 2*spring_length, lgrow = 1.5;

    filament * f;
    spring l; 
    bead a; 

    f      = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
   
    cout<<f->to_string();
    BOOST_CHECK_MESSAGE( f->get_nsprings() == nbead-1, "\ndifferent # of springs than expected");
    BOOST_CHECK_MESSAGE( f->get_nbeads() == nbead, "\ndifferent # of beads than expected");
    
    f->set_l0_max(l0max);
    f->grow(lgrow);

    BOOST_CHECK_MESSAGE( f->get_nsprings() == nbead, "\ndifferent # of springs than expected");
    BOOST_CHECK_MESSAGE( f->get_nbeads() == nbead + 1, "\ndifferent # of beads than expected");

    for (int i = 0; i < f->get_nsprings(); i++){
        cout<<"\nDEBUG: checking spring "<<i;
        BOOST_CHECK_MESSAGE( f->get_spring(i)->get_l0()==spring_length, "\nl0 is wrong");
        array<int, 2> check = {{i, i+1}};
        BOOST_CHECK_MESSAGE( f->get_spring(i)->get_aindex()==check, "\naindex is wrong");
    }

    delete f;
    
}

BOOST_AUTO_TEST_CASE( grow_fil_test )
{
    //Filament ENSEMBLE
    array<double, 2> fov = {{50,50}};
    array<int, 2> nq = {{2,2}}, state = {{0,0}}, findex = {{-1,-1}}, lindex = {{-1, -1}};
    vector<vector<double> > bead_sets;

    double dt = 1, temp = 0.004, vis = 0;
    string bc = "PERIODIC";
    
    double bead_rad = 0.5, spring_len = 1;
    double v0 = 0;
    
    double mstiff = 0.4, stretching = 0, bending = 0; //spring constants
    double kon = 10000, koff = 0, kend = 0, fstall = 3.85, rcut = 1;
    double frac_force = 10000000;
    
    vector<double> pos1={{0,0,bead_rad,0}},pos2={{1,0,bead_rad,0}},pos3={{2,0,bead_rad,0}},pos4={{3,0,bead_rad,0}};
    bead_sets.push_back(pos1);
    bead_sets.push_back(pos2);
    bead_sets.push_back(pos3);
    bead_sets.push_back(pos4);
    filament_ensemble * fe = new filament_ensemble(bead_sets, fov, nq, dt, temp, vis, spring_len, stretching, 1, bending, frac_force, bc,0,0);
    fe->quad_update_serial(); 
    filament * f = fe->get_filament(0);
    
    f->set_l0_max(2*spring_len);

    double mx1 = -0.05, mx2 = 0.5, mx3 = 1.5, my = 0.5, mang = pi/2, mlen = 1;
    
    motor * m1 = new motor(array<double, 3>{{mx1, my, mang}}, mlen, fe, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    motor * m2 = new motor(array<double, 3>{{mx2, my, mang}}, mlen, fe, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
    motor * m3 = new motor(array<double, 3>{{mx3, my, mang}}, mlen, fe, state, findex, lindex, fov, dt, temp, v0, mstiff, 1, kon, koff, kend, fstall, rcut, vis, bc);
   
    m1->attach(0);
    m2->attach(0);
    m3->attach(0);
   
    map<motor*, int>  mots0 = f->get_spring(0)->get_mots();
    map<motor*, int>  mots1 = f->get_spring(1)->get_mots();
    map<motor*, int>  mots2 = f->get_spring(2)->get_mots();


    BOOST_CHECK_MESSAGE(mots0.size() == 2, "\ninitially not two mots on spring 0");
    BOOST_CHECK_MESSAGE(mots1.size() == 1, "\nwrong number of motors initially on spring 1");
    BOOST_CHECK_MESSAGE(mots2.size() == 0, "\nwrong number of motors initially on spring 0");
    
    BOOST_CHECK_MESSAGE(mots0.at(m1) == 0, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(mots0.at(m2) == 0, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(mots1.at(m3) == 0, "\ninitialized wrong");
    
    BOOST_CHECK_MESSAGE(f->get_nsprings() == 3, "\nwrong number of springs initially");
    BOOST_CHECK_MESSAGE(f->get_nbeads() == 4, "\nwrong number of beads initially");

    // case 1: grows, but not enough to add a bead
    f->grow(0.5);
    
    BOOST_CHECK_MESSAGE(mots0.size() == 2, "\ninitially not two mots on spring 0");
    BOOST_CHECK_MESSAGE(mots1.size() == 1, "\nwrong number of motors initially on spring 1");
    BOOST_CHECK_MESSAGE(mots2.size() == 0, "\nwrong number of motors initially on spring 2");
    
    BOOST_CHECK_MESSAGE(mots0.at(m1) == 0, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(mots0.at(m2) == 0, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(mots1.at(m3) == 0, "\ninitialized wrong");
    
    BOOST_CHECK_MESSAGE(f->get_spring(0)->get_l0() == 1.5, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(f->get_spring(1)->get_l0() == spring_len, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(f->get_spring(2)->get_l0() == spring_len, "\ninitialized wrong");
    
    BOOST_CHECK_MESSAGE(f->get_nsprings() == 3, "\nwrong number of springs initially");
    BOOST_CHECK_MESSAGE(f->get_nbeads() == 4, "\nwrong number of beads initially");
    
    BOOST_CHECK_MESSAGE(m1->get_l_index()[0] == 0, "\nm1 on wrong spring");
    BOOST_CHECK_MESSAGE(m2->get_l_index()[0] == 0, "\nm2 on wrong spring");
    BOOST_CHECK_MESSAGE(m3->get_l_index()[0] == 1, "\nm3 on wrong spring");
    
    // case 1: grows, but enough to add a bead
    f->grow(0.5);
    
    BOOST_CHECK_MESSAGE(f->get_nsprings() == 4, "\nwrong number of springs");
    BOOST_CHECK_MESSAGE(f->get_nbeads() == 5, "\nwrong number of beads");
    
    mots0 = f->get_spring(0)->get_mots();
    mots1 = f->get_spring(1)->get_mots();
    mots2 = f->get_spring(2)->get_mots();
    map<motor*, int> mots3 = f->get_spring(3)->get_mots();
    
    BOOST_CHECK_MESSAGE(mots0.size() == 1, "\nafter spring add "+std::to_string(mots0.size())+ " mots on spring 0");
    BOOST_CHECK_MESSAGE(mots1.size() == 1, "\nafter spring add "+std::to_string(mots1.size())+ " mots on spring 1");
    BOOST_CHECK_MESSAGE(mots2.size() == 1, "\nafter spring add "+std::to_string(mots2.size())+ " mots on spring 2");
    BOOST_CHECK_MESSAGE(mots3.size() == 0, "\nafter spring add "+std::to_string(mots3.size())+ " mots on spring 3");
    
    BOOST_CHECK_MESSAGE(mots0.at(m1) == 0, "\nm1 not on spring 0");
    BOOST_CHECK_MESSAGE(mots1.at(m2) == 0, "\nm2 not on spring 1");
    BOOST_CHECK_MESSAGE(mots2.at(m3) == 0, "\nm3 not on spring 2");
    
    BOOST_CHECK_MESSAGE(f->get_spring(0)->get_l0() == spring_len, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(f->get_spring(1)->get_l0() == spring_len, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(f->get_spring(2)->get_l0() == spring_len, "\ninitialized wrong");
    BOOST_CHECK_MESSAGE(f->get_spring(3)->get_l0() == spring_len, "\ninitialized wrong");
   
    BOOST_CHECK_MESSAGE(m1->get_l_index()[0] == 0, "\nm1 on wrong spring");
    BOOST_CHECK_MESSAGE(m2->get_l_index()[0] == 1, "\nm2 on wrong spring");
    BOOST_CHECK_MESSAGE(m3->get_l_index()[0] == 2, "\nm3 on wrong spring");

    cout<<"\nDEBUG testing motor detachment";
    m2->detach_head(0,{{mx2,my}});
    mots1 = f->get_spring(1)->get_mots();
    BOOST_CHECK_MESSAGE(mots1.size() == 0, "\nafter spring detached, "+std::to_string(mots1.size())+ " mots on spring 1");
    
    m2->attach(0);
    mots1 = f->get_spring(1)->get_mots();
    BOOST_CHECK_MESSAGE(mots1.size() == 1, "\nafter spring reattached, "+std::to_string(mots1.size())+ " mots on spring 1");

    delete m1;
    delete m2;
    delete m3;
    delete fe;


}

// EOF
