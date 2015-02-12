#define BOOST_TEST_MODULE filament_test

#include "globals.h"
#include "Link.h"
#include "actin.h"
#include "filament.h"
#include <boost/test/unit_test.hpp>

/*
       
        ~filament();
    
        void set_shear(double g);

        void update(double t);
        
        
        vector<filament *> update_stretching();

        void update_shear(); 
        
        actin * get_rod(int i);

        vector<vector<vector<int> > > get_quadrants();
        
        string write_rods();
        
        string write_links();
        
        string to_string();
        
        bool operator==(const filament& that);    
    
        vector<actin *> get_rods(unsigned int first, unsigned int last);
        
        vector<filament *> fracture(int node);
        
        void update_forces(int index, double f1, double f2, double f3);
*/

BOOST_AUTO_TEST_CASE( constructors_test )
{
    int nrod = 10;
    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 1;
    double viscosity = 0.5;
    double stretching_stiffness = 100;
    double bending_stiffness = 1; 
    double fracture_force = 100;
    string bc="NONE";
    double startx = 0, starty = 0, startphi = 0;
    filament * f; 
    Link l; 
    actin a; 
   
    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    //Expect to have rods  at (0, 0), (2, 0), (4, 0), ..., (18, 0)
    //Expect to have links at (-1, 0), (1, 0), (3, 0), ...., (19, 0)
    
    for (int i = 0; i < nrod; i++){
        a = actin( 2*i, 0, 0, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_rod(i)), "\n" + f->get_rod(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, bending_stiffness, f, i-1, i);
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    l = Link( link_length, stretching_stiffness, bending_stiffness, f, nrod-1, -1);
    BOOST_CHECK_MESSAGE( l == *(f->get_link(nrod)), "\nLink : " + f->get_link(nrod)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );

    delete f;
    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    //Expect to have actin monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
    for (int i = 0; i < nrod; i++){
        a = actin( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_rod(i)), "\n" + f->get_rod(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, bending_stiffness, f, i-1, i);
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    l = Link( link_length, stretching_stiffness, bending_stiffness, f, nrod-1, -1);
    BOOST_CHECK_MESSAGE( l == *(f->get_link(nrod)), "\nLink : " + f->get_link(nrod)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );

    delete f;

} 

BOOST_AUTO_TEST_CASE( fracture_constructor )
{
    /*
    filament(vector<actin *> rodvec, double linkLength, double stretching_stiffness, double bending_stiffness, 
            double deltat, double temp, double fracture);
    */
    int nrod = 10;
    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 1;
    double viscosity = 0.5;
    double stretching_stiffness = 100;
    double bending_stiffness = 1; 
    double fracture_force = 100;
    double shear = 0;
    string bc = "NONE";
    filament * f; 
    Link l; 
    actin a;
    vector<actin *> rodvec;

    f = new filament(rodvec, link_length, stretching_stiffness, bending_stiffness, dt, temp, fracture_force, shear, bc);
    delete f; 
    
    for (int i = 0; i < nrod; i++){
        rodvec.push_back(new actin( 2*i, 0, 0, actin_length, xrange, yrange, xgrid, ygrid, viscosity ));
    }

    f = new filament(rodvec, link_length, stretching_stiffness, bending_stiffness, dt, temp, fracture_force, shear, bc);
    //Expect to have rods  at (0, 0), (2, 0), (4, 0), ..., (18, 0)
    //Expect to have links at (-1, 0), (1, 0), (3, 0), ...., (19, 0)
    
    for (int i = 0; i < nrod; i++){
        a = actin( 2*i, 0, 0, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_rod(i)), "\n" + f->get_rod(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, bending_stiffness, f, i-1, i);
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    
    l = Link( link_length, stretching_stiffness, bending_stiffness, f, nrod-1, -1);
    BOOST_CHECK_MESSAGE( l == *(f->get_link(nrod)), "\nLink : " + f->get_link(nrod)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    
    delete f;
} 

BOOST_AUTO_TEST_CASE( get_rods_test )
{
/*vector<actin *> get_rods(unsigned int first, unsigned int last);*/      
    
    int nrod = 10;
    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 1;
    double viscosity = 0.5;
    double stretching_stiffness = 100;
    double bending_stiffness = 1; 
    double fracture_force = 100;
    string bc = "NONE";
    double startx = 0, starty = 0, startphi = 0;
    filament * f; 
    Link l; 
    actin a; 

    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    vector<actin *> rods0to0  = f->get_rods(0,0);
    vector<actin *> rods0to3  = f->get_rods(0, 3);
    vector<actin *> rods5to7  = f->get_rods(5, 7);
    vector<actin *> rods9to13 = f->get_rods(9, 13);
    vector<actin *> rodsAll   = f->get_rods(0, nrod);
    
    BOOST_CHECK_EQUAL(rods0to0.size(), 0);
    BOOST_CHECK_EQUAL(rods0to3.size(), 3);
    BOOST_CHECK_EQUAL(rods5to7.size(), 2);
    BOOST_CHECK_EQUAL(rods9to13.size(), 1);
    BOOST_CHECK_EQUAL(rodsAll.size(), 10);


    //check that everythings the same with the original filament
    //Expect to have actin monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
    for (int i = 0; i < nrod; i++){
        a = actin( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_rod(i)), "\n" + f->get_rod(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, bending_stiffness, f, i-1, i);
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    l = Link( link_length, stretching_stiffness, bending_stiffness, f, nrod-1, -1);
    BOOST_CHECK_MESSAGE( l == *(f->get_link(nrod)), "\nLink : " + f->get_link(nrod)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );


    //check the rods arrays
    for (int i = 0; i < 3; i++){
        a = actin( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods0to3[i]), "\n" + rods0to3[i]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }

    for (int i = 5; i < 7; i++){
        a = actin( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods5to7[i-5]), "\n" + rods5to7[i-5]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 9; i < nrod; i++){
        a = actin( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(rods9to13[i-9]), "\n" + rods9to13[i-9]->to_string() + "\ndoes not equal\n" + a.to_string() );
    }
    
    for (int i = 0; i < nrod; i++){
        a = actin( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
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
    int nrod = 10;
    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 1;
    double viscosity = 0.5;
    double stretching_stiffness = 100;
    double bending_stiffness = 1; 
    double fracture_force = 100;
    string bc = "NONE";
    double startx = 0, starty = 0, startphi = 0;
    filament *f, *f1, *f2, *f2a, *f2b; 
    Link l; 
    actin a; 

    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);

    vector<filament *> newfils = f->fracture(5);

    f1 = new filament(startx, starty, startphi, 5, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    f2 = new filament(startx + 10, starty, startphi, 5, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    //check that everythings the same with the original filament
    //Expect to have actin monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
    for (int i = 0; i < nrod; i++){
        a = actin( 2*i, 0, 0, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_rod(i)), "\n" + f->get_rod(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, bending_stiffness, f, i-1, i);
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    l = Link( link_length, stretching_stiffness, bending_stiffness, f, nrod-1, -1);
    BOOST_CHECK_MESSAGE( l == *(f->get_link(nrod)), "\nLink : " + f->get_link(nrod)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );

    //check that the new filaments are the ones we expect
 
    BOOST_CHECK_MESSAGE(*(newfils[0]) == *f1, "\n" + newfils[0]->to_string() + "\ndoes not equal\n" + f1->to_string());
    BOOST_CHECK_MESSAGE(*(newfils[1]) == *f2, "\n" + newfils[1]->to_string() + "\ndoes not equal\n" + f2->to_string());

    f2a = new filament(startx + 10, starty, startphi, 4, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    f2b = new filament(startx + 18, starty, startphi, 1, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force, bc);
    vector<filament *> newfils2 = newfils[1]->fracture(4);
    BOOST_CHECK_MESSAGE(*(newfils2[0]) == *f2a, "\n" + newfils2[0]->to_string() + "\ndoes not equal\n" + f2a->to_string());
    BOOST_CHECK_MESSAGE(*(newfils2[1]) == *f2b, "\n" + newfils2[1]->to_string() + "\ndoes not equal\n" + f2b->to_string());
 
    // Try updating all new filaments
    f->update(1);
    newfils[0]->update(1);
    newfils[1]->update(1);

    delete f;
    delete f1;
    delete f2;
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
BOOST_AUTO_TEST_CASE( stretching_test_two_rods_parallel )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double viscosity = 0.5;
    double actin_length = 1;
    double link_length = 0.01;
    double stretching_stiffness = 1; 
    double bending_stiffness = 0.04; 
    double fovx = 50;
    double fovy = 50;
    int nx = 100;
    int  ny = 100;
    filament * f = new filament(); 
    int nsteps = 1000;
    //double angle = 0;

    actin * act1 = new actin(0, 0, 0, actin_length, fovx, fovy, nx, ny, viscosity);
    actin * act2 = new actin(2, 0, 0, actin_length, fovx, fovy, nx, ny, viscosity); 
    
    f->add_rod(act1, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act2, link_length, stretching_stiffness, bending_stiffness);

    cout<<f->to_string();
    
    double e0 = f->get_stretching_energy();
    double e1;
    vector<filament *> newfils;
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nStretching energy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching();
        if (newfils.size()>0)
            break;
        //f->update_bending();
        f->update(t);
        e1 = f->get_stretching_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nAngles not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    
    cout<<f->to_string();
    delete f;
    delete act1;
    delete act2;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
}

BOOST_AUTO_TEST_CASE( stretching_test_two_rods_perpendicular)
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double viscosity = 0.5;
    double actin_length = 1;
    double link_length = 0.01;
    double stretching_stiffness = 1; 
    double bending_stiffness = 0.04; 
    double fovx = 50;
    double fovy = 50;
    int nx = 100;
    int  ny = 100;
    filament * f = new filament(); 
    int nsteps = 1000;
    double angle = pi/2;

    actin * act1 = new actin(0, 0, 0, actin_length, fovx, fovy, nx, ny, viscosity);
    actin * act2 = new actin(1.04, 0.5, angle, actin_length, fovx, fovy, nx, ny, viscosity); 
    
    f->add_rod(act1, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act2, link_length, stretching_stiffness, bending_stiffness);

    
    double e0 = f->get_stretching_energy();
    double e1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nStretching energy at time "<<t<<" : "<<e0<<"\n";
        newfils = f->update_stretching();
        if (newfils.size()>0)
            break;
        f->update_bending();
        f->update(t);
        e1 = f->get_stretching_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nDistance not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    delete f;
    delete act1;
    delete act2;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
}

BOOST_AUTO_TEST_CASE( stretching_test_three_rods )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double viscosity = 0.5;
    double actin_length = 1;
    double link_length = 0.01;
    double stretching_stiffness = 1;
    double bending_stiffness = 0.04; 
    
    filament * f = new filament(); 
    int nsteps = 1000;

    actin * act1 = new actin(0, 0, 0, actin_length, 0, 0, 0, 0, viscosity);
    actin * act2 = new actin(1.5, 1.5, pi/2, actin_length, 0, 0, 0, 0, viscosity);
    actin * act3 = new actin(4, 2, pi, actin_length, 0, 0, 0, 0, viscosity);
    f->add_rod(act1, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act2, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act3, link_length, stretching_stiffness, bending_stiffness);
    
    double e0 = f->get_stretching_energy(), e1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
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

BOOST_AUTO_TEST_CASE( bending_test_two_rods_known_angle )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double viscosity = 0.5;
    double actin_length = 1;
    double link_length = 0.01;
    double stretching_stiffness = 10; 
    double bending_stiffness = 0.04; 
    double fovx = 50;
    double fovy = 50;
    int nx = 100;
    int  ny = 100;
    filament * f = new filament(); 
    int nsteps = 1000;
    double angle = pi/2;

    actin * act1 = new actin(0, 0, 0, actin_length, fovx, fovy, nx, ny, viscosity);
   
    actin * act2 = new actin(0.5 + link_length*cos(angle/2) + 0.5*actin_length*cos(angle), 0 + link_length*sin(angle/2) + 0.5*actin_length*sin(angle), 
            angle, actin_length, fovx, fovy, nx, ny, viscosity);
    f->add_rod(act1, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act2, link_length, stretching_stiffness, bending_stiffness);

    //cout<<f->to_string();
    
    double a0 = f->get_rod(1)->get_angle() - f->get_rod(0)->get_angle();
    double a1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching();
        if (newfils.size()>0)
            break;
        f->update_bending();
        f->update(t);
        a1 = f->get_rod(1)->get_angle() - f->get_rod(0)->get_angle();
        BOOST_CHECK_MESSAGE(fabs(a1) <= fabs(a0), "\nAngles not getting smaller at time " + to_string(t));
        a0 = a1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
    delete act1;
    delete act2;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
}    

BOOST_AUTO_TEST_CASE( bending_test_two_rods_known_negative_angle )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double viscosity = 0.5;
    double actin_length = 1;
    double link_length = 0.01;
    double stretching_stiffness = 1; 
    double bending_stiffness = 0.04; 
    double fovx = 50;
    double fovy = 50;
    int nx = 100;
    int  ny = 100;
    filament * f = new filament(); 
    int nsteps = 1000;

    actin * act1 = new actin(0, 0, 0, actin_length, fovx, fovy, nx, ny, viscosity);
    actin * act2 = new actin(0.5, -0.5, -pi/2, actin_length, fovx, fovy, nx, ny, viscosity);
    f->add_rod(act1, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act2, link_length, stretching_stiffness, bending_stiffness);
    
    double a0 = f->get_rod(1)->get_angle() - f->get_rod(0)->get_angle();
    double a1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        newfils = f->update_stretching();
        if (newfils.size()>0)
            break;
        f->update_bending();
        f->update(t);
        a1 = f->get_rod(1)->get_angle() - f->get_rod(0)->get_angle();
        BOOST_CHECK_MESSAGE(fabs(a1) <= fabs(a0), "\nAngles not getting smaller at time " + to_string(t));
        a0 = a1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
    delete act1;
    delete act2;
    if (newfils.size()>0){
        delete newfils[0];
        delete newfils[1];
    }
}    

BOOST_AUTO_TEST_CASE( bending_test_two_rods_random_angle )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    int nrod = 2;
    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 0.01;
    double viscosity = 0.5;
    double stretching_stiffness = 1;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    
    double startx = 0, starty = 0, startphi = 0;
    filament * f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, false, actin_length, link_length, stretching_stiffness, 
            bending_stiffness,fracture_force, bc);
    
    //cout<<f->to_string();
    int nsteps = 1000;

    double a0 = f->get_rod(1)->get_angle() - f->get_rod(0)->get_angle();
    double a1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        if (t%1==0)
            //cout<<"\nAngle at time "<<t<<" : "<<a0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
        a1 = f->get_rod(1)->get_angle() - f->get_rod(0)->get_angle();
        BOOST_CHECK_MESSAGE(fabs(a1) <= fabs(a0), "\nAngles not getting smaller at time " + to_string(t));
        a0 = a1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    
    delete f;
}    

BOOST_AUTO_TEST_CASE( bending_test_three_rods_known_angle )
{
    //Check if the angle between two rods is smaller after some small number of integration steps
    double viscosity = 0.5;
    double actin_length = 1;
    double link_length = 0.01;
    double stretching_stiffness = 1;
    double bending_stiffness = 0.04; 
    
    filament * f = new filament(); 
    int nsteps = 1000;

    actin * act1 = new actin(0, 0, 0, actin_length, 0, 0, 0, 0, viscosity);
    actin * act2 = new actin(0.5, 0.5, pi/2, actin_length, 0, 0, 0, 0, viscosity);
    actin * act3 = new actin(0, 1, pi, actin_length, 0, 0, 0, 0, viscosity);
    f->add_rod(act1, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act2, link_length, stretching_stiffness, bending_stiffness);
    f->add_rod(act3, link_length, stretching_stiffness, bending_stiffness);
    
    double e0 = f->get_bending_energy(), e1;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
            //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
        e1 = f->get_bending_energy();
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

BOOST_AUTO_TEST_CASE(bending_test_12_rods)
{

    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 0.01;
    double viscosity = 0.5;
    double stretching_stiffness = 1;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    int nrod = 12;
    filament * f; 
    //cout<<"\n"<<nrod<<" Rod Bending Test";
    double startx = 0, starty = 0, startphi = 0;
    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, false, actin_length, link_length, stretching_stiffness, 
            bending_stiffness,fracture_force, bc);
    
    double e0 = f->get_bending_energy(), e1;
    int nsteps = 1000;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
        e1 = f->get_bending_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nEnergy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    delete f; 
   
}

BOOST_AUTO_TEST_CASE(bending_test_12_rods_nonzero_temp)
{

    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-5;
    double temp = 0.004;
    double link_length = 0.01;
    double viscosity = 0.5;
    double stretching_stiffness = 1;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    int nrod = 12;
    filament * f; 
    //cout<<"\n"<<nrod<<" Rod Bending Test";
    double startx = 0, starty = 0, startphi = 0;
    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, false, actin_length, link_length, stretching_stiffness, 
            bending_stiffness,fracture_force, bc);
    
    double e0 = f->get_bending_energy(), e1;
    int nsteps = 100;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
        e1 = f->get_bending_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nEnergy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string();
    delete f; 
   
}

BOOST_AUTO_TEST_CASE(energy_test_2_rods)
{

    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-4;
    double temp = 0;
    double link_length = 0.01;
    double viscosity = 0.5;
    double stretching_stiffness = 1;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    int nrod = 2;
    filament * f; 
    //cout<<"\n"<<nrod<<" Rod Bending Test";
    double startx = 0, starty = 0, startphi = 0;
    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, false, actin_length, link_length, stretching_stiffness, 
            bending_stiffness,fracture_force, bc);
    
    double e0 = f->get_total_energy(), e1;
    int nsteps = 1000;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
        e1 = f->get_total_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nTotal Energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    delete f; 
   
}

BOOST_AUTO_TEST_CASE(energy_test_3_rods)
{

    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-4;
    double temp = 0;
    double link_length = 0.01;
    double viscosity = 0.5;
    double stretching_stiffness = 1;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    int nrod = 3;
    filament * f; 
    //cout<<"\n"<<nrod<<" Rod Bending Test";
    double startx = 0, starty = 0, startphi = 0;
    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, false, actin_length, link_length, stretching_stiffness, 
            bending_stiffness,fracture_force, bc);
    
    double e0 = f->get_total_energy(), e1;
    int nsteps = 1000;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nBending energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
        e1 = f->get_total_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nTotal Energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    delete f; 
   
}

BOOST_AUTO_TEST_CASE(energy_test_12_rods)
{

    double xrange = 50;
    double yrange = 50;
    int xgrid = (int)(2*xrange);
    int ygrid = (int)(2*yrange);
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0.00;
    double link_length = 0.01;
    double viscosity = 0.5;
    double stretching_stiffness = 0.4;
    double bending_stiffness = 0.04; 
    double fracture_force = 100000;
    string bc="REFLECTIVE";
    int nrod = 12;
    filament * f; 
    //cout<<"\n"<<nrod<<" Rod Bending Test";
    double startx = 0, starty = 0, startphi = 0;
    f = new filament(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, false, actin_length, link_length, stretching_stiffness, 
            bending_stiffness,fracture_force, bc);
    
    double e0 = f->get_total_energy(), e1;
    int nsteps = 100;
    vector<filament *> newfils;
    cout<<"\nINITIAL CONFIGURATION\n";
    cout<<f->to_string();
    for (int t = 0; t < nsteps; t++){
        //if (t%100==0)
        //cout<<"\nTotal energy at iteration "<<t<<" : "<<e0<<"\n";
        f->update_bending();
        newfils = f->update_stretching();
        f->update(t);
        e1 = f->get_total_energy();
        BOOST_CHECK_MESSAGE(fabs(e1) <= fabs(e0), "\nTotal Energy not getting smaller at time " + to_string(t));
        e0 = e1;
    }
    cout<<"\nFINAL CONFIGURATION\n";
    cout<<f->to_string(); 
    delete f; 
   
}
// EOF
