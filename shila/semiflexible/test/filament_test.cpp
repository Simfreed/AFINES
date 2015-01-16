#define BOOST_TEST_MODULE filament_test

#include "globals.h"
#include "Link.h"
#include "actin.h"
#include "filament.h"
#include <boost/test/unit_test.hpp>

/*
       
        filament();
        
        ~filament();
    
        void set_shear(double g);

        void update(double t);
        
        void update_bending();
        
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


BOOST_AUTO_TEST_CASE( add_polymer_test )
{

}
BOOST_AUTO_TEST_CASE( direction_test)
{
    /*
    actin * a = new actin(0.25, 0.5, 0.75, 1, 2, 2, 4, 4, 0.2);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( ae.get_direction(0)[0] , 0.731689 , tol); 
    BOOST_CHECK_CLOSE( ae.get_direction(0)[1] , 0.681639 , tol); 

    */
}

BOOST_AUTO_TEST_CASE( start_end_test)
{
/*
    actin * a = new actin(0, 0, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( ae.get_start(0)[0], -0.500000, tol); 
    BOOST_CHECK_CLOSE( ae.get_start(0)[1],  0.000000, tol); 
    BOOST_CHECK_CLOSE( ae.get_end(0)[0],  0.500000, tol); 
    BOOST_CHECK_CLOSE( ae.get_end(0)[1],  0.000000, tol); 
  */  
}

BOOST_AUTO_TEST_CASE( get_intpoint_test)
{
    /*
     * actin * a = new actin(0, 0, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double x = 0.1, y = 1;
    double tol = 0.001;
    double * intpoint = ae.get_intpoints(0,x,y);
    BOOST_CHECK_CLOSE( intpoint[0], 0.1, tol);
    BOOST_CHECK_CLOSE( intpoint[1], 0, tol);
    double ipx = intpoint[0];
    double ipy = intpoint[1];
    delete[] intpoint;
    BOOST_CHECK_CLOSE( ipx, 0.1, tol);
    BOOST_CHECK_CLOSE( ipy, 0, tol);
*/
}

BOOST_AUTO_TEST_CASE( update_bending_test )
{

    /*
     * double lnk_len = 1;
    double kl = 100;
    double kb = 2;
    string col = "yellow";
  */  
    /**************************
     * Case 1: 1 monomer      *
     * ***********************/
/*    actin * a0 = new actin(-1, 0, 0, 1, 0, 0, 0, 0, 0);
    link_ensemble * lks = new link_ensemble();
    actin_ensemble * ae = new actin_ensemble();
    ae->set_straight_filaments(true);
    ae->add_monomer(a0, 0);
    ae->connect_polymers(lks, lnk_len, kl, kb, col);
    ae->update_polymer_bending(0);

    double * frcs = ae->get_forces(0);
    BOOST_CHECK_EQUAL(frcs[0] , 0);
    BOOST_CHECK_EQUAL(frcs[1] , 0);
    BOOST_CHECK_EQUAL(frcs[2] , 0);

    delete lks;
    delete ae;
*/
    /**************************
     * Case 2: 2 monomers     *
     * ***********************/
/*    actin * a1 = new actin(-0.5, 0, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin(3.0/(2.0*sqrt(2.0)), 3.0/(2.0*sqrt(2.0)), pi/4, 1, 0, 0, 0, 0, 0);
    lks = new link_ensemble();
    ae = new actin_ensemble();
    ae->set_straight_filaments(true);
    
    ae->add_monomer(a1, 0);
    ae->add_monomer(a2, 0);
    ae->connect_polymers(lks, lnk_len, kl, kb, col);
    ae->update_polymer_bending(0);
    
    double tol = 0.001;
    
    //these are the first forces calculated in the update_polymer_bending function
    double fx = kb*(-1.5 + 3.0/(2.0*sqrt(2.0)));
    double fy = kb*(3.0/(2.0*sqrt(2.0)));
    frcs = ae->get_forces(0);
    BOOST_CHECK_CLOSE(frcs[0] , fx/2.0 , tol );
    BOOST_CHECK_CLOSE(frcs[1] , fy/2.0 , tol );
    BOOST_CHECK_CLOSE(frcs[2] , -0.25*(fy*(1.0/sqrt(2.0) + 1) - fx/sqrt(2.0)) , tol );

    frcs = ae->get_forces(1); 
    BOOST_CHECK_CLOSE(frcs[0] , sqrt(2.0)/4.0*(fx + fy) , tol );
    BOOST_CHECK_CLOSE(frcs[1] , sqrt(2.0)/4.0*(fy - fx) , tol );
    BOOST_CHECK_CLOSE(frcs[2] , 1.0/(2.0*sqrt(2.0))*(fy - fx) , tol );
    
    delete lks;
    delete ae;
*/
    /**************************
     * Case 3: n monomers    *
     * ***********************/
    // This would be really annoying to do on pen and paper, can't i just live with the fact
    // that the 2 monomer case passed?
}

// EOF
