#define BOOST_TEST_MODULE filament_test

#include "Link.h"
#include "globals.h"
#include <boost/test/unit_test.hpp>

/*
        filament(std::vector<actin *> rodvec, double linkLength, double stretching_stiffness, double bending_stiffness, 
                double deltat, double temp, double fracture);
       
        filament();
        
        ~filament();
    
        void set_shear(double g);

        void update(double t);
        
        void update_bending();
        
        std::vector<filament *> update_stretching();

        void update_shear(); 
        
        actin * get_rod(int i);

        std::vector<std::vector<std::vector<int> > > get_quadrants();
        
        std::string write_rods();
        
        std::string write_links();
        
        std::string to_string();
        
        bool operator==(const filament& that);    
    
        std::vector<actin *> get_rods(unsigned int first, unsigned int last);
        
        std::vector<filament *> fracture(int node);
        
        void update_forces(int index, double f1, double f2, double f3);
*/

BOOST_AUTO_TEST_CASE( constructors_test )
{
/*filament::filament(double startx, double starty, double startphi, int nrod, double fovx, double fovy, int nqx, int nqy, 
        double visc, double deltat, double temp, bool isStraight,
        double rodLength, double linkLength, double stretching_stiffness, double bending_stiffness,
        double frac_force)
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

    double actin_density = nmonomer * npolymer / (xrange * yrange);
    std::vector<double *> pos_set;
    double startx, starty, startphi;
    pos_set.push_back(pos2);
    filament f(startx, starty, startphi, nrod, xrange, yrange, xgrid, ygrid, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 
            bending_stiffness, fracture_force);
    
    //Expect to have rods  at (0, 0), (2, 0), (4, 0), ..., (18, 0)
    //Expect to have links at (-1, 0), (1, 0), (3, 0), ...., (19, 0)
    for (int i = 0; i < nmonomer; i++){
       actin a( 2*i, 0, 0, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
       BOOST_CHECK_MESSAGE( a == *(f.get_rod(i)), "\n" + f.get_rod(i)->to_string() + "does not equal\n" + a.to_string() );
       link l( link_length, stretching_stiffness, bending_stiffness, &f, i-1, i);
       BOOST_CHECK_MESSAGE( l == *(f.get_link(i)), "\n" + f.get_link(i)->to_string() + "does not equal\n" + l.to_string() );
    }
    link l( link_length, stretching_stiffness, bending_stiffness, &f, nrod-1, -1);
    BOOST_CHECK_MESSAGE( l == *(f.get_link(i)), "\n" + f.get_link(nrod)->to_string() + "does not equal\n" + l.to_string() );
    
   startx=-1; starty = 0; startphi= 3*pi/2;
    //Expect to have actin monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
   for (int i = 0; i < nmonomer; i++){
       actin a( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
       BOOST_CHECK_MESSAGE( a == *(f.get_rod(i)), "\n" + f.get_rod(i)->to_string() + "does not equal\n" + a.to_string() );
       link l( link_length, stretching_stiffness, bending_stiffness, &f, i-1, i);
       BOOST_CHECK_MESSAGE( l == *(f.get_link(i)), "\n" + f.get_link(i)->to_string() + "does not equal\n" + l.to_string() );
   }
   link l( link_length, stretching_stiffness, bending_stiffness, &f, nrod-1, -1);
   BOOST_CHECK_MESSAGE( l == *(f.get_link(i)), "\n" + f.get_link(nrod)->to_string() + "does not equal\n" + l.to_string() );
} 

BOOST_AUTO_TEST_CASE( fracture_constructor )
{
    
} 
BOOST_AUTO_TEST_CASE( friction_test )
{
}

BOOST_AUTO_TEST_CASE( force_test)
{
}

BOOST_AUTO_TEST_CASE( add_polymer_test )
{

}
BOOST_AUTO_TEST_CASE( direction_test)
{
    actin * a = new actin(0.25, 0.5, 0.75, 1, 2, 2, 4, 4, 0.2);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( ae.get_direction(0)[0] , 0.731689 , tol); 
    BOOST_CHECK_CLOSE( ae.get_direction(0)[1] , 0.681639 , tol); 
}

BOOST_AUTO_TEST_CASE( start_end_test)
{
    actin * a = new actin(0, 0, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( ae.get_start(0)[0], -0.500000, tol); 
    BOOST_CHECK_CLOSE( ae.get_start(0)[1],  0.000000, tol); 
    BOOST_CHECK_CLOSE( ae.get_end(0)[0],  0.500000, tol); 
    BOOST_CHECK_CLOSE( ae.get_end(0)[1],  0.000000, tol); 
    
}

BOOST_AUTO_TEST_CASE( get_intpoint_test)
{
    actin * a = new actin(0, 0, 0, 1, 0, 0, 0, 0, 0);
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

}

BOOST_AUTO_TEST_CASE( connect_polymers_test )
{
    double lnk_len = 1;
    double kl = 100;
    double kb = 2;
    std::string col = "yellow";
    
    /**************************
     * Case 1: 1 monomer      *
     * ***********************/
    actin * a0 = new actin(-1, 0, 0, 1, 0, 0, 0, 0, 0);
    link_ensemble * lks = new link_ensemble();
    actin_ensemble * ae = new actin_ensemble();
    ae->set_straight_filaments(true);
    ae->add_monomer(a0, 0);
    ae->connect_polymers(lks, lnk_len, kl, kb, col);
    
    Link * l0 = new Link(lnk_len, kl, kb, ae, -1,  0, col);
    Link * l1 = new Link(lnk_len, kl, kb, ae, 0, -1, col);
    
    BOOST_CHECK_EQUAL(lks->size(), 2);
    
    BOOST_CHECK_MESSAGE(*(lks->at(0)) == *l0, "\n" + lks->at(0)->to_string() + "does not equal\n" + l0->to_string());
    BOOST_CHECK_MESSAGE(*(lks->at(1)) == *l1, "\n" + lks->at(1)->to_string() + "does not equal\n" + l1->to_string());
    
    delete lks;
    delete ae;
    delete l0;
    delete l1;

    /**************************
     * Case 2: 2 monomers     *
     * ***********************/
    actin * a1 = new actin(-1, 0, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin(1, 0, 0, 1, 0, 0, 0, 0, 0);
    lks = new link_ensemble();
    ae = new actin_ensemble();
    
    ae->add_monomer(a1, 0);
    ae->add_monomer(a2, 0);
    ae->connect_polymers(lks, lnk_len, kl, kb, col);
    
    l0 = new Link(lnk_len, kl, kb, ae, -1, 0, col);
    l1 = new Link(lnk_len, kl, kb, ae, 0, 1, col);
    Link * l2 = new Link(lnk_len, kl, kb, ae, 1, -1, col);
    
    BOOST_CHECK_EQUAL(lks->size(), 3);
    
    BOOST_CHECK_MESSAGE(*(lks->at(0)) == *l0, "\n" + lks->at(0)->to_string() + "does not equal\n" + l0->to_string());
    BOOST_CHECK_MESSAGE(*(lks->at(1)) == *l1, "\n" + lks->at(1)->to_string() + "does not equal\n" + l1->to_string());
    BOOST_CHECK_MESSAGE(*(lks->at(2)) == *l2, "\n" + lks->at(2)->to_string() + "does not equal\n" + l2->to_string());
    
    delete lks;
    delete ae;
    delete l0;
    delete l1;
    delete l2;

    /**************************
     * Case 3: n monomers    *
     * ***********************/
    int nmonomer = 10;
    lks = new link_ensemble();
    ae = new actin_ensemble();
    ae->set_straight_filaments(true);

    for (int i = 0; i < nmonomer; i++){
        actin * a = new actin( i*2*nmonomer, 0, 0, 1, 0, 0, 0, 0, 0);
        ae->add_monomer(a, 0);
    }
    
    ae->connect_polymers(lks, lnk_len, kl, kb, col);

    for (int i = -1; i < nmonomer - 1; i++){
        l0 = new Link(lnk_len, kl, kb, ae, i, i+1, col);
        BOOST_CHECK_MESSAGE(*(lks->at(i+1)) == *l0, "\n" + lks->at(i+1)->to_string() + "does not equal\n" + l0->to_string());
        delete l0;
    }
    
    l0 = new Link(lnk_len, kl, kb, ae, nmonomer - 1, -1, col);
    BOOST_CHECK_MESSAGE(*(lks->at(nmonomer)) == *l0, "\n" + lks->at(nmonomer)->to_string() + "does not equal\n" + l0->to_string());
    
    delete lks;
    delete ae;
    delete l0;

}

BOOST_AUTO_TEST_CASE( update_bending_test )
{
    double lnk_len = 1;
    double kl = 100;
    double kb = 2;
    std::string col = "yellow";
    
    /**************************
     * Case 1: 1 monomer      *
     * ***********************/
    actin * a0 = new actin(-1, 0, 0, 1, 0, 0, 0, 0, 0);
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

    /**************************
     * Case 2: 2 monomers     *
     * ***********************/
    actin * a1 = new actin(-0.5, 0, 0, 1, 0, 0, 0, 0, 0);
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
    
    /**************************
     * Case 3: n monomers    *
     * ***********************/
    // This would be really annoying to do on pen and paper, can't i just live with the fact
    // that the 2 monomer case passed?
}

// EOF
