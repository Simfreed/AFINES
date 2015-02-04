#define BOOST_TEST_MODULE actin_ensemble_test

#include "actin_ensemble.h"
#include "link_ensemble.h"
#include "Link.h"
#include "globals.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    int nmonomer = 10;
    int npolymer = 2;
    double xrange = 50;
    double yrange = 50;
    double xgrid = 2*xrange;
    double ygrid = 2*yrange;
    double actin_length = 1;
    double dt = 1e-3;
    double temp = 0;
    double link_length = 1;
    double viscosity = 0.5;
    
    double actin_density = nmonomer * npolymer / (xrange * yrange);
    std::vector<double *> pos_set;
    double pos1[3] = {0,0,0};
    pos_set.push_back(pos1);
    double pos2[3] = {-1,0,3*pi/2};
    pos_set.push_back(pos2);
    actin_ensemble ae(actin_density,xrange,yrange,xgrid,ygrid,dt,temp,actin_length,viscosity,nmonomer,link_length, pos_set, -1);
    
    //Expect to have actin monomers at (0, 0), (2, 0), (4, 0), ..., (18, 0)
   for (int i = 0; i < nmonomer; i++){
       actin a( 2*i, 0, 0, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
       BOOST_CHECK_MESSAGE( a == *(ae.get_network()->at(i)), "\n" + ae.get_network()->at(i)->to_string() + "does not equal\n" + a.to_string() );
   }
    //Expect to have actin monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
   for (int i = 0; i < nmonomer; i++){
       actin a( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
       BOOST_CHECK_MESSAGE( a == *(ae.get_network()->at(i + nmonomer)), "\n" + ae.get_network()->at(i + nmonomer)->to_string() + "does not equal\n" + a.to_string() );
   }
} 

BOOST_AUTO_TEST_CASE( angle_correlations )
{
    int nmonomer = 10;
    int npolymer = 2;
    double xrange = 50;
    double yrange = 50;
    double dt = 1e-3; 
    double temp = 0;

    double xgrid = 2*xrange;
    double ygrid = 2*yrange;
    double actin_length = 1;
    double link_length = 1;
    double viscosity = 0.5;

    double actin_density = nmonomer * npolymer / (xrange * yrange);
    std::vector<double *> pos_set;
    double pos1[3] = {0,0,0};
    pos_set.push_back(pos1);
    double pos2[3] = {-1,0,3*pi/2};
    pos_set.push_back(pos2);
    
    actin_ensemble ae(actin_density,xrange,yrange,xgrid,ygrid,dt,temp,actin_length,viscosity,nmonomer,link_length, pos_set, -1);
    std::vector<double> ac0 = ae.get_angle_correlation(0);
    std::vector<double> ac1 = ae.get_angle_correlation(1);
    
    for (int i = 0; i < nmonomer; i++){
        BOOST_CHECK_MESSAGE(ac0[i]==1, "Polymer 0: ac["+std::to_string(i)+"] = "+std::to_string(ac0[i]) + " != 1");
        BOOST_CHECK_MESSAGE(ac1[i]==1, "Polymer 1: ac["+std::to_string(i)+"] = "+std::to_string(ac0[i]) + " != 1");
    }
    
    
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
