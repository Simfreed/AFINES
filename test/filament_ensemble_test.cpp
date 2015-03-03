#define BOOST_TEST_MODULE filament_ensemble_test

#include "filament_ensemble.h"
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
    double stretching = 100;
    double bending = .04;
    double frac_force = 150;
    string bc = "REFLECTIVE";

    
    double actin_density = nmonomer * npolymer / (xrange * yrange);
    std::vector<double *> pos_set;
    double pos1[3] = {0,0,0};
    pos_set.push_back(pos1);
    double pos2[3] = {-1,0,3*pi/2};
    pos_set.push_back(pos2);
    ATfilament_ensemble fe(actin_density, xrange, yrange, xgrid, ygrid, dt, temp,
        actin_length,  viscosity, nmonomer,  link_length, pos_set, stretching, bending, 
        frac_force, bc, -1);
    
    //Expect to have actin monomers at (0, 0), (2, 0), (4, 0), ..., (18, 0)
   for (int i = 0; i < nmonomer; i++){
       actin a( 2*i, 0, 0, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
       BOOST_CHECK_MESSAGE( a == *(fe.get_network()->at(0)->get_rod(i)), "\n" + fe.get_network()->at(0)->get_rod(i)->to_string() + "does not equal\n" + a.to_string() );
   }
    //Expect to have actin monomers at (-1, 0), (-1, -2), (-1, -4), ..., (-1, -18)
   for (int i = 0; i < nmonomer; i++){
       actin a( -1 , -2*i, 3*pi/2, actin_length, xrange, yrange, xgrid, ygrid, viscosity );
       BOOST_CHECK_MESSAGE( a == *(fe.get_network()->at(1)->get_rod(i)), "\n" + fe.get_network()->at(1)->get_rod(i)->to_string() + "does not equal\n" + a.to_string() );
   }
} 

// EOF
