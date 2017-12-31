#include "motor.h"
#include "filament_ensemble.h"
#include "globals.h"
#define BOOST_TEST_MODULE filament_ensemble_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{

    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2};
    vector<array<double, 3> > pos_sets;

    double dt = 1, temp = 0, vis = 0;
    string bc = "PERIODIC";
    double seed = -1;
    
    double bead_rad = 0.5, link_len = 1;
    
    double stretching = 0, bending = 0; //spring constants
    double frac_force = 0;
    
    double np = 3, nmon = 11, nmon_extra = 0, pextra_bead = 0;
    filament_ensemble * f = new filament_ensemble(np, nmon, nmon_extra, pextra_bead, fov, nq, dt, temp, 
            bead_rad, vis, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed);
  
    for (int i =0; i<np; i++){
        BOOST_CHECK_MESSAGE(f->get_filament(i)->get_nlinks() == nmon - 1, "\nfilament "<<i<<" has wrong number of links");
    }

//delete f;
   
    np = 100000; nmon = 2; nmon_extra=18; pextra_bead = 0.5; 
    f = new filament_ensemble(np, nmon, nmon_extra, pextra_bead, fov, nq, dt, temp, 
            bead_rad, vis, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed);
    
    int nm_tot = 0;
    vector<int> all_nms;
    for (int i =0; i<np; i++){
        nm_tot += f->get_filament(i)->get_nbeads();
    }

    double nm_mu = double(nm_tot) / double(np);
    double tol = 0.1;
    BOOST_CHECK_CLOSE( nm_mu,  nmon+nmon_extra*pextra_bead, tol);                   // 1 //

    double nl_tot_sq = 0, nl_diff = 0;
    for (int i =0; i<np; i++){
        nl_diff = f->get_filament(i)->get_nbeads() - nm_mu;
        nl_tot_sq += nl_diff*nl_diff;
    }
    BOOST_CHECK_CLOSE(double(nl_tot_sq)/double(np), nmon_extra*pextra_bead*(1.0-pextra_bead), tol); 

    delete f;
} 



/* Functions to test : 
        void attach(int hd);

        void brownian(double t, double gamma);

        void step_twoheads();

        void bead_update();

        void update_shape();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd);

        void update_pos_a_end(int hd, double pos);

*/
// EOF
