#define BOOST_TEST_MODULE filament_test

#include "globals.h"
#include "spring.h"
#include "bead.h"
#include "filament.h"
#include "filament_ensemble.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


BOOST_AUTO_TEST_CASE( constructors_test2 )
{

    array<double, 2> fov = {50,50};
    array<int, 2> nq = {2,2};
    vector<array<double, 3> > pos_sets;

    double dt = 1, temp = 0, vis = 0;
    string bc = "PERIODIC";
    double seed = -1;
    
    double actin_rad = 0.5, link_len = 1;
    
    double stretching = 0, bending = 0; //spring constants
    double frac_force = 0;
    
    double np = 3, nmon = 11, nmon_extra = 0, pextra_bead = 0;
    filament_ensemble * f = new filament_ensemble(np, nmon, nmon_extra, pextra_bead, fov, nq, dt, temp, 
            actin_rad, vis, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed);
  
    for (int i =0; i<np; i++){
        BOOST_CHECK_MESSAGE(f->get_filament(i)->get_nlinks() == nmon - 1, "\nfilament "<<i<<" has wrong number of links");
    }

    np = 100000; nmon = 2; nmon_extra=18; pextra_bead = 0.5; 
    f = new filament_ensemble(np, nmon, nmon_extra, pextra_bead, fov, nq, dt, temp, 
            actin_rad, vis, link_len, pos_sets, stretching, 1, bending, frac_force, bc, seed);
    
    int nm_tot = 0;
    vector<int> all_nms;
    for (int i =0; i<np; i++){
        nm_tot += f->get_filament(i)->get_nactins();
    }

    double nm_mu = double(nm_tot) / double(np);
    double tol = 0.1;
    BOOST_CHECK_CLOSE( nm_mu,  nmon+nmon_extra*pextra_bead, tol);                   // 1 //

    double nl_tot_sq = 0, nl_diff = 0;
    for (int i =0; i<np; i++){
        nl_diff = f->get_filament(i)->get_nactins() - nm_mu;
        nl_tot_sq += nl_diff*nl_diff;
    }
    BOOST_CHECK_CLOSE(double(nl_tot_sq)/double(np), nmon_extra*pextra_bead*(1.0-pextra_bead), tol); 

    delete f;
} 


BOOST_AUTO_TEST_CASE( constructors_test )
{
    double  xrange              = 50;
    double  yrange              = 50;
    int     xgrid               = (int)(2*xrange);
    int     ygrid               = (int)(2*yrange);
    
    int     nrod                = 10;
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
    spring l; 
    bead a; 
   
    f = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    //Expect to have beads  at (0, 0), (1, 0), (2, 0), ..., (9, 0)
    //Expect to have springs at (0.5, 0), (1.5, 0), (2.5, 0), ...., (8.5, 0)
    
    for (int i = 0; i < nrod - 1; i++){
        a = bead( i, 0, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_bead(i)), "\n" + f->get_bead(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = spring( spring_length, stretching_stiffness, 1, f, {i, i+1}, {{xrange, yrange}}, {{xgrid, ygrid}});
        BOOST_CHECK_MESSAGE( l == *(f->get_spring(i)), "\nspring : " + f->get_spring(i)->to_string() + "\ndoes not equal\nspring : " + l.to_string() );
    }
    a = bead( (nrod - 1), 0, bead_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_bead((nrod - 1))), "\n" + f->get_bead((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

    delete f;
    
    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament({{startx, starty, startphi}}, nrod, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    //Expect to have bead monomers at (-1, 0), (-1, -1), (-1, -2), ..., (-1, -9)
    for (int i = 0; i < nrod-1; i++){
        a = bead( -1 , -i, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_bead(i)), "\n" + f->get_bead(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = spring( spring_length, stretching_stiffness, 1, f, {i, i+1}, {{xrange, yrange}}, {{xgrid, ygrid}} );
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
    double  spring_length         = 1;
    double  bead_length        = spring_length/2;
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

    f = new filament(rodvec, {{xrange, yrange}}, {{xgrid, ygrid}}, spring_length, stretching_stiffness, 1, bending_stiffness, dt, temp, fracture_force, shear, bc);
    delete f; 
    
    for (int i = 0; i < nrod; i++){
        rodvec.push_back(new bead( i, 0, bead_length, viscosity ));
    }

    f = new filament(rodvec, {{xrange, yrange}}, {{xgrid, ygrid}}, spring_length, stretching_stiffness, 1, bending_stiffness, dt, temp, fracture_force, shear, bc);
    //Expect to have rods  at (0, 0), (2, 0), (4, 0), ..., (18, 0)
    //Expect to have springs at (-1, 0), (1, 0), (3, 0), ...., (19, 0)
    
    for (int i = 0; i < nrod-1; i++){
        a = bead( i, 0, bead_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_bead(i)), "\n" + f->get_bead(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = spring( spring_length, stretching_stiffness, 1, f, {i, i+1}, {{xrange, yrange}}, {{xgrid, ygrid}} );
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
    double  spring_length         = 1;
    double  bead_length        = spring_length/2;
    double  viscosity           = 0.5;
    double  stretching_stiffness= 100;
    double  bending_stiffness   = 1; 
    double  fracture_force      = 100;
    string  bc                  = "REFLECTIVE";
    
    double  startx = 0, starty = 0, startphi = 0;
    
    filament * f; 
    spring l; 
    bead a; 
   
    f = new filament({{startx, starty, startphi}}, nbead, {{xrange, yrange}}, {{xgrid, ygrid}}, 
            viscosity, dt, temp, true, bead_length, spring_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    vector<vector<array<int,2> > > quads, expected_quads;
    quads = f->get_quadrants();

    for(int i = 0; i < nbead - 1; i++){
        expected_quads.push_back(vector<array<int,2> >() );
        expected_quads[i].push_back({2*i+xgrid/2, ygrid/2});
        expected_quads[i].push_back({2*i+1+xgrid/2, ygrid/2});
        expected_quads[i].push_back({2*i+2+xgrid/2, ygrid/2});
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
BOOST_AUTO_TEST_CASE( filament_ensemble_update_spring_forces_from_quads ){ 

    //Check forces after one time-step to make sure they match expected values 
    double  viscosity               = 0.5;
    //double  bead_length            = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 10;
    double  fovy                    = 10;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 1.0; 

    vector<vector<double>> act = {{-0.5, 0.1, 0.25, 0}, {0.5, 0.1, 0.25, 0}, {0.0, 0.0, 0.25, 1}, {0.0, -1.0, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Test----------\n";
    cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);
	
    double f1_1_ex[2] = {0, 1.3634300623459383}; 
    double f2_1_ex[2] = {0, 1.3634300623459383}; 
    double f1_2_ex[2] = {0, -2.7268601246918767}; 
    double f2_2_ex[2] = {0, 0};
    array <double,2> f1_1; 
    array <double,2> f2_1; 
    array <double,2> f1_2; 
    array <double,2> f2_2;  
 
    network->quad_update_serial();  
    network->update_spring_forces_from_quads();

    //for(int f = 0; f < 2; f++) 
    //{ 
    //	network->update_spring_forces(f); 
    //} 
            
    f1_1 = network->get_filament(0)->get_bead(0)->get_force(); 
    f2_1 = network->get_filament(0)->get_bead(1)->get_force(); 
	
    f1_2 = network->get_filament(1)->get_bead(0)->get_force();
    f2_2 = network->get_filament(1)->get_bead(1)->get_force(); 
            
    for(int j = 0; j < 2; j++){
        
    	BOOST_CHECK_CLOSE(f1_1[j], f1_1_ex[j], 1e-10);
        BOOST_CHECK_CLOSE(f2_1[j], f2_1_ex[j], 1e-10);
	BOOST_CHECK_CLOSE(f1_2[j], f1_2_ex[j], 1e-10); 
   	BOOST_CHECK_CLOSE(f2_2[j], f2_2_ex[j], 1e-10); 

	cout << "f1_1: " << f1_1[j] << " f1_1_ex: " << f1_1_ex[j] << endl; 
	cout << "f2_1: " << f2_1[j] << " f2_1_ex: " << f2_1_ex[j] << endl; 
	cout << "f1_2: " << f1_2[j] << " f1_2_ex: " << f1_2_ex[j] << endl; 
	cout << "f2_2: " << f2_2[j] << " f2_2_ex: " << f2_2_ex[j] << endl; 
    }
} 

BOOST_AUTO_TEST_CASE( Filament_ensemble_update_test ){

    //Check positions to make sure they align with expected values after using update function 
    double  viscosity               = 0.5;
    //double  bead_length            = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 10;
    double  fovy                    = 10;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax		    = 0.25;  
    double  a                       = 1.0; 
    
    vector<vector<double>> act = {{-0.5, 0.1, 0.25, 0}, {0.5, 0.1, 0.25, 0}, {0.0, 0.0, 0.25, 1}, {0.0, -1.0, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Update_Test----------\n";
    cout<<"\nTwo filaments, each with two beads\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);
    double x_s[2][2]; 
    double y_s[2][2]; 
    double x_ev[2][2]={{-0.5, 0.5}, {0.0,0.0}}; 
    double y_ev[2][2]={{0.100578657690619854,0.100578657690619854}, {-0.001157315381239718,-1.0}}; 
    network->update(); 
    for(int i = 0; i < 2; i++){
     	for(int j = 0; j < 2; j++){
	    x_s[i][j]=network->get_filament(i)->get_bead(j)->get_xcm();
            y_s[i][j]=network->get_filament(i)->get_bead(j)->get_ycm(); 

	    cout<<"Simulated x value: "<< x_s[i][j]<<endl; 
	    cout<<"Expected x value: "<< x_ev[i][j]<<endl; 
	    cout<<"Simulated y value: "<< y_s[i][j]<<endl; 
	    cout<<"Expected y value: "<< y_ev[i][j]<<endl; 

            BOOST_CHECK_CLOSE(x_s[i][j], x_ev[i][j], 1e-1); 
	    BOOST_CHECK_CLOSE(y_s[i][j], y_ev[i][j], 1e-1); 
    	}

    }
}

BOOST_AUTO_TEST_CASE( Get_r_c_test ) 
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{-0.05, 0.353553, 0.25, 0}, {-0.807107, -0.353553, 0.25, 0}, {0.05, 0.353553, 0.25, 1}, {0.807107, -0.353553, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_R_C_Test----------\n";
    //cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);

    //array <double, 4> r_c;
    double delrx;  
    //network->update_delrx(); 
    delrx = network->get_delrx(); 
    double r1, r2, r3, r4, r1_s, r2_s, r3_s, r4_s; 

    cout << "delrx: " << delrx << endl; 

    //cout << " Problem 1" << endl; 

    r1_s=network->get_filament(0)->get_spring(0)->get_r_c(bc, delrx, 0.05, 0.353553); 
    r2_s=network->get_filament(0)->get_spring(0)->get_r_c(bc, delrx, 0.807107, -0.353553); 
    r3_s=network->get_filament(1)->get_spring(0)->get_r_c(bc, delrx, -0.05, 0.353553); 
    r4_s=network->get_filament(1)->get_spring(0)->get_r_c(bc, delrx, -0.807107, -0.353553);

    //cout << " Problem 2" << endl; 

    r1=0.1; 
    r2=1.1111402786768196; 
    r3=0.1; 
    r4=1.1111402786768196;  

    //cout << " Problem 3" << endl; 

    BOOST_CHECK_CLOSE(r1_s, r1, 1e-1);
    BOOST_CHECK_CLOSE(r2_s, r2, 1e-1);
    BOOST_CHECK_CLOSE(r3_s, r3, 1e-1);
    BOOST_CHECK_CLOSE(r4_s, r4, 1e-1);

    //cout << " Problem 4" << endl; 
    cout << "Simulated Value: " << r1_s << '\t' << "Expected Value: " << r1 << endl; 
    cout << "Simulated Value: " << r2_s << '\t' << "Expected Value: " << r2 << endl; 
    cout << "Simulated Value: " << r3_s << '\t' << "Expected Value: " << r3 << endl; 
    cout << "Simulated Value: " << r4_s << '\t' << "Expected Value: " << r4 << endl;
}
 
BOOST_AUTO_TEST_CASE( Get_r_c_test_2 )
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{-0.1, -0.5, 0.25, 0}, {-0.1, 0.5, 0.25, 0}, {0, -0.35, 0.25, 1}, {0.7071067811865475, 0.3571067811865475, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_R_C_Test_2----------\n";
    //cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);   

    array <double,4> r_c; 
    array <double,4> r_c_exp; 
    double delrx; 
    delrx = network->get_delrx(); 
    array <double,2> hx_1, hy_1, hx_2, hy_2; 

    r_c_exp[0] = 0.1; 
    r_c_exp[1] = 0.8071067811865475;
    r_c_exp[2] = 0.180278; 
    r_c_exp[3] = 0.671751; 
   
    hx_1 = network->get_filament(0)->get_spring(0)->get_hx(); 
    hy_1 = network->get_filament(0)->get_spring(0)->get_hy(); 
    hx_2 = network->get_filament(1)->get_spring(0)->get_hx(); 
    hy_2 = network->get_filament(1)->get_spring(0)->get_hy(); 			

    r_c[0] = network->get_filament(0)->get_spring(0)->get_r_c(bc, delrx, hx_2[0], hy_2[0]); 
    r_c[1] = network->get_filament(0)->get_spring(0)->get_r_c(bc, delrx, hx_2[1], hy_2[1]); 
    r_c[2] = network->get_filament(1)->get_spring(0)->get_r_c(bc, delrx, hx_1[0], hy_1[0]); 
    r_c[3] = network->get_filament(1)->get_spring(0)->get_r_c(bc, delrx, hx_1[1], hy_1[1]);    

    for(int i = 0; i < 4; i++) 
    { 
        BOOST_CHECK_CLOSE(r_c[i], r_c_exp[i], 1e-2); 
        cout << "Expected r_c: " << r_c_exp[i] << '\t' << "Simulated r_c: " << r_c[i] << endl; 
    } 
}

BOOST_AUTO_TEST_CASE( Get_r_c_test_3 )
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{0, -0.5, 0.25, 0}, {0, 0.5, 0.25, 0}, {-0.5, 0, 0.25, 1}, {0.5, 0, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_R_C_Test_3----------\n";
    //cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);

    array <double,4> r_c;
    array <double,4> r_c_exp;
    double delrx;
    delrx = network->get_delrx();
    array <double,2> hx_1, hy_1, hx_2, hy_2;
 
    r_c_exp[0] = 0.5; 
    r_c_exp[1] = 0.5; 
    r_c_exp[2] = 0.5; 
    r_c_exp[3] = 0.5; 

    hx_1 = network->get_filament(0)->get_spring(0)->get_hx(); 
    hy_1 = network->get_filament(0)->get_spring(0)->get_hy(); 
    hx_2 = network->get_filament(1)->get_spring(0)->get_hx(); 
    hy_2 = network->get_filament(1)->get_spring(0)->get_hy(); 

    r_c[0] = network->get_filament(0)->get_spring(0)->get_r_c(bc, delrx, hx_2[0], hy_2[0]);
    r_c[1] = network->get_filament(0)->get_spring(0)->get_r_c(bc, delrx, hx_2[1], hy_2[1]);
    r_c[2] = network->get_filament(1)->get_spring(0)->get_r_c(bc, delrx, hx_1[0], hy_1[0]);
    r_c[3] = network->get_filament(1)->get_spring(0)->get_r_c(bc, delrx, hx_1[1], hy_1[1]); 

    for(int i = 0; i < 4; i++) 
    { 
    	BOOST_CHECK_CLOSE(r_c[i], r_c_exp[i], 1e-3); 
        cout << "Expected r_c: " << r_c_exp[i] << '\t' << " Simulated r_c: " << r_c[i] << endl; 
    } 
}

BOOST_AUTO_TEST_CASE( get_line_intersect_test ) 
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{-0.5, 0.1, 0.25, 0}, {0.5, 0.1, 0.25, 0}, {-0.5, -0.1, 0.25, 1}, {0.5, -0.1, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Intersection_Test----------\n";
    //cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);

    bool test;  
    bool exp = false;
    double delrx; 
    delrx = network->get_delrx(); 

    test = network->get_filament(0)->get_spring(0)->get_line_intersect(bc, delrx, network->get_filament(1)->get_spring(0)); 

    BOOST_CHECK_MESSAGE(test == exp, "Expected return value is not equal to return value from get_line_intersect()");

    cout << "Expected Return Value: " << exp << '\t' << "Simulated Return Value " << test << endl; 
}

BOOST_AUTO_TEST_CASE( get_line_intersect_test_2 )
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{-0.5, 0, 0.25, 0}, {0.5, 0, 0.25, 0}, {0.35, -0.5, 0.25, 1}, {0.35, 0.5, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Intersection_Test_2----------\n";
    //cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);

    bool test;
    bool exp = true;
    double delrx;
    delrx = network->get_delrx();

    test = network->get_filament(0)->get_spring(0)->get_line_intersect(bc, delrx, network->get_filament(1)->get_spring(0));

    BOOST_CHECK_MESSAGE(test == exp, "Expected return value is not equal to return value from get_line_intersect()");

    cout << "Expected Return Value: " << exp << '\t' << "Simulated Return Value " << test << endl;
}

BOOST_AUTO_TEST_CASE( get_line_intersect_test_3 )
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{-0.5, -0.75, 0.25, 0}, {0.5, -0.75, 0.25, 0}, {-0.5, -0.25, 0.25, 1}, {0.5, 0.75, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Intersection_Test_3----------\n";
    //cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);

    bool test;
    bool exp = false;
    double delrx;
    delrx = network->get_delrx();

    test = network->get_filament(0)->get_spring(0)->get_line_intersect(bc, delrx, network->get_filament(1)->get_spring(0));

    BOOST_CHECK_MESSAGE(test == exp, "Expected return value is not equal to return value from get_line_intersect()");

    cout << "Expected Return Value: " << exp << '\t' << "Simulated Return Value " << test << endl;
}

BOOST_AUTO_TEST_CASE( get_line_intersect_test_4 )
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{-0.35355339059327373, -0.35355339059327373, 0.25, 0}, {0.35355339059327373, 0.35355339059327373, 0.25, 0}, {-0.1, 0.1, 0.25, 1}, {0.6071067811865475, -0.6071067811865475, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Intersection_Test_4----------\n";
    //cout<<"\nTwo filaments, each with two beads, update spring forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);

    bool test;
    bool exp = true;
    double delrx;
    delrx = network->get_delrx();

    test = network->get_filament(0)->get_spring(0)->get_line_intersect(bc, delrx, network->get_filament(1)->get_spring(0));

    BOOST_CHECK_MESSAGE(test == exp, "Expected return value is not equal to return value from get_line_intersect()");

    cout << "Expected return value: " << exp << '\t' << "Simulated return value: " << test << endl; 
}

BOOST_AUTO_TEST_CASE( test_update_potential )
{
    double  viscosity               = 0.5;
    double  spring_length             = 1;
    double  stretching_stiffness    = 1;
    double  bending_stiffness       = 0.04;
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 5;
    double  fovy                    = 5;
    int     nx                      = 20;
    int     ny                      = 20;
    double  ext                     = 5;
    string  bc                      = "PERIODIC";
    double  rmax                    = 0.25;
    double  a                       = 0.004;

    vector<vector<double>> act = {{0, -0.5, 0.25, 0}, {0, 0.5, 0.25, 0}, {-0.5, 0, 0.25, 1}, {0.5, 0, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Force_Test_----------\n";
    
    filament_ensemble * network;
    network = new filament_ensemble(act, {{fovx,fovy}}, {{nx,ny}}, dt, temp, viscosity, spring_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax, a);

    network->update(); 

    array <double,2> F1_s, F2_s, F3_s, F4_s; 
    array <double,2> F1_e, F2_e, F3_e, F4_e; 
    
    F1_s = network->get_filament(0)->get_bead(0)->get_force(); 
    F2_s = network->get_filament(0)->get_bead(1)->get_force(); 
    F3_s = network->get_filament(1)->get_bead(0)->get_force(); 
    F4_s = network->get_filament(1)->get_bead(1)->get_force(); 

    F1_e = {0.001414213562373095, 0.001414213562373095}; 
    F2_e = {0.001414213562373095, 0.001414213562373095};
    F3_e = {-0.001414213562373095, -0.001414213562373095};
    F4_e = {-0.001414213562373095, -0.001414213562373095}; 

    for(int i = 0; i < 2; i++) 
    { 
        BOOST_CHECK_MESSAGE(F1_s[i] == F1_e[i], "Simulated Value does not equal expected value"); 
        BOOST_CHECK_MESSAGE(F2_s[i] == F2_e[i], "Simulated Value does not equal expected value");
        BOOST_CHECK_MESSAGE(F3_s[i] == F3_e[i], "Simulated Value does not equal expected value");
        BOOST_CHECK_MESSAGE(F4_s[i] == F4_e[i], "Simulated Value does not equal expected value");

        cout << "Simulated force: " << F1_s[i] << '\t' << "Expected force: " << F1_e[i] << endl; 
	cout << "Simulated force: " << F2_s[i] << '\t' << "Expected force: " << F2_e[i] << endl;
      	cout << "Simulated force: " << F3_s[i] << '\t' << "Expected force: " << F3_e[i] << endl;
        cout << "Simulated force: " << F4_s[i] << '\t' << "Expected force: " << F4_e[i] << endl;
    }
} 
