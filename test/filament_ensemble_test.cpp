#define BOOST_TEST_MODULE filament_test

#include "globals.h"
#include "Link.h"
#include "actin.h"
#include "filament.h"
#include "filament_ensemble.h"
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
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    //Expect to have actins  at (0, 0), (1, 0), (2, 0), ..., (9, 0)
    //Expect to have links at (0.5, 0), (1.5, 0), (2.5, 0), ...., (8.5, 0)
    
    for (int i = 0; i < nrod - 1; i++){
        a = actin( i, 0, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_actin(i)), "\n" + f->get_actin(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, 1, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid});
        BOOST_CHECK_MESSAGE( l == *(f->get_link(i)), "\nLink : " + f->get_link(i)->to_string() + "\ndoes not equal\nLink : " + l.to_string() );
    }
    a = actin( (nrod - 1), 0, actin_length, viscosity );
    BOOST_CHECK_MESSAGE( a == *(f->get_actin((nrod - 1))), "\n" + f->get_actin((nrod - 1))->to_string() + "\ndoes not equal\n" + a.to_string() );

    delete f;
    
    startx=-1; starty = 0; startphi= 3*pi/2;

    f = new filament({startx, starty, startphi}, nrod, {xrange, yrange}, {xgrid, ygrid}, 
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);
    //Expect to have actin monomers at (-1, 0), (-1, -1), (-1, -2), ..., (-1, -9)
    for (int i = 0; i < nrod-1; i++){
        a = actin( -1 , -i, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_actin(i)), "\n" + f->get_actin(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, 1, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid} );
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

    f = new filament(rodvec, {xrange, yrange}, {xgrid, ygrid}, link_length, stretching_stiffness, 1, bending_stiffness, dt, temp, fracture_force, shear, bc);
    delete f; 
    
    for (int i = 0; i < nrod; i++){
        rodvec.push_back(new actin( i, 0, actin_length, viscosity ));
    }

    f = new filament(rodvec, {xrange, yrange}, {xgrid, ygrid}, link_length, stretching_stiffness, 1, bending_stiffness, dt, temp, fracture_force, shear, bc);
    //Expect to have rods  at (0, 0), (2, 0), (4, 0), ..., (18, 0)
    //Expect to have links at (-1, 0), (1, 0), (3, 0), ...., (19, 0)
    
    for (int i = 0; i < nrod-1; i++){
        a = actin( i, 0, actin_length, viscosity );
        BOOST_CHECK_MESSAGE( a == *(f->get_actin(i)), "\n" + f->get_actin(i)->to_string() + "\ndoes not equal\n" + a.to_string() );
        l = Link( link_length, stretching_stiffness, 1, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid} );
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
    
    int     nactin              = 5;
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
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 1, 
            bending_stiffness, fracture_force, bc);

    vector<vector<array<int,2> > > quads, expected_quads;
    quads = f->get_quadrants();

    for(int i = 0; i < nactin - 1; i++){
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
/*BOOST_AUTO_TEST_CASE( get_actins_test )
{
vector<actin *> get_actins(unsigned int first, unsigned int last);      
    
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
            viscosity, dt, temp, true, actin_length, link_length, stretching_stiffness, 1, 
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
        l = Link( link_length, stretching_stiffness, 1, f, {i, i+1}, {xrange, yrange}, {xgrid, ygrid} );
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
*/

BOOST_AUTO_TEST_CASE( Filament_ensemble_force_test_two_beads )
{
    //Check forces on actin beads after one time-step to make sure they match expected values 
    double  viscosity               = 0.5;
    //double  actin_length            = 0.5;
    double  link_length             = 1;
    double  stretching_stiffness    = 1; 
    double  bending_stiffness       = 0.04; 
    double  dt                      = 1e-3;
    double  temp                    = 0;
    double  frac                    = 1e10;
    double  fovx                    = 10;
    double  fovy                    = 10;
    int     nx                      = 20;
    int     ny                      = 20; 
    double  ext 		    = 5; 
    string  bc                      = "PERIODIC";  
    double  rmax 		    = 0.25; 
    
    vector<vector<double>> act = {{-0.5, 0.1, 0.25, 0}, {0.5, 0.1, 0.25, 0}, {-0.5, 0.0, 0.25, 1}, {0.5, 0.0, 0.25, 1}}; 

    cout<<"\n----------Filament_Ensemble_Test----------\n";
    cout<<"\nTwo filaments, each with two beads\n";

    filament_ensemble * network; 
 
    network = new filament_ensemble(act, {fovx,fovy}, {nx,ny}, dt, temp, viscosity, link_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax); 

//    actin * act1 = new actin(0, 0, actin_length, viscosity);
//    actin * act2 = new actin(2, 0, actin_length, viscosity); 
    
//    f->add_actin(act1, link_length, stretching_stiffness,1);
//    f->add_actin(act2, link_length, stretching_stiffness,1);

//    cout<<"\nfilament configuration (t = 0)\n"<<f->to_string();
    
    double e0_f1[2][2] = {{0, 4.8}, {0, 4.8}};
    double e0_f2[2][2] = {{0, -4.8}, {0, -4.8}}; 
    array<double,2> e1_f1;
    array<double,2> e1_f2;
    for (double f = 0; f < 2; f++){
        network->update_excluded_volume(f);
        for(int i = 0; i < network->get_filament(f)->get_nactins(); i++){
            if(f==1){e1_f2 = network->get_filament(f)->get_actin(i)->get_force();}
	    if(f==1){e1_f1 = network->get_filament(f-1)->get_actin(i)->get_force();}
	    for(int j = 0; j < 2; j++){
                if(f==1){
    		    BOOST_CHECK_CLOSE(e1_f1[j], e0_f1[i][j], 1e-10);// "e0_f1: "+ to_string(e0_f1[i][j]) + " e1_f1: " + to_string(e1_f1[j]));
		    BOOST_CHECK_CLOSE(e1_f2[j], e0_f2[i][j], 1e-10);// "e0_f2: "+ to_string(e0_f2[i][j]) + " e1_f2: " + to_string(e1_f2[j]));
                   
   		    cout << "f1: " << e1_f1[j] << " f1_ex: " << e0_f1[i][j] << endl; 
		    cout << "f2: " << e1_f2[j] << " f2_ex: " << e0_f2[i][j] << endl; 
		}
	    }        
        }
    }
}

BOOST_AUTO_TEST_CASE( Filament_enseble_update_link_forces_test ){ 

    //Check forces after one time-step to make sure they match expected values 
    double  viscosity               = 0.5;
    //double  actin_length            = 0.5;
    double  link_length             = 1;
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

    vector<vector<double>> act = {{-0.6, 0.5, 0.25, 0}, {-0.6, -0.5, 0.25, 0}, {-0.6, 0.0, 0.25, 1}, {0.4, 0.0, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Test----------\n";
    cout<<"\nTwo filaments, each with two beads, update link forces\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {fovx,fovy}, {nx,ny}, dt, temp, viscosity, link_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax);
	
    double f1_1_ex[2] = {0, 0}; 
    double f2_1_ex[2] = {0, 0}; 
    double f1_2_ex[2] = {0, 0}; 
    double f2_2_ex[2] = {0, 0};
    array <double,2> f1_1; 
    array <double,2> f2_1; 
    array <double,2> f1_2; 
    array <double,2> f2_2;  
    double net_sz = network->get_nfilaments(); 
    for(int f = 0; f < net_sz; f++){ network->update_link_forces(f);} 

    for(int f = 0; f < net_sz; f++){
            	
	if(f==0){ 
            f1_1 = network->get_filament(f)->get_actin(0)->get_force(); 
    	    f2_1 = network->get_filament(f)->get_actin(1)->get_force(); 
	} 
	if(f==1){
	    f1_2 = network->get_filament(f)->get_actin(0)->get_force(); 
	    f2_2 = network->get_filament(f)->get_actin(1)->get_force(); 
	}
    }
            
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
    //double  actin_length            = 0.5;
    double  link_length             = 1;
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
    
    vector<vector<double>> act = {{-0.6, 0.5, 0.25, 0}, {-0.6, -0.5, 0.25, 0}, {-0.5, 0.0, 0.25, 1}, {0.5, 0.0, 0.25, 1}};

    cout<<"\n----------Filament_Ensemble_Update_Test----------\n";
    cout<<"\nTwo filaments, each with two beads\n";

    filament_ensemble * network;

    network = new filament_ensemble(act, {fovx,fovy}, {nx,ny}, dt, temp, viscosity, link_length, stretching_stiffness, ext, bending_stiffness, frac, bc, rmax);
    double x_s[2][2]; 
    double y_s[2][2]; 
    double x_ev[2][2]={{-0.60101859, -0.60101859}, {-0.497962816728424,0.5}}; 
    double y_ev[2][2]={{0.5,-0.5}, {0,0}}; 
    network->update(); 
    for(int i = 0; i < 2; i++){
     	for(int j = 0; j < 2; j++){
	    x_s[i][j]=network->get_filament(i)->get_actin(j)->get_xcm();
            y_s[i][j]=network->get_filament(i)->get_actin(j)->get_ycm(); 

	    cout<<"Simulated x value: "<< x_s[i][j]<<endl; 
	    cout<<"Expected x value: "<< x_ev[i][j]<<endl; 
	    cout<<"Simulated y value: "<< y_s[i][j]<<endl; 
	    cout<<"Expected y value: "<< y_ev[i][j]<<endl; 

            BOOST_CHECK_CLOSE(x_s[i][j], x_ev[i][j], 1e-3); 
	    BOOST_CHECK_CLOSE(y_s[i][j], y_ev[i][j], 1e-3); 
    	}

    }
}

//
