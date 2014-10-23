#include "actin_ensemble.h"
#include "link_ensemble.h"
#include "motor_ensemble.h"
#include "globals.h"

// #define dt 0.0001 -- defined in globals.h

int main(int argc, char* argv[]){
    
    int seed=time(NULL);
    srand(seed);

/************************
 *      CONTROLS        *
 ***********************/
 
    // Space
    double xrange = 50.0;
    double yrange = 50.0;

    // Time 
    int count           = 0;
    double tinit        = 0.0;
    double tfinal       = 200;
    int print_dt     = 1;  // # of timesteps to print to file
    int stdout_dt    = 1; // # of timesteps to print to stdout
    double t            = tinit;
    
    // Environment
    double viscosity    = 0.5; //units?

    // Actin 
    double actin_length = 1; //(um) length of a monomer
    double npolymer     = 2; 

    // Links 
    double link_stretching_stiffness = 50.0; //pn / um
    double polymer_bending_modulus   = temperature * 10; //using kT * Lp for bending modulus, with Lp = 10 um
    std::string link_color           = "1"; //"blue";
    
    // Motors
    double motor_length=1.5;
    double motor_density=0.1;
    double motor_stiffness=50.0; //pN / um
    double vmotor=1.0;
    double m_kon=90.0;          
    double m_kend=5.0;
    double m_koff=1.0; 
   

    // Output
    std::string dir, afile, mfile, lfile;
    std::ofstream a_final, m_final, o_file;
	char numstr[21];
   
    std::vector<double *> actin_positions, motor_positions;
    double apos1[3] = {-5, 0.75, 0};
    double apos2[3] = {5, -0.75, pi};
    double mpos1[3] = {0, 0, pi/2};
    actin_positions.push_back(apos1);
    actin_positions.push_back(apos2);
    motor_positions.push_back(mpos1);
    /***********************
     * VARIABLES           *
     **********************/
    double nmonomer = 50;
//  double link_bending_stiffness = motor_stiffness/10;
    
    if (argc>1) {
        nmonomer                    =   atof(argv[1]);
        npolymer                    =   atof(argv[2]);
        actin_length                =   atof(argv[3]);
        motor_density               =   atof(argv[4]);
        tfinal                      =   atof(argv[5]);
        link_stretching_stiffness   =   atof(argv[6]);
        xrange                      =   atof(argv[7]);
        yrange                      =   atof(argv[8]);
        dir                         =        argv[9] ;
    }
    
    // DERIVED QUANTITIES :
    double xgrid  = 2*xrange;
    double ygrid  = 2*yrange;
    double actin_density = npolymer*nmonomer/(xrange*yrange);//0.65;
    double link_length               = actin_length/10; 
    double link_bending_stiffness    = polymer_bending_modulus * pow(1.0/actin_length,3);
    
    std::string output_file                         =   dir + "/data/output.txt";
    std::string actin_output                        =   dir + "/data/actin_final.txt";
    std::string myosin_output                       =   dir + "/data/myosin_final.txt";
    
    a_final.open(actin_output.c_str());
    m_final.open(myosin_output.c_str());
    o_file.open(output_file.c_str());
    
    o_file << " FILE: "                 << output_file     <<"\n";
    o_file << " Actin Density: "        << actin_density   << ", Actin Mean Length: "          << actin_length              << "\n";
    o_file << " Motor Density: "        << motor_density   << ", Motor Rest Length: "          << motor_length              << ", Motor Stiffness: "       << motor_stiffness        <<"\n";

    o_file << " Link Rest Length: "     << link_length     << ", Link Stretching Stiffness: "  << link_stretching_stiffness <<", Link Bending Stiffness: " << link_bending_stiffness <<"\n";
    o_file << " Simulation time: "      << tfinal - tinit  << ", dt: " << dt <<", Number of time steps between output files: "<< print_dt<<"\n";
    o_file.close();
    
    
    std::cout<<"Creating actin network..\n";
    actin_ensemble * net = new
        actin_ensemble(actin_density,xrange,yrange,xgrid,ygrid,actin_length,viscosity,nmonomer,link_length, actin_positions);
    std::cout<<"Creating link ensemble...\n";
    link_ensemble * lks = new link_ensemble();
    std::cout<<"Adding links to connect actin filament monomers...\n";
    net->connect_polymers( lks, link_length, link_stretching_stiffness, link_bending_stiffness, link_color );
    std::cout<<"Adding motors...\n";
    motor_ensemble * myosins = new motor_ensemble( motor_density, xrange, yrange, motor_length, 
                                             net, vmotor, motor_stiffness, m_kon, m_koff,
                                             m_kend, actin_length, viscosity, motor_positions);
    std::cout<<"Updating motors, filaments and crosslinks in the network..\n";
    
    while (t<=tfinal) {
        //print time count
		if (count%stdout_dt==0) {
			std::cout<<"Time counts: "<<count<<"\n";
		}

        net->update_bending();
        net->update();
        net->quad_update();

        //print to file
	    if (count%print_dt==0) {
	        
            sprintf(numstr, "%d", count/print_dt);
            
            afile = dir + "/txt_stack/afile" + numstr + ".txt";
            mfile = dir + "/txt_stack/mfile" + numstr + ".txt";
            lfile = dir + "/txt_stack/lfile" + numstr + ".txt";
			
            std::ofstream file_a, file_m, file_l;
			file_a.open(afile.c_str());
			file_m.open(mfile.c_str());
            file_l.open(lfile.c_str());
		    
            net->write(file_a);
            myosins->motor_write(file_m);
            lks->link_write(file_l);
			file_a.close();
			file_m.close();
            file_l.close();
            
		}
        //update network
        myosins->motor_walk();
        lks->link_walk(); 
        
        t+=dt;
		count++;
    }
    net->write(a_final);
    myosins->motor_write(m_final);
    
    //Delete all objects created
    std::cout<<"Here's where I think I delete things\n";
    
    delete lks;
    delete myosins;
    delete net;

    a_final.close();
    m_final.close();
    
    std::cout<<"\nTime counts: "<<count;
	std::cout<<"\nExecuted";
	std::cout<<"\n Done\n";
    return 0;
}
