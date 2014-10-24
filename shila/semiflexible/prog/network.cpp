#include "actin_ensemble.h"
#include "link_ensemble.h"
#include "motor_ensemble.h"
#include "globals.h"

#include <iostream>
#include <fstream> 
#include <iterator>
#include <array>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// #define dt 0.0001 -- defined in globals.h

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}


//main method
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
    double dt;
    int print_dt        = 1000;  // # of timesteps to print to file
    int stdout_dt       = 10000; // # of timesteps to print to stdout
    double t            = tinit;
    
    // Environment
    double viscosity           = 0.5; //um^2/s
    double temperature         = 0.004; //pN-um
    
    // Actin 
    double actin_length = 1; //(um) length of a ROD
    double npolymer     = 100; 
    double nmonomer = 10;

    // Links 
    double link_stretching_stiffness = 50.0; //pn / um
    std::string link_color           = "1"; //"blue";
    double polymer_bending_modulus   = 0.04; //using kT * Lp for bending modulus, assuming Lp = 10 um and kT = 0.004 um-pN
    
    // Motors
    double motor_length=0.5;
    double motor_density=0.1;
    double motor_stiffness=50.0; //pN / um
    double vmotor=1.0;
    double m_kon=90.0;          
    double m_kend=5.0;
    double m_koff=1.0; 
   
    std::string actin_pos_str, motor_pos_str;

    // Input
    std::string config_file;

    // Output
    std::string dir, afile, mfile, lfile;
    std::ofstream a_final, m_final, o_file;
	char numstr[21];
   
    /***********************
     * VARIABLES           *
     **********************/
    
    //Options allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("version, v", "print version std::string")
        ("help", "produce help message")
        ("config,c", po::value<std::string>(&config_file)->default_value("config/network.cgf"), "name of a configuration file")
        ;

    //Options allowed in a config file
    po::options_description config("Configuration");
    config.add_options()
        
        ("xrange", po::value<double>(&xrange)->default_value(50), "size of cell in horizontal direction (um)")
        ("yrange", po::value<double>(&yrange)->default_value(50), "size of cell in vertical direction (um)")
        
        ("dt", po::value<double>(&dt)->default_value(0.0001), "length of individual timestep in seconds")
        ("tfinal", po::value<double>(&tfinal)->default_value(100), "length of simulation in seconds")
        ("print_dt", po::value<int>(&print_dt)->default_value(1000), "number of timesteps between printing actin/link/motor positions to file")
        ("stdout_dt", po::value<int>(&stdout_dt)->default_value(10000), "number of timesteps between printing simulation progress to stdout")
        
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT that effects magnituded of Brownian component of simulation")
        
        ("nmonomer", po::value<double>(&nmonomer)->default_value(10), "number of monomers per filament")
        ("npolymer", po::value<double>(&npolymer)->default_value(100), "number of polymers in the network")
        ("actin_length", po::value<double>(&actin_length)->default_value(1), "Length of a single actin monomer")
        ("actin_pos_str", po::value<std::string> (&actin_pos_str), "Starting positions of actin polymers, commas delimit coordinates; spaces delimit positions")
        
        ("motor_density", po::value<double>(&motor_density)->default_value(0), "number of motors / area")
        ("motor_pos_str", po::value<std::string> (&motor_pos_str), "Starting positions of motors, commas delimit coordinates; spaces delimit positions")
        
        ("polymer_bending_modulus,kapb", po::value<double>(&polymer_bending_modulus)->default_value(0.04), "Bending modulus of a filament")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(100), "stiffness of link, pN/um")
        ("dir", po::value<std::string>(&dir)->default_value("out/test"), "output directory")
        ; 
    
    //Hidden options, will be allowed both on command line and 
    //in config file, but will not be shown to user
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file", po::value< std::vector<std::string> >(), "input file")
        ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);
    
    po::positional_options_description p;
    p.add("input-file", -1); ///wha in the world is this doing

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
    notify(vm);

    std::ifstream ifs(config_file.c_str());
    if (!ifs){
        std::cout<<"can not open config file: "<<config_file<<"\n";
        return 0;
    }
    else
    {
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
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
    o_file << " Motor unloaded speed: " << vmotor          << ", Motor binding rate: "         << m_kon                     <<"\n";
    o_file << " Motor unbinding rate: " << m_koff          << ", Motor end detachment rate: "  << m_kend                    <<", Viscosity: "              << viscosity              <<"\n";
    o_file << " Link Rest Length: "     << link_length     << ", Link Stretching Stiffness: "  << link_stretching_stiffness <<", Link Bending Stiffness: " << link_bending_stiffness <<"\n";
    o_file << " Simulation time: "      << tfinal - tinit  << ", dt: " << dt <<", Number of time steps between output files: "<< print_dt<<"\n";
    o_file.close();
    
    
    std::vector<double *> actin_position_ptrs = str2ptrvec(actin_pos_str, ";", ",");
    std::vector<double *> motor_position_ptrs = str2ptrvec(motor_pos_str, ";", ",");

    std::cout<<"Creating actin network..\n";
	actin_ensemble * net = new actin_ensemble(actin_density, xrange, yrange, xgrid, ygrid, dt, 
                                        temperature, actin_length, viscosity, nmonomer, link_length, 
                                        actin_position_ptrs, seed);
    std::cout<<"Creating link ensemble...\n";
    link_ensemble * lks = new link_ensemble();
    std::cout<<"Adding links to connect actin filament monomers...\n";
    net->connect_polymers( lks, link_length, link_stretching_stiffness, link_bending_stiffness, link_color );
    std::cout<<"Adding motors...\n";
    motor_ensemble * myosins = new motor_ensemble( motor_density, xrange, yrange, dt, temperature, 
                                             motor_length, net, vmotor, motor_stiffness, m_kon, m_koff,
                                             m_kend, actin_length, viscosity, motor_position_ptrs);
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
    
    int as = actin_position_ptrs.size(), ms = motor_position_ptrs.size();
    for (int i = 0; i < as; i++) delete [] actin_position_ptrs[i];
    for (int i = 0; i < ms; i++) delete [] motor_position_ptrs[i];

    a_final.close();
    m_final.close();
    
    std::cout<<"\nTime counts: "<<count;
	std::cout<<"\nExecuted";
	std::cout<<"\n Done\n";
    return 0;
}
