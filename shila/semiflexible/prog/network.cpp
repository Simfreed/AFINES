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

    /***********************
     * VARIABLES           *
     **********************/
    
    double xrange, yrange;                                                  //Space
    
    int    count = 0, print_dt, stdout_dt;                                  //Time 
    double tinit = 0.0, t = tinit, tfinal, dt;
    
    double viscosity, temperature;                                          //Environment
    
    double actin_length, npolymer, nmonomer;                                // Actin 
    std::string actin_pos_str;
    
    double link_length, polymer_bending_modulus, link_stretching_stiffness; // Links
    std::string link_color = "1"; //"blue"
    
    double a_motor_length=0.5, a_motor_density=0.1, a_motor_stiffness=50.0, // Active Motors (i.e., "myosin")
           a_motor_v=1.0, a_m_kon=90.0, a_m_kend=5.0, a_m_koff=1.0; 
    std::string a_motor_pos_str; 
    
    double p_motor_length=0.5, p_motor_density=0.1, p_motor_stiffness=50.0, // Passive Mtors (i.e., cross_linkers)
            p_motor_v=0, p_m_kon=90.0, p_m_kend=5.0, p_m_koff=1.0; 
    std::string p_motor_pos_str;
    
    std::string config_file;                                                // Input configuration
    
    std::string   dir,    afile,  amfile,  pmfile,  lfile;                  // Output
    std::ofstream o_file, file_a, file_am, file_pm, file_l;
   
    //Options allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("version, v", "print version std::string")
        ("help", "produce help message")
        ("config,c", po::value<std::string>(&config_file)->default_value("config/network.cfg"), "name of a configuration file")
        ;

    //Options allowed in a config file
    po::options_description config("Configuration");
    config.add_options()
        
        ("xrange", po::value<double>(&xrange)->default_value(50), "size of cell in horizontal direction (um)")
        ("yrange", po::value<double>(&yrange)->default_value(50), "size of cell in vertical direction (um)")
        
        ("dt", po::value<double>(&dt)->default_value(0.0001), "length of individual timestep in seconds")
        ("tfinal", po::value<double>(&tfinal)->default_value(1), "length of simulation in seconds")
        ("print_dt", po::value<int>(&print_dt)->default_value(10), "number of timesteps between printing actin/link/motor positions to file")
        ("stdout_dt", po::value<int>(&stdout_dt)->default_value(100), "number of timesteps between printing simulation progress to stdout")
       
        ("viscosity", po::value<double>(&viscosity)->default_value(0.5), "Implicity viscosity to determine friction [um^2 / s]")
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT [pN-um] that effects magnituded of Brownian component of simulation")
        
        ("nmonomer", po::value<double>(&nmonomer)->default_value(10), "number of monomers per filament")
        ("npolymer", po::value<double>(&npolymer)->default_value(3), "number of polymers in the network")
        ("actin_length", po::value<double>(&actin_length)->default_value(1), "Length of a single actin monomer")
        ("actin_pos_str", po::value<std::string> (&actin_pos_str), "Starting positions of actin polymers, commas delimit coordinates; spaces delimit positions")
        
        ("a_motor_density", po::value<double>(&a_motor_density)->default_value(0.001), "number of active motors / area")
        ("p_motor_density", po::value<double>(&p_motor_density)->default_value(0.001), "number of passive motors / area")
        ("a_motor_pos_str", po::value<std::string> (&a_motor_pos_str), "Starting positions of motors, commas delimit coordinates; spaces delimit positions")
        
        ("link_length", po::value<double>(&link_length)->default_value(actin_length/10), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.04), "Bending modulus of a filament")
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
    double link_bending_stiffness    = polymer_bending_modulus * pow(1.0/actin_length,3);
    
    std::vector<double *> actin_position_ptrs = str2ptrvec(actin_pos_str, ";", ",");
    std::vector<double *> a_motor_position_ptrs = str2ptrvec(a_motor_pos_str, ";", ",");
            
    
    afile  = dir + "/txt_stack/rods.txt";
    lfile  = dir + "/txt_stack/links.txt";
    amfile = dir + "/txt_stack/amotors.txt";
    pmfile = dir + "/txt_stack/pmotors.txt";
    file_a.open(afile.c_str());
    file_l.open(lfile.c_str());
    file_am.open(amfile.c_str());
    file_pm.open(pmfile.c_str());
		    

    // Create Network Objects
    std::cout<<"Creating actin network..\n";
	actin_ensemble * net = new actin_ensemble(actin_density, xrange, yrange, xgrid, ygrid, dt, 
                                        temperature, actin_length, viscosity, nmonomer, link_length, 
                                        actin_position_ptrs, seed);
    std::cout<<"Creating link ensemble...\n";
    link_ensemble * lks = new link_ensemble();
    std::cout<<"Adding links to connect actin filament monomers...\n";
    net->connect_polymers( lks, link_length, link_stretching_stiffness, link_bending_stiffness, link_color );
    std::cout<<"Adding active motors...\n";
    motor_ensemble * myosins = new motor_ensemble( a_motor_density, xrange, yrange, dt, temperature, 
                                             a_motor_length, net, a_motor_v, a_motor_stiffness, a_m_kon, a_m_koff,
                                             a_m_kend, actin_length, viscosity, a_motor_position_ptrs);
    std::cout<<"Adding passive motors (crosslinkers) ...\n";
    motor_ensemble * crosslks = new motor_ensemble( p_motor_density, xrange, yrange, dt, temperature, 
                                             p_motor_length, net, p_motor_v, p_motor_stiffness, p_m_kon, p_m_koff,
                                             p_m_kend, actin_length, viscosity, a_motor_position_ptrs);
    std::cout<<"Updating motors, filaments and crosslinks in the network..\n";
            
    std::string time_str = "t = 0\n";
    file_a << time_str;
    net->write(file_a);
    file_l << time_str;
    lks->link_write(file_l);
    file_am << time_str;
    myosins->motor_write(file_am);
    file_pm << time_str;
    crosslks->motor_write(file_pm);
   
    //Run the simulation
    while (t<=tfinal) {
        //print time count
		if (count%stdout_dt==0) {
			std::cout<<"Time counts: "<<count<<"\n";
		}

        net->update_bending();
        net->update();
        net->quad_update();
        
        //update network
        lks->link_walk(); 
        crosslks->motor_walk();
        myosins->motor_walk();
        
        t+=dt;
		count++;

        //print to file
	    if (count%print_dt==0) {
	        
            time_str = "t = "+std::to_string(t)+"\n";
            file_a << time_str;
            net->write(file_a);
            file_l << time_str;
            lks->link_write(file_l);
            file_am << time_str;
            myosins->motor_write(file_am);
            file_pm << time_str;
            crosslks->motor_write(file_pm);
            
		}
        
    }
    
    file_a.close();
    file_l.close();
    file_am.close();
    file_pm.close();
    
    //Delete all objects created
    std::cout<<"Here's where I think I delete things\n";
    
    delete lks;
    delete myosins;
    delete net;
    
    int as = actin_position_ptrs.size(), ms = a_motor_position_ptrs.size();
    for (int i = 0; i < as; i++) delete [] actin_position_ptrs[i];
    for (int i = 0; i < ms; i++) delete [] a_motor_position_ptrs[i];
    
    // Write the output configuration file
    std::string output_file                         =   dir + "/data/output.txt";
    o_file.open(output_file.c_str());
    o_file << " FILE: "                 << output_file     <<"\n";
    o_file << " Actin Density: "        << actin_density   << ", Actin Mean Length: "          << actin_length              << "\n";
    o_file << " Active Motor Density: "        << a_motor_density   << ", Active Motor Rest Length: "          << a_motor_length              << ", Active Motor Stiffness: "       << a_motor_stiffness        <<"\n";
    o_file << " Active Motor unloaded speed: " << a_motor_v          << ", Active Motor binding rate: "         << a_m_kon                     <<"\n";
    o_file << " Active Motor unbinding rate: " << a_m_koff          << ", Active Motor end detachment rate: "  << a_m_kend                    <<"\n";
    o_file << " Passive Motor Density: "        << p_motor_density   << ", Passive Motor Rest Length: "          << p_motor_length              << ", Passive Motor Stiffness: "       << p_motor_stiffness        <<"\n";
    o_file << " Passive Motor unloaded speed: " << p_motor_v          << ", Passive Motor binding rate: "         << p_m_kon                     <<"\n";
    o_file << " Passive Motor unbinding rate: " << p_m_koff          << ", Passive Motor end detachment rate: "  << p_m_kend                    <<"\n";
    o_file << " Link Rest Length: "     << link_length     << ", Link Stretching Stiffness: "  << link_stretching_stiffness <<", Link Bending Stiffness: " << link_bending_stiffness <<"\n";
    o_file << " Simulation time: "      << tfinal - tinit  << ", dt: " << dt <<", Number of time steps between output files: "<< print_dt<<", Viscosity: " << viscosity              <<"\n";
    o_file.close();
    
    
    std::cout<<"\nTime counts: "<<count;
	std::cout<<"\nExecuted";
	std::cout<<"\n Done\n";
    
    return 0;
}
