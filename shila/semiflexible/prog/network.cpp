#include "filament_ensemble.h"
#include "motor_ensemble.h"
#include "globals.h"

#include <iostream>
#include <fstream> 
#include <iterator>
#include <array>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}


//main method
int main(int argc, char* argv[]){
    

    /***********************
     * VARIABLES           *
     **********************/
    
    int seed;
    
    double xrange, yrange;                                                  //Space
    
    int    count = 0, nframes, nmsgs;                                  //Time 
    double tinit = 0.0, t = tinit, tfinal, dt;
    
    double viscosity, temperature;                                          //Environment
    string bnd_cnd;                                         //Allowed values: NONE, PERIODIC, REFLECTIVE

    double actin_length, npolymer, nmonomer;                                // Actin 
    string actin_pos_str;
    
    double link_length, polymer_bending_modulus, link_stretching_stiffness, fracture_force, bending_fracture_force; // Links
    double bending_correction_factor = 250;
    string link_color = "1"; //"blue"
    bool use_linear_bending;

    double a_motor_length=0.5, a_motor_v=1.0, a_motor_density, a_motor_stiffness, a_m_kon, a_m_kend, a_m_koff;// Active Motors (i.e., "myosin")
            
    string a_motor_pos_str; 
    
    double p_motor_length=0.5, p_motor_density, p_motor_stiffness, // Passive Mtors (i.e., cross_linkers)
            p_motor_v=0, p_m_kon, p_m_kend, p_m_koff; 
    string p_motor_pos_str;
    
    string config_file;                                                // Input configuration
    
    string   dir,    afile,  amfile,  pmfile,  lfile;                  // Output
    ofstream o_file, file_a, file_am, file_pm, file_l;

    double shear_rate;                                                      //External Force

    //Options allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("version, v", "print version string")
        ("help", "produce help message")
        ("config,c", po::value<string>(&config_file)->default_value("config/network.cfg"), "name of a configuration file")
        ;

    //Options allowed in a config file
    po::options_description config("Configuration");
    config.add_options()
        
        ("xrange", po::value<double>(&xrange)->default_value(50), "size of cell in horizontal direction (um)")
        ("yrange", po::value<double>(&yrange)->default_value(50), "size of cell in vertical direction (um)")
        
        ("dt", po::value<double>(&dt)->default_value(0.001), "length of individual timestep in seconds")
        ("tfinal", po::value<double>(&tfinal)->default_value(10), "length of simulation in seconds")
        ("nframes", po::value<int>(&nframes)->default_value(1000), "number of timesteps between printing actin/link/motor positions to file")
        ("nmsgs", po::value<int>(&nmsgs)->default_value(10000), "number of timesteps between printing simulation progress to stdout")
       
        ("viscosity", po::value<double>(&viscosity)->default_value(0.5), "Implicity viscosity to determine friction [um^2 / s]")
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT [pN-um] that effects magnituded of Brownian component of simulation")
        ("bnd_cnd,bc", po::value<string>(&bnd_cnd)->default_value("NONE"), "boundary conditions")
        
        ("nmonomer", po::value<double>(&nmonomer)->default_value(1), "number of monomers per filament")
        ("npolymer", po::value<double>(&npolymer)->default_value(3), "number of polymers in the network")
        ("actin_length", po::value<double>(&actin_length)->default_value(10), "Length of a single actin monomer")
        ("actin_pos_str", po::value<string> (&actin_pos_str)->default_value(""), "Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions")
        
        ("a_motor_density", po::value<double>(&a_motor_density)->default_value(0.001), "number of active motors / area")
        ("p_motor_density", po::value<double>(&p_motor_density)->default_value(0.001), "number of passive motors / area")
        ("a_motor_pos_str", po::value<string> (&a_motor_pos_str)->default_value(""), "Starting positions of motors, commas delimit coordinates; semicolons delimit positions")
        
        ("a_m_kon", po::value<double>(&a_m_kon)->default_value(90.0),"active motor on rate")
        ("a_m_koff", po::value<double>(&a_m_koff)->default_value(1),"active motor off rate")
        ("a_m_kend", po::value<double>(&a_m_kend)->default_value(5),"active motor off rate at filament end")
        ("a_motor_stiffness", po::value<double>(&a_motor_stiffness)->default_value(50),"active motor spring stiffness (pN/um)")
        
        ("p_m_kon", po::value<double>(&p_m_kon)->default_value(90),"passive motor on rate")
        ("p_m_koff", po::value<double>(&p_m_koff)->default_value(0),"passive motor off rate")
        ("p_m_kend", po::value<double>(&p_m_kend)->default_value(0),"passive motor off rate at filament end")
        ("p_motor_stiffness", po::value<double>(&p_motor_stiffness)->default_value(50),"passive motor spring stiffness (pN/um)")
        
        ("link_length", po::value<double>(&link_length)->default_value(0), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.04), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(1000000), "pN-- filament breaking point")
        ("bending_fracture_force", po::value<double>(&bending_fracture_force)->default_value(1000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(10), "stiffness of link, pN/um")
        ("use_linear_bending,linear", po::value<bool>(&use_linear_bending)->default_value(true),"option to send spring type of bending springs")
        ("shear_rate", po::value<double>(&shear_rate)->default_value(0), "shear rate in pN/(um*s)")
        
        ("dir", po::value<string>(&dir)->default_value("out/test"), "output directory")
        ("seed", po::value<int>(&seed)->default_value(time(NULL)), "Random number generator seed")
        ; 
    
    //Hidden options, will be allowed both on command line and 
    //in config file, but will not be shown to user
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file", po::value< vector<string> >(), "input file")
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

    ifstream ifs(config_file.c_str());
    if (!ifs){
        cout<<"can not open config file: "<<config_file<<"\n";
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
    
    if (polymer_bending_modulus < 0){ //This is a flag for using the temperature for the bending modulus
        polymer_bending_modulus = 10*temperature; // 10um * kT
    }

    double actin_density = npolymer*nmonomer/(xrange*yrange);//0.65;
    double link_bending_stiffness    = bending_correction_factor * polymer_bending_modulus * pow(1.0/actin_length,3);
    
    int n_bw_stdout = max(int((tfinal - tinit)/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal - tinit)/(dt*double(nframes))),1);

    vector<double *> actin_position_ptrs, a_motor_position_ptrs;
    if (actin_pos_str.size() > 0)
        actin_position_ptrs = str2ptrvec(actin_pos_str, ";", ",");
    if (a_motor_pos_str.size() > 0)
        a_motor_position_ptrs = str2ptrvec(a_motor_pos_str, ";", ",");
    
    srand(seed);
            
    
    afile  = dir + "/txt_stack/rods.txt";
    lfile  = dir + "/txt_stack/links.txt";
    amfile = dir + "/txt_stack/amotors.txt";
    pmfile = dir + "/txt_stack/pmotors.txt";
    file_a.open(afile.c_str());
    file_l.open(lfile.c_str());
    file_am.open(amfile.c_str());
    file_pm.open(pmfile.c_str());
		    

    // Create Network Objects
    cout<<"Creating actin network..\n";
	
    DLfilament_ensemble * net = new DLfilament_ensemble(actin_density, xrange, yrange, xgrid, ygrid, dt, 
                                        temperature, actin_length, viscosity, nmonomer, link_length, 
                                        actin_position_ptrs, 
                                        link_stretching_stiffness, link_bending_stiffness,
                                        fracture_force, bending_fracture_force, bnd_cnd, seed); 

    cout<<"Adding active motors...\n";
    motor_ensemble<DLfilament_ensemble> * myosins = new motor_ensemble<DLfilament_ensemble>( a_motor_density, xrange, yrange, dt, temperature, 
                                             a_motor_length, net, a_motor_v, a_motor_stiffness, a_m_kon, a_m_koff,
                                             a_m_kend, actin_length, viscosity, a_motor_position_ptrs);
    cout<<"Adding passive motors (crosslinkers) ...\n";
    motor_ensemble<DLfilament_ensemble> * crosslks = new motor_ensemble<DLfilament_ensemble>( p_motor_density, xrange, yrange, dt, temperature, 
                                             p_motor_length, net, p_motor_v, p_motor_stiffness, p_m_kon, p_m_koff,
                                             p_m_kend, actin_length, viscosity, a_motor_position_ptrs);
    cout<<"Updating motors, filaments and crosslinks in the network..\n";
            
    if (shear_rate != 0)
        net->set_shear_rate(shear_rate);
    
    if (use_linear_bending){
        net->set_bending_linear();
        cout<<"\nusing linear bending\n";
    }
    string time_str = "t = 0\n";
    file_a << time_str;
    net->write_rods(file_a);
    file_l << time_str;
    net->write_links(file_l);
    file_am << time_str;
    myosins->motor_write(file_am);
    file_pm << time_str;
    crosslks->motor_write(file_pm);
   
    //Run the simulation
    while (t<=tfinal) {
        //print time count
		if (count%n_bw_stdout==0) {
			cout<<"\nTime counts: "<<count;
		}

        //update network
        if (shear_rate != 0)
            net->update_shear();

        net->update_bending();
        net->update_stretching();
        net->update(t);
        net->quad_update();
        
        //update motors and cross linkers
        crosslks->motor_walk(t);
        myosins->motor_walk(t);
        
        t+=dt;
		count++;

        //print to file
	    if (count%n_bw_print==0) {
	        
            time_str = "t = "+to_string(t)+"\n";
            file_a << time_str;
            net->write_rods(file_a);
            file_l << time_str;
            net->write_links(file_l);
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
    cout<<"\nHere's where I think I delete things\n";
    
//    delete lks;
    delete myosins;
    delete crosslks;
    delete net;
   
    int as = actin_position_ptrs.size();
    for (int i = 0; i < as; i++) delete [] actin_position_ptrs[i];
    int ms = a_motor_position_ptrs.size();
    for (int i = 0; i < ms; i++) delete [] a_motor_position_ptrs[i];
    
    // Write the output configuration file
    string output_file                         =   dir + "/data/output.txt";
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
    o_file << " Simulation time: "      << tfinal - tinit  << ", dt: " << dt <<", dt between output files: "<< n_bw_print*dt<<", Viscosity: " << viscosity              <<"\n";
    o_file << " Boundary Conditions: " <<bnd_cnd<<"\n";
    o_file.close();
    
    
    cout<<"\nTime counts: "<<count;
	cout<<"\nExecuted";
	cout<<"\n Done\n";
    
    return 0;
}
