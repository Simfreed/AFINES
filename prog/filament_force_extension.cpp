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
    int xgrid, ygrid, grid_factor; 
    
    int    count = 0, nframes, nmsgs;                                  //Time 
    double tinit = 0.0, t = tinit, tfinal, dt;
    
    double viscosity, temperature;                                          //Environment
    string bnd_cnd;                                         //Allowed values: NONE, PERIODIC, REFLECTIVE

    double actin_length, npolymer, nmonomer;                                // Actin 
    string actin_pos_str;
    
    double link_length, polymer_bending_modulus, link_stretching_stiffness, fracture_force, bending_fracture_force; // Links
    string config_file, actin_in;
    
    string   dir,    afile,  lfile, thfile;                  // Output
    ofstream o_file, file_a, file_l, file_th;
    ios_base::openmode write_mode = ios_base::out;

    bool restart_actin;
    double pull, fene_pct;

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
        
        ("xrange", po::value<double>(&xrange)->default_value(100), "size of cell in horizontal direction (um)")
        ("yrange", po::value<double>(&yrange)->default_value(100), "size of cell in vertical direction (um)")
        ("grid_factor", po::value<int>(&grid_factor)->default_value(2), "number of grid boxes per um^2")
        
        ("dt", po::value<double>(&dt)->default_value(0.0001), "length of individual timestep in seconds")
        ("tfinal", po::value<double>(&tfinal)->default_value(10), "length of simulation in seconds")
        ("nframes", po::value<int>(&nframes)->default_value(1000), "number of timesteps between printing actin/link/motor positions to file")
        ("nmsgs", po::value<int>(&nmsgs)->default_value(10000), "number of timesteps between printing simulation progress to stdout")
       
        ("viscosity", po::value<double>(&viscosity)->default_value(0.001), "Dynamic viscosity to determine friction [mg / (um*s)]. At 20 C, is 0.001 for water")
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT [pN-um] that effects magnituded of Brownian component of simulation")
        ("bnd_cnd,bc", po::value<string>(&bnd_cnd)->default_value("PERIODIC"), "boundary conditions")
        ("pull", po::value<double>(&pull)->default_value(0.004), "amount of force to pull on ends of filament")
        ("fene_pct", po::value<double>(&fene_pct)->default_value(0.5), "pct of rest length of filament to allow outstretched until fene blowup")
        
        ("nmonomer", po::value<double>(&nmonomer)->default_value(51), "number of monomers per filament")
        ("npolymer", po::value<double>(&npolymer)->default_value(1), "number of polymers in the network")
        ("actin_length", po::value<double>(&actin_length)->default_value(0.5), "Length of a single actin monomer")
        ("actin_pos_str", po::value<string> (&actin_pos_str)->default_value(""), "Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions")
        
        ("link_length", po::value<double>(&link_length)->default_value(1), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.04), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(100000000), "pN-- filament breaking point")
        ("bending_fracture_force", po::value<double>(&bending_fracture_force)->default_value(1000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(10), "stiffness of link, pN/um")//probably should be about 70000 to correspond to actin
        
        ("actin_in", po::value<string>(&actin_in)->default_value(""), "input actin positions file")
        ("restart_actin", po::value<bool>(&restart_actin)->default_value(false), "if true, input actin positions file chosen by default")
        
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
   
    
    if (polymer_bending_modulus < 0){ //This is a flag for using the temperature for the bending modulus
        polymer_bending_modulus = 10*temperature; // 10um * kT
    }

    double actin_density = npolymer*nmonomer/(xrange*yrange);//0.65;
    double link_bending_stiffness    = polymer_bending_modulus/link_length;
    
    int n_bw_stdout = max(int((tfinal - tinit)/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal - tinit)/(dt*double(nframes))),1);

    // To Read positions from input strings in config file
    vector<array<double,3> > actin_position_arrs;
    if (actin_pos_str.size() > 0)
        actin_position_arrs   = str2arrvec(actin_pos_str, ":", ",");
    else 
        for(int i = 0; i < npolymer; i++)
            actin_position_arrs.push_back({0,0,0});

    // To Read positions from input files
    vector<vector<double> > actin_pos_vec;
    
    if (restart_actin)
        actin_in = dir + "/restart/in/actins.txt";

    if (actin_in.size() > 0)
        actin_pos_vec   = file2vecvec(actin_in, "\t");
    
    srand(seed);
            

    afile  = dir + "/txt_stack/actins.txt";
    lfile  = dir + "/txt_stack/links.txt";
    thfile = dir + "/data/thermo.txt";

    if (restart_actin) write_mode = ios_base::app;

    file_a.open(afile.c_str(), write_mode);
    file_l.open(lfile.c_str(), write_mode);
	file_th.open(thfile.c_str(), write_mode);


    // DERIVED QUANTITIES :
    xgrid = 0;
    ygrid = 0;

    // Create Network Objects
    cout<<"\nCreating actin network..";
    ATfilament_ensemble * net;
    if (actin_pos_vec.size() == 0){
        net = new ATfilament_ensemble(actin_density, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, actin_length, viscosity, nmonomer, link_length, 
                actin_position_arrs, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd, seed); 
    }else{
        net = new ATfilament_ensemble(actin_pos_vec, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, viscosity, link_length, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd); 
    }
   
    string time_str = "t = 0";


    file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
    net->write_actins(file_a);
    file_l << time_str<<"\tN = "<<to_string(net->get_nlinks());
    net->write_links(file_l);
    net->write_thermo(file_th);

    // Write the output configuration file
    string output_file                         =   dir + "/data/output.txt";
    o_file.open(output_file.c_str());
    o_file << " FILE: "                 << output_file     <<"\n";
    o_file << " Actin Density: "        << actin_density   << ", Actin Mean Length: "          << actin_length              << "\n";
    o_file << " Link Rest Length: "     << link_length     << ", Link Stretching Stiffness: "  << link_stretching_stiffness <<", Link Bending Stiffness: " << link_bending_stiffness <<"\n";
    o_file << " Simulation time: "      << tfinal - tinit  << ", dt: " << dt <<", dt between output files: "<< n_bw_print*dt<<", Viscosity: " << viscosity              <<"\n";
    o_file << " Boundary Conditions: " <<bnd_cnd<<"\n";
    o_file.close();
    
    /*vector<double> pulls;
    for (int i = -7; i <= -2; i++) 
        for (int j = 1; j <= 10; j++)
            pulls.push_back(j * pow(10, i));
    */
    while (t < tfinal) {

        if (count%n_bw_stdout==0) {
			cout<<"\nTime counts: "<<count;
            net->print_network_thermo();
            net->print_filament_lengths();
        }

        for (int i = 0; i < net->get_nfilaments(); i++){
            //net->get_filament(i)->pull_on_ends(pull);
            net->get_filament(i)->affine_pull(pull);
        }

        net->update();
        //clear the vector of fractured filaments
        //net->clear_broken();

        t+=dt;
		count++;

        //print to file
	    if (count%n_bw_print==0) {
	        
            time_str = "\nt = "+to_string(t);
            
            file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
            net->write_actins(file_a);
            
            file_l << time_str<<"\tN = "<<to_string(net->get_nlinks());
            net->write_links(file_l);
            
            file_th << time_str<<"\tN = "<<to_string(net->get_nlinks());
            net->write_thermo(file_th);
		}
        
    }

    file_a << "\n";
    file_l << "\n";
    file_th << "\n";

    file_a.close();
    file_l.close();
    file_th.close(); 
    //Delete all objects created
    cout<<"\nHere's where I think I delete things\n";
    
    delete net;
    
    
    
    cout<<"\nTime counts: "<<count;
	cout<<"\nExecuted";
	cout<<"\n Done\n";
    
    return 0;
}
