#include "filament_ensemble.h"
#include "motor_ensemble.h"
#include "globals.h"

#include <iostream>
#include <fstream> 
#include <iterator>
#include <array>
#include <boost/program_options.hpp>
#include <boost/any.hpp>
#include <typeinfo>

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
    
    int myseed;
    
    double xrange, yrange;                                                  //Space
    int xgrid, ygrid;
    double grid_factor; 
    
    int    count, nframes, nmsgs;                                  //Time 
    double t, tinit, tfinal, dt;
    
    double viscosity, temperature;                                          //Environment
    string bnd_cnd;                                         //Allowed values: NONE, PERIODIC, REFLECTIVE

    int npolymer, nmonomer, nmonomer_extra;                                                               // Actin 
    double actin_length, extra_bead_prob;                            
    string actin_pos_str;
    
    double link_length, polymer_bending_modulus, link_stretching_stiffness, fene_pct, fracture_force; // Links

    double a_motor_length, a_motor_v, a_motor_density, a_motor_stiffness, a_m_kon, a_m_kend, a_m_koff,
           a_m_stall, a_m_cut;// Active Motors (i.e., "myosin")
    string a_motor_pos_str; 
    
    string config_file, actin_in, a_motor_in;                                                // Input configuration
    
    string   dir, tdir, ddir,  amfile, pefile;                  // Output
    ofstream o_file, file_am, file_pe ;
    ios_base::openmode write_mode = ios_base::out;

    bool motor_intersect_flag, quad_off_flag, dead_head_flag;
    double a_linkage_prob;                                              
    int dead_head;

    bool restart;
    double restart_time;

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
        
        ("xrange", po::value<double>(&xrange)->default_value(10), "size of cell in horizontal direction (um)")
        ("yrange", po::value<double>(&yrange)->default_value(10), "size of cell in vertical direction (um)")
        ("grid_factor", po::value<double>(&grid_factor)->default_value(2.5), "number of grid boxes per um^2")
        
        ("dt", po::value<double>(&dt)->default_value(0.0001), "length of individual timestep in seconds")
        ("tinit", po::value<double>(&tinit)->default_value(0), "time that recording of simulation starts")
        ("tfinal", po::value<double>(&tfinal)->default_value(0.01), "length of simulation in seconds")
        ("nframes", po::value<int>(&nframes)->default_value(1000), "number of times between actin/link/motor positions to are printed to file")
        ("nmsgs", po::value<int>(&nmsgs)->default_value(10000), "number of times simulation progress is printed to stdout")
       
        ("viscosity", po::value<double>(&viscosity)->default_value(0.001), "Dynamic viscosity to determine friction [mg / (um*s)]. At 20 C, is 0.001 for water")
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT [pN-um] that effects magnituded of Brownian component of simulation")
        ("bnd_cnd,bc", po::value<string>(&bnd_cnd)->default_value("PERIODIC"), "boundary conditions")
        
        ("npolymer", po::value<int>(&npolymer)->default_value(3), "number of polymers in the network")
        ("nmonomer", po::value<int>(&nmonomer)->default_value(11), "number of monomers per filament (if extra_bead_prob != 0, then minimum #)")
        ("nmonomer_extra", po::value<int>(&nmonomer_extra)->default_value(0), "max # of monomers per filament")
        ("extra_bead_prob", po::value<double>(&extra_bead_prob)->default_value(0.5), "probability of adding an extra bead from nmonomer...nmonomer_extra")
        
        ("actin_length", po::value<double>(&actin_length)->default_value(0.5), "Length of a single actin monomer")
        ("actin_pos_str", po::value<string> (&actin_pos_str)->default_value(""), "Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions")
        
        ("a_motor_density", po::value<double>(&a_motor_density)->default_value(0.05), "number of active motors / um^2")
        ("a_motor_pos_str", po::value<string> (&a_motor_pos_str)->default_value(""), "Starting positions of motors, commas delimit coordinates; semicolons delimit positions")
        
        ("a_m_kon", po::value<double>(&a_m_kon)->default_value(1),"active motor on rate")
        ("a_m_koff", po::value<double>(&a_m_koff)->default_value(0.1),"active motor off rate")
        ("a_m_kend", po::value<double>(&a_m_kend)->default_value(0.1),"active motor off rate at filament end")
        ("a_motor_length", po::value<double>(&a_motor_length)->default_value(0.5),"active motor rest length (um)")
        ("a_motor_stiffness", po::value<double>(&a_motor_stiffness)->default_value(1),"active motor spring stiffness (pN/um)")
        ("a_motor_v", po::value<double>(&a_motor_v)->default_value(1),"active motor velocity (um/s)")
        
        ("a_m_stall", po::value<double>(&a_m_stall)->default_value(0.5),"force beyond which motors don't walk (pN)")
        ("a_m_cut", po::value<double>(&a_m_cut)->default_value(0.063),"cutoff distance for binding (um)")
        
        ("link_length", po::value<double>(&link_length)->default_value(1), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.068), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(100000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(1), "stiffness of link, pN/um")//probably should be about 70000 to correspond to actin
        ("fene_pct", po::value<double>(&fene_pct)->default_value(0.5), "pct of rest length of filament to allow outstretched until fene blowup")
        
        ("actin_in", po::value<string>(&actin_in)->default_value(""), "input actin positions file")
        ("a_motor_in", po::value<string>(&a_motor_in)->default_value(""), "input motor positions file")
        
        ("restart", po::value<bool>(&restart)->default_value(false), "if true, will restart simulation from last timestep recorded")
        ("restart_time", po::value<double>(&restart_time)->default_value(-1), "time to restart simulation from")

        ("dir", po::value<string>(&dir)->default_value("."), "output directory")
        ("myseed", po::value<int>(&myseed)->default_value(time(NULL)), "Random number generator myseed")
        
        ("motor_intersect_flag", po::value<bool>(&motor_intersect_flag)->default_value(false), "flag to put a motor at all filament intersections")
        ("a_linkage_prob", po::value<double>(&a_linkage_prob)->default_value(1), "If motor_intersect_flag, probability that two filaments that intersect will be motor-d")

        ("dead_head_flag", po::value<bool>(&dead_head_flag)->default_value(false), "flag to kill head <dead_head> of all motors")
        ("dead_head", po::value<int>(&dead_head)->default_value(0), "index of head to kill")
        
        ("quad_off_flag", po::value<bool>(&quad_off_flag)->default_value(false), "flag to turn off neighbor list updating")
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

    //double actin_density = double(npolymer*nmonomer)/(xrange*yrange);//0.65;
    //cout<<"\nDEBUG: actin_density = "<<actin_density; 
    double link_bending_stiffness    = polymer_bending_modulus / link_length;
    
    int n_bw_stdout = max(int((tfinal)/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal)/(dt*double(nframes))),1);
    int unprinted_count = int(double(tinit)/dt);

    tdir   = dir  + "/txt_stack";
    ddir   = dir  + "/data";
    fs::path dir1(tdir.c_str()), dir2(ddir.c_str());
    
    amfile = tdir + "/amotors.txt";
    pefile = ddir + "/pe.txt";
    
    if(fs::create_directory(dir1)) cerr<< "Directory Created: "<<amfile<<std::endl;
    if(fs::create_directory(dir2)) cerr<< "Directory Created: "<<pefile<<std::endl;
    
    // To Read positions from input strings in config file
    vector<array<double,3> > actin_position_arrs, a_motor_position_arrs;
    if (actin_pos_str.size() > 0)
        actin_position_arrs   = str2arrvec(actin_pos_str, ":", ",");
    if (a_motor_pos_str.size() > 0)
        a_motor_position_arrs = str2arrvec(a_motor_pos_str, ":", ",");
   
    // To get positions from input files: 
    vector<vector<double> > actin_pos_vec;
    vector<vector<double> > a_motor_pos_vec;
   
    if (actin_in.size() > 0)
        actin_pos_vec   = file2vecvec(actin_in, "\t");
    if (a_motor_in.size() > 0)
        a_motor_pos_vec = file2vecvec(a_motor_in, "\t");
    
    // To restart a whole trajectory from it's last full timestep : 
    if (restart){
        
        double tf_prev  = last_full_timestep(amfile);

        if (restart_time == -1 || restart_time > tf_prev)
            restart_time = tf_prev;

        cout<<"\nRestarting from t = "<<restart_time<<endl;

        double nprinted = restart_time / (dt*n_bw_print);
        a_motor_pos_vec = traj2vecvec(amfile, "\t ", restart_time);
        
        //   copy whole file into temp
        //   while hasn't reached tf in temp file:
        //      write from copy into afile
        write_first_tsteps(amfile, restart_time);
        write_first_nlines( pefile, (int) nprinted);
        
        tinit       = restart_time;
        write_mode  = ios_base::app;
    }
    
    
    set_seed(myseed);
    

    
    //const char* path = _filePath.c_str();
    
    file_am.open(amfile.c_str(), write_mode);
	file_pe.open(pefile.c_str(), write_mode);

    // DERIVED QUANTITIES :
    if(a_motor_density == 0 && a_motor_pos_vec.size() == 0 && !motor_intersect_flag){
        xgrid = 1;
        ygrid = 1;
    }
    else{
        xgrid  = (int) round(grid_factor*xrange);
        ygrid  = (int) round(grid_factor*yrange);
    }

    // Create Network Objects
    cout<<"\nCreating actin network..";
    filament_ensemble * net;
    if (actin_pos_vec.size() == 0 && actin_in.size() == 0){
        net = new filament_ensemble(npolymer, nmonomer, nmonomer_extra, extra_bead_prob, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, actin_length, viscosity, link_length, 
                actin_position_arrs, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd, myseed); 
    }else{
        net = new filament_ensemble(actin_pos_vec, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, viscosity, link_length, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd); 
    }
   
    if (motor_intersect_flag) a_motor_pos_vec = net->link_link_intersections(a_motor_length, a_linkage_prob); 
    if (quad_off_flag) 
        net->turn_quads_off();
    else
        net->quad_update_serial();
    cout<<"\nAdding active motors...";
    motor_ensemble * myosins;
    
    if (a_motor_pos_vec.size() == 0 && a_motor_in.size() == 0)
        myosins = new motor_ensemble( a_motor_density, {xrange, yrange}, dt, temperature, 
                a_motor_length, net, a_motor_v, a_motor_stiffness, fene_pct, a_m_kon, a_m_koff,
                a_m_kend, a_m_stall, a_m_cut, viscosity, a_motor_position_arrs, bnd_cnd);
    else
        myosins = new motor_ensemble( a_motor_pos_vec, {xrange, yrange}, dt, temperature, 
                a_motor_length, net, a_motor_v, a_motor_stiffness, fene_pct, a_m_kon, a_m_koff,
                a_m_kend, a_m_stall, a_m_cut, viscosity, bnd_cnd);
    if (dead_head_flag) myosins->kill_heads(dead_head);

    // Write the full configuration file
    string output_file  =   dir + "/data/config_full.cfg";
    o_file.open(output_file.c_str());
    
    boost::any val;
    for(po::variables_map::const_iterator it=vm.begin(); it!=vm.end(); ++it){
        
        if (it->first == "config") continue;

        val=it->second.value();
        
        if(typeid(bool) == val.type())
            o_file << it->first <<"="<< boost::any_cast<bool>(val) <<endl;
        else if(typeid(int) == val.type())
            o_file << it->first <<"="<< boost::any_cast<int>(val) <<endl;
        else if(typeid(double) == val.type())
            o_file << it->first <<"="<< boost::any_cast<double>(val) <<endl;
        else if(typeid(string) == val.type())
            o_file << it->first <<"="<< boost::any_cast<string>(val) <<endl;
    }
   
    // Run the simulation
    cout<<"\nUpdating motors in the network..";
    string time_str; 
    count=0;
    t = tinit;

    while (t <= tfinal) {
        
        //print to file
	    if (t+dt/100 >= tinit && (count-unprinted_count)%n_bw_print==0) {
	        
            if (t>tinit) time_str ="\n";
            time_str += "t = "+to_string(t);
            
            file_am << time_str<<"\tN = "<<to_string(myosins->get_nmotors());
            myosins->motor_write(file_am);
            
            file_pe << net->get_stretching_energy()<<"\t"<<net->get_bending_energy()<<"\t"<<
                myosins->get_potential_energy()<<endl;
            file_am<<std::flush;
            file_pe<<std::flush;
		}
        
        if (count%n_bw_stdout==0) {
			cout<<"\nTime counts: "<<count;
            myosins->print_ensemble_thermo();
        }

        //update motors
        myosins->motor_walk(t);
        
        t+=dt;
		count++;

    }

    file_am << "\n";

    file_am.close();
    file_pe.close(); 
    
    //Delete all objects created
    cout<<"\nHere's where I think I delete things\n";
    
//    delete lks;
    delete myosins;
    
    cout<<"\nTime counts: "<<count;
	cout<<"\nExecuted";
	cout<<"\n Done\n";
    
    return 0;
}
