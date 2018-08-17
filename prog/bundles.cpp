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
    
    int    count, nframes, nmsgs, quad_update_period;
    double t, tinit, tfinal, dt;
    
    double viscosity, temperature;                                          //Environment
    string bnd_cnd;                                         //Allowed values: NONE, PERIODIC, REFLECTIVE

    double actin_length, npolymer, nmonomer;                                // Actin 
    string actin_pos_str;
    
    double link_length, polymer_bending_modulus, link_stretching_stiffness, fene_pct, fracture_force; // Links

    double spacer1_length, spacer1_v=0, spacer1_density, spacer1_stiffness, s1_kon, s1_kend, s1_koff, s1_koff2,
           s1_stall, s1_cut, s1_bend, s1_ang;// Active Motors (i.e., "myosin")
    string spacer1_pos_str; 
    
    double spacer2_length, spacer2_density, spacer2_stiffness, // Passive Mtors (i.e., cross_linkers)
            spacer2_v=0, s2_kon, s2_kend, s2_koff, s2_koff2, s2_stall, s2_cut, s2_bend, s2_ang; 
    string spacer2_pos_str;
    
    string config_file, actin_in, spacer1_in, spacer2_in;                                                // Input configuration
    
    string   dir, tdir, ddir,  afile,  amfile,  pmfile,  lfile, thfile, pefile;                  // Output
    ofstream o_file, file_a, file_am, file_pm, file_l, file_th, file_pe;
    ios_base::openmode write_mode = ios_base::out;

    bool s2_intersect_flag, s1_intersect_flag, s1_dead_head_flag, s2_dead_head_flag, static_cl_flag, quad_off_flag, write_unbound_flag;
    double s2_linkage_prob, s1_linkage_prob;                                              
    int s1_dead_head, s2_dead_head;

    double rmax; 
    double kexv; 

    bool restart;
    double restart_time;
    
    double kgrow, lgrow, l0min, l0max;
    int nlink_max;

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
        ("grid_factor", po::value<double>(&grid_factor)->default_value(2), "number of grid boxes per um^2")
        
        ("dt", po::value<double>(&dt)->default_value(0.0001), "length of individual timestep in seconds")
        ("tinit", po::value<double>(&tinit)->default_value(0), "time that recording of simulation starts")
        ("tfinal", po::value<double>(&tfinal)->default_value(0.01), "length of simulation in seconds")
        ("nframes", po::value<int>(&nframes)->default_value(1000), "number of times between actin/link/motor positions to are printed to file")
        ("nmsgs", po::value<int>(&nmsgs)->default_value(10000), "number of times simulation progress is printed to stdout")
       
        ("viscosity", po::value<double>(&viscosity)->default_value(0.001), "Dynamic viscosity to determine friction [mg / (um*s)]. At 20 C, is 0.001 for water")
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT [pN-um] that effects magnituded of Brownian component of simulation")
        ("bnd_cnd,bc", po::value<string>(&bnd_cnd)->default_value("PERIODIC"), "boundary conditions")
        
        ("nmonomer", po::value<double>(&nmonomer)->default_value(11), "number of monomers per filament")
        ("npolymer", po::value<double>(&npolymer)->default_value(3), "number of polymers in the network")
        ("actin_length", po::value<double>(&actin_length)->default_value(0.5), "Length of a single actin monomer")
        ("actin_pos_str", po::value<string> (&actin_pos_str)->default_value(""), "Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions")
        
        ("spacer1_density", po::value<double>(&spacer1_density)->default_value(0.05), "number of xlink 1s / um^2")
        ("spacer2_density", po::value<double>(&spacer2_density)->default_value(0.05), "number of xlink 2s / um^2")
        ("spacer1_pos_str", po::value<string> (&spacer1_pos_str)->default_value(""), "Starting positions of motors, commas delimit coordinates; semicolons delimit positions")
        ("spacer2_pos_str", po::value<string> (&spacer2_pos_str)->default_value(""), "Starting positions of crosslinks, commas delimit coordinates; semicolons delimit positions")
        
        ("s1_kon", po::value<double>(&s1_kon)->default_value(2),"xlink 1 on rate")
        ("s1_koff", po::value<double>(&s1_koff)->default_value(0.1),"xlink 1 off rate")
        ("s1_kend", po::value<double>(&s1_kend)->default_value(0.1),"xlink 1 off rate at filament end")
        ("s1_koff2", po::value<double>(&s1_koff2)->default_value(0.1),"xlink 1 off rate when both heads bound")
        ("spacer1_length", po::value<double>(&spacer1_length)->default_value(0.15),"filamin rest length (um)")
        ("spacer1_stiffness", po::value<double>(&spacer1_stiffness)->default_value(1),"xlink 1 spring stiffness (pN/um)")
        ("spacer1_v", po::value<double>(&spacer1_v)->default_value(0),"xlink 1 velocity (um/s)")
        
        ("s1_stall", po::value<double>(&s1_stall)->default_value(0.5),"force beyond which xlinks don't walk (pN)")
        ("s1_cut", po::value<double>(&s1_cut)->default_value(0.063),"cutoff distance for binding (um)")
        
        ("s1_bend", po::value<double>(&s1_bend)->default_value(0.04),"bending force constant of spacer  (pN) (10kT/um by default)")
        ("s1_ang", po::value<double>(&s1_ang)->default_value(pi/2),"equilibrium angle of spacer1-filament system")
        
        ("s2_kon", po::value<double>(&s2_kon)->default_value(2),"xlink 2 on rate")
        ("s2_koff", po::value<double>(&s2_koff)->default_value(0.1),"xlink 2 off rate")
        ("s2_kend", po::value<double>(&s2_kend)->default_value(0.1),"xlink 2 off rate at filament end")
        ("s2_koff2", po::value<double>(&s2_koff2)->default_value(0.1),"xlink 2 off rate when both heads bound")
        ("spacer2_length", po::value<double>(&spacer2_length)->default_value(0.03),"xlink 2 rest length (um) (default: filamin)")
        ("spacer2_stiffness", po::value<double>(&spacer2_stiffness)->default_value(1),"xlink 2 spring stiffness (pN/um)")
       
        ("s2_stall", po::value<double>(&s2_stall)->default_value(0.5),"force beyond which xlinks don't walk (pN)")
        ("s2_cut", po::value<double>(&s2_cut)->default_value(0.063),"cutoff distance for binding (um)")
        
        ("s2_bend", po::value<double>(&s2_bend)->default_value(0.04),"bending force constant of spacer 2  (pN) (10kT/um by default)")
        ("s2_ang", po::value<double>(&s2_ang)->default_value(pi/2),"equilibrium angle of spacer2-filament system")

        ("link_length", po::value<double>(&link_length)->default_value(1), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.068), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(100000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(1), "stiffness of link, pN/um")//probably should be about 70000 to correspond to actin
        ("fene_pct", po::value<double>(&fene_pct)->default_value(0.5), "pct of rest length of filament to allow outstretched until fene blowup")
        
        ("actin_in", po::value<string>(&actin_in)->default_value(""), "input actin positions file")
        ("spacer1_in", po::value<string>(&spacer1_in)->default_value(""), "input motor positions file")
        ("spacer2_in", po::value<string>(&spacer2_in)->default_value(""), "input crosslinker positions file")
        
        ("restart", po::value<bool>(&restart)->default_value(false), "if true, will restart simulation from last timestep recorded")
        ("restart_time", po::value<double>(&restart_time)->default_value(-1), "time to restart simulation from")
        
        ("dir", po::value<string>(&dir)->default_value("."), "output directory")
        ("myseed", po::value<int>(&myseed)->default_value(time(NULL)), "Random number generator myseed")
        
        ("s2_intersect_flag", po::value<bool>(&s2_intersect_flag)->default_value(false), "flag to put spacer 2 at all filament intersections")
        ("s1_intersect_flag", po::value<bool>(&s1_intersect_flag)->default_value(false), "flag to put spacer 1 at all filament intersections")
        ("s2_linkage_prob", po::value<double>(&s2_linkage_prob)->default_value(1), "If s2_intersect_flag, probability that two filaments that intersect will be linked")
        ("s1_linkage_prob", po::value<double>(&s1_linkage_prob)->default_value(1), "If s1_intersect_flag, probability that two filaments that intersect will be motor-d")

        ("s1_dead_head_flag", po::value<bool>(&s1_dead_head_flag)->default_value(false), "flag to kill head <s1_dead_head> of all motors")
        ("s1_dead_head", po::value<int>(&s1_dead_head)->default_value(0), "index of head to kill")
        
        ("s2_dead_head_flag", po::value<bool>(&s2_dead_head_flag)->default_value(false), "flag to kill head <s1_dead_head> of all crosslinks")
        ("s2_dead_head", po::value<int>(&s2_dead_head)->default_value(0), "index of head to kill")
      	("rmax", po::value<double>(&rmax)->default_value(0.25), "cutoff distance for interactions between actins beads and filaments")

        ("kexv", po::value<double>(&kexv)->default_value(1.0), "parameter of exv force calculation") 
        
        ("static_cl_flag", po::value<bool>(&static_cl_flag)->default_value(false), "flag to indicate compeletely static xlinks; i.e, no walking, no detachment")
        ("quad_off_flag", po::value<bool>(&quad_off_flag)->default_value(false), "flag to turn off neighbor list updating")
        ("quad_update_period", po::value<int>(&quad_update_period)->default_value(1), "number of timesteps between actin/link/motor position updates to update quadrants")
        
        ("write_unbound_flag", po::value<bool>(&write_unbound_flag)->default_value(false), "if true, all spacer positions written, if false only doubly bound written")
        
        //Options for filament growth
        ("kgrow", po::value<double>(&kgrow)->default_value(0), "rate of filament growth")
        ("lgrow", po::value<double>(&lgrow)->default_value(0), "additional length of filament upon growth")
        ("l0min", po::value<double>(&l0min)->default_value(0), "minimum length a link can shrink to before disappearing")
        ("l0max", po::value<double>(&l0max)->default_value(0), "maximum length a link can grow to before breaking into two links")
        ("nlink_max", po::value<int>(&nlink_max)->default_value(25), "maximum number of links allowed on filament")
        
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
    cout<<"\nDEBUG: actin_density = "<<actin_density; 
    double link_bending_stiffness    = polymer_bending_modulus / link_length;
    
    int n_bw_stdout = max(int((tfinal)/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal)/(dt*double(nframes))),1);
    int unprinted_count = int(double(tinit)/dt);

    tdir   = dir  + "/txt_stack";
    ddir   = dir  + "/data";
    fs::path dir1(tdir.c_str()), dir2(ddir.c_str());
   

    string fname_add = write_unbound_flag? "" : "_bound";

    afile  = tdir + "/actins.txt";
    lfile  = tdir + "/links.txt";
    amfile = tdir + "/spacers1"+fname_add+".txt";
    pmfile = tdir + "/spacers2"+fname_add+".txt";
    thfile = ddir + "/filament_e.txt";
    pefile = ddir + "/pe.txt";
    
    if(fs::create_directory(dir1)) cerr<< "Directory Created: "<<afile<<std::endl;
    if(fs::create_directory(dir2)) cerr<< "Directory Created: "<<thfile<<std::endl;
    
    // To Read positions from input strings in config file
    vector<array<double,3> > actin_position_arrs, spacer1_position_arrs, spacer2_position_arrs;
    if (actin_pos_str.size() > 0)
        actin_position_arrs   = str2arrvec(actin_pos_str, ":", ",");
    if (spacer1_pos_str.size() > 0)
        spacer1_position_arrs = str2arrvec(spacer1_pos_str, ":", ",");
    if (spacer2_pos_str.size() > 0)
        spacer2_position_arrs = str2arrvec(spacer2_pos_str, ":", ",");
   

    // To Read positions from input files
    vector<vector<double> > actin_pos_vec;
    vector<vector<double> > spacer1_pos_vec, spacer2_pos_vec;
   
    if (actin_in.size() > 0)
        actin_pos_vec   = file2vecvec(actin_in, "\t");
    if (spacer1_in.size() > 0)
        spacer1_pos_vec = file2vecvec(spacer1_in, "\t");
    if (spacer2_in.size() > 0)
        spacer2_pos_vec = file2vecvec(spacer2_in, "\t");
    
    // To restart a whole trajectory from it's last full timestep : 
    if (restart){
        
        double tf_prev  = min(last_full_timestep(afile), last_full_timestep(lfile));
        if (spacer1_density > 0)
            tf_prev = min(tf_prev, last_full_timestep(amfile));
        if (spacer2_density > 0)
            tf_prev = min(tf_prev, last_full_timestep(pmfile));

        if (restart_time == -1 || restart_time > tf_prev)
            restart_time = tf_prev;

        cout<<"\nRestarting from t = "<<restart_time<<endl;

        double nprinted = restart_time / (dt*n_bw_print);

        actin_pos_vec   = traj2vecvec(afile, "\t ", restart_time);
        spacer1_pos_vec = traj2vecvec(amfile, "\t ", restart_time);
        spacer2_pos_vec = traj2vecvec(pmfile, "\t ", restart_time);
        
        // for actins, links, amotors, pmotors: 
        // do: 
        //   copy whole file into temp
        //   while hasn't reached tf in temp file:
        //      write from copy into afile
        write_first_tsteps(afile,  restart_time);
        write_first_tsteps(lfile,  restart_time);
        write_first_tsteps(amfile, restart_time);
        write_first_tsteps(pmfile, restart_time);
        
        write_first_tsteps(thfile, restart_time);
        write_first_nlines( pefile, (int) nprinted);

    
        tinit       = restart_time;
        write_mode  = ios_base::app;
    }
    
    
    set_seed(myseed);
    
    file_a.open(afile.c_str(), write_mode);
    file_l.open(lfile.c_str(), write_mode);
    file_am.open(amfile.c_str(), write_mode);
    file_pm.open(pmfile.c_str(), write_mode);
	file_th.open(thfile.c_str(), write_mode);
	file_pe.open(pefile.c_str(), write_mode);


    // DERIVED QUANTITIES :
    if(spacer1_density == 0 && spacer1_pos_vec.size() == 0 &&
            spacer2_density==0 && spacer2_pos_vec.size() == 0 &&
            !s2_intersect_flag && !s1_intersect_flag){
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
        net = new filament_ensemble(actin_density, {{xrange, yrange}}, {{xgrid, ygrid}}, dt, 
                temperature, actin_length, viscosity, nmonomer, link_length, 
                actin_position_arrs, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd, myseed, rmax, kexv); 
    }else{
        net = new filament_ensemble(actin_pos_vec, {{xrange, yrange}}, {{xgrid, ygrid}}, dt, 
                temperature, viscosity, link_length, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd, rmax, kexv); 
    }
  
    net->set_growing(kgrow, lgrow, l0min, l0max, nlink_max);

    if (s2_intersect_flag) spacer2_pos_vec = net->spring_spring_intersections(spacer2_length, s2_linkage_prob); 
    if (s1_intersect_flag) spacer1_pos_vec = net->spring_spring_intersections(spacer1_length, s1_linkage_prob); 
    if (quad_off_flag) net->turn_quads_off();

    cout<<"\nAdding xlink 1s...";
    spacer_ensemble * big_xlinks;
    
    if (spacer1_pos_vec.size() == 0)
        big_xlinks = new spacer_ensemble( spacer1_density, {xrange, yrange}, dt, temperature, 
                spacer1_length, net, spacer1_v, spacer1_stiffness, fene_pct, s1_kon, s1_koff,
                s1_kend, s1_stall, s1_cut, viscosity, spacer1_position_arrs, bnd_cnd);
    else
        big_xlinks = new spacer_ensemble( spacer1_pos_vec, {xrange, yrange}, dt, temperature, 
                spacer1_length, net, spacer1_v, spacer1_stiffness, fene_pct, s1_kon, s1_koff,
                s1_kend, s1_stall, s1_cut, viscosity, bnd_cnd);
    big_xlinks->set_bending(s1_bend, s1_ang);
    big_xlinks->set_binding_two(s1_kon, s1_koff2, s1_koff2);

    if (s1_dead_head_flag) big_xlinks->kill_heads(s1_dead_head);

    cout<<"Adding xlink 2s (crosslinkers) ...\n";
    spacer_ensemble * small_xlinks; 
    
    if(spacer2_pos_vec.size() == 0)
        small_xlinks = new spacer_ensemble( spacer2_density, {xrange, yrange}, dt, temperature, 
                spacer2_length, net, spacer2_v, spacer2_stiffness, fene_pct, s2_kon, s2_kend,
                s2_kend, s2_stall, s2_cut, viscosity, spacer2_position_arrs, bnd_cnd);
    else
        small_xlinks = new spacer_ensemble( spacer2_pos_vec, {xrange, yrange}, dt, temperature, 
                spacer2_length, net, spacer2_v, spacer2_stiffness, fene_pct, s2_kon, s2_kend,
                s2_kend, s2_stall, s2_cut, viscosity, bnd_cnd);
    small_xlinks->set_bending(s2_bend, s2_ang);
    small_xlinks->set_binding_two(s2_kon, s2_koff2, s2_koff2);
    if (s2_dead_head_flag) small_xlinks->kill_heads(s2_dead_head);

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
    cout<<"\nUpdating motors, filaments and crosslinks in the network..";
    string time_str; 
    count=0;
    t = tinit;

    //Run the simulation
    while (t < tfinal) {
        
        //print to file
	    if (t+dt/100 >= tinit && (count-unprinted_count)%n_bw_print==0) {
	        
            if (t>tinit) time_str ="\n";
            time_str += "t = "+to_string(t);
            
//            file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
//            net->write_actins(file_a);
            
            file_l << time_str<<"\tN = "<<to_string(net->get_nsprings());
            net->write_springs(file_l);
            
            file_am << time_str<<"\tN = "<<to_string(big_xlinks->get_nmotors());
            
            if (write_unbound_flag)
                big_xlinks->motor_write(file_am);
            else
                big_xlinks->motor_write_doubly_bound(file_am);
            
            file_pm << time_str<<"\tN = "<<to_string(small_xlinks->get_nmotors());
            
            if (write_unbound_flag)
                small_xlinks->motor_write(file_pm);
            else
                small_xlinks->motor_write_doubly_bound(file_pm);
            
            file_th << time_str<<"\tN = "<<to_string(net->get_nsprings());
            net->write_thermo(file_th);

            file_pe << net->get_stretching_energy()<<"\t"<<net->get_bending_energy()<<"\t"<<
                big_xlinks->get_potential_energy()<<"\t"<<small_xlinks->get_potential_energy()<<endl;
            
            file_a<<std::flush;
            file_l<<std::flush;
            file_am<<std::flush;
            file_pm<<std::flush;
            file_th<<std::flush;
            file_pe<<std::flush;

		}
        
        //print time count

        if (count%n_bw_stdout==0) {
			cout<<"\nTime counts: "<<count;
		    //net->print_filament_thermo();
            net->print_network_thermo();
            small_xlinks->print_ensemble_thermo();
            big_xlinks->print_ensemble_thermo();
        }

        //update network
        net->update();//updates all forces, velocities and positions of filaments

        if ( ! quad_off_flag && count % quad_update_period == 0)
            net->quad_update_serial();
        //update cross linkers
        if (static_cl_flag)
            small_xlinks->motor_update();
        else
            small_xlinks->motor_walk(t);
       
        //update motors
        big_xlinks->motor_walk(t);
        
        //clear the vector of fractured filaments
        net->clear_broken();

        
        t+=dt;
		count++;

        
    }

    file_a << "\n";
    file_l << "\n";
    file_am << "\n";
    file_pm << "\n";
    file_th << "\n";

    file_a.close();
    file_l.close();
    file_am.close();
    file_pm.close();
    file_th.close(); 
    file_pe.close(); 
    //Delete all objects created
    cout<<"\nHere's where I think I delete things\n";
    
//    delete lks;
    delete big_xlinks;
    delete small_xlinks;
    delete net;
    
    
    
    cout<<"\nTime counts: "<<count;
	cout<<"\nExecuted";
	cout<<"\n Done\n";
    
    return 0;
}
