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
    
    int    count, nframes, nmsgs, n_bw_shear, quad_update_period;            //Time 
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
    
    double p_motor_length, p_motor_density, p_motor_stiffness, // Passive Mtors (i.e., cross_linkers)
            p_motor_v, p_m_kon, p_m_kend, p_m_koff, p_m_stall, p_m_cut;
    string p_motor_pos_str;
    
    string config_file, actin_in, a_motor_in, p_motor_in;                                                // Input configuration
    
    string   dir, tdir, ddir,  afile,  amfile,  pmfile,  lfile, thfile, pefile, kefile;                  // Output
    ofstream o_file, file_a, file_am, file_pm, file_l, file_th, file_pe, file_ke;
    ios_base::openmode write_mode = ios_base::out;

    double strain_pct, time_of_strain, pre_strain, d_strain_pct, d_strain_amp;                                      //External Force
    double d_strain, prev_d_strain, d_strain_freq, time_of_dstrain;
    bool link_intersect_flag, motor_intersect_flag, dead_head_flag, p_dead_head_flag, static_cl_flag, quad_off_flag;
    bool diff_strain_flag, osc_strain_flag;
    double p_linkage_prob, a_linkage_prob;                                              
    int dead_head, p_dead_head;
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
        ("n_bw_shear", po::value<int>(&n_bw_shear)->default_value(1000000000), "number of timesteps between subsequent shears")
       
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
        ("p_motor_density", po::value<double>(&p_motor_density)->default_value(0.05), "number of passive motors / um^2")
        ("a_motor_pos_str", po::value<string> (&a_motor_pos_str)->default_value(""), "Starting positions of motors, commas delimit coordinates; semicolons delimit positions")
        ("p_motor_pos_str", po::value<string> (&p_motor_pos_str)->default_value(""), "Starting positions of crosslinks, commas delimit coordinates; semicolons delimit positions")
        
        ("a_m_kon", po::value<double>(&a_m_kon)->default_value(1),"active motor on rate")
        ("a_m_koff", po::value<double>(&a_m_koff)->default_value(0.1),"active motor off rate")
        ("a_m_kend", po::value<double>(&a_m_kend)->default_value(0.1),"active motor off rate at filament end")
        ("a_motor_length", po::value<double>(&a_motor_length)->default_value(0.4),"active motor rest length (um)")
        ("a_motor_stiffness", po::value<double>(&a_motor_stiffness)->default_value(1),"active motor spring stiffness (pN/um)")
        ("a_motor_v", po::value<double>(&a_motor_v)->default_value(1),"active motor velocity (um/s)")
        
        ("a_m_stall", po::value<double>(&a_m_stall)->default_value(0.5),"force beyond which motors don't walk (pN)")
        ("a_m_cut", po::value<double>(&a_m_cut)->default_value(0.063),"cutoff distance for binding (um)")
        
        ("p_m_kon", po::value<double>(&p_m_kon)->default_value(1),"passive motor on rate")
        ("p_m_koff", po::value<double>(&p_m_koff)->default_value(0.1),"passive motor off rate")
        ("p_m_kend", po::value<double>(&p_m_kend)->default_value(0.1),"passive motor off rate at filament end")
        ("p_motor_length", po::value<double>(&p_motor_length)->default_value(0.150),"passive motor rest length (um) (default: filamin)")
        ("p_motor_stiffness", po::value<double>(&p_motor_stiffness)->default_value(1),"passive motor spring stiffness (pN/um)")
        ("p_motor_v", po::value<double>(&p_motor_v)->default_value(0),"passive motor velocity (um/s)")
       
        ("p_m_stall", po::value<double>(&p_m_stall)->default_value(0),"force beyond which xlinks don't walk (pN)")
        ("p_m_cut", po::value<double>(&p_m_cut)->default_value(0.063),"cutoff distance for binding (um)")


        ("link_length", po::value<double>(&link_length)->default_value(1), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.068), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(100000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(1), "stiffness of link, pN/um")
        ("fene_pct", po::value<double>(&fene_pct)->default_value(0.5), "pct of rest length of filament to allow outstretched until fene blowup")
        
        ("strain_pct", po::value<double>(&strain_pct)->default_value(0), "pct that the boundarys get sheared")
        ("time_of_strain", po::value<double>(&time_of_strain)->default_value(0), "time at which the step strain occurs")

        ("d_strain_freq", po::value<double>(&d_strain_freq)->default_value(1), "differential strain frequency")
        ("d_strain_pct", po::value<double>(&d_strain_pct)->default_value(0), "differential strain amplitude")
        ("time_of_dstrain", po::value<double>(&time_of_dstrain)->default_value(10000), "time when differential strain starts")
        
        ("actin_in", po::value<string>(&actin_in)->default_value(""), "input actin positions file")
        ("a_motor_in", po::value<string>(&a_motor_in)->default_value(""), "input motor positions file")
        ("p_motor_in", po::value<string>(&p_motor_in)->default_value(""), "input crosslinker positions file")
        
        ("restart", po::value<bool>(&restart)->default_value(false), "if true, will restart simulation from last timestep recorded")
        ("restart_time", po::value<double>(&restart_time)->default_value(-1), "time to restart simulation from")

        ("dir", po::value<string>(&dir)->default_value("."), "output directory")
        ("myseed", po::value<int>(&myseed)->default_value(time(NULL)), "Random number generator myseed")
        
        ("link_intersect_flag", po::value<bool>(&link_intersect_flag)->default_value(false), "flag to put a cross link at all filament intersections")
        ("motor_intersect_flag", po::value<bool>(&motor_intersect_flag)->default_value(false), "flag to put a motor at all filament intersections")
        ("p_linkage_prob", po::value<double>(&p_linkage_prob)->default_value(1), "If link_intersect_flag, probability that two filaments that intersect will be linked")
        ("a_linkage_prob", po::value<double>(&a_linkage_prob)->default_value(1), "If motor_intersect_flag, probability that two filaments that intersect will be motor-d")

        ("dead_head_flag", po::value<bool>(&dead_head_flag)->default_value(false), "flag to kill head <dead_head> of all motors")
        ("dead_head", po::value<int>(&dead_head)->default_value(0), "index of head to kill")
        
        ("p_dead_head_flag", po::value<bool>(&p_dead_head_flag)->default_value(false), "flag to kill head <dead_head> of all crosslinks")
        ("p_dead_head", po::value<int>(&p_dead_head)->default_value(0), "index of head to kill")

      	("rmax", po::value<double>(&rmax)->default_value(0.25), "cutoff distance for interactions between actins beads and filaments")

	("kexv", po::value<double>(&kexv)->default_value(1.0), "parameter of exv force calculation") 
        
        ("static_cl_flag", po::value<bool>(&static_cl_flag)->default_value(false), "flag to indicate compeletely static xlinks; i.e, no walking, no detachment")
        ("quad_off_flag", po::value<bool>(&quad_off_flag)->default_value(false), "flag to turn off neighbor list updating")
        ("quad_update_period", po::value<int>(&quad_update_period)->default_value(1), "number of timesteps between actin/link/motor position updates to update quadrants")
        
        ("diff_strain_flag", po::value<bool>(&diff_strain_flag)->default_value(false), "flag to turn on linear differential strain")
        ("osc_strain_flag", po::value<bool>(&osc_strain_flag)->default_value(false), "flag to turn on oscillatory differential strain")
        
        //Options for filament growth
        ("kgrow", po::value<double>(&kgrow)->default_value(0), "rate of filament growth")
        ("lgrow", po::value<double>(&lgrow)->default_value(0), "additional length of filament upon growth")
        ("l0min", po::value<double>(&l0min)->default_value(0), "minimum length a link can shrink to before disappearing")
        ("l0max", po::value<double>(&l0max)->default_value(2), "maximum length a link can grow to before breaking into two links")
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

    double actin_density = double(npolymer*nmonomer)/(xrange*yrange);//0.65;
    cout<<"\nDEBUG: actin_density = "<<actin_density; 
    double link_bending_stiffness    = polymer_bending_modulus / link_length;
    
    int n_bw_stdout = max(int((tfinal)/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal)/(dt*double(nframes))),1);
    int unprinted_count = int(double(tinit)/dt);

    tdir   = dir  + "/txt_stack";
    ddir   = dir  + "/data";
    fs::path dir1(tdir.c_str()), dir2(ddir.c_str());
    
    afile  = tdir + "/actins.txt";
    lfile  = tdir + "/links.txt";
    amfile = tdir + "/amotors.txt";
    pmfile = tdir + "/pmotors.txt";
    thfile = ddir + "/filament_e.txt";
    pefile = ddir + "/pe.txt"; 
    kefile = ddir + "/ke.txt"; 
    
    if(fs::create_directory(dir1)) cerr<< "Directory Created: "<<afile<<std::endl;
    if(fs::create_directory(dir2)) cerr<< "Directory Created: "<<thfile<<std::endl;
    
    // To Read positions from input strings in config file
    vector<array<double,3> > actin_position_arrs, a_motor_position_arrs, p_motor_position_arrs;
    if (actin_pos_str.size() > 0)
        actin_position_arrs   = str2arrvec(actin_pos_str, ":", ",");
    if (a_motor_pos_str.size() > 0)
        a_motor_position_arrs = str2arrvec(a_motor_pos_str, ":", ",");
    if (p_motor_pos_str.size() > 0)
        p_motor_position_arrs = str2arrvec(p_motor_pos_str, ":", ",");
   

    // To get positions from input files: 
    vector<vector<double> > actin_pos_vec;
    vector<vector<double> > a_motor_pos_vec, p_motor_pos_vec;
   
    if (actin_in.size() > 0)
        actin_pos_vec   = file2vecvec(actin_in, "\t");
    if (a_motor_in.size() > 0)
        a_motor_pos_vec = file2vecvec(a_motor_in, "\t");
    if (p_motor_in.size() > 0)
        p_motor_pos_vec = file2vecvec(p_motor_in, "\t");
    
    // To restart a whole trajectory from it's last full timestep : 
    if (restart){
        
        double tf_prev  = min(last_full_timestep(afile), last_full_timestep(lfile));
        if (a_motor_density > 0)
            tf_prev = min(tf_prev, last_full_timestep(amfile));
        if (p_motor_density > 0)
            tf_prev = min(tf_prev, last_full_timestep(pmfile));

        if (restart_time == -1 || restart_time > tf_prev)
            restart_time = tf_prev;

        cout<<"\nRestarting from t = "<<restart_time<<endl;

        double nprinted = restart_time / (dt*n_bw_print);

        actin_pos_vec   = traj2vecvec(afile, "\t ", restart_time);
        a_motor_pos_vec = traj2vecvec(amfile, "\t ", restart_time);
        p_motor_pos_vec = traj2vecvec(pmfile, "\t ", restart_time);
        
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
    	write_first_nlines( kefile, (int) nprinted);
    
        tinit       = restart_time;
        write_mode  = ios_base::app;
    }
    
    
    set_seed(myseed);
    

    
    //const char* path = _filePath.c_str();
    
    file_a.open(afile.c_str(), write_mode);
    file_l.open(lfile.c_str(), write_mode);
    file_am.open(amfile.c_str(), write_mode);
    file_pm.open(pmfile.c_str(), write_mode);
	file_th.open(thfile.c_str(), write_mode);
	file_pe.open(pefile.c_str(), write_mode); 
	file_ke.open(kefile.c_str(), write_mode);


    // DERIVED QUANTITIES :
    if(a_motor_density == 0 && a_motor_pos_vec.size() == 0 &&
            p_motor_density==0 && p_motor_pos_vec.size() == 0 &&
            !link_intersect_flag && !motor_intersect_flag){
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
    if (link_intersect_flag) p_motor_pos_vec = net->spring_spring_intersections(p_motor_length, p_linkage_prob); 
    if (motor_intersect_flag) a_motor_pos_vec = net->spring_spring_intersections(a_motor_length, a_linkage_prob); 
    
    if (quad_off_flag) net->turn_quads_off();

    cout<<"\nAdding active motors...";
    xlink_ensemble * myosins;
    
    if (a_motor_pos_vec.size() == 0 && a_motor_in.size() == 0)
        myosins = new xlink_ensemble( a_motor_density, {{xrange, yrange}}, dt, temperature, 
                a_motor_length, net, a_motor_v, a_motor_stiffness, fene_pct, a_m_kon, a_m_koff,
                a_m_kend, a_m_stall, a_m_cut, viscosity, a_motor_position_arrs, bnd_cnd);
    else
        myosins = new xlink_ensemble( a_motor_pos_vec, {{xrange, yrange}}, dt, temperature, 
                a_motor_length, net, a_motor_v, a_motor_stiffness, fene_pct, a_m_kon, a_m_koff,
                a_m_kend, a_m_stall, a_m_cut, viscosity, bnd_cnd);
    if (dead_head_flag) myosins->kill_heads(dead_head);

    cout<<"Adding passive motors (crosslinkers) ...\n";
    xlink_ensemble * crosslks; 
    
    if(p_motor_pos_vec.size() == 0 && p_motor_in.size() == 0)
        crosslks = new xlink_ensemble( p_motor_density, {{xrange, yrange}}, dt, temperature, 
                p_motor_length, net, p_motor_v, p_motor_stiffness, fene_pct, p_m_kon, p_m_koff,
                p_m_kend, p_m_stall, p_m_cut, viscosity, p_motor_position_arrs, bnd_cnd);
    else
        crosslks = new xlink_ensemble( p_motor_pos_vec, {{xrange, yrange}}, dt, temperature, 
                p_motor_length, net, p_motor_v, p_motor_stiffness, fene_pct, p_m_kon, p_m_koff,
                p_m_kend, p_m_stall, p_m_cut, viscosity, bnd_cnd);
    if (p_dead_head_flag) crosslks->kill_heads(p_dead_head);

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
    pre_strain = strain_pct * xrange;
    d_strain_amp = d_strain_pct * xrange;
    prev_d_strain = 0;

    cout<<"\nDEBUG: time of pre_strain = "<<time_of_strain;
    //Run the simulation
    if (time_of_strain == 0){
        cout<<"\nDEBUG: t = "<<t<<"; adding pre_strain of "<<pre_strain<<" um here";
        net->update_delrx( pre_strain );
        net->update_shear();
    }
    while (t <= tfinal) {
        
        //print to file
	    if (t+dt/100 >= tinit && (count-unprinted_count)%n_bw_print==0) {
	        
            if (t>tinit) time_str ="\n";
            time_str += "t = "+to_string(t);
            
            file_a << time_str<<"\tN = "<<to_string(net->get_nbeads());
            net->write_beads(file_a);

            file_l << time_str<<"\tN = "<<to_string(net->get_nsprings());
            net->write_springs(file_l);

            file_am << time_str<<"\tN = "<<to_string(myosins->get_nmotors());
            myosins->motor_write(file_am);
            
            file_pm << time_str<<"\tN = "<<to_string(crosslks->get_nmotors());
            crosslks->motor_write(file_pm);

            file_th << time_str<<"\tN = "<<to_string(net->get_nfilaments());
            net->write_thermo(file_th);

            file_pe <<net->get_stretching_energy()<<"\t"<<net->get_bending_energy()<<"\t"<<net->get_exv_energy()<<"\t"<<myosins->get_potential_energy()<<"\t"<<crosslks->get_potential_energy()<<endl;

	    file_ke <<net->get_kinetic_energy_vir()<<"\t"<<myosins->get_kinetic_energy()<<"\t"<<crosslks->get_kinetic_energy()<<endl;  
            
            file_a<<std::flush;
            file_l<<std::flush;
            file_am<<std::flush;
            file_pm<<std::flush;
            file_th<<std::flush;
            file_pe<<std::flush;
            file_ke<<std::flush;
            
		}
        
        //print time count
        if (time_of_strain!=0 && close(t, time_of_strain, dt/(10*time_of_strain))){
            //Perform the shear here
            cout<<"\nDEBUG: t = "<<t<<"; adding pre_strain of "<<pre_strain<<" um here";
            net->update_delrx( pre_strain );
            net->update_shear();
        }

        if (osc_strain_flag && t >= time_of_dstrain){
            d_strain = d_strain_amp * sin(2*pi*d_strain_freq * ( t - time_of_dstrain) );
            net->update_delrx( pre_strain + d_strain );
            net->update_d_strain( d_strain - prev_d_strain );
            cout<<"\nDEBUG: t = "<<t<<"; adding d_strain of "<<(d_strain-prev_d_strain)<<" um here; total strain = "<<(pre_strain+d_strain);
            prev_d_strain = d_strain;
        }
        
        if (diff_strain_flag && t >= time_of_dstrain && count%n_bw_shear==0){
            net->update_delrx( pre_strain + d_strain_amp );
            net->update_d_strain( d_strain_amp );
            pre_strain += d_strain_amp;
            cout<<"\nDEBUG: t = "<<t<<"; adding d_strain of "<<d_strain_amp<<" um here; total strain = "<<pre_strain<<" um";
        }

        if (count%n_bw_stdout==0) {
			cout<<"\nTime counts: "<<count;
		    //net->print_filament_thermo();
            net->print_network_thermo();
            crosslks->print_ensemble_thermo();
            myosins->print_ensemble_thermo();
        }

        //update network
        net->update();//updates all forces, velocities and positions of filaments

        if ( ! quad_off_flag && count % quad_update_period == 0)
            net->quad_update_serial();

        //update cross linkers
        if (static_cl_flag)
            crosslks->motor_update();
        else
            crosslks->motor_walk(t);
       
        //update motors
        myosins->motor_walk(t);
        
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
    file_ke.close();

    //Delete all objects created
    cout<<"\nHere's where I think I delete things\n";
    
//    delete lks;
    delete myosins;
    delete crosslks;
    delete net;
    
    
    
    cout<<"\nTime counts: "<<count;
	cout<<"\nExecuted";
	cout<<"\n Done\n";
    
    return 0;
}
