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
    
    int myseed;
    
    double xrange, yrange;                                                  //Space
    int xgrid, ygrid;
    double grid_factor; 
    
    int    count, nframes, nmsgs;
    double t, tinit, tfinal, dt;
    
    double viscosity, temperature;                                          //Environment
    string bnd_cnd;                                         //Allowed values: NONE, PERIODIC, REFLECTIVE

    double actin_length, npolymer, nmonomer;                                // Actin 
    string actin_pos_str;
    
    double link_length, polymer_bending_modulus, link_stretching_stiffness, fene_pct, fracture_force; // Links

    double spacer1_length, spacer1_v=0, spacer1_density, spacer1_stiffness, s1_kon, s1_kend, s1_koff,
           s1_stall, s1_break, s1_bind;// Active Motors (i.e., "myosin")
    string spacer1_pos_str; 
    
    double spacer2_length, spacer2_density, spacer2_stiffness, // Passive Mtors (i.e., cross_linkers)
            spacer2_v=0, s2_kon, s2_kend, s2_koff, s2_stall, s2_break, s2_bind; 
    string spacer2_pos_str;
    
    string config_file, actin_in, spacer1_in, spacer2_in;                                                // Input configuration
    
    string   dir,    afile,  amfile,  pmfile,  lfile, thfile, pefile;                  // Output
    ofstream o_file, file_a, file_am, file_pm, file_l, file_th, file_pe;
    ios_base::openmode write_mode = ios_base::out;

    bool link_intersect_flag, motor_intersect_flag, dead_head_flag, p_dead_head_flag, static_cl_flag, quad_off_flag;
    double p_linkage_prob, a_linkage_prob;                                              
    int dead_head, p_dead_head;

    bool restart_actin, restart_spacer1, restart_spacer2;

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
        
        ("s1_kon", po::value<double>(&s1_kon)->default_value(100),"xlink 1 on rate")
        ("s1_koff", po::value<double>(&s1_koff)->default_value(20),"xlink 1 off rate")
        ("s1_kend", po::value<double>(&s1_kend)->default_value(20),"xlink 1 off rate at filament end")
        ("spacer1_length", po::value<double>(&spacer1_length)->default_value(0.15),"filamin rest length (um)")
        ("spacer1_stiffness", po::value<double>(&spacer1_stiffness)->default_value(1),"xlink 1 spring stiffness (pN/um)")
        ("spacer1_v", po::value<double>(&spacer1_v)->default_value(0),"xlink 1 velocity (um/s)")
        
        ("s2_kon", po::value<double>(&s2_kon)->default_value(100),"xlink 2 on rate")
        ("s2_koff", po::value<double>(&s2_koff)->default_value(20),"xlink 2 off rate")
        ("s2_kend", po::value<double>(&s2_kend)->default_value(20),"xlink 2 off rate at filament end")
        ("spacer2_length", po::value<double>(&spacer2_length)->default_value(0.03),"xlink 2 rest length (um) (default: filamin)")
        ("spacer2_stiffness", po::value<double>(&spacer2_stiffness)->default_value(1),"xlink 2 spring stiffness (pN/um)")
       
        ("s2_stall", po::value<double>(&s2_stall)->default_value(0),"force beyond which xlinks don't walk (pN)")
        ("s2_break", po::value<double>(&s2_break)->default_value(10),"force constant for xlink detachment (related to rupture force F_r, P(detach | F_r) -> 1) (pN)")
        ("s2_bind", po::value<double>(&s2_bind)->default_value(0.04),"binding energy of xlink (pN-um) (10kT by default)")

        ("s1_stall", po::value<double>(&s1_stall)->default_value(10),"force beyond which motors don't walk (pN)")
        ("s1_break", po::value<double>(&s1_break)->default_value(10),"force constant for motor detachment (pN)")
        ("s1_bind", po::value<double>(&s1_bind)->default_value(0.04),"binding energy of motor (pN-um) (10kT by default)")

        ("link_length", po::value<double>(&link_length)->default_value(1), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.04), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(100000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(1), "stiffness of link, pN/um")//probably should be about 70000 to correspond to actin
        ("fene_pct", po::value<double>(&fene_pct)->default_value(0.5), "pct of rest length of filament to allow outstretched until fene blowup")
        
        ("actin_in", po::value<string>(&actin_in)->default_value(""), "input actin positions file")
        ("spacer1_in", po::value<string>(&spacer1_in)->default_value(""), "input motor positions file")
        ("spacer2_in", po::value<string>(&spacer2_in)->default_value(""), "input crosslinker positions file")
        
        ("restart_actin", po::value<bool>(&restart_actin)->default_value(false), "if true, input actin positions file chosen by default")
        ("restart_spacer1", po::value<bool>(&restart_spacer1)->default_value(false), "if true, input motor positions file chosen by default")
        ("restart_spacer2", po::value<bool>(&restart_spacer2)->default_value(false), "if true, input crosslinker positions file chosen by default")
        
        ("dir", po::value<string>(&dir)->default_value("out/test"), "output directory")
        ("myseed", po::value<int>(&myseed)->default_value(time(NULL)), "Random number generator myseed")
        
        ("link_intersect_flag", po::value<bool>(&link_intersect_flag)->default_value(false), "flag to put a cross link at all filament intersections")
        ("motor_intersect_flag", po::value<bool>(&motor_intersect_flag)->default_value(false), "flag to put a motor at all filament intersections")
        ("p_linkage_prob", po::value<double>(&p_linkage_prob)->default_value(1), "If link_intersect_flag, probability that two filaments that intersect will be linked")
        ("a_linkage_prob", po::value<double>(&a_linkage_prob)->default_value(1), "If motor_intersect_flag, probability that two filaments that intersect will be motor-d")

        ("dead_head_flag", po::value<bool>(&dead_head_flag)->default_value(false), "flag to kill head <dead_head> of all motors")
        ("dead_head", po::value<int>(&dead_head)->default_value(0), "index of head to kill")
        
        ("p_dead_head_flag", po::value<bool>(&p_dead_head_flag)->default_value(false), "flag to kill head <dead_head> of all crosslinks")
        ("p_dead_head", po::value<int>(&p_dead_head)->default_value(0), "index of head to kill")
        
        ("static_cl_flag", po::value<bool>(&static_cl_flag)->default_value(false), "flag to indicate compeletely static xlinks; i.e, no walking, no detachment")
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

    double actin_density = npolymer*nmonomer/(xrange*yrange);//0.65;
    cout<<"\nDEBUG: actin_density = "<<actin_density; 
    double link_bending_stiffness    = polymer_bending_modulus / link_length;
    
    int n_bw_stdout = max(int(tfinal/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal - tinit)/(dt*double(nframes))),1);
    int unprinted_count = int(double(tinit)/dt);

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
    
    if (restart_actin)
        actin_in = dir + "/restart/in/actins.txt";
    if (restart_spacer1)
        spacer1_in = dir + "/restart/in/amotors.txt";
    if (restart_spacer2)
        spacer2_in = dir + "/restart/in/pmotors.txt";
    
    if (actin_in.size() > 0)
        actin_pos_vec   = file2vecvec(actin_in, "\t");
    if (spacer1_in.size() > 0)
        spacer1_pos_vec = file2vecvec(spacer1_in, "\t");
    if (spacer2_in.size() > 0)
        spacer2_pos_vec = file2vecvec(spacer2_in, "\t");
    
    set_seed(myseed);
    
    afile  = dir + "/txt_stack/actins.txt";
    lfile  = dir + "/txt_stack/links.txt";
    amfile = dir + "/txt_stack/amotors.txt";
    pmfile = dir + "/txt_stack/pmotors.txt";
    thfile = dir + "/data/thermo.txt";
    pefile = dir + "/data/pe.txt";

    if (restart_actin || restart_spacer1 || restart_spacer2) write_mode = ios_base::app;

    file_a.open(afile.c_str(), write_mode);
    file_l.open(lfile.c_str(), write_mode);
    file_am.open(amfile.c_str(), write_mode);
    file_pm.open(pmfile.c_str(), write_mode);
	file_th.open(thfile.c_str(), write_mode);
	file_pe.open(pefile.c_str(), write_mode);


    // DERIVED QUANTITIES :
    if(spacer1_density == 0 && spacer1_pos_vec.size() == 0 &&
            spacer2_density==0 && spacer2_pos_vec.size() == 0 &&
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
    if (actin_pos_vec.size() == 0){
        net = new filament_ensemble(actin_density, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, actin_length, viscosity, nmonomer, link_length, 
                actin_position_arrs, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd, myseed); 
    }else{
        net = new filament_ensemble(actin_pos_vec, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, viscosity, link_length, 
                link_stretching_stiffness, fene_pct, link_bending_stiffness,
                fracture_force, bnd_cnd); 
    }
   
    if (link_intersect_flag) spacer2_pos_vec = net->link_link_intersections(spacer2_length, p_linkage_prob); 
    if (motor_intersect_flag) spacer1_pos_vec = net->link_link_intersections(spacer1_length, a_linkage_prob); 
    if (quad_off_flag) net->turn_quads_off();

    cout<<"\nAdding xlink 1s...";
    spacer_ensemble * big_xlinks;
    
    if (spacer1_pos_vec.size() == 0)
        big_xlinks = new spacer_ensemble( spacer1_density, {xrange, yrange}, dt, temperature, 
                spacer1_length, net, spacer1_v, spacer1_stiffness, fene_pct, s1_kon, s1_koff,
                s1_kend, s1_stall, s1_break, s1_bind, viscosity, spacer1_position_arrs, bnd_cnd);
    else
        big_xlinks = new spacer_ensemble( spacer1_pos_vec, {xrange, yrange}, dt, temperature, 
                spacer1_length, net, spacer1_v, spacer1_stiffness, fene_pct, s1_kon, s1_koff,
                s1_kend, s1_stall, s1_break, s1_bind, viscosity, bnd_cnd);
    big_xlinks->set_bending(0.04, pi/2);
    if (dead_head_flag) big_xlinks->kill_heads(dead_head);

    cout<<"Adding xlink 2s (crosslinkers) ...\n";
    spacer_ensemble * small_xlinks; 
    
    if(spacer2_pos_vec.size() == 0)
        small_xlinks = new spacer_ensemble( spacer2_density, {xrange, yrange}, dt, temperature, 
                spacer2_length, net, spacer2_v, spacer2_stiffness, fene_pct, s2_kon, s2_kend,
                s2_kend, s2_stall, s2_break, s2_bind, viscosity, spacer2_position_arrs, bnd_cnd);
    else
        small_xlinks = new spacer_ensemble( spacer2_pos_vec, {xrange, yrange}, dt, temperature, 
                spacer2_length, net, spacer2_v, spacer2_stiffness, fene_pct, s2_kon, s2_kend,
                s2_kend, s2_stall, s2_break, s2_bind, viscosity, bnd_cnd);
    small_xlinks->set_bending(0.04, pi/2);
    if (p_dead_head_flag) small_xlinks->kill_heads(p_dead_head);

    // Write the output configuration file
    string output_file                         =   dir + "/data/output.txt";
    o_file.open(output_file.c_str());
    o_file << " FILE: "                 << output_file     <<"\n";
    o_file << " Actin Density: "        << actin_density   << ", Actin Mean Length: "          << actin_length              << "\n";
    o_file << " Active Motor Density: "        << spacer1_density   << ", Active Motor Rest Length: "          << spacer1_length              << ", Active Motor Stiffness: "       << spacer1_stiffness        <<"\n";
    o_file << " Active Motor unloaded speed: " << spacer1_v          << ", Active Motor binding rate: "         << s1_kon                     <<"\n";
    o_file << " Active Motor unbinding rate: " << s1_koff          << ", Active Motor end detachment rate: "  << s1_kend                    <<"\n";
    o_file << " Passive Motor Density: "        << spacer2_density   << ", Passive Motor Rest Length: "          << spacer2_length              << ", Passive Motor Stiffness: "       << spacer2_stiffness        <<"\n";
    o_file << " Passive Motor unloaded speed: " << spacer2_v          << ", Passive Motor binding rate: "         << s2_kon                     <<"\n";
    o_file << " Passive Motor unbinding rate: " << s2_koff          << ", Passive Motor end detachment rate: "  << s2_kend                    <<"\n";
    o_file << " Link Rest Length: "     << link_length     << ", Link Stretching Stiffness: "  << link_stretching_stiffness <<", Link Bending Stiffness: " << link_bending_stiffness <<"\n";
    o_file << " Simulation time: "      << tfinal - tinit  << ", dt: " << dt <<", dt between output files: "<< n_bw_print*dt<<", Viscosity: " << viscosity              <<"\n";
    o_file << " Boundary Conditions: " <<bnd_cnd<<"\n";
    o_file.close();
    
    cout<<"\nUpdating motors, filaments and crosslinks in the network..";
    string time_str; 
    count=0;
    t = 0;

    //Run the simulation
    while (t < tfinal) {
        
        //print to file
	    if (t+dt/100 >= tinit && (count-unprinted_count)%n_bw_print==0) {
	        
            if (t>tinit) time_str ="\n";
            time_str += "t = "+to_string(t);
            
            file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
            net->write_actins(file_a);
            
            file_l << time_str<<"\tN = "<<to_string(net->get_nlinks());
            net->write_links(file_l);
            
            file_am << time_str<<"\tN = "<<to_string(big_xlinks->get_nmotors());
            big_xlinks->motor_write(file_am);
            
            file_pm << time_str<<"\tN = "<<to_string(small_xlinks->get_nmotors());
            small_xlinks->motor_write(file_pm);
            
            file_th << time_str<<"\tN = "<<to_string(net->get_nlinks());
            net->write_thermo(file_th);

            file_pe << net->get_stretching_energy()<<"\t"<<net->get_bending_energy()<<"\t"<<
                big_xlinks->get_potential_energy()<<"\t"<<small_xlinks->get_potential_energy()<<endl;
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
