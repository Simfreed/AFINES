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
    bool use_linear_bending;

    double a_motor_length, a_motor_v, a_motor_density, a_motor_stiffness, a_m_kon, a_m_kend, a_m_koff;// Active Motors (i.e., "myosin")
    string a_motor_pos_str; 
    
    double p_motor_length, p_motor_density, p_motor_stiffness, // Passive Mtors (i.e., cross_linkers)
            p_motor_v=0, p_m_kon, p_m_kend, p_m_koff; 
    string p_motor_pos_str;
    
    string config_file, filament_type, actin_in, a_motor_in, p_motor_in;                                                // Input configuration
    
    string   dir,    afile,  amfile,  pmfile,  lfile, thfile;                  // Output
    ofstream o_file, file_a, file_am, file_pm, file_l, file_th;

    double shear_rate, shear_freq, shear_stop, strain_pct;                     //External Force
    
    bool link_intersect_flag;

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
        ("grid_factor", po::value<int>(&grid_factor)->default_value(2), "number of grid boxes per um^2")
        
        
        ("dt", po::value<double>(&dt)->default_value(0.001), "length of individual timestep in seconds")
        ("tfinal", po::value<double>(&tfinal)->default_value(10), "length of simulation in seconds")
        ("nframes", po::value<int>(&nframes)->default_value(1000), "number of timesteps between printing actin/link/motor positions to file")
        ("nmsgs", po::value<int>(&nmsgs)->default_value(10000), "number of timesteps between printing simulation progress to stdout")
       
        ("viscosity", po::value<double>(&viscosity)->default_value(1), "Dynamic viscosity to determine friction [mg / (um*s)]. At 20 C, is 1 for water")
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT [pN-um] that effects magnituded of Brownian component of simulation")
        ("bnd_cnd,bc", po::value<string>(&bnd_cnd)->default_value("REFLECTIVE"), "boundary conditions")
        
        ("nmonomer", po::value<double>(&nmonomer)->default_value(1), "number of monomers per filament")
        ("npolymer", po::value<double>(&npolymer)->default_value(3), "number of polymers in the network")
        ("actin_length", po::value<double>(&actin_length)->default_value(0.5), "Length of a single actin monomer")
        ("actin_pos_str", po::value<string> (&actin_pos_str)->default_value(""), "Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions")
        
        ("a_motor_density", po::value<double>(&a_motor_density)->default_value(0.001), "number of active motors / um^2")
        ("p_motor_density", po::value<double>(&p_motor_density)->default_value(0.001), "number of passive motors / um^2")
        ("a_motor_pos_str", po::value<string> (&a_motor_pos_str)->default_value(""), "Starting positions of motors, commas delimit coordinates; semicolons delimit positions")
        ("p_motor_pos_str", po::value<string> (&p_motor_pos_str)->default_value(""), "Starting positions of crosslinks, commas delimit coordinates; semicolons delimit positions")
        
        ("a_m_kon", po::value<double>(&a_m_kon)->default_value(90.0),"active motor on rate")
        ("a_m_koff", po::value<double>(&a_m_koff)->default_value(1),"active motor off rate")
        ("a_m_kend", po::value<double>(&a_m_kend)->default_value(5),"active motor off rate at filament end")
        ("a_motor_length", po::value<double>(&a_motor_length)->default_value(0.5),"active motor rest length (um)")
        ("a_motor_stiffness", po::value<double>(&a_motor_stiffness)->default_value(50),"active motor spring stiffness (pN/um)")
        ("a_motor_v", po::value<double>(&a_motor_v)->default_value(1),"active motor velocity (um/s)")
        
        ("p_m_kon", po::value<double>(&p_m_kon)->default_value(90),"passive motor on rate")
        ("p_m_koff", po::value<double>(&p_m_koff)->default_value(0.01),"passive motor off rate")
        ("p_m_kend", po::value<double>(&p_m_kend)->default_value(0.01),"passive motor off rate at filament end")
        ("p_motor_length", po::value<double>(&p_motor_length)->default_value(0.5),"passive motor rest length (um)")
        ("p_motor_stiffness", po::value<double>(&p_motor_stiffness)->default_value(50),"passive motor spring stiffness (pN/um)")
        
        ("link_length", po::value<double>(&link_length)->default_value(1), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.04), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(100000000), "pN-- filament breaking point")
        ("bending_fracture_force", po::value<double>(&bending_fracture_force)->default_value(1000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(1), "stiffness of link, pN/um")//probably should be about 70000 to correspond to actin
        ("use_linear_bending,linear", po::value<bool>(&use_linear_bending)->default_value(false),"option to send spring type of bending springs")
        
        ("shear_rate", po::value<double>(&shear_rate)->default_value(0), "shear rate in pN/um")
        ("strain_pct", po::value<double>(&strain_pct)->default_value(0), "pct that the boundarys get sheared")
        ("shear_freq", po::value<double>(&shear_freq)->default_value(1), "Frequency of shearing in # of timesteps (i.e., shear_freq = 0.33 ==> shear every third time step)")
        ("shear_stop", po::value<double>(&shear_stop)->default_value(1e12), "amount of time to shear for [s]")
        
        ("actin_in", po::value<string>(&actin_in)->default_value(""), "input actin positions file")
        ("a_motor_in", po::value<string>(&a_motor_in)->default_value(""), "input motor positions file")
        ("p_motor_in", po::value<string>(&p_motor_in)->default_value(""), "input crosslinker positions file")
        
        ("dir", po::value<string>(&dir)->default_value("out/test"), "output directory")
        ("seed", po::value<int>(&seed)->default_value(time(NULL)), "Random number generator seed")
        ("filament_type", po::value<string>(&filament_type)->default_value("BAOAB"), "type of filament / integrator to use.\
                                                                            Options: (1) BD = Ermak Yeh Brownian Dynamics\
                                                                                     (2) BAOAB = Charlie's Overdamped BAOAB\
                                                                                     (3) LLF = Langevin Leap Frog")
        
        ("link_intersect_flag", po::value<bool>(&link_intersect_flag)->default_value(false), "output directory")
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
    if(a_motor_density == 0 && p_motor_density==0){
        xgrid = 0;
        ygrid = 0;
    }
    else{
        xgrid  = (int) grid_factor*xrange;
        ygrid  = (int) grid_factor*yrange;
    }

    if (polymer_bending_modulus < 0){ //This is a flag for using the temperature for the bending modulus
        polymer_bending_modulus = 10*temperature; // 10um * kT
    }

    double actin_density = npolymer*nmonomer/(xrange*yrange);//0.65;
    
    double link_bending_stiffness    = polymer_bending_modulus * pow(1.0/link_length , 1);
    
    int n_bw_stdout = max(int((tfinal - tinit)/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal - tinit)/(dt*double(nframes))),1);

    // To Read positions from input strings in config file
    vector<array<double,3> > actin_position_arrs, a_motor_position_arrs, p_motor_position_arrs;
    if (actin_pos_str.size() > 0)
        actin_position_arrs   = str2arrvec(actin_pos_str, ":", ",");
    if (a_motor_pos_str.size() > 0)
        a_motor_position_arrs = str2arrvec(a_motor_pos_str, ":", ",");
    if (p_motor_pos_str.size() > 0)
        p_motor_position_arrs = str2arrvec(p_motor_pos_str, ":", ",");
    
    // To Read positions from input files
    vector<vector<double> > actin_pos_vec;
    vector<vector<double> > a_motor_pos_vec, p_motor_pos_vec;
    if (actin_in.size() > 0)
        actin_pos_vec   = file2vecvec(actin_in, "\t");
    if (a_motor_in.size() > 0)
        a_motor_pos_vec = file2vecvec(a_motor_in, "\t");
    if (p_motor_in.size() > 0)
        p_motor_pos_vec = file2vecvec(p_motor_in, "\t");
    
    srand(seed);
            
    
    afile  = dir + "/txt_stack/actins.txt";
    lfile  = dir + "/txt_stack/links.txt";
    amfile = dir + "/txt_stack/amotors.txt";
    pmfile = dir + "/txt_stack/pmotors.txt";
    thfile = dir + "/data/thermo.txt";

    file_a.open(afile.c_str());
    file_l.open(lfile.c_str());
    file_am.open(amfile.c_str());
    file_pm.open(pmfile.c_str());
	file_th.open(thfile.c_str());


    // Create Network Objects
    cout<<"\nCreating actin network..";
    ATfilament_ensemble * net;
    if (actin_pos_vec.size() == 0){
        net = new ATfilament_ensemble(actin_density, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, actin_length, viscosity, nmonomer, link_length, 
                actin_position_arrs, 
                link_stretching_stiffness, link_bending_stiffness,
                fracture_force, bnd_cnd, seed); 
    }else{
        net = new ATfilament_ensemble(actin_pos_vec, {xrange, yrange}, {xgrid, ygrid}, dt, 
                temperature, viscosity, link_length, 
                link_stretching_stiffness, link_bending_stiffness,
                fracture_force, bnd_cnd); 
    }
   
    if (link_intersect_flag) p_motor_pos_vec = net->get_intersections(p_motor_length); 

    cout<<"\nAdding active motors...";
    motor_ensemble<ATfilament_ensemble> * myosins;
    
    if (a_motor_pos_vec.size() == 0)
        myosins = new motor_ensemble<ATfilament_ensemble>( a_motor_density, {xrange, yrange}, dt, temperature, 
                a_motor_length, net, a_motor_v, a_motor_stiffness, a_m_kon, a_m_koff,
                a_m_kend, actin_length, viscosity, a_motor_position_arrs, bnd_cnd);
    else
        myosins = new motor_ensemble<ATfilament_ensemble>( a_motor_pos_vec, {xrange, yrange}, dt, temperature, 
                a_motor_length, net, a_motor_v, a_motor_stiffness, a_m_kon, a_m_koff,
                a_m_kend, actin_length, viscosity, bnd_cnd);

    cout<<"Adding passive motors (crosslinkers) ...\n";
    motor_ensemble<ATfilament_ensemble> * crosslks; 
    
    if(p_motor_pos_vec.size() == 0)
        crosslks = new motor_ensemble<ATfilament_ensemble>( p_motor_density, {xrange, yrange}, dt, temperature, 
                p_motor_length, net, p_motor_v, p_motor_stiffness, p_m_kon, p_m_kend,
                p_m_kend, actin_length, viscosity, p_motor_position_arrs, bnd_cnd);
    else
        crosslks = new motor_ensemble<ATfilament_ensemble>( p_motor_pos_vec, {xrange, yrange}, dt, temperature, 
                p_motor_length, net, p_motor_v, p_motor_stiffness, p_m_kon, p_m_kend,
                p_m_kend, actin_length, viscosity, bnd_cnd);

    cout<<"\nUpdating motors, filaments and crosslinks in the network..";

    double shear_dt = dt/shear_freq;
    if (strain_pct != 0){
        // Based on completing a shear of strain_pct * xrange within the simulation, 
        // By shearing at 1 time step and then allowing the system to relax for t = shear_dt 
        // strain_dist = strain_pct * xrange = shear_rate * actin_mobility * boundary_height * tfinal / 2
        double strain_dist = strain_pct * yrange;
        shear_rate = 8 * pi * actin_length * viscosity * strain_dist / (yrange * tfinal * shear_freq);
        //shear_rate = 16. * pi * actin_length * viscosity * xrange * strain_pct/ (yrange * tfinal);
        
        net->set_shear_rate(shear_rate);
        net->set_shear_dt(shear_dt);
        net->set_shear_stop(shear_stop);
    }

    string time_str = "t = 0";

    file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
    net->write_actins(file_a);
    file_l << time_str<<"\tN = "<<to_string(net->get_nlinks());
    net->write_links(file_l);
    file_am << time_str<<"\tN = "<<to_string(myosins->get_nmotors());
    myosins->motor_write(file_am);
    file_pm << time_str<<"\tN = "<<to_string(crosslks->get_nmotors());
    crosslks->motor_write(file_pm);
    file_th << time_str;
    net->write_thermo(file_th);

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
    
    //Run the simulation
    while (t<tfinal) {
        //print time count
		if (count%n_bw_stdout==0) {
			cout<<"\nTime counts: "<<count;
		    //net->print_filament_thermo();
            net->print_network_thermo();
        }

        //update network
        net->update();//updates all forces, velocities and positions of filaments

        //update motors and cross linkers
        crosslks->motor_walk(t);
        myosins->motor_walk(t);
        
        //clear the vector of fractured filaments
        net->clear_broken();

        t+=dt;
		count++;

        //print to file
	    if (count%n_bw_print==0) {
	        
            time_str = "\nt = "+to_string(t);
            
            file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
            net->write_actins(file_a);
            
            file_l << time_str<<"\tN = "<<to_string(net->get_nlinks());
            net->write_links(file_l);
            
            file_am << time_str<<"\tN = "<<to_string(myosins->get_nmotors());
            myosins->motor_write(file_am);
            
            file_pm << time_str<<"\tN = "<<to_string(crosslks->get_nmotors());
            crosslks->motor_write(file_pm);
            
            file_th << time_str<<"\tN = "<<to_string(net->get_nlinks());
            net->write_thermo(file_th);
		}
        
    }
    
    file_a.close();
    file_l.close();
    file_am.close();
    file_pm.close();
    file_th.close(); 
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
