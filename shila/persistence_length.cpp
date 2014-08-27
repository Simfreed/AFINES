#include "actin_myosin_flexible.cpp"
#define tinit 0.0
#define tfinal 10 
//#define dt 0.001 -- defined previously
#define print_dt 100

int main(int argc, char* argv[]){
    
    int seed=time(NULL);
    srand(seed);
    
    std::string output_file="output.txt";
    std::string actin_output="actin_final.txt";
    std::string myosin_output="myosin_final.txt";
    std::string persistence_length_output = "angle_correlations.txt"; 
    // CONTROLS :
    double motor_length=0.5;
    double motor_density=0.5;
    double motor_stiffness=50.0; //pN / um
    double vmotor=1.0;
    double m_kon=90.0;          
    double m_kend=5.0;
    double m_koff=1.0; 
    double viscosity=0.5;
    

    // VARIABLES :
    double npolymer = 1, nmonomer = 100;
    double actin_length=30/nmonomer; //length of a monomer
    
    // DERIVED QUANTITIES :
    
    double actin_density= npolymer*nmonomer/(xrange*yrange);//0.65;
    double link_length = actin_length/10; 
    double link_stiffness = motor_stiffness/2;
    
    std::string link_color = "2"; //"blue";
    
    if (argc>1) {
        nmonomer = std::atof(argv[1]);
        link_stiffness = std::atof(argv[2]);
    }
    
    
    std::ofstream a_final, m_final, p_final;
    a_final.open(actin_output.c_str());
    m_final.open(myosin_output.c_str());
    p_final.open(persistence_length_output.c_str());
    std::ofstream o_file;
    o_file.open(output_file.c_str());
    o_file << " FILE: " << output_file <<"\n"<< "Actin Density: " << actin_density  << ", Actin Mean Length: " << actin_length << "\n";
    o_file << " Motor Density: " << motor_density << ", Motor Rest Length: " << motor_length << ", Motor Stiffness: " << motor_stiffness<<"\n";
    o_file << ", Motor unloaded speed: " << vmotor << ", Motor binding rate: " << m_kon <<"\n"<<"Motor unbinding rate: " << m_koff << ", Motor end detachment rate: " << m_kend<<", Viscosity: "<<viscosity<<"\n";
    o_file << " Link Rest Length: "<< link_length <<", Link Stiffness: " << link_stiffness <<"\n";
    o_file << " Simulation time: " << tfinal - tinit << ", dt: " << dt <<", Number of time steps between output files: "<< print_dt<<"\n";
    o_file.close();
    
    
    int count=0;
    int timesteps=ceil(tfinal/dt);
	char afile[timesteps+2];
	char mfile[timesteps+2];
    char tfile[timesteps+2];
    
    
    std::cout<<"\nCreating actin network..";
	actin_ensemble net=actin_ensemble(actin_density,xrange,yrange,xgrid,ygrid,actin_length,viscosity,link_length);
    std::cout<<"\nAdding motors..";
    motor_ensemble myosins=motor_ensemble(motor_density,xrange,yrange,motor_length,&net,vmotor,motor_stiffness,m_kon,m_koff,m_kend,actin_length,viscosity);
    std::cout<<"\nAdding links to connect actin filament monomers...";
    net.connect_polymers( &myosins, link_length, link_stiffness, link_color );

    double t=tinit;
    std::cout<<"\nUpdating motors, filaments and crosslinks in the network..";
    
    std::map<double, std::map<int, std::map<double, double> > > angle_correlations_time; //maps time --> {polymer --> {length -> angle correlation} }
    std::map<int, std::map<double, double> > angle_correlations_sum;
    std::map<int, std::map<double, double> >::iterator it;
    std::map<double, double>::iterator it1;
        

    while (t<=tfinal) {

        //print time count
		if (count%1000==0) {
			std::cout<<"\nTime counts: "<<count;
		}
        
        //update network
        net.update();
        net.quad_update();
        //print to file
		angle_correlations_time[t] = net.get_all_angle_correlations();
		/*if (count%print_dt==0) {
		
            sprintf (afile, "afile%d.txt", count/print_dt);
			sprintf (mfile, "mfile%d.txt", count/print_dt);
            sprintf (tfile, "tfile%d.txt", count/print_dt);
			std::ofstream file_a, file_m, file_t;
			file_a.open(afile);
			file_m.open(mfile);
            file_t.open(tfile);
			net.write(file_a);
			myosins.motor_write(file_m);
            myosins.motor_tension(file_t);
			file_a.close();
			file_m.close();
            file_t.close();
            
		}
        */
        angle_correlations_time[ t ] = net.get_all_angle_correlations();
        for (it = angle_correlations_time[t].begin(); it != angle_correlations_time[t].end(); ++it){
            // it->first is the polymer index
            // it->second is a map of length --> correlation
            for (it1 = it->second.begin(); it1 != it->second.end(); ++it1){
                angle_correlations_sum[it->first][it1->first] += it1->second;
            }
        }
        
        //myosins.reshape();
        myosins.motor_walk();
        //
        t+=dt;
		count++;
    }
    net.write(a_final);
    myosins.motor_write(m_final);
    
    //write the correlation file:
    //format of correlation file (assuming all polymers have the same monomer length, basically)
    
    // distance     correlation_polymer1        correlation polymer_2       correlation_polymer_3       ....
    
    std::vector<std::string> lines;
    lines.push_back("header");
    int line_count = 0, ntimesteps = (int)((tfinal - tinit)/dt);
//    for (int i = 0; i < angle_correlations_sum.size(); i++){
    for (it = angle_correlations_sum.begin(); it != angle_correlations_sum.end(); ++it){
        // it->first is the polymer index
        // it->second is a map of length --> correlation
        lines[0] += "\t" + std::to_string(it->first);
        for (it1 = it->second.begin(); it1 != it->second.end(); ++it1){
            lines.push_back(std::to_string(it1->first) + "\t" + std::to_string(it1->second/ntimesteps));
        }
    }
    
    for (int l = 0; l<lines.size(); l++){
        p_final<<lines[l]<<"\n";
    }
    a_final.close();
    m_final.close();
    p_final.close(); 

    std::cout<<"\nTime counts: "<<count;
	std::cout<<"\nExecuted";
	std::cout<<"\n Done\n";
    return 0;
}

