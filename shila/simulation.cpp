#include "actin_myosin_flexible.cpp"
#define tinit 0.0
#define tfinal 10 
//#define dt 0.001 -- defined previously
#define print_dt 100
int main(int argc, char* argv[])
{
    int seed=time(NULL);
    srand(seed);
    
    double npolymer_desired = 1, nmonomer_desired = 10;
    std::string mainpath="./";
    std::string output_file="output.txt";
    std::string actin_output="actin_final.txt";
    std::string myosin_output="myosin_final.txt";
    double actin_length=3.0; //length of a monomer
    double actin_density= npolymer_desired*nmonomer_desired/(xrange*yrange);//0.65;
    double motor_length=0.5;
    double motor_density=0.5;
    double motor_stiffness=50.0;
    double vmotor=1.0;
    double m_kon=90.0;          
    double m_kend=5.0;
    double m_koff=1.0; 
    double viscosity=0.5;
   
    if (argc>1) {
        mainpath=argv[1];
        output_file=argv[2];
        actin_length=std::atof(argv[3]);
        actin_density=std::atof(argv[4]);
        motor_length=std::atof(argv[5]);
        motor_density=std::atof(argv[6]);
        motor_stiffness=std::atof(argv[7]);
        vmotor=std::atof(argv[8]);
        m_kon=std::atof(argv[9]);
        m_kend=std::atof(argv[10]);
        m_koff=std::atof(argv[11]);
        actin_output=argv[12];
        myosin_output=argv[13];
        viscosity=std::atof(argv[14]);
    }
    double link_length = actin_length/10; 
    double link_stiffness = motor_stiffness/2;
    
    std::ofstream a_final, m_final;
    a_final.open(actin_output.c_str());
    m_final.open(myosin_output.c_str());
    
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
    std::cout<<"\nAdding links to connect motors...";
    std::map<int, std::vector<double> > * link_map_ptr = net.get_link_map();
    std::map<int, std::vector<int> > * mono_map_ptr = net.get_mono_map();
    int monomer_index;
    std::vector<double> motor_coords;
    // loop through each actin polymer
    std::string link_color = "2"; //"blue";
    for (int i = 0; i < mono_map_ptr->size(); i++){
//        std::cout<<"DEBUG: i:"<<i<<";\tmono_map_ptr->at(i).size(); "<<mono_map_ptr->at(i).size()<<"\n";        
        for (int j = 0; j < mono_map_ptr->at(i).size(); j++){
//            for(std::vector<int>::iterator it=mono_map_ptr->at(i).begin(); it<mono_map_ptr->at(i).end() - 1; it++)
            monomer_index = mono_map_ptr->at(i).at(j);
            motor_coords = link_map_ptr->at(monomer_index);
//            std::cout<<"\nDEBUG: Adding a link to the "<<i<<"th polymer at the "<<monomer_index<<"th monomer\n";            
//            std::cout<<"DEBUG: motor_coords[0] : "<<motor_coords[0]<<";motor_coords[1] : "<<motor_coords[1]<<";motor_coords[2] : "<<motor_coords[2]<<"\n";
            myosins.add_motor( motor( motor_coords[0], motor_coords[1], motor_coords[2], 
                       link_length, &net, 1, 1, monomer_index, monomer_index+1, xrange, yrange, 0, link_stiffness,
                       0, 0, 0, actin_length, viscosity, link_color) );
        }
    }
    double t=tinit;
    std::cout<<"\nUpdating motors, filaments and crosslinks in the network..";
    
    std::map<double, std::map<int, std::map<double, double> > > angle_correlations_time; //maps time --> {polymer --> {length -> angle correlation} }
    while (t<=tfinal) {
        //print time count
		if (count%1000==0) {
			std::cout<<"\nTime counts: "<<count;
		}
        
        
        //update network
        
        net.update();
        net.quad_update();
        //print to file
		
		if (count%print_dt==0) {
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
        angle_correlations_time[ t ] = net.get_all_angle_correlations();
        //myosins.reshape();
        myosins.motor_walk();
        //
        t+=dt;
		count++;
    }
    net.write(a_final);
    myosins.motor_write(m_final);
    a_final.close();
    m_final.close();
    
    std::cout<<"\nTime counts: "<<count;
	std::cout<<"\nExecuted";
	std::cout<<"\n Done\n";
    
    return 0;
}

