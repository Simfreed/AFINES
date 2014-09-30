#include "actin_ensemble.h"
#include "link_ensemble.h"
#include "motor_ensemble.h"
#include "globals.h"

#define xrange 100.0
#define yrange 100.0
#define tinit 0.0
#define tfinal 10 
// #define dt 0.0001 -- defined previously
#define print_dt 100
#define stdout_dt 1000

int main(int argc, char* argv[]){
    
    double xgrid=2*xrange;
    double ygrid=2*yrange;

    int seed=time(NULL);
    srand(seed);

/************************
 *      CONTROLS        *
 ***********************/
 
    // Motors
    double motor_length=0.5;
    double motor_density=0;
    double motor_stiffness=50.0; //pN / um
    double vmotor=1.0;
    double m_kon=90.0;          
    double m_kend=5.0;
    double m_koff=1.0; 
   
    // Actin 
    double actin_length=1; //(um) length of a monomer
    double npolymer = 1; 

    // Links 
    double link_stretching_stiffness = motor_stiffness;
    double link_length = actin_length/10; 
    double polymer_bending_modulus = temperature * 10; //using kT * Lp for bending modulus, with Lp = 10 um
    double link_bending_stiffness = polymer_bending_modulus * pow(1.0/actin_length,3);
    std::string link_color = "1"; //"blue";
    
    // Environment
    double viscosity=0.5;

    // Time 
    int count=0;
    double t=tinit;
    double ntimesteps = (tfinal - tinit)/dt;

    // Output
    std::string dir, afile, mfile, lfile;
    std::ofstream a_final, m_final, pl, pl_fourier, o_file;
	char numstr[21];
   
    
    /***********************
     * VARIABLES           *
     **********************/
    double nmonomer = 100;
//  double link_bending_stiffness = motor_stiffness/10;
    
    if (argc>1) {
        nmonomer                    =   atof(argv[1]);
        link_stretching_stiffness   =   atof(argv[2]);
        dir                         =        argv[3] ;
    }
    
    // DERIVED QUANTITIES :
    double actin_density = npolymer*nmonomer/(xrange*yrange);//0.65;
    std::string output_file                         =   dir + "/data/output.txt";
    std::string actin_output                        =   dir + "/data/actin_final.txt";
    std::string myosin_output                       =   dir + "/data/myosin_final.txt";
    std::string persistence_length_output           =   dir + "/data/angle_correlations.txt"; 
    std::string persistence_length_fourier_output   =   dir + "/data/fourier_modes.txt";
    
    // Fourier Modes
    int n_modes = nmonomer - 1;
    std::map<int, std::vector<double> > fm;
/*    for(int i = 0; i < n_modes; i ++){
        fm[i] = std::vector<double>();
    }
*/    // Angle corrleations
    std::vector<double> angle_correlations;
    std::vector<double> angle_correlations_sum;

    
    a_final.open(actin_output.c_str());
    m_final.open(myosin_output.c_str());
    pl.open(persistence_length_output.c_str());
    pl_fourier.open(persistence_length_fourier_output.c_str());
    o_file.open(output_file.c_str());
    
    o_file << " FILE: "                 << output_file     <<"\n";
    o_file << " Actin Density: "        << actin_density   << ", Actin Mean Length: "          << actin_length              << "\n";
    o_file << " Motor Density: "        << motor_density   << ", Motor Rest Length: "          << motor_length              << ", Motor Stiffness: "       << motor_stiffness        <<"\n";
    o_file << " Motor unloaded speed: " << vmotor          << ", Motor binding rate: "         << m_kon                     <<"\n";
    o_file << " Motor unbinding rate: " << m_koff          << ", Motor end detachment rate: "  << m_kend                    <<", Viscosity: "              << viscosity              <<"\n";
    o_file << " Link Rest Length: "     << link_length     << ", Link Stretching Stiffness: "  << link_stretching_stiffness <<", Link Bending Stiffness: " << link_bending_stiffness <<"\n";
    o_file << " Simulation time: "      << tfinal - tinit  << ", dt: " << dt <<", Number of time steps between output files: "<< print_dt<<"\n";
    o_file.close();
    
    
    std::cout<<"Creating actin network..\n";
	actin_ensemble * net = new actin_ensemble(actin_density,xrange,yrange,xgrid,ygrid,actin_length,viscosity,nmonomer,link_length);
    std::cout<<"Creating link ensemble...\n";
    link_ensemble * lks = new link_ensemble();
    std::cout<<"Adding links to connect actin filament monomers...\n";
    net->connect_polymers( lks, link_length, link_stretching_stiffness, link_bending_stiffness, link_color );
    std::cout<<"Adding motors...\n";
    motor_ensemble * myosins = new motor_ensemble( motor_density, xrange, yrange, motor_length, 
                                             net, vmotor, motor_stiffness, m_kon, m_koff,
                                             m_kend, actin_length, viscosity);
    std::cout<<"Updating motors, filaments and crosslinks in the network..\n";
    
    while (t<=tfinal) {
        //print time count
		if (count%stdout_dt==0) {
			std::cout<<"Time counts: "<<count<<"\n";
		}

        net->update_bending();
        net->update();
        net->quad_update();

        //print to file
	    if (count%print_dt==0) {
	        
            sprintf(numstr, "%d", count/print_dt);
            
            afile = dir + "/txt_stack/afile" + numstr + ".txt";
            mfile = dir + "/txt_stack/mfile" + numstr + ".txt";
            lfile = dir + "/txt_stack/lfile" + numstr + ".txt";
			
            std::ofstream file_a, file_m, file_l;
			file_a.open(afile.c_str());
			file_m.open(mfile.c_str());
            file_l.open(lfile.c_str());
		    
            net->write(file_a);
            myosins->motor_write(file_m);
            lks->link_write(file_l);
			file_a.close();
			file_m.close();
            file_l.close();
            
		}
        //update network
        myosins->motor_walk();
        lks->link_walk(); 
        
        angle_correlations = net->get_angle_correlation(0); //assume one polymer
        angle_correlations_sum = sum_vecs(angle_correlations_sum, angle_correlations);
        
        for (int n = 1; n <= n_modes; n++){
            //assuming only one polymer
            fm[n].push_back(net->get_fourier_mode(n, 0));
        }
        
        t+=dt;
		count++;
    }
    net->write(a_final);
    myosins->motor_write(m_final);
    
    //write the correlation file:
    //format of correlation file (assuming all polymers have the same monomer length, basically)
    
    // distance     correlation_polymer1        correlation polymer_2       correlation_polymer_3       ....
    
    for (unsigned int i = 1; i < angle_correlations_sum.size(); i++){
        pl<< i * actin_length << "\t" << angle_correlations_sum[i]/ntimesteps<<"\n";
    }
    
    // write the fourier mode file 
    for (int n = 1; n<= n_modes; n++){
        pl_fourier<< ( 1/(n*n*pi*pi) )<<"\t"<<mode_var(fm[n], 0)<<"\n"; 
    }
   
    //Delete all objects created
    std::cout<<"Here's where I think I delete things\n";
    
    delete lks;
    delete myosins;
    delete net;

    a_final.close();
    m_final.close();
    pl.close(); 
    pl_fourier.close();
    std::cout<<"\nTime counts: "<<count;
	std::cout<<"\nExecuted";
	std::cout<<"\n Done\n";
    return 0;
}
