//#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std; 

int main()
{
        int beads = 5;
	int max_fil = 5; 
	int min_fil = 2; 
	int num_motors = 50; 
        double spacing_x = 1.25;
	double spacing_y = 1.0;
	double start_point_x;
  	double start_point_y;  
        int filaments = beads - 1;  
	double x_1; 
	double y_1; 
	double x_2; 
	double y_2;
	int num_fil; 
	char xi [100]; 
     	char xf [100]; 
	char yi [100];
	char yf [100];
	double motor_l = -1.0; 

	//ofstream myfile;
        //myfile.open("spacers2_5.txt");

 	for(int a = min_fil; a < max_fil + 1; a++)
  	{
		char file [100];
		num_fil = a; 
		sprintf(file, "spacers2_%d.txt", num_fil);  
	
                ofstream myfile; 
		myfile.open(file); 
          
		cout << "I work!" << endl; 
  
                int box = filaments * (num_fil - 1);  	
		double div = num_motors / box; 
  		double fil_div = spacing_x / div; 	//length to seperate filaments  
		double fil_row = div * filaments; 

		cout << "Spacers per 'box': " << div << endl; 
		cout << "Spacers per Row: " << fil_row << endl; 
		cout << "Divisions: " << fil_div << endl; 

		double start_point_x; 
		double start_point_y; 

		if(beads % 2 != 0) start_point_x = (filaments / 2) * spacing_x;	//centered at 0 
                else start_point_x = ((filaments - 1) / 2) * spacing_x + (0.5 * spacing_x);  
                if(num_fil % 2 != 0) start_point_y = ((num_fil - 1) / 2) * spacing_y; 
	        else start_point_y = ((num_fil / 2) - 1) * spacing_y + (0.5 * spacing_y);	
	
	 	//Starting with the top filamenr 
	 	//
       		for(int b = 0; b < num_fil - 1; b++)
		{
			y_1 = start_point_y;
	                y_2 = motor_l;
			y_1 = y_1 - b * spacing_y; 
			y_2 = y_2 - b * spacing_y;  

			cout << "I still work!" << endl; 
                        cout << "y_1: " << y_1 << endl; 
			cout << "y_2: " << y_2 << endl;  

			for(int c = 0; c < fil_row; c++)
			{
				x_1 = start_point_x; 
	                        x_2 = 0.0; 
				x_1 = x_1 - c * fil_div;    

				cout << "Start point x: " << start_point_x; 
                              	cout << "x_1: " << x_1 << endl; 
				cout << "x_2: " << x_2 << endl; 
               			cout << "Done!" << endl; 
				//cout << "===============================" << endl; 

                                //sprintf(xi,"%d", x_1); 
			 	//sprintf(xf,"%d", x_2);
				//sprintf(yi,"%d", y_1);
				//sprintf(yf,"%d", y_2);
 
				//char* xi = (char*)(&x_1); 
				//char* yi = (char*)(&y_1); 
				//char* xf = (char*)(&x_2); 
				//char* yf = (char*)(&y_2); 

                                //cout << "xi: " <<  xi << endl;  
 
  				char tab = '\t'; 
                                char newline = '\n'; 
  			
				//cout << "xi: " << x_1 << "yi: " << y_1 << endl; 
  				//cout << "xf: " << x_2 << "yf: " << y_2 << endl; 
	
				myfile << x_1 << tab << y_1 << tab << x_2 << tab << y_2 << tab << "-1" << tab << "-1" << tab << "-1" << tab << "-1" << newline;
			}

			cout << "===========================" << endl; 
		}

		myfile.close(); 
	} 

	return 0; 
}
