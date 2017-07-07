//#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std; 

int main()
{
        int beads = 5;
	int max_fil = 5; 
	int min_fil = 4; 
	int num_motors = 36; 
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

	ofstream myfile;
        myfile.open("spacers2_5.txt");

 	for(int a = min_fil; a < max_fil; a++)
  	{
		char file [100];
		num_fil = a + 1; 
		sprintf(file, "spacers2_%d.txt", num_fil);  
	
                //ofstream myfile; 
		//myfile.open("spacers2_5.txt"); 
          
		cout << "I work!" << endl; 
  
                int box = filaments * (num_fil - 1);  	
		double div = num_motors / box; 
  		double fil_div = spacing_x / div; 	//length to seperate filaments  
		double fil_row = div * filaments; 

		int start_point_x; 
		int start_point_y; 

		if(beads % 2 != 0) start_point_x = (filaments / 2) * spacing_x;	//centered at 0 
                else start_point_x = ((filaments - 1) / 2) * spacing_x + (0.5 * spacing_x);  
                if(num_fil % 2 != 0) ((num_fil - 1) / 2) * spacing_y; 
	        else start_point_y = ((num_fil / 2) - 1) * spacing_y + (0.5 * spacing_y);		
	 	//Starting with the top filamenr 
	 	//
       		for(int b = 0; b < filaments; b++)
		{
			y_1 = start_point_y;
	                y_2 = start_point_y - (spacing_y * num_fil);
			y_1 -= b * spacing_y; 
			y_2 -= b * spacing_y;  

			cout << "I still work!" << endl; 
                        cout << y_1 << endl;  

			for(int c = 0; c < num_motors; c++)
			{
				x_1 = start_point_x;
				cout << "x_1: " << x_1 << endl; 
	                        x_2 = start_point_x - (filaments * spacing_x);	
				x_1 -= c * div;  
 				cout << "X_1 new: " << x_1 << endl; 
				x_2 -= c * div; 
                            
            			cout << "Done!" << endl; 
                                sprintf(xi,"%d", x_1); 
			 	sprintf(xf,"%d", x_2);
				sprintf(yi,"%d", y_1);
				sprintf(yf,"%d", y_2);
  
  				char tab = '\t'; 
                                  			
				//cout << "xi: " << x_1 << "yi: " << y_1 << endl; 
  				//cout << "xf: " << x_2 << "yf: " << y_2 << endl; 
	
				myfile << xi << tab << yi << tab << xf << tab << yf << tab << "-1" << tab << "-1" << tab << "-1" << tab << "-1" << endl;
			}
		}

		myfile.close(); 
	} 

	return 0; 
}
