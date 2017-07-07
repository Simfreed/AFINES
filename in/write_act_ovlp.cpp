#include <iostream> 
#include <fstream>
using namespace std; 

int main()
{
	int num_bead = 6; 
	int num_fil = 2;  
	double r = 0.25; 
	double spacing_x = 0.5; 
  	double spacing_y = 0.25; 
	char tab = '\t'; 
	char newline = '\n'; 

	double x; 
	double y; 
	double starting_x; 
	double starting_y; 
	double bead; 
	int fil; 
	char file[100]; 

	sprintf(file, "actins_f%i_b%i.txt", num_fil, num_bead); 
	ofstream myfile; 
	myfile.open(file); 

	for(int fil_index = 0; fil_index < num_fil ; fil_index ++)
	{
		fil = fil_index; 
                
		if(num_bead % 2 == 0) starting_x = ((num_bead/2 - 1) *spacing_x) + (0.5 * spacing_x); 
		else starting_x = ((num_bead - 1) / 2) * spacing_x; 
		
		if(num_fil % 2 == 0) starting_y = ((num_fil/2 - 1) * spacing_y) + (0.5 * spacing_y); 
		else starting_y = ((num_fil - 1) / 2) * spacing_y; 

		y = starting_y - fil * spacing_y; 

		for(int bead_index = 0; bead_index < num_bead; bead_index++) 
		{
			x = starting_x - bead_index * spacing_x; 

			myfile << x  << tab << y << tab << r << tab << fil << newline; 
		}
	}
	
	myfile.close(); 
	return 0; 
}
