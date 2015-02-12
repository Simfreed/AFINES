/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "iostream"
#include "iomanip"
#include <math.h>
#include "fstream"
#include "string"
#include "stdint.h"
#include <stdlib.h>
#include "vector"
#include <map>
#include "time.h"
/* distances in microns, time in seconds, forces in pN */

#define pi 3.14159265358979323
#define eps 0.001


//These two things should probably be inputs to various classes
#define dt 0.0001
#define temperature 0//0.004

class motor_ensemble; //forward declaration of motor ensemble because I need it in actin_ensemble
class link_ensemble; 

/*generic functions to be used below*/
double rng(double start, double end)
{
	return start+(end-start)*((double)rand()/(RAND_MAX));
}

int pr(int num)
{
	if (num==0) {
		return 1;	
	}
	else {
		return 0;
	}
}

double rng_exp(double mean)
{
    double u;
    u=rand() / (RAND_MAX + 1.);;
    return  -mean*log(u);
}


double rng_n(double mean, double var)
{
    static double U, V;
    static int phase = 0;
    double Z;

    if(phase == 0) {
        U = (rand() + 1.) / (RAND_MAX + 2.);
        V = rand() / (RAND_MAX + 1.);
        Z = sqrt(-2 * log(U)) * sin(2 * pi * V);
    } else
        Z = sqrt(-2 * log(U)) * cos(2 * pi * V);

    phase = 1 - phase;

    return mean+var*Z;
}

int event(double rate, double timestep)
{
    if (rng(0,1.0)<rate*timestep) {
        return 1;
    }
    else
        return 0;
}

double dis_points(double x1, double y1, double x2, double y2)
{
    double dis=sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    return dis;
}

double velocity(double vel0, double force, double fstall)
{
    if (force>=fstall) {
        return 0;
    }
    else if (force>0 && force<fstall) {
        return vel0*(1-(fabs(force)/fstall));
    }
    else if (force<=0 && force>=-fstall) {
        return vel0*(1+(fabs(force)/fstall));
    }
    else{
        return 2*vel0;
    }
}

double cross(double ax, double ay, double bx, double by)
{
    return ax*by-bx*ay;
}

double dot(double x1, double y1, double x2, double y2)
{
    return x1*x2+y1*y2;
}

double mean(std::vector<double> vals)
{
    double sum = 0;
    for (int i = 0; i < vals.size(); i++){
        sum+= vals[i];
    }
    return sum / vals.size();

}

double var(std::vector<double> vals)
{
    double m = mean(vals), sum = 0;
    for (int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

double mode_var(std::vector<double> vals, double m)
{
    double sum = 0;
    for (int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

std::vector<double> sum_vecs(std::vector<double> v1, std::vector<double> v2)
{
    std::vector<double> s;
    if (v1.empty())
        s = v2;
    else if( v2.empty())
        s = v1;
    else if (v1.size() != v2.size())
        return s;
    else{
        for (int i = 0; i < v1.size(); i++){
            s.push_back(v1[i] + v2[i]);
        }
    }
    return s;
}
//actin filament class
class actin
{
public:
actin(double xcm, double ycm, double angle, double len, double fovx, double fovy, int nx, int ny, double vis)
{
    x=xcm;
    y=ycm;
    phi=angle;
    ld=len;
    diameter=ld/40;
    start[0]=x-ld*0.5*cos(phi);
    start[1]=y-ld*0.5*sin(phi);
    end[0]=x+ld*0.5*cos(phi);
    end[1]=y+ld*0.5*sin(phi);
    a_vis=vis;
    //unit vector
        e[0]=cos(phi);
        e[1]=sin(phi);
        //unit normal
        n[0] = -e[1];
        n[1] = e[0];
		//motor-induced forces
		forces[0]=0; //along the filament
        forces[1]=0; //perpendicular to the filament
        forces[2]=0; //torque
        
        //quadrant numbers crossed by the actin in x-direction
        quad.clear();
        tmp.clear();
		int lower_limit, upper_limit, index;
        if(start[0] <= end[0])
        {
            lower_limit = int(floor(start[0]/fovx*nx));
            if(lower_limit > 0){lower_limit--;};
            upper_limit = int(ceil(end[0]/fovx*nx));
            if(upper_limit < nx-1){upper_limit++;};
            
            for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);};
        }
        else
        {
            lower_limit = int(floor(end[0]/fovx*nx));
            if(lower_limit > 0){lower_limit--;};
            upper_limit = int(ceil(start[0]/fovx*nx));
            if(upper_limit < nx-1){upper_limit++;};
            for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);};
        };
        quad.push_back(tmp);
        
        //quadrant numbers crossed by the actin in y-direction
        tmp.clear();
        if(start[1] <= end[1])
        {
            lower_limit = int(floor(start[1]/fovy*ny));
            if(lower_limit > 0){lower_limit--;};
            upper_limit = int(ceil(end[1]/fovy*ny));
            if(upper_limit < ny-1){upper_limit++;};
            for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);}
        }
        else
        {
            lower_limit = int(floor(end[1]/fovy*ny));
            if(lower_limit > 0){lower_limit--;};
            upper_limit = int(ceil(start[1]/fovy*ny));
            if(upper_limit < ny-1){upper_limit++;};
            for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);}
        };
        quad.push_back(tmp);
	}
    ~actin(){};
    
    //shortest(perpendicular) distance between an arbitray point and the filament
    double get_distance(double xp, double yp)
    {
        double l2=pow(dis_points(start[0],start[1],end[0],end[1]),2);
        if (l2==0) {
            return dis_points(xp,yp,start[0],start[1]);
        }
        double tp=dot(xp-start[0],yp-start[1],end[0]-start[0],end[1]-start[1])/l2;
        if (tp<0) {
            return dis_points(xp,yp,start[0],start[1]);
        }
        else if(tp>1.0){
            return dis_points(xp,yp,end[0],end[1]);
        }
        else{
            double px=start[0]+tp*(end[0]-start[0]);
            double py=start[1]+tp*(end[1]-start[1]);
            return dis_points(xp,yp,px,py);
        }
    }
    
    double* get_intpoint(double xp, double yp)
    {
        double* points;
        double coordinates[2];
        double l2=pow(dis_points(start[0],start[1],end[0],end[1]),2);
        if (l2==0) {
            coordinates[0]=start[0];
            coordinates[1]=start[1];
        }
        double tp=dot(xp-start[0],yp-start[1],end[0]-start[0],end[1]-start[1])/l2;
        if (tp<0) {
            coordinates[0]=start[0];
            coordinates[1]=start[1];
        }
        else if(tp>1.0){
            coordinates[0]=end[0];
            coordinates[1]=end[1];
        }
        else{
            coordinates[0]=start[0]+tp*(end[0]-start[0]);
            coordinates[1]=start[1]+tp*(end[1]-start[1]);
        }
        points=coordinates;
        return points;
    }
    
	double get_int_angle(double xp, double yp)
    {
        double angle;
        double xcor,ycor;
        double slope=(end[1]-start[1])/(end[0]-start[0]);
        double yintercept=y-slope*x;
        xcor=(slope*yp + xp - slope*yintercept)/(slope*slope + 1);
        ycor=(slope*slope*yp + slope*xp + yintercept)/(1 + slope*slope);
        angle=atan2((ycor-yp),(xcor-xp));
        return angle;
    }
	
    double*  get_direction()
	{
		return e;
	}
    
    double get_length()
    {
        return ld;
    }
    
    double* get_forces()
    {
        double *fpr;
        fpr=forces;
        return fpr;
    }
	
	void update_force(double f1, double f2, double f3)
	{
		forces[0]+=f1;
		forces[1]+=f2;
		forces[2]+=f3;
	}
    
    double* get_friction()
    {
        double fric[3];
        double *fcr;
        fric[0]=2*pi*a_vis*ld/log(ld/diameter);
        fric[1]=2*fric[0];
        fric[2]=fric[0]*pow(ld,2)/6;
        fcr=fric;
        return fcr;
    }
    
	double* getpos()
	{
		double pos[3];
		double *ptr;
		pos[0]=x*e[0]+y*e[1];
		pos[1]=x*n[0]+y*n[1];
		pos[2]=phi;
		ptr=pos;
		return ptr;
	}
    
    double* getposcm()
	{
		double poscm[3];
		double *ptrs;
		poscm[0]=x;
		poscm[1]=y;
		poscm[2]=phi;
		ptrs=poscm;
		return ptrs;
	}
    
    double* getendpts()
    {
        double endpts[4];
        double *pts;
        endpts[0]=start[0];
        endpts[1]=start[1];
        endpts[2]=end[0];
        endpts[3]=end[1];
        pts=endpts;
        return pts;
    }
    
    std::vector<std::vector<int> > get_quadrants()
    { 
        return quad; 
    }
	
private:
	double x,y,phi,ld, start[2], end[2], e[2], n[2], forces[3];
    double xperp, xpar, friction_perp, friction_par, diameter, a_vis;
    std::vector<std::vector<int> > quad; //vector of two vectors(x and y quadrants) of integers
    std::vector<int> tmp;
};


//actin network class

class actin_ensemble
{
public:
    actin_ensemble(double density, double fovx, double fovy, int nx, int ny, double len, double vis, int nmonomer, double link_len)
    {
        view=0.9;
        fov[0]=fovx;
        fov[1]=fovy;
        nq[0]=nx;
        nq[1]=ny;
        rho=density;
        av_vel=0;
        visc=vis;
//        nmonomer_min = 100; //hard coded number of min/max monomers per filament
//        nmonomer_max = nmonomer_min;
//        int nmonomer = (nmonomer_max + nmonomer_min)/2;
        npolymer=int(ceil(density*fov[0]*fov[1]) / nmonomer);
        ld=len;//rng_n(len,1.0);
        link_ld = link_len;
        std::cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
        std::cout<<"DEBUG: Number of monomers per filament:"<<nmonomer<<"\n"; 
        std::cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
        std::vector<double> link_map_value; // xcm, ycm, angle
        double * rod_end_pts;
        double  xcm, ycm, theta;
        for (int i=0; i<npolymer; i++) {
            theta=3*pi/2;//rng(0,2*pi);
            //nmonomer = (int) rng(nmonomer_min, nmonomer_max);
            //the start of the polymer: 
//            network.push_back(actin(rng(-0.5*(view*fovx-ld),0.5*(view*fovx-ld)), rng(-0.5*(view*fovy-ld),0.5*(view*fovy-ld)),
//                        theta,ld,fov[0],fov[1],nq[0],nq[1],visc));
            std::cout<<"WARNING: STARTING ACTIN FILAMENT POSITION CHOSEN DETERMINISTICALLY\n";
            network.push_back(actin(0,0,theta,ld,fov[0],fov[1],nq[0],nq[1],visc));
            //Add the quadrants of the first rod
            std::vector<std::vector<int> > tmp_quads=network.back().get_quadrants();
            for (int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
                for (int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                    quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(network.size()-1);
                }
            }
            // add monomers to the polymer
            mono_map[i] = empty_vector; 
            for (int j=0; j<nmonomer-1; j++) {
                
                // Calculate a link (represented as a dead motor) at the end of the rod
                ltheta = theta; //rng(0, 2*pi);
                rod_end_pts = network.back().getendpts();

                lx=rod_end_pts[2] + 0.5*link_ld*cos(ltheta);
                ly=rod_end_pts[3] + 0.5*link_ld*sin(ltheta);

                // Calculate the Next rod on the actin polymer--  continues from the static motor
                //theta = rng(0,2*pi);
                xcm = lx + link_ld*0.5*cos(ltheta) + ld*0.5*cos(theta);
                ycm = ly + link_ld*0.5*sin(ltheta) + ld*0.5*sin(theta);
                
                // Check that this monomer is in the field of view, otherwise start a new polymer:
                if ( xcm > (0.5*(view*fovx - ld)) || xcm < (-0.5*(view*fovx - ld)) 
                        || ycm > (0.5*(view*fovy - ld)) || ycm < (-0.5*(view*fovy - ld)) )
                {
                    std::cout<<"\nDEBUG:"<<j+1<<"th monomer of "<<i<<"th polymer outside field of view; stopped building polymer\n";
                    break;
                }else{
                    
                    // Load the link map: {monomer index} --> {link_xcm, link_ycm, link_angle} 
                    link_map_value.clear();
                    link_map_value.push_back(lx);
                    link_map_value.push_back(ly);
                    link_map_value.push_back(ltheta);
                    link_map[ network.size()-1 ] = link_map_value;
                    mono_map[ i ].push_back(network.size() -1);
                    
                    // Add the actin monomer
                    network.push_back( actin(xcm, ycm, theta, ld, fov[0], fov[1], nq[0], nq[1], visc) );
                    
                    // Add it's quadrants:
                    std::vector<std::vector<int> > tmp_quads=network.back().get_quadrants();
                    for (int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
                        for (int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                            quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(network.size()-1);
                        }
                    }
                } 
                
            }
        }
        
        
    }
    ~actin_ensemble(){};
    
    
    void quad_update()
    {
        quad_fils.clear();
        for (int i=0; i<network.size(); i++) {
            std::vector<std::vector<int> > tmp_quads=network[i].get_quadrants();
            for (int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
                for (int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                    quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(i);
                }
            }
        }
        
        
    }
    
    std::vector<actin>* get_network()
	{
		return &network;
	}
    
    //given motor quadrant, return the indices and the distances to all actin filaments in the neighboring quadrants 
    std::map<int,double> get_dist(double x, double y)
    {
        int nxx=int(floor(x/fov[0]*nq[0]));
        int nyy=int(floor(y/fov[1]*nq[1]));
        t_map.clear();
        if(!quad_fils[nxx][nyy].empty())
        {
            for(std::vector<int>::iterator it=quad_fils[nxx][nyy].begin(); it<quad_fils[nxx][nyy].end(); it++)
            {   
                t_map[*it] = network[*it].get_distance(x,y);
            }
        }
        return t_map;
    }
    
    double* get_direction(int index)
	{
		return network[index].get_direction();
	}
	
    
    double* get_intpoints(int index, double xp, double yp)
    {
        return network[index].get_intpoint(xp,yp);
    }
	
	double get_int_direction(int index, double xp, double yp)
	{
		return network[index].get_int_angle(xp,yp);
	}
    
    double* get_position(int index)
    {
        return network[index].getposcm();
    }
    
    double get_alength(int index)
    {
        return network[index].get_length();
    }
    
    double* get_ends(int index)
    {
        return network[index].getendpts();
    }
    
    void update()
    {
        av_vel=0;
        for (int i=0; i<network.size(); i++) {
            vpar=(network[i].get_forces()[0])/network[i].get_friction()[0]  + sqrt(6*temperature/(dt*network[i].get_friction()[0]))*rng_n(0,1);
            vperp=(network[i].get_forces()[1])/network[i].get_friction()[1] + sqrt(6*temperature/(dt*network[i].get_friction()[1]))*rng_n(0,1);
            vx=vpar*cos(network[i].getpos()[2])-vperp*sin(network[i].getpos()[2]);
            vy=vpar*sin(network[i].getpos()[2])+vperp*cos(network[i].getpos()[2]);
            omega=network[i].get_forces()[2]/network[i].get_friction()[2] + sqrt(6*temperature/(dt*network[i].get_friction()[2]))*rng_n(0,1);
            
            alength=network[i].get_length();
            
            xnew=network[i].getposcm()[0]+dt*vx;
            ynew=network[i].getposcm()[1]+dt*vy;
            phinew=network[i].getpos()[2]+dt*omega;
            
            
            a_ends[0]=xnew-alength*0.5*cos(phinew);
            a_ends[1]=ynew-alength*0.5*sin(phinew);
            a_ends[2]=xnew+alength*0.5*cos(phinew);
            a_ends[3]=ynew+alength*0.5*sin(phinew);
            
            if (a_ends[0]<=-fov[0]*0.5 || a_ends[0]>=fov[0]*0.5 || a_ends[2]<=-fov[0]*0.5 || a_ends[2]>=fov[0]*0.5)
            {
                vx=-vx;//xnew=network[i].getposcm()[0]-dt*vx;//  
                omega=-omega;//phinew=network[i].getpos()[2]-dt*omega;//omega=-omega;
                
            }
            if (a_ends[1]<=-fov[1]*0.5 || a_ends[1]>=fov[1]*0.5 || a_ends[3]<=-fov[1]*0.5 || a_ends[3]>=fov[1]*0.5)
            {
                vy=-vy;//ynew=network[i].getposcm()[1]-dt*vy;
                omega=-omega;//phinew=network[i].getpos()[2]-dt*omega;
            }
			
            xnew=network[i].getposcm()[0]+dt*vx;
            ynew=network[i].getposcm()[1]+dt*vy;
            phinew=network[i].getpos()[2]+dt*omega;
			network.at(i)=actin(xnew,ynew,phinew,alength,fov[0],fov[1],nq[0],nq[1],visc);
        }
        
        
    }
    
    void update_forces(int index, double f1, double f2, double f3)
	{
		network[index].update_force(f1,f2,f3);
	}
    
    
    void write(std::ofstream& fout)
    {
        for (int i=0; i<network.size(); i++) {
            fout<<network.at(i).getendpts()[0]<<"\t"<<network.at(i).getendpts()[1]<<"\t"<<network.at(i).getendpts()[2]-network.at(i).getendpts()[0]<<"\t"<<network.at(i).getendpts()[3]-network.at(i).getendpts()[1]<<"\n";
        } 
    }
    /*
    wikipedia says < cos(theta) > = Exp[-L/P]
    where L is the length of the polymer and P is the persistence length
    and theta is the angle between subsequent rods
    This function will return a map of x --> <cos theta>
    where x is a distance and <cos theta> is the correlation between cosines of the angles between subsequent
    rods. 
    */ 
    std::vector<double> get_angle_correlation(int polymer_index)
    {
        double sum = 0, phi1, phi2;
        std::vector<double> corr;
        std::vector<int> monos = mono_map[polymer_index]; 
        for(int i = 0; i < monos.size() - 1; i++){
            
            phi1 = network.at(monos[i]).getposcm()[2]; 
            phi2 = network.at(monos[i+1]).getposcm()[2]; 
            
            sum += cos(phi2 - phi1);
            
            corr.push_back(sum/(i+1));
        
        }
        return corr;

    }
    
    std::map<int, std::vector<double> > get_all_angle_correlations()
    {
        std::map<int, std::vector<double> > corrs;
        
        for(int i = 0; i < mono_map.size(); i++){
            corrs[i] = this->get_angle_correlation(i);
        }

        return corrs;
    }
    
    /* Based on derivation in :
     * Flexural Rigidity of Microtubules and Actin Filaments etc. 
     * Gittes, Micky, Nettelton Howard
     * Journal of Cell Biology; Feb 1993
     */
    double get_fourier_mode(int n, int polymer_index){
        
        std::vector<int> mons = mono_map[polymer_index];
        int s = mons.size();
        double L = s * ld, sum = 0, phi, sk;
        
        for (int i = 0; i < s; i++){

            phi = network[ mons[ i ] ].getposcm()[2];
            sk = i * ld + ld/2;
            sum += phi * ld * cos( n * pi * sk / L);

        }

        return sqrt(2.0/L)*sum;
    
    }
    
    std::map<int, std::vector<double> > * get_link_map(){
        return &link_map;
    }
    
    std::map<int, std::vector<int> > * get_mono_map(){
        return &mono_map;
    }
    
void connect_polymers(motor_ensemble * mots, double link_length, double link_stiffness, std::string link_color);
void connect_polymers(link_ensemble * links, 
        double link_length, double link_stiffness, std::string link_color,
        double b_link_stiffness, std::string b_link_color);

private:
    double fov[2], rho, ld, xnew, ynew, phinew, vpar, omega, vperp, vx, vy, alength, view, a_ends[4], av_vel, visc;
    double link_ld, lx, ly, ltheta;
    int npolymer, nmonomer_max, nmonomer_min, nq[2], xn, yn, qxcm, qycm;
    std::vector<actin> network;
    std::vector<int> empty_vector;
    std::map<int, std::vector<double> > link_map; //Maps {monomer} --> {link_xcm, link_ycm, link_angle} 
    std::map<int, std::vector<int> > mono_map; //Maps {filament index} --> {list of monomer indices}
    std::map<int, std::map<int,std::vector<int> > > quad_fils;
    std::map<int,double> t_map;
    
};


//motor class
class motor 
{
public:
    motor(double mx, double my, double mang, double mlen, actin_ensemble* network, int state0, int state1, int aindex0, int aindex1, double fovx, double fovy, double v0, double stiffness, double ron, double roff, double rend, double actin_len, double vis, std::string col)
    {
        vs=v0;//rng_n(v0,0.4);//rng(v0-0.3,v0+0.3);
        dm=0.25;//actin_len/10;
        mk=stiffness;//rng(10,100);
        fmax=mk*dm*2;//rng(1,20);
        mld=mlen;
        kon=ron;
        koff=roff;
        kend=rend;
        mphi=mang;
        hx[0]=mx-0.5*mld*cos(mphi);
        hy[0]=my-0.5*mld*sin(mphi);
        hx[1]=hx[0]+mld*cos(mphi);
        hy[1]=hy[0]+mld*sin(mphi);
        mobility=log(10)/(4*pi*vis*mld);
		state[0]=state0;
        state[1]=state1;
        aindex[0]=aindex0;//   actin index for head in state[0]
        aindex[1]=aindex1;// actin index for head in state[1]
        actin_network=network;
//		pos_actin[0]=0; // I don't think this variable get's ACCESSED anywhere 
        pos_a_end[0]=0; //distance from pointy end -- by default 0
        pos_a_end[1]=0;

        if (state0){
            pos_a_end[0] = fmin(dis_points(hx[0],hy[0],actin_network->get_ends(aindex[0])[0],actin_network->get_ends(aindex[0])[1]),
                                dis_points(hx[0],hy[0],actin_network->get_ends(aindex[0])[2],actin_network->get_ends(aindex[0])[3]));
		    std::cout<<"DEBUG: distance of head 0 from pointy end of actin monomer "<<aindex[0]<<" : "<<pos_a_end[0]<<"\n";
        }
        if (state1){
            pos_a_end[1] = fmin(dis_points(hx[1],hy[1],actin_network->get_ends(aindex[1])[0],actin_network->get_ends(aindex[1])[1]),
                                dis_points(hx[1],hy[1],actin_network->get_ends(aindex[1])[2],actin_network->get_ends(aindex[1])[3]));
		    std::cout<<"DEBUG: distance of head 1 from pointy end of actin monomer "<<aindex[1]<<" : "<<pos_a_end[0]<<"\n";
        }
//        pos_actin[1]=0;
		fov[0]=fovx;
		fov[1]=fovy;
        color = col; 
        
    }
    ~motor(){};
    
    //return motor state with a given head number
    int* get_states() 
    {
        int* sptr;
        sptr=state;
        return sptr;
    }
    
     
    double* get_heads()
    {
        double h[4];
        double *gh;
        h[0]=hx[0];
        h[1]=hy[0];
        h[2]=hx[1];
        h[3]=hy[1];
        gh=h;
        return gh;
    }
    
    std::string get_color()
    {
        return color;
    }

    double tension()
    {
        double lf=dis_points(hx[0],hy[0],hx[1],hy[1]);
        return mk*(lf-mld)/fmax;
    }
    
    //check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
	void attach(int hd)
	{
        dist.clear();
        dist=actin_network->get_dist(hx[hd],hy[hd]);
        if(!dist.empty()){
            for (std::map<int,double>::iterator it=dist.begin(); it!=dist.end(); ++it)
            { 
                if (it->second <= dm && aindex[pr(hd)]!=it->first) {
                    onrate=kon*exp(-((it->second)*(it->second))/(dm*dm));
                    if (event(onrate,dt)==1) {
                        //update state
                        state[hd]=1;
                        aindex[hd]=it->first;
                        //update head position
                        hx[hd]=actin_network->get_intpoints(it->first,hx[hd],hy[hd])[0];
                        hy[hd]=actin_network->get_intpoints(it->first,hx[hd],hy[hd])[1];
                        if (state[pr(hd)]==0) {
                            hx[pr(hd)]=hx[hd]+pow(-1,hd)*mld*cos(mphi);
                            hy[pr(hd)]=hy[hd]+pow(-1,hd)*mld*sin(mphi);
                        }
                        else {
                            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
                        }
                        
//                        pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
                        pos_a_end[hd]=dis_points(hx[hd],hy[hd],actin_network->get_ends(aindex[hd])[2],actin_network->get_ends(aindex[hd])[3]);
                        break;
                    }
                }
            }
        }	
	} 
	
	//perform brownian motion if head unattached
	void brownian()
	{
		if (state[0]==0 && state[1]==0) {
            
			xm[0]=hx[0]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hx[0]-hx[1]+mld*cos(mphi));
			xm[1]=hx[1]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) + mk*dt*(hx[0]-hx[1]+mld*cos(mphi));
			ym[0]=hy[0]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
			ym[1]=hy[1]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) + mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
            reflect(xm[0],xm[1],ym[0],ym[1]);
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
		}
		else if (state[0]==0 || state[1]==0) {
			int hd=state[0];
            
			xm[hd]=hx[hd]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hx[hd]-hx[pr(hd)]+pow(-1,hd)*mld*cos(mphi));
			ym[hd]=hy[hd]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hy[hd]-hy[pr(hd)]+pow(-1,hd)*mld*sin(mphi));
            xm[pr(hd)]=hx[pr(hd)];
            ym[pr(hd)]=hy[pr(hd)];
            reflect(xm[0],xm[1],ym[0],ym[1]);
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
		}
		else {
			return;
		}
        
	}
    
    //stepping and detachment kinetics of a single bound head 
	void step_onehead(int hd)
	{
        
		if (event(koff,dt)==1) {
			state[hd]=0;
			aindex[hd]=-1;
//			pos_actin[hd]=0;
            pos_a_end[hd]=0;
            hx[hd]=hx[pr(hd)]-pow(-1,hd)*mld*cos(mphi);
			hy[hd]=hy[pr(hd)]-pow(-1,hd)*mld*sin(mphi);
			
		}
		else
		{
            pos_temp=pos_a_end[hd]+dt*vs;
			move_end_detach(hd,vs,pos_temp);
            
			
		}
	}
	
	//stepping and detachment kinetics of a motor with both heads attached
    void step_twoheads()
    {
		stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
        std::cout<<"DEBUG: step_twoheads: color = "<< color<< "\tstretch = "<<stretch<<"\n";
        //fm = vec(Fm).(-vec(u)) 
		fm[0]=mk*((hx[0]-hx[1]+mld*cos(mphi))*actin_network->get_direction(aindex[0])[0] + (hy[0]-hy[1]+mld*sin(mphi))*actin_network->get_direction(aindex[0])[1]);
        vm[0]=velocity(vs,fm[0],fmax);
		fm[1]=mk*(-(hx[0]-hx[1]+mld*cos(mphi))*actin_network->get_direction(aindex[1])[0] - (hy[0]-hy[1]+mld*sin(mphi))*actin_network->get_direction(aindex[1])[1]);
        
		vm[1]=velocity(vs,fm[1],fmax);
		/*impose force-dependent bell's law on detachment rates*/
		offrate[0]=koff*exp(fabs(fm[0])/fmax);
		offrate[1]=koff*exp(fabs(fm[1])/fmax);
        
		if (event(offrate[0],dt)==1) {
			state[0]=0;
			aindex[0]=-1;
			hx[0]=hx[1]-mld*cos(mphi);
			hy[0]=hy[1]-mld*sin(mphi);
//			pos_actin[0]=0;
            pos_a_end[0]=0;
            pos_temp=pos_a_end[1]+dt*vs;
			move_end_detach(1,vs,pos_temp);
		}
		else if (event(offrate[1],dt)==1) {
            state[1]=0;
            aindex[1]=-1;
            hx[1]=hx[0]+mld*cos(mphi);
            hy[1]=hy[0]+mld*sin(mphi);
//            pos_actin[1]=0;
            pos_a_end[1]=0;
            pos_temp=pos_a_end[0]+dt*vs;
			move_end_detach(0,vs,pos_temp);
        }
        
        else {
            
            pos_temp=pos_a_end[0]+dt*vm[0];
            move_end_detach(0,vm[0],pos_temp);
            pos_temp=pos_a_end[1]+dt*vm[1];
            move_end_detach(1,vm[1],pos_temp); 
        }
	}
	
	
    
    
    
    
    void actin_update()
    {
        if (state[0]==1 && state[1]==1) {
			stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
            std::cout<<"DEBUG: actin_update: color = "<< color<< "\tstretch = "<<stretch<<"\n";
            forcex[0]=-mk*(hx[0]-hx[1]+mld*cos(mphi));
            forcex[1]=-forcex[0];
            forcey[0]=-mk*(hy[0]-hy[1]+mld*sin(mphi));
            forcey[1]=-forcey[0];
            force_par[0]=forcex[0]*actin_network->get_direction(aindex[0])[0] + forcey[0]*actin_network->get_direction(aindex[0])[1];
            force_perp[0]=-forcex[0]*actin_network->get_direction(aindex[0])[1] + forcey[0]*actin_network->get_direction(aindex[0])[0];
            force_par[1]=forcex[1]*actin_network->get_direction(aindex[1])[0] + forcey[1]*actin_network->get_direction(aindex[1])[1];
            force_perp[1]=-forcex[1]*actin_network->get_direction(aindex[1])[1] + forcey[1]*actin_network->get_direction(aindex[1])[0];
            
            torque[0]=cross(hx[0]-actin_network->get_position(aindex[0])[0],hy[0]-actin_network->get_position(aindex[0])[1],forcex[0],forcey[0]);
            torque[1]=cross(hx[1]-actin_network->get_position(aindex[1])[0],hy[1]-actin_network->get_position(aindex[1])[1],forcex[1],forcey[1]);
            actin_network->update_forces(aindex[0],force_par[0],force_perp[0],torque[0]);
            actin_network->update_forces(aindex[1],force_par[1],force_perp[1],torque[1]);
        }
        else
            return;
        
    }
    
    void update_shape()
    {
        if (state[0]==1 && state[1]==1) {
            hx[0]=actin_network->get_ends(aindex[0])[2]-pos_a_end[0]*actin_network->get_direction(aindex[0])[0];
            hy[0]=actin_network->get_ends(aindex[0])[3]-pos_a_end[0]*actin_network->get_direction(aindex[0])[1];
//            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_position(aindex[0])[0],actin_network->get_position(aindex[0])[1]);
            hx[1]=actin_network->get_ends(aindex[1])[2]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
            hy[1]=actin_network->get_ends(aindex[1])[3]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
//            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_position(aindex[1])[0],actin_network->get_position(aindex[1])[1]);
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
        }
        else if(state[0]==1 && state[1]==0)
        {
            hx[0]=actin_network->get_ends(aindex[0])[2]-pos_a_end[0]*actin_network->get_direction(aindex[0])[0];
            hy[0]=actin_network->get_ends(aindex[0])[3]-pos_a_end[0]*actin_network->get_direction(aindex[0])[1];
//            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_position(aindex[0])[0],actin_network->get_position(aindex[0])[1]);
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
        }
        else if(state[0]==0 && state[1]==1)
        {
            hx[1]=actin_network->get_ends(aindex[1])[2]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
            hy[1]=actin_network->get_ends(aindex[1])[3]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
//            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_position(aindex[1])[0],actin_network->get_position(aindex[1])[1]);
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
        }
        else
        {
            return;
        }
        
        
    }
    
private:
    double hx[2],hy[2], xm[2], ym[2], mphi,mld,mobility,fm[2],vm[2],offrate[2], onrate, stretch,forcex[2],forcey[2],torque[2], force_par[2],force_perp[2],vs, pos_temp;
    int state[2], aindex[2];
    double dm,fmax,mk, kon, koff, kend;
    std::map<int, double> dist;
    std::string color;
    actin_ensemble *actin_network;
	double /*pos_actin[2],*/ pos_a_end[2], fov[2];
    
    
    inline void move_end_detach(int hd, double speed, double pos)
    {
        stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
        std::cout<<"DEBUG: move_end_detach, before, hd = "<<hd<<" : color = "<< color<< "\tstretch = "<<stretch<<"\n";
        if (pos>=actin_network->get_alength(aindex[hd])) { //not sure why only "greater than or equal"
            if (event(kend,dt)==1) {
                std::cout<<"The new myosin position of head "<<hd<<" is OFF the actin filament AND detaching\n";
                state[hd]=0;
                aindex[hd]=-1;
//                pos_actin[hd]=0;
                pos_a_end[hd]=0;
                hx[hd]=hx[pr(hd)]-pow(-1,hd)*mld*cos(mphi);
                hy[hd]=hy[pr(hd)]-pow(-1,hd)*mld*sin(mphi);
            }
            else {
                std::cout<<"The new myosin position of head "<<hd<<" is OFF the actin filament BUT not detaching\n";
                hx[hd]=actin_network->get_ends(aindex[hd])[2]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
                hy[hd]=actin_network->get_ends(aindex[hd])[3]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
                mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
//                pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
            }
        }
        else {
            std::cout<<"The new myosin position of head "<<hd<<" is "<<pos<<" away from the pointy end of the actin filament \n"; //still on the actin filament";
            pos_a_end[hd]=pos;
            hx[hd]=actin_network->get_ends(aindex[hd])[2]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
            hy[hd]=actin_network->get_ends(aindex[hd])[3]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
//            pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
        }
        stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
        std::cout<<"DEBUG: move_end_detach, after, hd = "<<hd<<" : color = "<< color<< "\tstretch = "<<stretch<<"\n";
        
    }
    
    inline void reflect(double x1, double x2, double y1, double y2)
    {
        if (-fov[0]*0.5<x1<fov[0]*0.5 && -fov[0]*0.5<x2<fov[0]*0.5 && -fov[1]*0.5<y1<fov[1]*0.5 && -fov[1]*0.5<y2<fov[1]*0.5) {
            hx[0]=x1;
            hx[1]=x2;
            hy[0]=y1;
            hy[1]=y2;
        }
        else if (x1>=fov[0]*0.5 || x1<=-fov[0]*0.5)
        {
            hx[1]=x2;
            hy[0]=y1;
            hy[1]=y2;   
        }
        else if (x2>=fov[0]*0.5 || x2<=-fov[0]*0.5)
        {
            hx[0]=x1;
            hy[0]=y1;
            hy[1]=y2;
        }
        else if(y1>=fov[1]*0.5 || y1<=fov[1]*0.5)
        {
            hx[0]=x1;
            hx[1]=x2;
            hy[1]=y2;
        }
        else{
            hx[0]=x1;
            hx[1]=x2;
            hy[0]=y1;
        }
    }
    
    
};


//motor ensemble class
class motor_ensemble
{
public:
    motor_ensemble(double mdensity, double fovx, double fovy, double mlen, actin_ensemble* network, double v0, double stiffness, double ron, double roff, double rend, double actin_len, double vis)
    {
        fov[0]=fovx;
        fov[1]=fovy;
        mrho=mdensity;
        mld=mlen;
        nm=int(ceil(mrho*fov[0]*fov[1]));
        std::cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";
        a_network=network;
        alpha=0.8;
        color = "0.5";//"green"; 
        for (int i=1; i<=nm; i++) {
            motorx=rng(-0.5*(fovx*alpha-mld),0.5*(fovx*alpha-mld));
            motory=rng(-0.5*(fovy*alpha-mld),0.5*(fovy*alpha-mld));
            mang=rng(0,2*pi);
            n_motors.push_back(motor(motorx,motory,mang,mld,a_network,0,0,-1,-1,fov[0],fov[1],v0,stiffness,ron,roff,rend,actin_len,vis,color));
        }
    }
    ~motor_ensemble(){};
    
    void motor_walk()
    {
        
        for (int i=0; i<n_motors.size(); i++) {
            
            s[0]=n_motors[i].get_states()[0];
            s[1]=n_motors[i].get_states()[1];
			
            
            if (s[0]==0 && s[1]==0) {
                n_motors[i].attach(0);
                n_motors[i].attach(1);
				n_motors[i].brownian();
            }
            else if (s[0]==0 && s[1]==1) {
                n_motors[i].attach(0);
				n_motors[i].brownian();
                n_motors[i].step_onehead(1);
            }
            else if (s[0]==1 && s[1]==0) {
                n_motors[i].attach(1);
				n_motors[i].brownian();
                n_motors[i].step_onehead(0);
            }
            else {
                n_motors[i].step_twoheads();
            }
            
            n_motors[i].actin_update();
        }
        
    }
    
    void reshape()
    {
        for (int i=0; i<n_motors.size(); i++) {
            n_motors[i].update_shape();
        }
    }
    
    
    
    void motor_write(std::ofstream& fout)
    {
        for (int i=0; i<n_motors.size(); i++) {
            double stretch=dis_points(n_motors[i].get_heads()[0],n_motors[i].get_heads()[1],n_motors[i].get_heads()[2],n_motors[i].get_heads()[3])-mld;
         /*   if (stretch>3*0.25) {
                continue;
            }
            else{
         */   fout<<n_motors[i].get_heads()[0]<<"\t"<<n_motors[i].get_heads()[1]<<"\t"<<n_motors[i].get_heads()[2]-n_motors[i].get_heads()[0]<<"\t"<<n_motors[i].get_heads()[3]-n_motors[i].get_heads()[1]<<"\t"<<n_motors[i].get_color()<<"\n";
           //}
        } 
    }
    
    void motor_tension(std::ofstream& fout)
    {
        for (int i=0; i<n_motors.size(); i++) {
            fout<<n_motors[i].tension()<<"\n";
        }
    }
    
    void add_motor(motor m)
    {
        n_motors.push_back(m);
    }

private:
    double fov[2], mrho, mld, mang, motorx, motory, mphi, mcor[3],alpha;
    int nm, s[2],a[2];
    actin_ensemble *a_network;
    std::vector<motor> n_motors;  
    std::string color;
};

//Link class
class Link
{
public:
    Link(double len, double stiffness, actin_ensemble* network, 
            int aindex0, int aindex1, std::string col)
    {
        lk              =   stiffness;
        ld              =   len;
        aindex[0]       =   aindex0;
        aindex[1]       =   aindex1;
        actin_network   =   network;
		
        // Set the coordinates of the heads:
        this->step();
        color           =   col; 
        
    }
    ~Link(){}
    
     
    double* get_heads()
    {
        double h[4];
        double *gh;
        h[0]=hx[0];
        h[1]=hy[0];
        h[2]=hx[1];
        h[3]=hy[1];
        gh=h;
        return gh;
    }
    
    std::string get_color()
    {
        return color;
    }

	// stepping kinetics
    void step()
    {
    
        //CONVENTION: head 0 will be connected to the POINTY end of a filament
        //            head 1 will be connected to the BARBED end of a filament
        hx[0]=actin_network->get_ends(aindex[0])[2];
        hy[0]=actin_network->get_ends(aindex[0])[3];
        hx[1]=actin_network->get_ends(aindex[0])[0];
        hy[1]=actin_network->get_ends(aindex[0])[1];
        
        phi=atan2(hy[1]-hy[0],hx[1]-hx[0]);
        
	}
	
    void actin_update()
    {
        stretch         =   dis_points(hx[0],hy[0],hx[1],hy[1])-ld;
        std::cout<<"DEBUG: actin_update: color = "<< color<< "\tstretch = "<<stretch<<"\n";
    
        forcex[0]       =   lk * stretch * cos(phi); //-lk*(hx[0]-hx[1]+ld*cos(phi)); <-- old formula, probably equivalent
        forcex[1]       =   -forcex[0];
        forcey[0]       =   lk * stretch * sin(phi); //-lk*(hy[0]-hy[1]+ld*sin(phi)); <-- old formula, probably equivalent
        forcey[1]       =   -forcey[0];
        
        force_par[0]    =   forcex[0]*actin_network->get_direction(aindex[0])[0] + forcey[0]*actin_network->get_direction(aindex[0])[1];
        force_perp[0]   =   -forcex[0]*actin_network->get_direction(aindex[0])[1] + forcey[0]*actin_network->get_direction(aindex[0])[0];
        force_par[1]    =   forcex[1]*actin_network->get_direction(aindex[1])[0] + forcey[1]*actin_network->get_direction(aindex[1])[1];
        force_perp[1]   =   -forcex[1]*actin_network->get_direction(aindex[1])[1] + forcey[1]*actin_network->get_direction(aindex[1])[0];

        torque[0]       =   cross(hx[0]-actin_network->get_position(aindex[0])[0],hy[0]-actin_network->get_position(aindex[0])[1],forcex[0],forcey[0]);
        torque[1]       =   cross(hx[1]-actin_network->get_position(aindex[1])[0],hy[1]-actin_network->get_position(aindex[1])[1],forcex[1],forcey[1]);
        
        actin_network->update_forces(aindex[0],force_par[0],force_perp[0],torque[0]);
        actin_network->update_forces(aindex[1],force_par[1],force_perp[1],torque[1]);
    }
    
private:
    double hx[2],hy[2], phi, ld, stretch, forcex[2], forcey[2], torque[2], force_par[2],force_perp[2], lk;
    int aindex[2];
    std::string color;
    actin_ensemble *actin_network;

};

//link ensemble class
class link_ensemble
{
public:
    
    link_ensemble(){};
    ~link_ensemble(){};
    
    void link_walk()
    {
        for (int i=0; i<links.size(); i++) {
            
            links[i].step();
            links[i].actin_update();
        }
    }
    
    void link_write(std::ofstream& fout)
    {
        for (int i=0; i<links.size(); i++) {
            fout<<links[i].get_heads()[0]<<"\t"<<links[i].get_heads()[1]<<"\t"
                <<links[i].get_heads()[2]-links[i].get_heads()[0]<<"\t"<<
                links[i].get_heads()[3]-links[i].get_heads()[1]<<"\t"<<
                links[i].get_color()<<"\n";
        } 
    }
    
    void add_link(Link l)
    {
        links.push_back(l);
    }

private:
    std::vector<Link> links;
};

// Adds dead motors between monomers of an actin polymer to act as springs.
// Needs a motor_ensemble object to work, first 
void actin_ensemble::connect_polymers(motor_ensemble * mots, double link_length, double link_stiffness, 
        std::string link_color){
    int monomer_index;
    std::vector<double> motor_coords;
    for (int i = 0; i < mono_map.size(); i++){
        for (int j = 0; j < mono_map[i].size(); j++){
            monomer_index = mono_map[i][j];
            motor_coords = link_map[monomer_index];
            mots->add_motor( motor( motor_coords[0], motor_coords[1], motor_coords[2], 
                        link_length, this, 1, 1, monomer_index, monomer_index+1, fov[0], fov[1], 0, link_stiffness,
                        0, 0, 0, ld, visc, link_color) );
        }
    }
}
// Adds dead motors between monomers of an actin polymer to act as springs.
// Add additional dead motors between center of masses of monomers to account for bending energy
// Needs a motor_ensemble object to work, first 
void actin_ensemble::connect_polymers(link_ensemble * links, 
        double link_length, double link_stiffness, std::string link_color,
        double b_link_stiffness, std::string b_link_color){
    
    int mono1, mono2;
    double x1, x2, y1, y2, ang, link_length2, b_link_length1, b_link_length2;
    std::vector<double> motor_coords;
    for (int i = 0; i < mono_map.size(); i++){
        for (int j = 0; j < mono_map[i].size(); j++){
            
            // Join monomers within the polymer
            mono1 = mono_map[i][j];
            mono2 = mono1 + 1; //mono_map[i][j+1];
            
            link_length2 = 10; //dis_points(x1, x2, y1, y2) + 100 * eps;
            std::cout<<"DEBUG: motor angle = "<<ang<<"\n";
            links->add_link( Link( link_length, link_stiffness, this, mono1, mono2,  link_color) );
            
            // Add bending energy links
            x1      = network[mono1].getposcm()[0];
            y1      = network[mono1].getposcm()[1];
            x2      = network[mono2].getposcm()[0];
            y2      = network[mono2].getposcm()[1];
            b_link_length1 = network[mono1].get_length()/2 + network[mono2].get_length()/2 + link_length;
            b_link_length2 = dis_points(x1,y1,x2,y2);
//            std::cout<<"DEBUG: b_link_length2-b_link_length1 : "<<b_link_length2-b_link_length1<<"\n";

//            mots->add_motor( motor( (x2+x1)/2, (y2+y1)/2, motor_coords[2],//atan2(y2-y1,x2-x1), 
//                        b_link_length1, this, 1, 1, mono1, mono2, fov[0], fov[1], 0, b_link_stiffness,
//                        0, 0, 0, ld, visc, b_link_color) ); 
        
        }
    }
}








