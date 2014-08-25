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
#include "math.h"
#include "fstream"
#include "string"
#include "stdint.h"
#include "vector"
#include <map>
#include "time.h"

/* distances in microns, time in seconds, forces in pN */

#define pi 3.14159265358979323


#define xrange 50.0
#define yrange 50.0
#define xgrid 100.0
#define ygrid 100.0
#define tinit 0.0
#define tfinal 10 
#define dt 0.001
#define print_dt 100


#define temperature 0.004

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
        angle=atan((ycor-yp)/(xcor-xp));
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
    actin_ensemble(double density, double fovx, double fovy, int nx, int ny, double len, double vis, double link_len)
    {
        view=0.9;
        fov[0]=fovx;
        fov[1]=fovy;
        nq[0]=nx;
        nq[1]=ny;
        rho=density;
        av_vel=0;
        visc=vis;
        nmonomer_min = 10; //hard coded number of min/max monomers per filament
        nmonomer_max = 100;
        int nmonomer = (nmonomer_max - nmonomer_min)/2;
        npolymer=int(ceil(density*fov[0]*fov[1]) / nmonomer);
        ld=len;//rng_n(len,1.0);
        link_ld = link_len;
        std::cout<<"\nDEBUG: Number of filament:"<<npolymer<<"\n";
        
        std::vector<double> link_map_value; // xcm, ycm, angle
        double * rod_end_pts;
        double  xcm, ycm, theta;
        for (int i=0; i<npolymer; i++) {
            theta=rng(0,2*pi);
            nmonomer = (int) rng(nmonomer_min, nmonomer_max);
            //the start of the polymer: 
            network.push_back(actin(rng(-0.5*(view*fovx-ld),0.5*(view*fovx-ld)), rng(-0.5*(view*fovy-ld),0.5*(view*fovy-ld)),
                        theta,ld,fov[0],fov[1],nq[0],nq[1],visc));
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
//                    std::cout<<"\nDEBUG:"<<j+1<<"th monomer of "<<i<<"th polymer outside field of view; stopped building polymer\n";
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
   
    std::map<double, double> get_angle_correlation(int polymer_index)
    {
        //wikipedia says < cos(theta) > = Exp[-L/P]
        //where L is the length of the polymer and P is the persistence length
        //and theta is the angle between subsequent rods
        // This function will return a map of x --> <cos theta>
        // where x is a distance and <cos theta> is the correlation between cosines of the angles between subsequent
        // rods. 
        double phi1, phi2;
        double sum = 0;
        std::map<double, double> corr;
        for(int i = 0; i < mono_map[polymer_index].size() - 1; i++){
            phi1 = network.at(mono_map[polymer_index][i]).getpos()[2];
            phi2 = network.at(mono_map[polymer_index][i + 1]).getpos()[2];
            sum += cos(phi2 - phi1);
            corr[ (i+1) * alength] = sum/(i+1);
        }
        return corr;

    }
    
    std::map<int, std::map<double, double> > get_all_angle_correlations()
    {
        std::map<int, std::map<double, double> > corrs;
        
        for(int i = 0; i < mono_map.size(); i++){
            corrs[i] = this->get_angle_correlation(i);
        }

        return corrs;
    }

    std::map<int, std::vector<double> > * get_link_map(){
        return &link_map;
    }
    
    std::map<int, std::vector<int> > * get_mono_map(){
        return &mono_map;
    }

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
		pos_actin[0]=0; // I don't think this variable get's ACCESSED anywhere 
        pos_a_end[0]=0; //distance from pointy end -- by default 0
        pos_a_end[1]=0;
        if (state0)
            pos_a_end[0] = dis_points(hx[0],hy[0],actin_network->get_ends(aindex[0])[2],actin_network->get_ends(aindex[0])[3]);
		if (state1)
            pos_a_end[1] = dis_points(hx[1],hy[1],actin_network->get_ends(aindex[1])[2],actin_network->get_ends(aindex[1])[3]);

        pos_actin[1]=0;
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
                            mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
                        }
                        
                        pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
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
            mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
		}
		else if (state[0]==0 || state[1]==0) {
			int hd=state[0];
            
			xm[hd]=hx[hd]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hx[hd]-hx[pr(hd)]+pow(-1,hd)*mld*cos(mphi));
			ym[hd]=hy[hd]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hy[hd]-hy[pr(hd)]+pow(-1,hd)*mld*sin(mphi));
            xm[pr(hd)]=hx[pr(hd)];
            ym[pr(hd)]=hy[pr(hd)];
            reflect(xm[0],xm[1],ym[0],ym[1]);
            mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
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
			pos_actin[hd]=0;
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
			pos_actin[0]=0;
            pos_a_end[0]=0;
            pos_temp=pos_a_end[1]+dt*vs;
			move_end_detach(1,vs,pos_temp);
		}
		else if (event(offrate[1],dt)==1) {
            state[1]=0;
            aindex[1]=-1;
            hx[1]=hx[0]+mld*cos(mphi);
            hy[1]=hy[0]+mld*sin(mphi);
            pos_actin[1]=0;
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
            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_position(aindex[0])[0],actin_network->get_position(aindex[0])[1]);
            hx[1]=actin_network->get_ends(aindex[1])[2]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
            hy[1]=actin_network->get_ends(aindex[1])[3]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_position(aindex[1])[0],actin_network->get_position(aindex[1])[1]);
            mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
        }
        else if(state[0]==1 && state[1]==0)
        {
            hx[0]=actin_network->get_ends(aindex[0])[2]-pos_a_end[0]*actin_network->get_direction(aindex[0])[0];
            hy[0]=actin_network->get_ends(aindex[0])[3]-pos_a_end[0]*actin_network->get_direction(aindex[0])[1];
            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_position(aindex[0])[0],actin_network->get_position(aindex[0])[1]);
            mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
        }
        else if(state[0]==0 && state[1]==1)
        {
            hx[1]=actin_network->get_ends(aindex[1])[2]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
            hy[1]=actin_network->get_ends(aindex[1])[3]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_position(aindex[1])[0],actin_network->get_position(aindex[1])[1]);
            mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
        }
        else
        {
            return;
        }
        
        
    }
    
private:
    double hx[2],hy[2],vhx[2], xm[2], ym[2], vhy[2], mphi,mld,mobility,fm[2],vm[2],offrate[2], onrate, stretch,forcex[2],forcey[2],torque[2], force_par[2],force_perp[2],vs, pos_temp;
    int state[2], aindex[2];
    double dm,fmax,mk, kon, koff, kend;
    std::map<int, double> dist;
    std::string color;
    actin_ensemble *actin_network;
	double pos_actin[2], pos_a_end[2], fov[2];
    
    
    inline void move_end_detach(int hd, double speed, double pos)
    {
        if (pos>=actin_network->get_alength(aindex[hd])) { //not sure why only "greater than or equal"
            if (event(kend,dt)==1) {
                //std::cout<<"\nThe new myosin position of head "<<hd<<" is OFF the actin filament AND detaching";
                state[hd]=0;
                aindex[hd]=-1;
                pos_actin[hd]=0;
                pos_a_end[hd]=0;
                hx[hd]=hx[pr(hd)]-pow(-1,hd)*mld*cos(mphi);
                hy[hd]=hy[pr(hd)]-pow(-1,hd)*mld*sin(mphi);
            }
            else {
                //std::cout<<"\nThe new myosin position of head "<<hd<<" is OFF the actin filament BUT not detaching";
                hx[hd]=actin_network->get_ends(aindex[hd])[2]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
                hy[hd]=actin_network->get_ends(aindex[hd])[3]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
                mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
                pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
            }
        }
        else {
            //std::cout<<"\nThe new myosin position of head "<<hd<<" is "<<pos<<" away from the pointy end of the actin filament "; //still on the actin filament";
            pos_a_end[hd]=pos;
            hx[hd]=actin_network->get_ends(aindex[hd])[2]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
            hy[hd]=actin_network->get_ends(aindex[hd])[3]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
            mphi=atan((hy[1]-hy[0])/(hx[1]-hx[0]));
            pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
        }
        
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
        color = "1";//"green"; 
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







