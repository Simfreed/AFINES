#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include "math.h"
#include <list>

class thing{
    private: 
        std::vector<int *> container;
    public:
        void append(int * i){
            container.push_back(i);
            delete i;
        }
        int geti(int index){
            return *(container[index]);
        }
        void clear(){
            int s = container.size();
            std::cout<<"container.size() = "<<s<<"\n";
            for (int i = 0; i<s; i++)
                container.pop_back();
        }
};

class thingmanipulator{
    public:
        void putInSquares(thing * t){
            //int * j;
            for(int i=0; i<10; i++){
                //int * j = new int(i*i);
                t->append(new int(i*i));
                //delete j; 
                //free(j);
            }
        }
        /*double* getsomedata()
        {
            double dat[3];
            double *datptr;
            dat[0] = 333;
            dat[1] = 111;
            dat[2] = 222;
            datptr = dat;
            return dat;
        }*/

};
int main()
{

    int * j;
    thing p;
    thingmanipulator tm;
    /*   thing * pptr = &p;
  
    for(int i=0; i<10; i++){
        j = new int(i*i);
        pptr->append(j);
    }
    */
    tm.putInSquares(&p);
    for(int i=0; i<10; i++){
        std::cout<<"thing["<<i<<"] = "<<p.geti(i)<<"\n";
    }
/*    std::cout<<"dat[0] = "<<tm.getsomedata()[0]<<"\n";
    std::cout<<"dat[1] = "<<tm.getsomedata()[1]<<"\n";
    std::cout<<"dat[2] = "<<tm.getsomedata()[2]<<"\n";
  */  
//    p.clear();
    
    return 0; 
};
