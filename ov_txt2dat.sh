#!/bin/bash 

for value in {1,2,5,10,20,50,100,200,500,1000}
do 
    links="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/ov/$value links" 
    actins="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/ov/$value actins"
    amotors="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/ov/$value amotors"
    pmotors="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/ov/$value pmotors"

    $links
    $actins
    $amotors
    $pmotors 

done 
