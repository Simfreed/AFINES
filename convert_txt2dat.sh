#!/bin/bash 

for value in {1,2,5,10,20,50,100,200,500}
do 
    links="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/exv_h/$value links" 
    actins="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/exv_h/$value actins"
    amotors="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/exv_h/$value amotors"
    pmotors="python versatile_framework_paper/python/txt2dat.py /project/weare-dinner/ibunge/exv_h/$value pmotors"

    $links
    $actins
    $amotors
    $pmotors

done 
