#!/bin/bash
# Shell options
shopt -s extglob

###########################
# Begin work section      #
###########################

make
kl=100 # pn/um
nm=10
main_dr="persistence_length/kl${kl}_nm$nm"

png_dr="$main_dr/png_stack"
txt_dr="$main_dr/txt_stack"
data_dr="$main_dr/data"

echo "deleting and recreating directory $main_dr"
rm -r $main_dr
mkdir $main_dr
mkdir "$png_dr"
mkdir "$txt_dr"
mkdir "$data_dr"

echo "./per $nm $kl $main_dr"
./per $nm $kl $main_dr

pushd $txt_dr
echo "gnuplotting"
gnuplot "../../../plotfiles.gp"
echo "making movie"
ffmpeg -r 1/0.0075 -i time%01d.png -c:v libx264 -r 30 -pix_fmt yuv420p kl${kl}_nm${nm}.mp4
popd
mv $txt_dr/*.png $png_dr
mv $txt_dr/*.mp4 $data_dr

