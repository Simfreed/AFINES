#!/bin/bash
# Shell options
shopt -s extglob
module load gnuplot
module load ffmpeg

#SBATCH --job-name=network-small
#SBATCH --output=log/network-small.out
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --account=pi-dinner
#SBATCH --cpus-per-task=1

###########################
# Begin work section      #
###########################

make network
exc="bin/nt"
plotter="plot/network.gp"

nm=50
np=10
len=0.2
rho_mot=0.2
tf=200
kl=100 # pn/um
xr=50
yr=50
main_dr="out/network/np10_nm50"

png_dr="$main_dr/png_stack"
txt_dr="$main_dr/txt_stack"
data_dr="$main_dr/data"

echo "deleting and recreating directory $main_dr"
rm -r $main_dr
mkdir $main_dr
mkdir "$png_dr"
mkdir "$txt_dr"
mkdir "$data_dr"

echo "./$exc $nm $np $len $rho_mot $tf $kl $xr $yr $main_dr"
./$exc $nm $np $len $rho_mot $tf $kl $xr $yr $main_dr

pushd $txt_dr
echo "gnuplotting"
gnuplot "../../../../$plotter"
echo "making movie"
ffmpeg -r 1/0.0075 -i time%01d.png -c:v libx264 -r 30 -pix_fmt yuv420p network.mp4
popd
mv $txt_dr/*.png $png_dr
mv $txt_dr/*.mp4 $data_dr

