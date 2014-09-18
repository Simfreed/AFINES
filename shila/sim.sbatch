#!/bin/bash
# Shell options
shopt -s extglob
module load gnuplot
module load ffmpeg

# set the number of nodes to 2
#SBATCH --ntasks=1
#SBATCH --job-name=actomyosin_sim
#SBATCH --output=actomyosin_sim.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --account=pi-dinner
#SBATCH --cpus-per-task=1

###########################
# Begin work section      #
###########################

mf="makefile.sim"
exc="sim"
plotter="plot_sim.gp"
make -f $mf
kl=100 # pn/um
kb=1 # pn
nm=10
main_dr="simulation/kl${kl}_kb${kb}_nm$nm"

png_dr="$main_dr/png_stack"
txt_dr="$main_dr/txt_stack"
data_dr="$main_dr/data"

echo "deleting and recreating directory $main_dr"
rm -r $main_dr
mkdir $main_dr
mkdir "$png_dr"
mkdir "$txt_dr"
mkdir "$data_dr"

./$exc $nm $kl $kb $main_dr

pushd $txt_dr
echo "gnuplotting"
gnuplot "../../../$plotter"
echo "making movie"
ffmpeg -r 1/0.0075 -i time%01d.png -c:v libx264 -r 30 -pix_fmt yuv420p kl${kl}_nm${nm}.mp4
popd
mv $txt_dr/*.png $png_dr
mv $txt_dr/*.mp4 $data_dr

