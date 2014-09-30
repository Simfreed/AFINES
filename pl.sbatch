#!/bin/bash
# Shell options
shopt -s extglob
module load gnuplot
module load ffmpeg

# set the number of nodes to 2
#SBATCH --ntasks=1
#SBATCH --job-name=persistence_length
#SBATCH --output=log/pl.out
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --account=pi-dinner
#SBATCH --cpus-per-task=1

###########################
# Begin work section      #
###########################

make persistence_length
exc="bin/pl"
plotter="plot/pl.gp"
kl=100 # pn/um
nm=50
main_dr="out/pl/kl${kl}_nm$nm"

png_dr="$main_dr/png_stack"
txt_dr="$main_dr/txt_stack"
data_dr="$main_dr/data"

echo "deleting and recreating directory $main_dr"
rm -r $main_dr
mkdir $main_dr
mkdir "$png_dr"
mkdir "$txt_dr"
mkdir "$data_dr"

echo "./$exc $nm $kl $main_dr"
./$exc $nm $kl $main_dr

pushd $txt_dr
echo "gnuplotting"
gnuplot "../../../../$plotter"
echo "making movie"
ffmpeg -r 1/0.0075 -i time%01d.png -c:v libx264 -r 30 -pix_fmt yuv420p kl${kl}_nm${nm}.mp4
popd
mv $txt_dr/*.png $png_dr
mv $txt_dr/*.mp4 $data_dr

