#!/bin/bash


cd $BLUR/blur

#rm -rf frame_file
#mkdir frame_file

#for f in *.png; do magick "$f" "frame_file/${f%.png}.png"; done

ffmpeg -framerate 10 -pattern_type glob -i "$BLUR/blur/frame_file/*-vof.png" -c:v libx264 -pix_fmt yuv420p vof.mp4 -loglevel debug

cd $BLUR/blur

rm -rf velo_file
mkdir velo_file

#while read -r p; do
#    magick "$p" "velo_file/png-${p%.png}.png"
#done < velo.txt

#ffmpeg -framerate 15 -f concat -i "$BLUR/blur/velo.txt" -c:v libx264 -pix_fmt yuv420p velo.mp4 -loglevel debug


cd ../

