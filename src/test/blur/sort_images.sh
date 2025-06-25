#!/bin/bash



cd $BLUR/blur/blur-$1


frame_rate=10;
image_types="velo vof grid vort"
echo $image_types;

rm -r image-0-0
rm -r image-0-1
rm -r image-0-2
rm -r image-0-3

rm -r image-1-0
rm -r iamge-1-1
rm -r image-1-2
rm -r image-1-3

rm -r image-2-0
rm -r image-2-1
rm -r image-2-2
rm -r image-2-3


#for f in 
# 0-i - regex - 0-*-

	for (( i=0; i<=3; i++)); do
		for type in $image_types; do
			while IFS= read -r -d $'\0' file; do
				local_name=$(basename "$file");
				IFS=- read iteration image_num level rest <<< "$local_name"
				padded=$(printf "%04d" "$image_num")
				echo $padded
				echo "image-$i/$type/$iteration-$image_num-$padded-$rest"
				mkdir -p "image-$i/$type"
				magick "$file" "image-$i/$type/$iteration-$image_num-$padded-$rest"	
				echo "Found file $file";
			
			done < <(find $BLUR/blur/blur-$1 -name "$i-*-$type.png" -print0)
			echo "separate";
			rm "video-$i-$type-$1.mp4";
			ffmpeg -framerate $frame_rate -pattern_type glob -i "$BLUR/blur/image-$i/$type/*-$type.png" -c:v libx264 -pix_fmt yuv420p "video-$i-$type-$1.mp4" -loglevel debug
		done
	done




cd ../
