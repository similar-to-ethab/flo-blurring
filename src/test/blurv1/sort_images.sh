#!/bin/bash



cd $BLUR/blur


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
for (( j=0; j<=2; j++)); do
	for (( i=0; i<=3; i++)); do
		for type in $image_types; do
			while IFS= read -r -d $'\0' file; do
				local_name=$(basename "$file");
				IFS=- read blank_num i_num image_num rest <<< "$local_name"
				padded=$(printf "%04d" "$image_num")
				echo $padded
				echo "image-$j-$i/$type/$blank_num-$i_num-$padded-$rest"
				mkdir -p "image-$j-$i/$type"
				magick "$file" "image-$j-$i/$type/$blank_num-$i_num-$padded-$rest"	
				echo "Found file $file";
			
			done < <(find $BLUR/blur -name "$j-$i-*-$type.png" -print0)
			echo "separate";
			rm "video-$j-$i-$type.mp4";
			ffmpeg -framerate $frame_rate -pattern_type glob -i "$BLUR/blur/image-$j-$i/$type/*-$type.png" -c:v libx264 -pix_fmt yuv420p "video-$j-$i-$type.mp4" -loglevel debug
		done
	done
done



cd ../
