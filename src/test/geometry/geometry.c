// #include "embed.h"
#include "saint-venant.h"
//#include "run.h"
#include "view.h"
#define LEVEL 6 
#define L0 1.

//use rectange to make a rectange on the board
void rectange (coord center, double x_size, double y_size) {
	scalar cs[]; //the actual boundary
    face vector fs[]; //where it can flow i think?
}

int main(){
	//initialize grid
	size (L0);
	init_grid(1 << LEVEL);

	//start simulation
	run();
}

event init (t=0){
	//define characteristics of rectangle here
	//need center, x size, y size
	coord center = {0.5,0.5};
	double x_size = 1;
	double y_size = 1;

	//call rectange here	
	//rectange(center, x_size, y_size);
	//figure out how to visualize iti
}

event movie (t=1){
	clear();
	draw_vof("f");
	box();
	save("image.ppm");
}



