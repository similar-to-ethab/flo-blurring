scalar contact_angle;

#include "ghost-cell-IBM/src/ibm-gcm.h"
#include "ghost-cell-IBM/src/my-centered.h"
#include "ghost-cell-IBM/src/ibm-gcm-events.h"
//#include "ghost-cell-IBM/src/ibm-gcm-vof.h"
#include "ghost-cell-IBM/src/my-two-phase.h"
#include "ghost-cell-IBM/src/my-tension.h"
#include "ghost-cell-IBM/src/contact-ibm3D.h"
#include "view.h"
#include "tag.h"

//for make_blur (file management)
#include <sys/types.h>   // sometimes needed first on BSD/macOS
#include <sys/stat.h>    // struct stat, stat(), mkdir()
#include <stdio.h>       // sprintf()



#define LEVEL 10

int num = 0; //video file 1 so its not annoying as fuck
int j = 0;
int k = 0;
double uemax = 0.5;
double t_end = 2.5;

const double H = 1;  // gas nozzle
const double R = 4;  // liquid inlet
const double rho_l = 828.0730897;
const double rho_g = 1;

const double SIGMA = 1.;


const double l = 50;				// * H
const double d = R;					// * H
const double B = R;					// * H
const double L = 20.;				// * H
const double vertical_width = H; 	// =H

const coord center = {9,0};

double MUl = 0.001442976;
double MUg = 0.076097998;
double Ul = 0.652787709;
double Ug[7] = {135.1389837, 168.9237297, 202.7084756, 236.4932215, 270.2779675, 304.0627134, 337.8474593};
double theta0[1] = {70};

const double photo_time_step = 0.025;

///////Velocity boundary conditions ===

// Left boundary
u.n[left] = dirichlet(fabs(y) < (d/2 + 0.1) ? Ul : Ug[j]);
u.t[left] = neumann(0);

// Right boundary
u.n[right] = neumann(0);
u.t[right] = neumann(0);

// Top boundary
u.n[top] = neumann(0);
u.t[top] = neumann(0);

// Bottom boundary
u.n[bottom] = neumann(0);
u.t[bottom] = neumann(0);

// Embedded boundaries (e.g., obstacles)
u_x_ibm_dirichlet(0)
u_y_ibm_dirichlet(0)

//Pressure boundary conditions

// Right boundary
p[right] = dirichlet(0);
pf[right] = dirichlet(0); // Filtered pressure

// Top boundary
p[top] = dirichlet(0);
pf[top] = dirichlet(0); // Filtered pressure
// Bottom boundary
p[bottom] = dirichlet(0);
pf[bottom] = dirichlet(0);

///////f boundary conditions

// Left boundary
f[left] = dirichlet(fabs(y) < (d/2 + 0.1) ? 1 : 0);

// Right boundary
f[right] = dirichlet(0);

int make_blur()
{
    char path[64];
    sprintf(path, "blur-%d", LEVEL);

    struct stat st;
        mkdir(path, 0700);
    return 0;
}

double rectangle (double x, double  y, coord center, double w, double h) {
	double Px_Plus = x - w/2. - center.x;
	double Px_Minus = x + w/2. - center.x;
	double Px = max (Px_Plus, -Px_Minus);
	double Py_Plus = y  - h/2. - center.y;
	double Py_Minus = y + h/2. - center.y;
	double Py = max (Py_Plus, -Py_Minus);

	return max (Px, Py);
}


// Compute the signed distance from point (x,y) to the line segment from point 'a' to point 'b'
static double sd_segment(double x, double y, coord a, coord b) {
    double vx = b.x - a.x;
    double vy = b.y - a.y;
    double wx = x - a.x;
    double wy = y - a.y;

    double c1 = vx * wx + vy * wy;
    if (c1 <= 0.0)
        return sqrt(wx * wx + wy * wy);

    double c2 = vx * vx + vy * vy;
    if (c2 <= c1) {
        double dx = x - b.x;
        double dy = y - b.y;
        return sqrt(dx * dx + dy * dy);
    }

    double t = c1 / c2;
    double projx = a.x + t * vx;
    double projy = a.y + t * vy;
    double dx = x - projx;
    double dy = y - projy;
    return sqrt(dx * dx + dy * dy);
}

// Helper function to compute the “signed area” (2D cross product) for three points.
double triangle_sign(coord p, coord q, coord r) {
    return (p.x - r.x) * (q.y - r.y) - (q.x - r.x) * (p.y - r.y);
}

// Test whether a point p is inside the triangle defined by vertices a, b, and c.
int point_in_triangle(coord p, coord a, coord b, coord c) {
    double d1 = triangle_sign(p, a, b);
    double d2 = triangle_sign(p, b, c);
    double d3 = triangle_sign(p, c, a);

    // The point is inside if it is on the same side of each edge.
    int has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    int has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

// SDF for a triangle defined by the vertices 'a', 'b', and 'c'.
// Returns the distance from (x,y) to the triangle: negative inside, positive outside.
double triangle(double x, double y, coord a, coord b, coord c) {
    // Compute the distance from the point to each edge of the triangle.
    double d_ab = sd_segment(x, y, a, b);
    double d_bc = sd_segment(x, y, b, c);
    double d_ca = sd_segment(x, y, c, a);

    // picks smallest of these distances.
    double d = d_ab;
    if (d_bc < d) d = d_bc;
    if (d_ca < d) d = d_ca;

    // If the point lies inside the triangle, return -d (negative distance).
    coord p = { x, y };
    if (point_in_triangle(p, a, b, c))
        return -d;
    else
        return d;
}


double geometry (double x, double y) {

	double right_edge_full = H;

	double channel_dist = d;

	coord left_side = {center.x - l/2,0};
	coord right_side = {center.x + l/2,0};

	coord vert_coord = {left_side.x + L + vertical_width/2,0};

	coord channel_position_top = {left_side.x + L/2 - 50, channel_dist + d/2 + B/2};
	coord channel_position_bot = {left_side.x + L/2 - 50,-channel_dist - d/2 - B/2};

	coord nozzle_end = {left_side.x + L + vertical_width + right_edge_full, 0};

	coord right_box_coord = {(right_side.x + nozzle_end.x)/2 + 5,center.y};
	double right_box = rectangle (x,y,right_box_coord, right_side.x - nozzle_end.x + 0.5 + 10, 200.5);

	double hor_line = rectangle (x,y, center, 200, d);

	double vert_line = rectangle (x,y, vert_coord, H, 2*channel_dist + d +2*B); 					// the | in the +
	double cross = min (hor_line, vert_line); 								// this creates a cross
																			//
	double upper_channel = rectangle (x,y, channel_position_top, L + 100, B ); 				// this is the upperchannel where the flow geos through
	double lower_channel = rectangle (x,y, channel_position_bot, L + 100, B ); 				// this is the lowerchannel where the flow goes through
	double channels = min(upper_channel, lower_channel);

	double pre_triangles = min(cross, min(right_box, channels));

	double outer_triangle_height = 3;
	double outer_triangle_width = 3;

	coord outer_top_triangle_coords[] = {
		{nozzle_end.x - d/2 + outer_triangle_width, center.y + d/2 + outer_triangle_height - H},
		{nozzle_end.x - d/2, 										center.y + d/2 - H},
		{nozzle_end.x - d/2 + outer_triangle_width,										center.y + d/2 - H}
	};

	coord outer_bot_triangle_coords[] = {
		{outer_top_triangle_coords[0].x, center.y - d/2 - outer_triangle_height + H},
		{outer_top_triangle_coords[1].x, center.y - d/2 + H},
		{outer_top_triangle_coords[2].x, center.y - d/2 + H}
	};

	double outer_triangles = min(
					triangle(x,y, outer_top_triangle_coords[0], outer_top_triangle_coords[1], outer_top_triangle_coords[2]),
					triangle(x,y, outer_bot_triangle_coords[0], outer_bot_triangle_coords[1], outer_bot_triangle_coords[2]));


	coord inner_top_triangle_coords[] = {
		{vert_coord.x - vertical_width/2 + 1/2*vertical_width, center.y + d/2 + H },//channel_position_top.y - B/2},
		{vert_coord.x - vertical_width/2 + 1/2*vertical_width, center.y + d + 2*H + H},
		{vert_coord.x - vertical_width/2 - channel_dist, center.y + d +2*H + H}
	};

	coord inner_bot_triangle_coords[] = {
		{inner_top_triangle_coords[0].x, -inner_top_triangle_coords[0].y},
		{inner_top_triangle_coords[1].x, -inner_top_triangle_coords[1].y},
		{inner_top_triangle_coords[2].x, -inner_top_triangle_coords[2].y}
	};
	double inner_triangles = min (
					triangle(x,y, inner_top_triangle_coords[2],inner_top_triangle_coords[0],inner_top_triangle_coords[1]),
					triangle(x,y, inner_bot_triangle_coords[0],inner_bot_triangle_coords[1],inner_bot_triangle_coords[2]));

	double triangles = min (inner_triangles, outer_triangles);
	return -min(pre_triangles, triangles);
}


// Define solid/nozzle shape
static void solid_domain (scalar c, face vector cf) {
  solid(ibm, ibmf, geometry(x,y));
}

// Define which domain parts are liquid vs gas
static void lq_domain (scalar c) {
	vertex scalar phi[];
	face vector u[];
	foreach_vertex () {
    	if (x < (center.x - l/2 + L - vertical_width)  && (fabs(y) <= d/2 + 0.1)){
				phi[]=1.;

		}
		else {
				phi[]=-1.;

		}
			// if (x < -H && fabs(y) < R)
  		/*if (x < -H)
    		phi[] = 1.;
    	else
    		phi[] = -1.;
			*/
  	}
	boundary ({phi});
  	fractions (phi, c);
	draw_vof("ibm");
	c.refine = c.prolongation = fraction_refine;

}


int main (int argc, char * argv[]) {
	size (l);
  init_grid (1 << (LEVEL - 4));
  origin (-L*H, -(l*H/2));

  rho1 = rho_l;
  rho2 = rho_g;
  f.sigma = SIGMA;

  const scalar c[] = theta0[0]*pi/180;
  contact_angle = c;

  TOLERANCE = 1e-4 [*];

	mu1 = MUl;
	mu2 = MUg;

    const scalar cc[] = theta0[0]*pi/180;
    contact_angle = cc;

	for (j = 0; j <= 6; j++) {
		run();
		num=0;
  	}

}

event init (t = 0) {
	make_blur();
  // mask(y > 5*H ? top: y < -5*H ? bottom : none);

  // Refine interfaces until they reach the maximum level
  astats ss;
  int ic = 0;
  do {
    ic++;
    lq_domain (f);
    solid_domain (ibm, ibmf);
    ss = adapt_wavelet ({ibm, f}, (double[]) {1.e-30, 1.e-30},
		        maxlevel = LEVEL , minlevel = LEVEL - 2);
  } while ((ss.nf || ss.nc) && ic < 1000);
  refine (ibm[] > 0 && ibm[] < 1 && level < LEVEL);

  // Reinitalize the solid and fluid/gas fields
  lq_domain (f);
  solid_domain (ibm, ibmf);
//	fractions_cleanup (ibm, ibmf);
  // Initialize velocity
  foreach() {
    foreach_dimension()
      u.x[] = 0;
  }
}

event logfile (i++; t <= t_end) {
  foreach_face()
    if (!fm.x[])
	  uf.x[] = 0;
		fflush(stdout);
  fprintf (stderr, "%d %g %d %d\n", i, t, j, k);

  char name[80];
  sprintf (name, "blur-%d/%d-%d-xh2.txt", LEVEL, k, j);
  FILE * fp2 = fopen (name, "a");
  foreach_point(-2*H, 0)
    fprintf (fp2, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp2);

  sprintf (name, "blur-%d/%d-%d-xh1.txt", LEVEL, k, j);
  FILE * fp1 = fopen (name, "a");
  foreach_point(-1*H, 0)
    fprintf (fp1, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp1);

  sprintf (name, "blur-%d/%d-%d-xh0.txt", LEVEL, k, j);
  FILE * fp0 = fopen (name, "a");
  foreach_point(0, 0)
    fprintf (fp0, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp0);
}

event interface (t += 0.3125; t <= t_end) {
  char name[80];
  sprintf (name, "blur-%d/%d-%d-int-%g.txt", LEVEL, k, j, t);
  FILE * fpi = fopen(name, "w");
  output_facets (f, fp = fpi);
  fclose (fpi);
}



event snapshot (t += photo_time_step) {
  scalar omega[];
  vorticity (u, omega);
  foreach()
    omega[] *= H/Ug[j];

  char name[80];



  sprintf (name, "blur-%d/%d-%d-%d-vort.png", LEVEL, j, num, LEVEL);
  FILE * fp1 = fopen (name, "w");
  view (fov=18.9,tx=-0.188,width = 1200, height = 1200);
  clear();
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f", filled = 1, fc = {1,0,0}, lw = 5);
  squares ("omega", max = 10, min = -10);
  save (fp = fp1);
  fclose (fp1);

  scalar umag[];
  foreach()
    umag[] = f[] * sqrt(sq(u.x[]) + sq(u.y[])) / Ug[j];

  sprintf (name, "blur-%d/%d-%d-%d-velo.png", LEVEL, j, num, LEVEL);
  FILE * fp3 = fopen (name, "w");
  clear();
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f");
  squares ("umag", min = 0);
  save (fp = fp3);
  fclose (fp3);

  sprintf (name, "blur-%d/%d-%d-%d-grid.png", LEVEL, j, num, LEVEL);
  FILE * fp4 = fopen (name, "w");
  view (fov=19.5,tx=-0.188, width = 1200, height = 1200);
  clear();
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f");
  cells();
  save (fp = fp4);
  fclose (fp4);
}

event movie (t += 0.025; t <= t_end) {

  char name[80];
  sprintf (name, "blur-%d/%d-%d-%d-vof.png", LEVEL, j, num, LEVEL); // im taking out t so i can do fmpeg easier
  FILE * fp1 = fopen (name, "w");
  view (fov=18.9, tx=-0.188, width = 1200, height = 1200);
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f");
  squares ("f", min = 0, max = 1);
  save (fp = fp1);
  fclose (fp1);
  num++;
}

event adapt (i++) {
  adapt_wavelet ({ibm,f,u}, (double[]){1.e-15,1.e-4,uemax,uemax},
  maxlevel = LEVEL, minlevel = LEVEL - 4);
  unrefine (x > 24*H && level > 3);
}

event droplets (i+=10) {
	remove_droplets (f, minsize = 3);
}
/*
event dumps(i += 50; t <= t_end){
	scalar * dumplist = {u,p,f,ibm};
	char name[80];
	sprintf(name, "dumps/dump-%d-%d-%g", j,k,t);
	dump (file = name, list = dumplist);

}
*/
