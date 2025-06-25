#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "./contact-embed.h"
#include "view.h"
#include "tag.h"

#define LEVEL 12
#define RHOR 80

int num = 0; //video file 1 so its not annoying as fuck
int j = 3;
int k = 0;
double uemax = 0.5;
double t_end = 2.5 ;

const double H = 1 ;  // gas nozzle
const double R = 3;  // liquid inlet
const double RHO1 = 1;
const double SIGMA = 1.;


const double l = 50;				// * H
const double d = 5.9;					// * H
const double B = 2.9;					// * H
const double L = 14.9;				// * H
const double vertical_width = 0.9; 	// =H

const coord center = {9,0};

double MUl[4] = {2.e-2, 5.e-3, 5.e-3, 5.e-3}; //5.e-3
double MUg[4] = {1.e-1, 2.5e-2, 2.5e-2, 2.5e-2}; //2.5e-2
double Ul[4] = {2, 2, 1, 0.5}; //0.5
double Ug[4] = {40, 160, 160, 160/3}; //160
double theta0[3] = {70, 30, 90}; //70



u.n[left] = dirichlet (fabs(y) < (d/2 + 0.1)? Ul[j]:Ug[j]);
//u.n[left] = dirichlet (160);

p[left] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

u.t[embed] = dirichlet (0);
u.n[embed] = dirichlet (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);

u.n[top] = dirichlet (0);
u.n[bottom] = dirichlet (0);

f[right] = dirichlet(0);
f[left] = dirichlet (fabs(y) < (d/2 + 0.1)? 1: 0);
 
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

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
	double right_edge_end = H/2;

	double channel_dist = 3*H;

	coord left_side = {center.x - l/2,0};
	coord right_side = {center.x + l/2,0};

	coord vert_coord = {left_side.x + L + vertical_width/2,0};

	coord channel_position_top = {left_side.x + L/2, channel_dist + d/2 + B/2};
	coord channel_position_bot = {left_side.x + L/2,-channel_dist - d/2 - B/2};

	coord nozzle_end = {left_side.x + L + vertical_width + right_edge_full, 0};
	
	coord right_box_coord = {(right_side.x + nozzle_end.x)/2,center.y};
	double right_box = rectangle (x,y,right_box_coord, right_side.x - nozzle_end.x + 0.5, 200.5);

	double hor_line = rectangle (x,y, center, 200, d);

	double vert_line = rectangle (x,y, vert_coord, H * 0.9, 2*channel_dist + d +2*B); 					// the | in the +
	double cross = min (hor_line, vert_line); 								// this creates a cross
																			//
	double upper_channel = rectangle (x,y, channel_position_top, L, B ); 				// this is the upperchannel where the flow geos through
	double lower_channel = rectangle (x,y, channel_position_bot, L, B ); 				// this is the lowerchannel where the flow goes through
	double channels = min(upper_channel, lower_channel);

	double pre_triangles = min(cross, min(right_box, channels));

	double outer_triangle_height = 3;
	double outer_triangle_width = 3;

	coord outer_top_triangle_coords[] = {
		{nozzle_end.x - H/2 + outer_triangle_width, center.y + d/2 + outer_triangle_height},
		{nozzle_end.x - H/2, 										center.y + d/2},
		{nozzle_end.x - H/2 + outer_triangle_width,										center.y + d/2}
	};

	coord outer_bot_triangle_coords[] = {
		{outer_top_triangle_coords[0].x, center.y - d/2 - outer_triangle_height},
		{outer_top_triangle_coords[1].x, center.y - d/2},
		{outer_top_triangle_coords[2].x, center.y - d/2}
	};
	
	double outer_triangles = min(
					triangle(x,y, outer_top_triangle_coords[0], outer_top_triangle_coords[1], outer_top_triangle_coords[2]),
					triangle(x,y, outer_bot_triangle_coords[0], outer_bot_triangle_coords[1], outer_bot_triangle_coords[2]));


	coord inner_top_triangle_coords[] = {
		{vert_coord.x - vertical_width/2 + 3/2*vertical_width, center.y + d/2},//channel_position_top.y - B/2},
		{vert_coord.x - vertical_width/2 + 3/2*vertical_width, center.y + d/2 + 3*H},
		{vert_coord.x - vertical_width/2 - 2*H, center.y + d/2 +3*H}
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
	//return triangles;
	return -min(pre_triangles, triangles);
}


// Define solid/nozzle shape
static void solid_domain (scalar c, face vector cf) {
  solid(cs, fs, geometry(x,y));
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
	draw_vof("cs");
	c.refine = c.prolongation = fraction_refine;

}


int main (int argc, char * argv[]) {
	size (l);
  init_grid (1 << (LEVEL - 4));
  origin (-L*H, -(l*H/2));

  rho1 = RHO1;
  rho2 = rho1 / RHOR;
  f.sigma = SIGMA;

  const scalar c[] = theta0[0]*pi/180;
  contact_angle = c;

  TOLERANCE = 1e-4 [*];

  //for (j = 0; j <= 3; j++) {
  int j = 3;
  mu1 = MUl[j];
  mu2 = MUg[j];

    //run();
	//num=0;
  //}

  for (k = 0; k <= 2; k++) {
    const scalar cc[] = theta0[k]*pi/180;
    contact_angle = cc;
    run();
	num=0;
  }

}

event init (t = 0) {
  // mask(y > 5*H ? top: y < -5*H ? bottom : none);

  // Refine interfaces until they reach the maximum level
  astats ss;
  int ic = 0;
  do {
    ic++;
    lq_domain (f);
    solid_domain (cs, fs);
    ss = adapt_wavelet ({cs, f}, (double[]) {1.e-30, 1.e-30},
		        maxlevel = LEVEL , minlevel = LEVEL - 2);
  } while ((ss.nf || ss.nc) && ic < 1000);
  refine (cs[] > 0 && cs[] < 1 && level < LEVEL);

  // Reinitalize the solid and fluid/gas fields
  lq_domain (f);
  solid_domain (cs, fs);
//	fractions_cleanup (cs, fs);
  // Initialize velocity
  foreach() {
    foreach_dimension()
      u.x[] = 0;
  }
}

event logfile (i++; t <= t_end) {
	fflush(stdout);
  fprintf (stderr, "%d %g %d %d\n", i, t, j, k);

  char name[80];
  sprintf (name, "%d-%d-xh2.txt", k, j);
  FILE * fp2 = fopen (name, "a");
  foreach_point(-2*H, 0)
    fprintf (fp2, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp2);

  sprintf (name, "%d-%d-xh1.txt", k, j);
  FILE * fp1 = fopen (name, "a");
  foreach_point(-1*H, 0)
    fprintf (fp1, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp1);

  sprintf (name, "%d-%d-xh0.txt", k, j);
  FILE * fp0 = fopen (name, "a");
  foreach_point(0, 0)
    fprintf (fp0, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp0);
}

event interface (t += 0.3125; t <= t_end) {
  char name[80];
  sprintf (name, "%d-%d-int-%g.txt", k, j, t);
  FILE * fpi = fopen(name, "w");
  output_facets (f, fp = fpi);
  fclose (fpi);
}

event snapshot (t += 0.025) {
  scalar omega[];
  vorticity (u, omega);
  foreach()
    omega[] *= H/Ug[j];

  char name[80];
  sprintf (name, "%d-%d-%d-vort.png", k, j, num);
  FILE * fp1 = fopen (name, "w");
  view (fov=18.9,tx=-0.188,width = 1200, height = 1200);
  clear();
  draw_vof ("cs", "fs", filled = -1);
  draw_vof ("f", filled = 1, fc = {1,0,0}, lw = 5);
  squares ("omega", max = 10, min = -10);
  save (fp = fp1);
  fclose (fp1);

  scalar umag[];
  foreach()
    umag[] = f[] * sqrt(sq(u.x[]) + sq(u.y[])) / Ug[j];

  sprintf (name, "%d-%d-%d-velo.png", k, j, num);
  FILE * fp3 = fopen (name, "w");
  clear();
  draw_vof ("cs", "fs", filled = -1);
  draw_vof ("f");
  squares ("umag", min = 0);
  save (fp = fp3);
  fclose (fp3);

  sprintf (name, "%d-%d-%d-grid.png", k, j, num);
  FILE * fp4 = fopen (name, "w");
  view (fov=19.5,tx=-0.188, width = 1200, height = 1200);
  clear();
  draw_vof ("cs", "fs", filled = -1);
  draw_vof ("f");
  cells();
  save (fp = fp4);
  fclose (fp4);
}

event movie (t += 0.025; t <= t_end) {

  char name[80];
  sprintf (name, "%d-%d-%d-vof.png",k, j, num); // im taking out t so i can do fmpeg easier
  FILE * fp1 = fopen (name, "w");
  view (fov=18.9, tx=-0.188, width = 1200, height = 1200);
  draw_vof ("cs", "fs", filled = -1);
  draw_vof ("f");
  squares ("f", min = 0, max = 1);
  save (fp = fp1);
  fclose (fp1);
  num++;
}

event adapt (i++) {
  adapt_wavelet ({cs,f,u}, (double[]){1.e-15,1.e-4,uemax,uemax},
  maxlevel = LEVEL, minlevel = LEVEL - 4);
  unrefine (x > 26*H && level > 3);
}

event droplets (i+=10) {
	remove_droplets (f, minsize = 2);
}

//event end (t = t_end) {
//	system ("./sort_images.sh");
//}
