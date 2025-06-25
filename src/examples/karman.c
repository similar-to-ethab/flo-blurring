/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

An example of 2D viscous flow around a simple solid boundary. Fluid is
injected to the left of a channel bounded by solid walls with a slip
boundary condition. A passive tracer is injected in the bottom half of
the inlet.

![Animation of the vorticity field.](karman/vort.mp4)(loop)

![Animation of the tracer field.](karman/f.mp4)(loop)

We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer *f*. */

#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];


/**
The domain is eight units long, centered vertically. */

//creates a rectangular plane
double rectangle (double x, double y, coord center, double x_size, double y_size) {

	double Px_Plus = x - x_size/2. - center.x;
	double Px_Minus = x + x_size/2. - center.x;
	double Px = max (Px_Plus, -Px_Minus);

	double Py_Plus = y  - y_size/2. - center.y;
	double Py_Minus = y + y_size/2. - center.y;
	double Py = max (Py_Plus, -Py_Minus);

	double rect = max(Px,Py);

	return rect;
}

double geometry (double x, double y) {
	// def coord variables here
	coord center = {3.5,0};

	coord top_left= {-0.5,4};
	coord top_right = {7.5,4};

	coord bot_left = {-0.5,-4};
	coord bot_right = {7.5,-4};

	// define size of fluid channels here
	double main_channel_size = 2;
	double small_channel_size = 1;

	// big rectangle that covers entire grid
	double rectMax = rectangle (x,y, center, 8,8);

	// rectangles that cover each axis at a size of 2
	double rectXaxis = rectangle (x,y, center, 8, main_channel_size);
		// the y-axis will be different because I need to account for the smaller channels so there is no "overlap"
	double yHeight = 2 * (center.y + (top_left.y+center.y)/2 + small_channel_size/2);
	double rectYaxis = rectangle (x,y, center, main_channel_size,yHeight);

	// combines the two axis rectangles into one, shaped like a '+'
	double rectAxis = min(rectXaxis,rectYaxis);

	// takes the difference between the big one and the '+'
	double rectCorners = max(rectMax,-rectAxis);

	// add in the smaller channels
	coord top_channel_center = {(top_left.x + center.x)/2, (top_left.y + center.y)/2};
	double top_channel = rectangle (x,y, top_channel_center, 4 /*+ main_channel_size/2*/, small_channel_size);

	coord bot_channel_center = {(bot_left.x + center.x)/2, (bot_left.y +center.y)/2};
	double bot_channel = rectangle (x,y, bot_channel_center, 4 /*+ main_channel_size/2*/, small_channel_size);
	// combine them
	double small_channel = min (top_channel, bot_channel);
	// adding it all together
	return max(rectCorners,-small_channel);
}

int main()
{
  L0 = 8. [1];
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;

  /**
  When using bview we can interactively control the Reynolds number
  and maximum level of refinement. */

  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);
  run();
}

/**
We set a constant viscosity based on the Reynolds number, the cylinder
diameter $D$ and the inflow velocity $U0$. */

double D = 0.125, U0 = 1.;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

/**
The fluid is injected on the left boundary with velocity $U0$. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{



	solid (cs, fs, geometry(x,y));
	FILE * jp = stdout;
  	int num = 0;
	foreach(){
		if (num%16==0){
			fprintf(jp,"%g %g %g\n",x,y, geometry(x,y));

		}
	num++;
  }
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the vorticity and tracer fields... */

event movies (i += 4; t <= 20.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-4},{7.5,4}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-0.5,-4},{7.5,4}},
	      linear = false, min = 0, max = 1, mask = m);
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
