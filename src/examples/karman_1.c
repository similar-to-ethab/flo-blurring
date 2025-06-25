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

int main()
{
	

	L0 = 10. [1];
  	origin (-2, -2);
  	N = 512;
  	mu = muv;

  	/**
  	When using bview we can interactively control the Reynolds number
	and maximum level of refinement. */
  
  	display_control (Reynolds, 10, 1000);
  	display_control (maxlevel, 6, 12);
  
 	 run(); 
}

double rectangle (double x, double y, coord center, double width, double height) {

    double dx = fabs(x - center.x) - width / 2;
    double dy = fabs(y - center.y) - height / 2;

    // Max of dx and dy gives the signed distance
    return max(dx, dy);	
}


double geometry () {
	coord center = {0,0};
	double r1 = rectangle (1,1, center, 2.0, 1.0);
	double r2 = rectangle (0.5,0.5, center, 1.0, 0.5);

	double c = max(r1, -r2);
	return c;
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
/*
u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
*/
event init (t = 0)
{

  /**
  The domain is the intersection of a channel of width unity and a
  circle of diameter 0.125. */
  
  //solid (cs, fs, geometry());
  //geometry()
//solid(cs, fs, geometry());
solid (cs, fs, intersection (intersection (0.5 - 0.5, 0.5 + 0.5), sqrt(sq(x) + sq(y)) - 0.125/2.));
  /**
  We set the initial velocity field. */
  
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

event movies (i += 4; t <= 1.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  //foreach()
 //   m[] = cs[] - 0.5;
  output_ppm (omega, file = "vort.mp4", box = {{-4.0,-4.0},{4.0,4.0}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-4,-4},{4,4}},
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
