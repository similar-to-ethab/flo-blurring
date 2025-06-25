#include "ibm/ibm-gcm1.h"
#include "ibm/my-centered.h"
#include "ibm/ibm-gcm-events.h"
#include "ibm/my-two-phase.h"
#include "ibm/my-tension.h"
#include "ibm/contact-ibm3D.h"
#include "tag.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#include "profiling.h"

#define cylinder(a, b, r) (sq(a) + sq(b) - sq(r)) // or cone if x = r
#define cone(a, b, r, angle) (((sq(a) + sq(b))/sq(tan(angle))) - sq(r))

// Geometry constants
const double H = 1.;       // gas nozzle outlet
const double B = 3.;       // gas nozzle inlet (3*H)
const double R = 3.;       // liquid inlet radius (3*H)
const double W = 15.;      // width of liquid nozzle (5*R)
const double thetai  = 45*pi/180.; // chamfer angle
const double thetai2 = 45*pi/180.; // outlet angle

int maxlevel = 8;
int minlevel = 5;
double tend = 2.5001;
double tout = 0.0125;	// interval for outputting images and dumps
double uemax = 1., semax = 1e-2, femax = 1e-1;

// dimensionless properties normalized by sigma, H, and rhol
double rhor = 80;    // density ratio
double rhol = 1;     // liquid density
double sigma = 1.;   // surface tension coefficient
double MUl = 5e-3;   // liquid viscosity
double MUg = 2.5e-2; // gas viscosity
double Ul = 0.5;     // liquid velocity
double Ug = 160;     // gas velocity
double theta0 = 60;  // liquid contact angle


// Boundary Conditions
u.n[left] = dirichlet (distance(y, z) <= 1.1*R? Ul: distance(y, z) <= 2.1*R+B? Ug: 0);
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
f[left]	  = dirichlet (distance(y, z) <= 1.1*R? 1: 0);

p[right]  = dirichlet (0);
pf[right] = dirichlet (0);
f[right]  = dirichlet (0);

u_x_ibm_dirichlet(0) // no slip  walls for solid boundary
u_y_ibm_dirichlet(0)
u_z_ibm_dirichlet(0)

	
// create the solid nozzle using C.S.G & the immersed boundary method
// inside solid   phi < 0
// outside solid  phi > 0
// on interface   phi = 0
void solid_domain (scalar s, face vector sf) {
  double chamfer = 2*H/tan(thetai);
  double offset  = 2*R/tan(thetai) - chamfer;
  double offset2 = 1.5*H - R/tan(thetai2);
  double tshell  = H;

  vertex scalar phi[];
  foreach_vertex() {
    if (x < -chamfer) // 1.
      phi[] = max(-min(max(cylinder(y, z, 2*R+B), -cylinder(y, z, 2*R)), cylinder(y, z, R)), cylinder(y, z, 2*R+B+tshell));
    else if (x >= -chamfer && x < 0) // 2.
      phi[] = max(-min(max(cylinder(y, z, 2*R+B), -(cone(y, z, x - offset, thetai))), cylinder(y, z, R)), cylinder(y, z, 2*R+B+tshell));
    else if (x >= 0 && x < H) // 3.
      phi[] = max(-cylinder(y, z, 2*R+B), cylinder(y, z, 2*R+B+tshell));
    else if (x >= H && x <= 2*H) // 4.
      phi[] = max(max(-cylinder(y, z, R), -(cone(y, z, x - offset2, thetai2))),cylinder(y, z, 2*R+B+tshell)) ;
    else
      phi[] = 1;
  }
  boundary({phi});
  fractions (phi, s, sf);
  fractions_cleanup (s, sf);
}


int main (int argc, char * argv[]) {
  double l0 = 50*H;
  
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    minlevel = atoi (argv[2]);
  if (argc > 3)
    tend = atof (argv[3]);
  if (argc > 4)
    rhor = atof (argv[4]);
  if (argc > 5)
    rhol = atof (argv[5]);
  if (argc > 6)
    sigma = atof (argv[6]);
  if (argc > 7)
    MUl = atof (argv[7]);
  if (argc > 8)
    MUg = atof (argv[8]);
  if (argc > 9)
    Ul = atof (argv[9]);
  if (argc > 10)
    Ug = atof (argv[10]);
  if (argc > 11)
    theta0 = atof (argv[11]);
  if (argc > 12)
    uemax = atof (argv[12]);
  if (argc > 13)
    tout  = atof (argv[13]);
  if (argc > 14)
    semax = atof (argv[14]);
  if (argc > 15)
    femax = atof (argv[15]);
  if (argc > 16)
    l0    = atof (argv[16]);

  size (l0);
  init_grid (1 << (minlevel));
  origin (-W, -(l0/2), -(l0/2));
  
  rho1 = rhol;
  rho2 = rho1 / rhor;
  f.sigma = sigma;

  mu1 = MUl;
  mu2 = MUg;

  const scalar c[] = theta0*pi/180;
  contact_angle = c;

  TOLERANCE = 1e-3;

  run();
}


event init (t = 0) {
  if (!restore("dump-1.275")) {
    scalar ibm1[], f0[];
    face vector ibmf1[]; 

    astats ss;
    int ic = 0;
    do {
      ic++;
      solid_domain(ibm1, ibmf1);
      fraction (f0, x <= -H && distance (y, z) < 1.1*R? 2*H - x: 0);
      ss = adapt_wavelet ({ibm1, f0}, (double[]) {1.e-30, 1.e-30}, // nrelax decreases with lower tolerances for some reason
      	        maxlevel = maxlevel , minlevel = minlevel); 
      fprintf(stderr, "%d refined = %d coarsened = %d\n", ic, ss.nf, ss.nc);
    } while ((ss.nf || ss.nc) && ic < 20);

    solid_domain(ibm, ibmf);
    boundary ({ibm, ibmf});

    fraction (f, x <= -H && distance (y, z) < 3.1*H? 2*H - x: 0);
    boundary ({f});
  }
  else {
    solid_domain (ibm, ibmf);
    boundary (all);
  }
}


event logfile (i++; t <= tend) {
  foreach_face()
    if(!fm.x[] && uf.x[])
      uf.x[] = 0;
 
  double maxu_x = 0, maxu_y = 0, maxu_z = 0;
  foreach_face(reduction(max:maxu_x) reduction(max:maxu_y) reduction(max:maxu_z)) {
    if (fm.x[] && fabs(uf.x[]) > fabs(maxu_x))
      maxu_x = uf.x[];
  }

  fprintf (stderr, "%d %g %d %d %d %d %d %d %g %g %g\n", 
		     i, t, mgpf.i, mgpf.nrelax, mgu.i, 
		     mgu.nrelax, mgp.i, mgp.nrelax, maxu_x, maxu_y, maxu_z);
  fflush(stderr);
}


// create the fake nozzle for clearer visualization
void fake_solid_domain (scalar s, face vector sf) {
  double chamfer = 2*H/tan(thetai);
  double offset  = 2*R/tan(thetai) - chamfer;
  double tshell  = H;

  vertex scalar phi[];
  foreach_vertex() {
    if (x < -chamfer) // 1.
      phi[] = max(-min(max(cylinder(y, z, 2*R+B), -cylinder(y, z, 2*R)), cylinder(y, z, R)), cylinder(y, z, 2*R+B+tshell));
    else if (x >= -chamfer && x < 0) // 2.
      phi[] = max(-min(max(cylinder(y, z, 2*R+B), -(cone(y, z, x - offset, thetai))), cylinder(y, z, R)), cylinder(y, z, 2*R+B+tshell));
    else
      phi[] = 1;
  }
  boundary({phi});
  fractions (phi, s, sf);
  fractions_cleanup (s, sf);
}


#if 1
event movie (t += tout; t <= tend) {
//event movie (i++; t <= tend) {

  scalar ibmFake[];
  face vector ibmfFake[]; 
  fake_solid_domain (ibmFake, ibmfFake);
  
  scalar f1[], f2[], f3[], umag[];

  double umin = HUGE, umax = -HUGE;
  foreach(reduction(max:umax) reduction(min:umin)) {
    f1[] = f[]*ibm[];
    f2[] = y > 0? 0: f[]*ibm[];
    f3[] = z > 0? 0: f[]*ibm[];
    umag[] = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
    if (umag[] > umax)
      umax = umag[];
    if (umag[] < umin)
      umin = umag[];
  }

  char name[80];
  sprintf (name, "img/fvof-%g.png", t);
  view (fov = 15, camera = "iso", tx = -0.15, bg = {0.2969,0.3984,0.5977}, width = 2000, height = 2000);
  clear();
  draw_vof ("f", fc = {1, 0.522, 0.439});
  draw_vof ("ibmFake", "ibmFakef");
  save (name);

  sprintf (name, "img/f1vof-%g.png", t);
  clear();
  draw_vof ("f1", fc = {1, 0.522, 0.439});
  save (name);

  sprintf (name, "img/f2vof-%g.png", t);
  clear();
  draw_vof ("f2", fc = {1, 0.522, 0.439});
  save (name);

  sprintf (name, "img/f3vof-%g.png", t);
  clear();
  draw_vof ("f3", fc = {1, 0.522, 0.439});
  save (name);

  sprintf (name, "img/velo-%g.png", t);
  view (fov = 20, camera = "front", tx=-0.1, bg = {1,1,1}, width = 2000, height = 2000);
  clear();
  squares ("umag", cbar = true, mid = true, min = umin, max = umax);
  save (name);

  sprintf (name, "img/f-%g.png", t);
  clear();
  squares ("f");
  save (name);

  sprintf (name, "img/cells-%g.png", t);
  clear();
  squares ("ibm");
  cells();
  save (name);

  sprintf (name, "img/f3front-%g.png", t);
  view (fov = 7, camera = "front", tx=0.1, bg = {1,1,1}, width = 2000, height = 2000);
  clear();
  draw_vof ("f3", fc = {1, 0.522, 0.439});
  save (name);

  sprintf (name, "img/f2top-%g.png", t);
  view (fov = 7, camera = "top", tx=0.1, bg = {1,1,1}, width = 2000, height = 2000);
  clear();
  draw_vof ("f2", fc = {1, 0.522, 0.439});
  save (name);
}
#endif

event dumps (t += tout; t <= tend) {
  scalar * dumplist = {u,p,f,ibm};
  char name[80];
  sprintf(name, "dumps/dump-%g", t);
  dump (file = name, list = dumplist);
}


// remove small droplets
event small_droplets (i += 50) {
  remove_droplets (f);
}


// count the number of droplets. taken from "atomisation.c"
event count_droplets (t += tout) {
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);

  double v[n];
  coord b[n];

  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;

  foreach (serial)
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
        b[j].x += dv()*f[]*p.x;
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  static FILE * fpd = fopen ("drops", "w");
  int ni = 0;
  for (int j = 0; j < n; j++) {
    if (v[j] > 0) {
        ni++;
        fprintf (fpd, "%d %g %d %g %g %g %d %d\n", i, t,
	                      j, v[j], b[j].x/v[j], b[j].y/v[j], n, ni);
    }
  }
  fflush (fpd);

  double maxv = -HUGE;
  int maxtag = 0, nreal = 0;
  for (int j = 0; j < n; j++) {
    if (v[j] > 0 && v[j] < 1e10)
      nreal++;
    if (v[j] > maxv && v[j] < 1e10) {
      maxv = v[j];
      maxtag = j;
    }
  }

  double area = 0;
  foreach(reduction(+:area))
    if (m[] == maxtag && on_interface(f)) {
      coord nf = interface_normal (point, f), pc;
      double alpha = plane_alpha (f[], nf);
      area += pow(Delta, dimension - 1)*plane_area_center (nf, alpha, &pc)*ibm[];
    }
  static FILE * fp = fopen ("surface_area", "w");
  if (i == 0) fprintf (fp, "i t area #drops\n");
  fprintf (fp, "%d %g %g %d\n", i, t, area, nreal); // was n instead of nreal by mistake
  fflush (fp);
}


event mixing_stats (i += 25)
{
    double totalVolume = pi*sq(R)*W;
    double liquidVolume = 0, area = 0;
    foreach(reduction(+:liquidVolume) reduction(+:area)) {
      if (x < 0 && fabs(y) < 1.1*R && fabs(z) < 1.1*R)
        liquidVolume += f[]*dv();
      if (on_interface(f)) {
        coord nf = interface_normal (point, f), pc;
        double alpha = plane_alpha (f[], nf);
        area += pow(Delta, dimension - 1)*plane_area_center (nf, alpha, &pc)*ibm[];
      }
    }

    norm unorm = normf (u.x);
    norm vnorm = normf (u.y);
    norm wnorm = normf (u.z);

    static FILE * fp = fopen ("mixing_stats", "w");
    if (i == 0) fprintf (fp, "i t lv lvFraction totalArea urms vrms wrms\n");
    fprintf (fp, "%d %g %g %g %g %g %g %g\n", i, t, liquidVolume, liquidVolume/totalVolume, 
                  area, unorm.rms, vnorm.rms, wnorm.rms);
    fflush (fp);
}

event adapt (i+=2) {
  scalar ibm1[], f1[];
  vector u1[];
  foreach(serial) {
    ibm1[] = ibm[];
    f1[]   = x > 28*H? 0: f[];
    u1.x[] = x > 28*H? u.x[-1]: u.x[];
    u1.y[] = x > 28*H? u.y[-1]: u.y[];
    u1.z[] = x > 28*H? u.z[-1]: u.z[];
  }

  fprintf(stderr, "starting adapt_wavelet\n");
  fflush(stderr);
  astats ad = adapt_wavelet ({ibm1,f1,u1}, (double[]){semax,femax,uemax,uemax,uemax}, 
 	  maxlevel = maxlevel, minlevel = minlevel);
  fprintf(stderr, "refined = %d coarsend = %d, starting unrefine\n", ad.nf, ad.nc);
  fflush(stderr);
#if 0
  foreach() {
    if (x > 28.1*H && fabs(y) > 1.5*R && fabs(z) > 1.5*R) {
      foreach_dimension()
        u.x[] = 0;
    }
    if (x > 28.1*H)
	f[] = 0;
  }
  boundary ({f, u});
#endif

  int nc = 0;
  unrefine2 (x > 28*H && level > minlevel+1, nc);
  fprintf (stderr, "unrefined %d cells. starting unrefine 3\n", nc);
#if 0
    fprintf (stderr, "starting unrefine 1\n");
    static const int too_fine = 1 << user;
    foreach_cell() {
      if (is_leaf(cell))
	continue;
      if (is_local(cell) && (x > 28*H && level > minlevel+3))
	cell.flags |= too_fine;
    }
    fprintf (stderr, "starting unrefine 2\n");
    int nc = 0;
    for (int _l = depth(); _l >= 0; _l--) {
      foreach_cell() {
	if (is_leaf(cell))
	  continue;
	if (level == _l) {
	  if (is_local(cell) && (cell.flags & too_fine)) {
	    coarsen_cell (point, all);
	    cell.flags &= ~too_fine;
	    nc++;
	  }
	  continue;
	}
      }
      mpi_boundary_coarsen (_l, too_fine);
    }

  mpi_all_reduce (nc, MPI_INT, MPI_SUM);
  fprintf (stderr, "unrefined %d cells. starting unrefine 3\n", nc);
  //mpi_boundary_update ({cm, fm, u, p});
  mpi_boundary_update(all);
  fprintf(stderr, "unrefined done\n");
  fflush(stderr);
#endif
}

