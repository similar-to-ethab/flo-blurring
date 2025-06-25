#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"
#include "view.h"
#include "output_stl.h"

double cubeF (double x, double y, double z, coord center, double size)
{
		double P1_Plus = x - size/2. + center.x;
  double P1_Minus = x + size/2. + center.x;
  double P1 = max (P1_Plus, -P1_Minus);

    double P2_Plus = y - size/2. + center.y;
  double P2_Minus = y + size/2. + center.y;
  double P2 = max (P2_Plus, -P2_Minus);

  double P3_Plus = z - size/2. + center.z;
  double P3_Minus = z + size/2. + center.z;
  double P3 = max (P3_Plus, -P3_Minus);

    double c = max (P1, max (P2, P3));

  return c;
}

double sphere (double x, double y, double z, coord center, double radius)
{
  return (sq(x - center.x) + sq (y - center.y) + sq (z - center.z)
          - sq (radius));
}

double geometry (double x, double y, double z)
{
  coord center = {0,0,0};
    double s = sphere (x, y, z, center, 0.25);
  double c = cubeF (x, y, z, center, 0.38);
    double sIc = max (s, c);

	  double cylinder1 = sq(x) + sq(y) - sq(0.12);
  double cylinder2 = sq(z) + sq(y) - sq(0.12);
  double cylinder3 = sq(x) + sq(z) - sq(0.12);

    double cylinderInter = min (cylinder1, cylinder2);
  double cylinderUnion = min (cylinderInter, cylinder3);

    return max(sIc, -cylinderUnion);
}

int main() {
		  size (1 [0]);
  origin (-0.5, -0.5, -0.5);
    init_grid (1<<7);
	  scalar f[];
	    int iteration = 0;
  do
    fraction (f, geometry (x, y, z));
  while (adapt_wavelet({f}, (double []){0.2}, maxlevel = 9, 2).nf &&
         iteration++ <= 10);
    view (fov = 13.0359, quat = {-0.353553,0.353553,0.146447,0.853553},
        bg = {1,1,1}, width = 600, height = 600);
  draw_vof ("f");
  draw_vof ("f", edges = true);
  save ("vof.png");
  stl_output_binary (f, "csg.stl");
}
