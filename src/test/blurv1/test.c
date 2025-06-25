#include <math.h>
#include <stdio.h>

// Define min macro if not already provided.
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

// Define a coordinate structure.
typedef struct {
    double x;
    double y;
} coord;

// Returns the distance from point (px,py) to the segment from point a to b.
double sdSegment(double x, double y, coord a, coord b) {
    double vx = b.x - a.x;
    double vy = b.y - a.y;
    double wx = x - a.x;
    double wy = y - a.y;
    
    double t = (vx * wx + vy * wy) / (vx * vx + vy * vy);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    double projx = a.x + t * vx;
    double projy = a.y + t * vy;
    return hypot(x - projx, y - projy);
}

// Computes a 2D cross product (edge function) for the segment from a to b relative to point (px,py).
double edgeFunction(coord a, coord b, double x, double y) {
    return (x - a.x) * (b.y - a.y) - (y - a.y) * (b.x - a.x);
}

// Signed Distance Function for a triangle.
// Parameters:
//   x, y     : Evaluation point.
//   corner   : The first vertex (A) of the triangle.
//   w        : Horizontal width of the triangle.
//   m        : Slope factor which, multiplied by w, gives the vertical offset.
//   val_sign : Orientation flag (+1 or -1) to flip the vertical direction.
// The triangle vertices are defined as:
//   A = corner
//   B = (corner.x + w, corner.y)
//   C = (corner.x + w, corner.y + m * w * val_sign)
double sdfTriangle(double x, double y, coord corner, double w, double m, int val_sign) {
    // Define triangle vertices.
    coord A = corner;
    coord B = { corner.x + w, corner.y };
    coord C = { corner.x + w, corner.y + m * w * val_sign };

    // Compute distances to each of the triangle's edges.
    double dAB = sdSegment(x, y, A, B);
    double dBC = sdSegment(x, y, B, C);
    double dCA = sdSegment(x, y, C, A);
    
    // Determine the minimum distance.
    double d = min(dAB, min(dBC, dCA));

    // Compute the area using the edge function to determine orientation.
    double area = edgeFunction(A, B, C.x, C.y);
    int inside = 0;
    if (area > 0.0) {
        if (edgeFunction(A, B, x, y) >= 0.0 &&
            edgeFunction(B, C, x, y) >= 0.0 &&
            edgeFunction(C, A, x, y) >= 0.0)
            inside = 1;
    } else {
        if (edgeFunction(A, B, x, y) <= 0.0 &&
            edgeFunction(B, C, x, y) <= 0.0 &&
            edgeFunction(C, A, x, y) <= 0.0)
            inside = 1;
    }
    
    // Return negative distance if inside the triangle.
    return inside ? -d : d;
}

// Main function to test the sdfTriangle.
int main() {
    // Define the triangle.
    coord corner = { 0.0, 0.0 };
    double w = 1.0;
    double m = 0.5;
    int val_sign = 1;
    
    // Test point.
    double x = 0.7, y = -0.3;
    
    // Compute the SDF.
    double distance = sdfTriangle(x, y, corner, w, m, val_sign);
    
    printf("The signed distance to the triangle is: %f\n", distance);
    return 0;
}

