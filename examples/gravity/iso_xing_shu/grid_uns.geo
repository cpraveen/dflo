Mesh.RecombineAll=1;           //recombine all defined surfaces
Mesh.Algorithm=8;              //delquad mesher
Mesh.RecombinationAlgorithm=1; //blossom

r = 1.0; // Radius of the circle

n  = 700;
lc = 2*Pi*r/n; // Characteristic length

xc = 0.0;  // x coordinate of circle centre
yc = 0.0;  // y coordinate of circle centre

Point(1)  = {xc, yc, 0, lc};
 
Point(2)  = {xc + r*Cos(Pi/4), yc + r*Sin(Pi/4), 0, lc};
Point(3)  = {xc + r*Cos(3*Pi/4), yc + r*Sin(3*Pi/4), 0, lc};
Point(4)  = {xc + r*Cos(5*Pi/4), yc + r*Sin(5*Pi/4), 0, lc};
Point(5)  = {xc + r*Cos(7*Pi/4), yc + r*Sin(7*Pi/4), 0, lc};

Circle(1) = {5, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {3, 1, 4};
Circle(4) = {4, 1, 5};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Line(1) = {1,2,3,4};    // boundary

Physical Surface(100) = {1};
