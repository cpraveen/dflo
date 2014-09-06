Mesh.RecombineAll=1;           //recombine all defined surfaces
Mesh.Algorithm=8;              //delquad mesher
Mesh.RecombinationAlgorithm=1; //blossom

xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = 1.0;

// y length is 3 times x length
nx = 11;
ny = nx;
cl = 1.0/(nx-1);

Point(1) = {xmin, ymin, 0, 0.5*cl};
Point(2) = {xmax, ymin, 0, 0.6*cl};
Point(3) = {xmax, ymax, 0, 0.7*cl};
Point(4) = {xmin, ymax, 0, 0.8*cl};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Line(1) = {1}; // bottom
Physical Line(2) = {2}; // right
Physical Line(3) = {3}; // top
Physical Line(4) = {4}; // left
Physical Surface(10) = {1};
