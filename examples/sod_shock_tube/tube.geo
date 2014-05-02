Mesh.RecombineAll=1;
Mesh.RecombinationAlgorithm=1; //blossom

Lx = 1.0;

nx = 101; // number of points along x
ny = 11;  // number of points along y

dx = Lx/(nx-1);
Ly = dx*(ny-1);

Point(1) = {0,  0,  0};
Point(2) = {Lx, 0,  0};
Point(3) = {Lx, Ly, 0};
Point(4) = {0,  Ly, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};

Transfinite Line{1,3} = nx;
Transfinite Line{2,4} = ny;

Physical Surface(100) = {1};

Physical Line(0) = {1,3}; // side walls
Physical Line(1) = {2};   // outlet
Physical Line(2) = {4};   // inlet
