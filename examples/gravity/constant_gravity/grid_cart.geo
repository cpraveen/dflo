xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = 1.0;

// y length is 3 times x length
nx = 101;
ny = nx;

Point(1) = {xmin, ymin, 0};
Point(2) = {xmax, ymin, 0};
Point(3) = {xmax, ymax, 0};
Point(4) = {xmin, ymax, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};

Recombine Surface(1);
Transfinite Line{1,3} = nx;
Transfinite Line{2,4} = ny;

Physical Line(1) = {1}; // bottom
Physical Line(2) = {2}; // right
Physical Line(3) = {3}; // top
Physical Line(4) = {4}; // left
Physical Surface(10) = {1};
