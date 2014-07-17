Mesh.RecombineAll=1;
Mesh.RecombinationAlgorithm=1; //blossom

Lx = 4.0;
Ly = 1.0;
x0 = 1.0/6.0;

ny = 101;
dy = Ly/(ny-1);

Printf("h = %e", dy);

n1 = Ceil(x0/dy);
n2 = Ceil((Lx-x0)/dy);

xmin = x0 - n1*dy;
xmax = x0 + n2*dy;

Printf("xmin = %e", xmin);
Printf("xmax = %e", xmax);

Point(1) = {xmin,  0,  0};
Point(2) = {x0  ,  0,  0};
Point(3) = {xmax, 0,  0};
Point(4) = {xmax, Ly, 0};
Point(5) = {x0  ,  Ly, 0};
Point(6) = {xmin,  Ly, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,5};

Line Loop(1) = {1,7,5,6};
Plane Surface(1) = {1};
Transfinite Surface(1) = {1,2,5,6};

Line Loop(2) = {2,3,4,-7};
Plane Surface(2) = {2};
Transfinite Surface(2) = {2,3,4,5};

Transfinite Line{1,5} = n1 + 1;
Transfinite Line{2,4} = n2 + 1;
Transfinite Line{3,6,7} = ny;

Physical Surface(100) = {1,2};

Physical Line(0) = {1};   // bottom  (x < x0)
Physical Line(1) = {2};   // bottom  (x > x0)
Physical Line(2) = {3};   // right
Physical Line(3) = {4,5}; // top
Physical Line(4) = {6};   // left
