H  = 3.0;
L1 = 1.0;
L2 = 4.0;
theta = 9.5; // in degrees

n1 = 10;
n2 = 30;
n3 = 20;

// convert to radians
theta = Pi * theta / 180;

Point(1) = {0, 0, 0};
Point(2) = {L1, 0, 0};
Point(3) = {L1 + L2, Tan(theta)*L2, 0};
Point(4) = {L1+L2, H, 0};
Point(5) = {L1,    H, 0};
Point(6) = {0,     H, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,5};

Line Loop(1) = {1,7,5,6};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,5,6};

Line Loop(2) = {2,3,4,-7};
Ruled Surface(2) = {2};
Transfinite Surface(2) = {2,3,4,5};

Recombine Surface(1);
Recombine Surface(2);

Transfinite Line {1,5} = n1;
Transfinite Line {2,4} = n2;
Transfinite Line {6,7,3} = n3;

Physical Line(1) = {1,2,4,5}; // solid walls, top and bottom
Physical Line(2) = {6};       // inflow
Physical Line(3) = {3};       // outflow

Physical Surface(4) = {1,2};
