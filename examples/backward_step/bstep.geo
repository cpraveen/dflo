// Mesh for backward facing step

Mesh.RecombineAll=1;
Mesh.RecombinationAlgorithm=1; //blossom

L1 = 1;    // Length of the horizontal part of the step
L2 = 12;   // length of domain after the step
h = 6;     // Height of vertical part of the step
H = 11;    // Total domain height 

n  = 16;
cl = 1.0/n;  // Mesh size

n1 = L1*n + 1;
n2 = L2*n + 1;
n3 = (H-h)*n + 1;
n4 = h*n + 1;

Point(1) = {0, h, 0};
Point(2) = {L1, h, 0};
Point(3) = {L1, 0, 0};
Point(4) = {L1+L2, 0, 0};
Point(5) = {L1+L2, h, 0};
Point(6) = {L1+L2, H, 0};
Point(7) = {L1, H, 0};
Point(8) = {0, H, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 5};
Line(10) = {2, 7};

Line Loop(1) = {1, 10, 7, 8};
Plane Surface(1) = {1};
Transfinite Surface(1) = {1,2,7,8};

Line Loop(2) = {2, 3, 4, -9};
Plane Surface(2) = {2};
Transfinite Surface(2) = {2,3,4,5};

Line Loop(3) = {5, 6, -10, 9};
Plane Surface(3) = {3};
Transfinite Surface(3) = {5,6,7,2};

Transfinite Line{1, -7} = n1;
Transfinite Line{3, 9, -6} = n2;
Transfinite Line{-8, 10, 5} = n3;
Transfinite Line{-2,4} = n4;

Physical Line(1) = {8};        		// inflow
Physical Line(2) = {1,2};  		// wall
Physical Line(3) = {3,4,5,6,7};         // outlet

Physical Surface(100) = {1,2,3};
