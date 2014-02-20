Mesh.RecombineAll=1;
Mesh.RecombinationAlgorithm=1; //blossom

L1=0.6;
L2=2.4;
h1=0.2;
h2=0.8;

cl=0.01;

n1=40;
n2=20;
n3=40;
n4=40;

n1 = L1/cl;
n2 = h1/cl;
n3 = h2/cl;
n4 = L2/cl;

Point(1) = {0, 0, 0};
Point(2) = {L1, 0, 0};
Point(3) = {L1, h1, 0};
Point(4) = {L1+L2, h1, 0};
Point(5) = {L1+L2, h1+h2, 0};
Point(6) = {L1, h1+h2, 0};
Point(7) = {0, h1+h2, 0};
Point(8) = {0, h1, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {3,8};
Line(10) = {6,3};

Line Loop(1) = {1,2,9,8};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,8};

Line Loop(2) = {-9,-10,6,7};
Ruled Surface(2) = {2};
Transfinite Surface(2) = {8,3,6,7};

Line Loop(3) = {3,4,5,10};
Ruled Surface(3) = {3};
Transfinite Surface(3) = {3,4,5,6};

Transfinite Line{1,9,6} = n1;
Transfinite Line{2,8} = n2;
Transfinite Line{4,7,10} = n3;
Transfinite Line{3,5} = n4;

Physical Surface(100000) = {1,2,3};

Physical Line(1) = {7,8};        // in flow
Physical Line(2) = {1,2,3,5,6};  // bottom and top
Physical Line(3) = {4};          // outlet
