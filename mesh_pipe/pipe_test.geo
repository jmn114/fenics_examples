D = 0.1;
size = 1e-1;

Point (1) = {2*D, 0.5*D + D/2, 0, size};
Point (2) = {2*D + D/2, 0.5*D + D/2, 0, size};
Point (3) = {2*D, 0.5*D + D, 0, size};
Point (4) = {2*D - D/2, 0.5*D + D/2, 0, size};
Point (5) = {2*D, 0.5*D, 0, size};
Point (6) = {0, 0, 0, size};
Point (7) = {5*D, 0, 0, size};
Point (8) = {5*D, 2*D, 0, size};
Point (9) = {0, 2*D, 0, size};

Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};
Line (5) = {6, 7};
Line (6) = {7, 8};
Line (7) = {8, 9};
Line (8) = {9, 6};

Physical Line (1) = {8}; //left
Physical Line (2) = {7}; //top
Physical Line (3) = {5}; //bottom
Physical Line (4) = {6}; //right
Physical Line (5) = {1, 2, 3, 4}; //circle

Line Loop (11) = {1, 2, 3, 4};
Line Loop (12) = {-5, -6, -7, -8};
Plane Surface (13) = {12, 11};
Physical Surface (14) = {13};

//Field[1]=Attractor;
//Field[1].EdgesList={3,4,5};

//Field[2]=MathEval;
//Field[2].F="0.005+0.5*F1";
//Background Field=2;

View "comments" {T2(10, 15, 0){ StrCat("File created on ", Today) };};
Coherence;
