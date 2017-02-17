D = 0.1;
size = 3e-2;

Point (1) = {10*D, 0.05*D + D/2, 0, size/10};
Point (2) = {10*D + D/2, 0.05*D + D/2, 0, size/10};
Point (3) = {10*D, 0.05*D + D, 0, size/10};
Point (4) = {10*D - D/2, 0.05*D + D/2, 0, size/10};
Point (5) = {10*D, 0.05*D, 0, size/10};
Point (6) = {0, 0, 0, size};
Point (7) = {10*D - D/2, 0, 0, size/10};
Point (8) = {10*D + D/2, 0, 0, size/10};
Point (9) = {30*D, 0, 0, size};
Point (10) = {30*D, 4*D, 0, size};
Point (11) = {10*D + D/2, 4*D, 0, size/10}; 
Point (12) = {10*D - D/2, 4*D, 0, size/10};
Point (13) = {0, 4*D, 0, size};

Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};
Line (5) = {6, 13};
Line (6) = {13, 12};
Line (7) = {12, 11};
Line (8) = {11, 10};
Line (9) = {6, 7};
Line (10) = {7,8};
Line (11) = {8, 9};
Line (12) = {9, 10};

Physical Line (1) = {5}; //left
Physical Line (2) = {6,7,8}; //top
Physical Line (3) = {9,10,11}; //bottom
Physical Line (4) = {12}; //right
Physical Line (5) = {1, 2, 3, 4}; //circle

Line Loop (13) = {1, 2, 3, 4};
Line Loop (14) = {-5, -6, -7, -8, 9, 10, 11, 12};
Plane Surface (15) = {14, 13};
Physical Surface (16) = {15};
/*
Line Loop (13) = {1, 2, 3, 4};
Line Loop (14) = {-5, -6, -7, -8, 9, 10, 11, 12};
Plane Surface (15) = {14, 13};

Physical Volume("internal") = {1};
  Extrude {0, 0, 0.01} {
   Surface{15};
   Layers{1};
   Recombine;
  }
Physical Surface("inlet") = {32};
Physical Surface("outlet") = {48};
Physical Surface("top") = {52,56,60};
Physical Surface("bottom") = {36,40,44};
Physical Surface("cylinder") = {64,68,72,76};
Physical Surface("sides") = {15,77};
*/
Field[1]=Attractor;
Field[1].EdgesList={1,2,3,4,7,10};

Field[2]=MathEval;
Field[2].F="0.005+0.15*F1";
Background Field=2;

View "comments" {T2(10, 15, 0){ StrCat("File created on ", Today) };};
