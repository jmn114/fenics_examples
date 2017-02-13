Point (1) = {-0.5, 0, 0, 0.005};
Point (2) = {-0.5, 0.1, 0, 0.005};
Line (1) = {1, 2};

//Extrude to make box
Extrude {1.0,0,0} {
  Line{1}; Layers{201};
}
// top
Physical Line(1) = {4};
// bottom
Physical Line(2) = {3};
// ends
Physical Line(3) = {2,1};
Physical Surface(1) = {5};
