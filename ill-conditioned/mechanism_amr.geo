// This code was created by pygmsh vunknown.
p0 = newp;
size = 1.0;
size2 = 0.1;
Point(p0) = {0.0, 0.0, 0.0, size2};
p1 = newp;
Point(p1) = {120, 0.0, 0.0, size};
p2 = newp;
Point(p2) = {120, 50, 0.0, size2};
p3 = newp;
Point(p3) = {120.0, 60.0, 0.0, size2};
p8 = newp;
Point(p8) = {60.0, 60.0, 0.0, size};
p4 = newp;
Point(p4) = {0.0, 60.0, 0.0, size2};
p5 = newp;
Point(p5) = {0.0, 50.0, 0.0, size2};
p6 = newp;
Point(p6) = {0.0, 2.0, 0.0, size2};
p9 = newp;
Point(p9) = {0.0, 25.0, 0.0, size};
p10 = newp;
Point(p10) = {60.0, 0.0, 0.0, size2};
//+
Line(1) = {1, 10};
//+
Line(2) = {10, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 9};
//+
Line(9) = {9, 8};
//+
Line(10) = {8, 1};
//+
Curve Loop(1) = {6, 7, 8, 9, 10, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {10};
//+
Physical Curve(2) = {7};
//+
Physical Curve(3) = {4};
//+
Physical Curve(4) = {5, 6};
//+
Physical Surface(1) = {1};
