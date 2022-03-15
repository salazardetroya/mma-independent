// This code was created by pygmsh vunknown.
size = 1.0;
p0 = newp;
Point(p0) = {0.0, 0.0, 0.0, size};
p1 = newp;
Point(p1) = {100, 0.0, 0.0, size};
p2 = newp;
Point(p2) = {100, 40, 0.0, size};
p3 = newp;
Point(p3) = {60.0, 40.0, 0.0, size};
p4 = newp;
Point(p4) = {60.0, 60.0, 0.0, size};
p5 = newp;
Point(p5) = {40.0, 60.0, 0.0, size};//+
p6 = newp;
Point(p6) = {40.0, 100.0, 0.0, size};//+//+
p7 = newp;
Point(p7) = {0.0, 100.0, 0.0, size};//+//+
p8 = newp;
Point(p8) = {0.0, 60.0, 0.0, size};//+//+
p9 = newp;
Point(p9) = {0.0, 40.0, 0.0, size};//+//+//+//+
p10 = newp;
Point(p10) = {40.0, 40.0, 0.0, size * 0.01};
p11 = newp;
Point(p11) = {95.0, 40.0, 0.0, size};//+
Line(1) = {8, 9};
//+
Line(2) = {9, 10};
//+
Line(3) = {10, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {2, 3};
//+
Line(6) = {3, 12};
//+
Line(7) = {12, 4};
//+
Line(8) = {4, 5};
//+
Line(9) = {5, 6};
//+
Line(10) = {6, 7};
//+
Line(11) = {7, 8};
//+
Line(12) = {9, 6};
//+
Line(13) = {6, 11};
//+
Line(14) = {11, 4};
//+
Line(15) = {11, 11};
//+
Line(16) = {10, 11};
//+
Curve Loop(1) = {11, 1, 12, 10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, 13, -16, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {13, 14, 8, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, -14, -16, 3, 4, 5, 6};
//+
Plane Surface(4) = {4};

Physical Curve(1) = {11};
//+
Physical Curve(2) = {6};
//+
Physical Surface(1) = {3};
//+
Physical Surface(0) = {1, 2, 4};
