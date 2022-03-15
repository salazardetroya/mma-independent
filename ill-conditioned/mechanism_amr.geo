// This code was created by pygmsh vunknown.
p0 = newp;
size = 0.5;
Point(p0) = {0.0, 0.0, 0.0, size};
p1 = newp;
Point(p1) = {120, 0.0, 0.0, size};
p2 = newp;
Point(p2) = {120, 50, 0.0, size};
p3 = newp;
Point(p3) = {120.0, 60.0, 0.0, size};
p4 = newp;
Point(p4) = {0.0, 60.0, 0.0, size};
p5 = newp;
Point(p5) = {0.0, 50.0, 0.0, size};
p6 = newp;
Point(p6) = {0.0, 2.0, 0.0, size};
p7 = newp;
Point(p7) = {22, 45.0, 0.0, size};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 1};
//+
Curve Loop(1) = {4, 5, 6, 7, 1, 2, 3};
Plane Surface(1) = {1};

Field[1] = Distance;
Field[1].PointsList = {p7};
Field[1].CurvesList = {1};
Field[1].NumPointsPerCurve = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.2;
Field[2].SizeMax = 5.0;
Field[2].DistMin = 10.0;
Field[2].DistMax = 20.0;
Background Field = 2;

//+
Physical Curve(1) = {7};
//+
Physical Curve(2) = {5};
//+
Physical Curve(3) = {3};
//+
Physical Curve(4) = {4};
//+
Physical Surface(1) = {1};
