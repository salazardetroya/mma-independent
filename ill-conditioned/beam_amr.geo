// This code was created by pygmsh vunknown.
p0 = newp;
Point(p0) = {0.0, 0.0, 0.0, 1.0};
p1 = newp;
Point(p1) = {33.333333333333336, 0.0, 0.0, 1.0};
p2 = newp;
Point(p2) = {66.66666666666667, 0.0, 0.0, 1.0};
p3 = newp;
Point(p3) = {100.0, 0.0, 0.0, 1.0};
p4 = newp;
Point(p4) = {100.0, 16.0, 0.0, 0.08};
p5 = newp;
Point(p5) = {100.0, 24.0, 0.0, 0.08};
p6 = newp;
Point(p6) = {100.0, 40.0, 0.0, 1.0};
p7 = newp;
Point(p7) = {66.66666666666667, 40.0, 0.0, 1.0};
p8 = newp;
Point(p8) = {33.333333333333336, 40.0, 0.0, 1.0};
p9 = newp;
Point(p9) = {0.0, 40.0, 0.0, 1.0};
l0 = newl;
Line(l0) = {p0, p1};
l1 = newl;
Line(l1) = {p1, p2};
l2 = newl;
Line(l2) = {p2, p3};
l3 = newl;
Line(l3) = {p3, p4};
l4 = newl;
Line(l4) = {p4, p5};
l5 = newl;
Line(l5) = {p5, p6};
l6 = newl;
Line(l6) = {p6, p7};
l7 = newl;
Line(l7) = {p7, p8};
l8 = newl;
Line(l8) = {p8, p9};
l9 = newl;
Line(l9) = {p9, p0};
ll0 = newll;
Line Loop(ll0) = {l0, l1, l2, l3, l4, l5, l6, l7, l8, l9};
s0 = news;
Plane Surface(s0) = {ll0};
Physical Surface(1) = {s0};
Physical Line(4) = {l4};
Physical Line(3) = {l9};

Field[1] = Distance;
Field[1].PointsList = {p7};
Field[1].CurvesList = {l3, l4, l5, l9};
Field[1].NumPointsPerCurve = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.2;
Field[2].SizeMax = 5.0;
Field[2].DistMin = 10.0;
Field[2].DistMax = 20.0;
Background Field = 2;
