(* ::Package:: *)

TorqueAttitudeTrackingControllerTime[
{n1_,n2_,n3_,\[Omega]1_,\[Omega]2_,\[Omega]3_},
{{ndd0x_,ndd0y_,ndd0z_},{ndd1x_,ndd1y_,ndd1z_},{ndd2x_,ndd2y_,ndd2z_}},
ControllerGains_
]:=Module[{

	(*state*)
	x = {n1,n2,n3,\[Omega]1,\[Omega]2,\[Omega]3},
	n = {n1,n2,n3},
	\[Omega] = {\[Omega]1,\[Omega]2,\[Omega]3},

	ndd0 = {ndd0x,ndd0y,ndd0z},
	ndd1 = {ndd1x,ndd1y,ndd1z},
	ndd2 = {ndd2x,ndd2y,ndd2z},

	(*part 1*)
	\[Omega]nd,\[Omega]cl,
	(*part 3*)
	\[Tau]\[Omega]nd,\[Tau]\[Omega]cl,\[Tau]cl,

(*gains*)
kp =kp/.ControllerGains,
kd =kd/.ControllerGains,
kVp=kVp/.ControllerGains,
kVd=kVd/.ControllerGains,

X0,V0,W0,
X,V,W,

(*defining some useful functions*)
(*Orthogonal Projection Operator, in R3*)
OP = IdentityMatrix[3] - KroneckerProduct[#, #]&,
(*Skew symmetric matrix in R3*)
Skew = {{0,-#[[3]],#[[2]]},{#[[3]],0,-#[[1]]},{-#[[2]],#[[1]],0}}&

},

(*part 1*)
\[Omega]nd = Skew[ndd0].ndd1;
(*backstepping term*)
\[Omega]cl = \[Omega]nd - kp Skew[ndd0].n;

X0 = Skew[\[Omega]].n;
V0 = kVp(1-n.ndd0);
W0 = -kVp kp (Skew[n].ndd0).(Skew[n].ndd0);

(*part 2*)
\[Tau]\[Omega]nd = Skew[ndd0].ndd2;
\[Tau]\[Omega]cl = \[Tau]\[Omega]nd - kp Skew[ndd1].n - kp Skew[ndd0].Skew[\[Omega]].n;
\[Tau]cl = OP[n].(\[Tau]\[Omega]cl + (\[Omega]cl.n)Skew[n].\[Omega]cl - kd(\[Omega] - \[Omega]cl) + kVp/kVd Skew[n].ndd0);

X = Join[X0, OP[n].\[Tau]cl];
V = V0 + kVd/2 (Skew[n].(\[Omega]-\[Omega]cl)).(Skew[n].(\[Omega]-\[Omega]cl));
W = W0 - kVd kd (Skew[n].(\[Omega]-\[Omega]cl)).(Skew[n].(\[Omega]-\[Omega]cl));

{\[Tau]cl,V,W,X}
]

