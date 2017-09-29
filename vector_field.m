(* ::Package:: *)

(* ::Chapter:: *)
(*Vector field*)


UAVBarVectorField[
(*1st input: position of system*)
{p1_,p2_,p3_,(*manipulator's linear position: in inertial frame*)
r11_,r12_,r13_,r21_,r22_,r23_,r31_,r32_,r33_,(*manipulator's attitude*)
P1_,P2_,P3_,(*uav's linear position: in inertial frame*)
R11_,R12_,R13_,R21_,R22_,R23_,R31_,R32_,R33_(*uav's attitude*)
},
(*2nd input: velocity of system*)
{
v1_,v2_,v3_,(*manipulator's linear velocity: in inertial frame*)
\[Omega]1_,\[Omega]2_,\[Omega]3_,(*manipulator's body-framed angular velocity*)
V1_,V2_,V3_,(*uav's linear velocity: in inertial frame*)
\[CapitalOmega]1_,\[CapitalOmega]2_,\[CapitalOmega]3_(*uav's body-framed angular velocity*)
},
(*3rd input: input to system*)
{
U_(*uav's thrust input*),
\[Tau]1_,\[Tau]2_,\[Tau]3_(*uav's body-frame torque input*),
\[Tau]m1_,\[Tau]m2_(*torque inputs on manipulator*)
},
(*4th input: physical constants*)
{
g_(*gravity*),
(**)
(*propeties of rigid body 1: manipulator*)
m_(*mass*),
(*moment of inertia: must be symmetric*)
jxy_,jzz_,
(**)
(*propeties of rigid body 2: UAV*)
M_(*mass*),
(*moment of inertia: must be symmetric*)
Jxx_,Jyy_,Jzz_,
(**)
(*propeties related to restriction between rigid bodies*)
l_(*distance between center's of mass*)
}
]:=Module[{

(*state components' names: these are dummy names: are other names would work*)

(*RIDIG BODY 1*)
(*manipulator's pose*)
p={p1,p2,p3},(*manipulator's linear position: in inertial frame*)
r={{r11,r12,r13},{r21,r22,r23},{r31,r32,r33}},(*manipulator's attitude*)
(*manipulator's twist*)
v={v1,v2,v3},(*manipulator's linear velocity: in inertial frame*)
\[Omega]={\[Omega]1,\[Omega]2,\[Omega]3}, (*manipulator's body-framed angular velocity*)

(*moment of inertia*)(*in bar jxx=jyy*)
j=DiagonalMatrix[{jxy,jxy,jzz}],

(*RIDIG BODY 2*)
(*uav's pose*)
P={P1,P2,P3},(*uav's linear position: in inertial frame*)
R={{R11,R12,R13},{R21,R22,R23},{R31,R32,R33}},(*uav's attitude*)
(*uav's twist*)
V={V1,V2,V3},(*uav's linear velocity: in inertial frame*)
\[CapitalOmega]={\[CapitalOmega]1,\[CapitalOmega]2,\[CapitalOmega]3}, (*uav's body-framed angular velocity*)

(*moment of inertia*)
J=DiagonalMatrix[{Jxx,Jyy,Jzz}],

(*input components: these are dummy names*)
u={
U(*uav's thrust input*),
\[Tau]1,\[Tau]2,\[Tau]3(*uav's body-frame torque input*),
\[Tau]m1,\[Tau]m2(*torque inputs on manipulator*)
},
(*canonical basis vectors*)
e1={1,0,0},
e2={0,1,0},
e3={0,0,1},

(*skew matrix*)
Skew={{0,-#3,#2},{#3,0,-#1},{-#2,#1,0}}&,
(*orthogonal projection*)
OP=IdentityMatrix[3]-KroneckerProduct[{#1,#2,#3},{#1,#2,#3}]&,

xp,xv,
n,T,X

},

(*position of system*)
xp=Join[p,Flatten[r],P,Flatten[R]];
(*velocity of system*)
xv=Join[v,\[Omega],V,\[CapitalOmega]];

(*cable unit vector*)
n=(P-p)/l;
(*tension on cable*)
T=1/(1+jxy/l^2 (m+M)/(m M)) (jxy/l^2 IdentityMatrix[3]+(m M)/(m+M) KroneckerProduct[e3,e3]).(U/M Transpose[r].R.e3 +l/jxy (Skew@@e3).{\[Tau]m1,\[Tau]m2,0}+l (e3 \[Omega].(OP@@e3).\[Omega]-jzz/jxy \[Omega].e3 (OP@@e3).\[Omega]));

(*vector field*)
X=
Join[(*KINEMATIC EQUATIONS*)
v,(*linear velocity of manipulator*)
Flatten[r.(Skew@@\[Omega])],
V,(*linear velocity of manipulator*)
Flatten[R.(Skew@@\[CapitalOmega])],
(*DYNAMIC EQUATIONS*)
1/m (r.T-m g e3),(*linear acceleration of manipulator*)
Inverse[j].(-(Skew@@\[Omega]).j.\[Omega] +l (Skew@@e3).T+{\[Tau]m1,\[Tau]m2,0}),(*angular acceleration of manipulator*)
1/M (U R.e3-r.T-M g e3),(*linear acceleration of uav*)
Inverse[J].(-(Skew@@\[CapitalOmega]).J.\[CapitalOmega] +{\[Tau]1,\[Tau]2,\[Tau]3}-Transpose[R].r.{\[Tau]m1,\[Tau]m2,0})(*angular acceleration of uav*)
];
{X,T}
]


UAVBarState[l_(*cable length*)]:=Module[{
Skew={{0,-#3,#2},{#3,0,-#1},{-#2,#1,0}}&,
SkewInv={-#6,#3,-#2}&,
e3={0,0,1},

(*any x in State set can be generated from the mapping below*)
Rx={{1,0,0},{0,Cos[#1],-Sin[#1]},{0,Sin[#1],Cos[#1]}}&,
Ry={{Cos[#1],0,Sin[#1]},{0,1,0},{-Sin[#1],0,Cos[#1]}}&,
Rz={{Cos[#1],-Sin[#1],0},{Sin[#1],Cos[#1],0},{0,0,1}}&,

(*RIDIG BODY 1*)
(*manipulator's pose*)
p={p1,p2,p3},(*manipulator's linear position: in inertial frame*)
r={{r11,r12,r13},{r21,r22,r23},{r31,r32,r33}},(*manipulator's attitude*)
(*manipulator's twist*)
v={v1,v2,v3},(*manipulator's linear velocity: in inertial frame*)
\[Omega]={\[Omega]1,\[Omega]2,\[Omega]3}, (*manipulator's body-framed angular velocity*)

(*RIDIG BODY 2*)
(*uav's pose*)
P={P1,P2,P3},(*uav's linear position: in inertial frame*)
R={{R11,R12,R13},{R21,R22,R23},{R31,R32,R33}},(*uav's attitude*)
(*uav's twist*)
V={V1,V2,V3},(*uav's linear velocity: in inertial frame*)
\[CapitalOmega]={\[CapitalOmega]1,\[CapitalOmega]2,\[CapitalOmega]3}, (*uav's body-framed angular velocity*)

x,xGeneric,u,xRandom
},
Module[{
Rot=Rz[#1].Ry[#2].Rx[#3]&
},
Module[{
DRot=Simplify[Evaluate[D[Rot[#1,#2,#3],{{#1,#2,#3}}]].{#4,#5,#6}]&,
\[Omega]Rot=Simplify[SkewInv@@Flatten[Transpose[Rot[#1,#2,#3]].Evaluate[D[Rot[#1,#2,#3],{{#1,#2,#3}}]].{#4,#5,#6}]]&
},

x=Join[p,Flatten[r],P,Flatten[R],v,\[Omega],V,\[CapitalOmega]];

xGeneric=Join[
{p1,p2,p3},
Flatten[Rot[\[Psi],\[Theta],\[Phi]]],
{p1,p2,p3}+l Rot[\[Psi],\[Theta],\[Phi]].e3,
Flatten[Rot[\[CapitalPsi],\[CapitalTheta],\[CapitalPhi]]],
{v1,v2,v3},
\[Omega]Rot[\[Psi],\[Theta],\[Phi],\[Omega]\[Psi],\[Omega]\[Theta],\[Omega]\[Phi]],
{v1,v2,v3}+l DRot[\[Psi],\[Theta],\[Phi],\[Omega]\[Psi],\[Omega]\[Theta],\[Omega]\[Phi]].e3,
\[Omega]Rot[\[CapitalPsi],\[CapitalTheta],\[CapitalPhi],\[CapitalOmega]\[CapitalPsi],\[CapitalOmega]\[CapitalTheta],\[CapitalOmega]\[CapitalPhi]]
];

xRandom=Join[
{p1,p2,p3},
Flatten[Rot[\[Psi],\[Theta],\[Phi]]],
{p1,p2,p3}+l Rot[\[Psi],\[Theta],\[Phi]].e3,
Flatten[Rot[\[CapitalPsi],\[CapitalTheta],\[CapitalPhi]]],
{v1,v2,v3},
\[Omega]Rot[\[Psi],\[Theta],\[Phi],\[Omega]\[Psi],\[Omega]\[Theta],\[Omega]\[Phi]],
{v1,v2,v3}+l DRot[\[Psi],\[Theta],\[Phi],\[Omega]\[Psi],\[Omega]\[Theta],\[Omega]\[Phi]].e3,
\[Omega]Rot[\[CapitalPsi],\[CapitalTheta],\[CapitalPhi],\[CapitalOmega]\[CapitalPsi],\[CapitalOmega]\[CapitalTheta],\[CapitalOmega]\[CapitalPhi]]
]/.Thread[{p1,p2,p3,\[Psi],\[Theta],\[Phi],\[CapitalPsi],\[CapitalTheta],\[CapitalPhi],v1,v2,v3,\[Omega]\[Psi],\[Omega]\[Theta],\[Omega]\[Phi],\[CapitalOmega]\[CapitalPsi],\[CapitalOmega]\[CapitalTheta],\[CapitalOmega]\[CapitalPhi]}-> RandomReal[{-1,1},18]]&;

(*input components: these are dummy names*)
u={
U(*uav's thrust input*),
\[Tau]1,\[Tau]2,\[Tau]3(*uav's body-frame torque input*),
\[Tau]m1,\[Tau]m2(*torque inputs on manipulator*)
};

{x,xGeneric,u,{p,r,P,R,v,\[Omega],V,\[CapitalOmega]},xRandom}
]
]
]

