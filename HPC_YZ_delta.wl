(* ::Package:: *)

(* ::Subsection:: *)
(*Setup:*)


ClearAll["`*"]


(* Arguements 1, 2 ,3 and 4 correspond to math -nopromt- script codename *) 
(* Every input after that is for code configuration *)


LPA=0;


inputs = {5, 6, 7, 8, 9, 10, 11}
{Ncores, Ns, dim, uor, yor, zsor} = Rationalize[Chop[$CommandLine[[#]]]]& /@ inputs
LaunchKernels[Ncores];
dims=$CommandLine[[7]];


(* Code below is for local implementation *)
(*ClearAll["`*"]
SetDirectory["/Users/gregoryjohnson/Desktop/Code/Projects/FRG_main"]
Ns=2; dim=3;
uor=5; yor=1; zsor=2;
dims=3;
Ncores=8;*)


(* ::Subsubsection::Closed:: *)
(*Generate File Names:*)


(* ::Input:: *)
(*(* Theory and Truncation String *)*)


TheoryString = StringJoin[ ToString /@ {"N_", N\[Phi]s, "_d_", dims, "_a_", as, "_Uor_", uor, "_Yor_", yor, "_Zsor_", zsor}];
TheoryString = StringReplace[ TheoryString, {"."->"p", "-"->"m"}];
TheoryString = StringJoin[ TheoryString, ".wl"];


FileNameImportRGeqns = StringJoin[{"RGeqns_",TheoryString}];
FileNameImport\[Phi]RGeqns = StringJoin[{"phi_RGeqns_",TheoryString}];
FileNameImportFPdata = StringJoin[{"FPdata_",TheoryString}];
FileNameOutputCryWF = StringJoin[{"WFCriticalDataYLE_",ToString[Ns],TheoryString}];


(* ::Subsubsection::Closed:: *)
(*Import FRG:*)


Print["Just before FRG Import"]


Import[FileNameJoin[{Directory[],"Packages/ONvars_YZ.wl"}]]
Import[FileNameJoin[{Directory[], "RGeqns/ON_YLE_flow/", FileNameImportRGeqns}]];
Import[FileNameJoin[{Directory[],"FPData/ON_YLE_flow/", FileNameImportFPdata}]];


FPIntNap5=Cases[FPIntNap5,{_Integer,__}] //NumericalSort;


apmsINTdata=Table[{i,1/2},{i,1,10}]
(* \[Eta]pms vs \[CapitalDelta]pms *)
(* Turn for \[CapitalDelta]pms - remove line for \[Eta]pms *)
(*apmsINTdata=apms\[CapitalDelta]Data*)
(* Turn for Hcpms - remove line for \[Eta]pms *)
(*apmsINTdata=HcPMSData*)


(* ::Subsection::Closed:: *)
(*Solvers:*)


Print["Just before Solvers"]


(* ::Subsubsection::Closed:: *)
(*IC Finder from FP:*)


(* ::Item:: *)
(*Construct fixed point solver. Use \[Alpha]=1/2 solution for initial conditions to recompute the solution for whatever \[Alpha] value the code may need.*)


(* ::Item:: *)
(*Detail the initial condition function*)


(* ::Item:: *)
(*Do i use fixed mass or fixed renormalized mass*)


WFSolver[uo_, zo_, yo_?NumericQ, Ms_, as_, wp_, prec_ ]:=Module[{(*Guess,*)VarInit,VarS}, 
VarInit=Drop[Table[{VarList[[i]],Guess[[i]]},{i,1,Length[VarList]}],{2}];
VarS=FindRoot[Table[(\[Beta]FYLE2[index,dim,VarList, \[Eta]Fys2[dim,VarList,Ms,as],Ms,as]/. Subscript[u, 1]->-2 \[Rho] Subscript[u, 2])==0,{index,1,Length[VarList]-1}], VarInit,
WorkingPrecision->wp,PrecisionGoal->prec,AccuracyGoal->prec,MaxIterations->1 10^3];
Guess=VarList /. VarS;
VarS
]


InitialConditionsFunc[\[Delta]t_?NumericQ,m2_?NumericQ, Ms_(*, as_*),wp_,prec_]:=
Module[{EqPrep,Rho0,InitialC,Initial,rhofp,EqPrepZs,EqPrepY,R,PotentialNearYLE,InitialZs, InitialY, VarS, as},

(* Finding fixed point PMS solution for N at given truncation using PMS alpha values*)
Guess=SetPrecision[VarList /.FPIntNap5[[Ms,3]], 100]; 
as=SetPrecision[apmsINTdata[[Ms,2]],150];
VarS=WFSolver[uor, yor, zsor, Ms, as, 100, 100];

(* Add in thermal perturbation and re-expand functions about field value with fixed mass of m2 *)
(* Return the set of initial conditions for the flows *)
PotentialNearYLE=Sum[((R-\[Rho])^i Subscript[u, i])/i!,{i,1,uor}]//.Subscript[u, 1]->-2 \[Rho] Subscript[u, 2]/.VarS//Expand;
rhofp=\[Rho]/.VarS;
EqPrep=PotentialNearYLE + \[Delta]t R;
EqPrepY=Sum[Subscript[y, i](R-\[Rho])^i/(i!),{i,0,yor}] /.VarS;
EqPrepZs=Sum[Subscript[zs, i](R-\[Rho])^i/(i!),{i,0,zsor}] /.VarS;
Rho0=R/.NSolve[D[EqPrep,R]+2 R D[EqPrep,{R,2}]==m2 && rhofp/5<R<rhofp*5, R, Reals, WorkingPrecision->wp][[1]];
InitialC=Table[Subscript[u, i][0]==D[EqPrep,{R,i}]/.R->Rho0,{i,2,uor}];
InitialY=Table[Subscript[y, i][0]==D[EqPrepY,{R,i}]/.R->Rho0,{i,0,yor}];
InitialZs=Table[Subscript[zs, i][0]==D[EqPrepZs,{R,i}]/.R->Rho0,{i,1,zsor}];
Initial=If[LPA==1,
	Flatten[Join[{\[Rho][0]==Rho0}, InitialC ]],
	Flatten[Join[{\[Rho][0]==Rho0}, InitialC, InitialY, InitialZs]]];
Initial
]


(* ::Subsubsection::Closed:: *)
(*YLE Solver 2 WF:*)


YLESolver2[t0_, rin_, m_, dim_, wp_, prec_, Ms_]:=Module[{d,Vars, Varsoft,ReducedVars,ReducedVarsoft,SetVariables, BSOE,IV,tm, as},

(* t0 is IR cutoff, rin is thermal perturbation, m is mass perturbation m=Subscript[m, s]^2, dim is dimension, wp and prec are precisions, Ms is Number of field components *)

Off[NIntegrate::precw]; Off[NDSolve::precw];
Vars=VarList;
Varsoft=Table[VarList[[i]][t],{i,1,Length[VarList]}] ;
SetVariables=Table[Vars[[i]]-> Varsoft[[i]],{i,1,Length[Vars]}];
d=dim;

ReducedVars=Drop[Vars,{2}];
ReducedVarsoft=Drop[Varsoft,{2}];

BSOE=Join[Table[D[ ReducedVarsoft[[i]],t]==(\[Beta]FYLE2[i,d, Vars, \[Eta]Fys2[d,Vars,Ms,as],Ms,as] /. SetVariables),{i,1,Length[ReducedVarsoft]}],
{Z'[t]==-(\[Eta]Fys2[d,Vars,Ms,as]/. SetVariables) Z[t]}] /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]};

If[Ms==1, If[LPA==0, BSOE=Drop[BSOE,{uor+1,uor+1+yor}]]]; (* At N=1, flows for Y decople. So, drops these flows as they are unneeded *)

(*If[LPA==1,ReducedVars=ReducedVars, 
ReducedVars=Drop[ReducedVars,{-1-(yor)}];
BSOE=Drop[BSOE,{-1-(yor+1)}];];
ReducedVarsoft=Table[ReducedVars[[i]][t],{i,1,Length[ReducedVars]}];*)

tm=0;
as=SetPrecision[apmsINTdata[[Ms,2]],150];
IV=Join[InitialConditionsFunc[rin,m, Ms, 100, 100],{Z[tm]==1}];
Broken=NDSolve[ {BSOE, IV}, 
Join[ReducedVars,{Z}], {t, tm, t0} ,WorkingPrecision->wp,PrecisionGoal-> prec, AccuracyGoal->prec, Method->{"StiffnessSwitching" ,"EquationSimplification"->"Solve"} ,
MaxSteps-> 2 10^3
]
];


(* ::Subsection::Closed:: *)
(*Computing critical data related to delta scaling form:*)


(*wp=80; prec=wp/2;*)
wp=60; prec=wp/2;
$MinPrecision=prec;
ParallelEvaluate[$MinPrecision=prec];


(*ParallelEvaluate[Off[FindRoot::precw]];*)


Off[Precision::precsm,FittedModel::precw,NonlinearModelFit::precw,FindRoot::precw, NDSolve::ndsz]
ParallelEvaluate[Off[Precision::precsm,FittedModel::precw,NonlinearModelFit::precw,FindRoot::precw,NDSolve::ndsz]];


Off[N::precsm,NumericalMath`LowerPrecision::precsm]


(* ::Subsubsection::Closed:: *)
(*Finding \[Delta]:  *)


Print["Just before delta"]


t0=-40;


dt0=0; si=-20; sf=-16; ss=(sf-si)/((*2*) Ncores-1);


{time,\[Delta]\[Chi]Data}=ParallelTable[{Log10[Sqrt[2\[Rho][t0]] ((10^s Exp[-2t0])(*/Z[t0]*) -2 \[Rho][t0] Subscript[u, 2][t0])Sqrt[Z[t0]]Exp[t0 (dim+2)/2]],
Log10[Sqrt[(2 \[Rho][t0]Exp[t0]^(dim-2))/Z[t0]]],
Log10[ Sqrt[(10^s) /Z[t0]]], Log10[10^s]}
/. YLESolver2[t0, dt0, 10^s,dim,wp,prec, Ns][[1]],{s,si,sf,ss}] //AbsoluteTiming;
Print["\[Delta] time=", time]


\[Delta]Data=Table[{\[Delta]\[Chi]Data[[i,1]],\[Delta]\[Chi]Data[[i,2]]},{i,1,Length[\[Delta]\[Chi]Data]}];
\[Xi]Data=Table[{\[Delta]\[Chi]Data[[i,1]],\[Delta]\[Chi]Data[[i,3]]},{i,1,Length[\[Delta]\[Chi]Data]}];
\[Chi]Data=Table[{\[Delta]\[Chi]Data[[i,1]],\[Delta]\[Chi]Data[[i,4]]},{i,1,Length[\[Delta]\[Chi]Data]}];


Clear[\[Delta]1,H0]
FitD2=NonlinearModelFit[\[Delta]Data,1/\[Delta]1*x + Log10[Bc0],{\[Delta]1,Bc0},x];
\[Delta]=\[Delta]1/.FitD2["BestFitParameters"];
Bc=Bc0 /.FitD2["BestFitParameters"];


FitD2=NonlinearModelFit[\[Delta]Data,1/\[Delta]1*x + Log10[Bc0]+Log10[1+d 10^(t x)],{{\[Delta]1,\[Delta]},{Bc0,Bc},{d,0},{t,3/10}},x, AccuracyGoal->prec,PrecisionGoal->prec, WorkingPrecision->wp, MaxIterations->10000];
Show[ListPlot[\[Delta]Data],Plot[FitD2[x],{x,-40,0}]];
\[Delta]=\[Delta]1/.FitD2["BestFitParameters"];
Bc=Bc0 /.FitD2["BestFitParameters"];
\[Theta]\[Delta]=t/.FitD2["BestFitParameters"];
a\[Delta]=d /.FitD2["BestFitParameters"];


FitD2=NonlinearModelFit[\[Delta]Data,1/\[Delta]1*x + Log10[Bc0]+Log10[1+d 10^(t x)+d2 10^( 2 t x)],{{\[Delta]1,\[Delta]},{Bc0,Bc},{d,0},{d2,0},{t,3/10}},
x, AccuracyGoal->prec,PrecisionGoal->prec, WorkingPrecision->wp, MaxIterations->10000];
FitD2["ParameterConfidenceIntervalTable"]
\[Delta]=\[Delta]1/.FitD2["BestFitParameters"];
Bc=Bc0 /.FitD2["BestFitParameters"];
\[Theta]\[Delta]=t/.FitD2["BestFitParameters"];
a\[Delta]=d /.FitD2["BestFitParameters"];
a2\[Delta]=d2 /.FitD2["BestFitParameters"];


(* ::Subsection:: *)
(*Saving:*)


Save[FileNameJoin[{Directory[], "CritData_delta/", FileNameOutputCryWF}], {\[Delta]Data,\[Xi]Data,\[Chi]Data}]


Save[FileNameJoin[{Directory[], "CritData_delta/", FileNameOutputCryWF}], {\[Delta],Bc,\[Theta]\[Delta],a\[Delta],a2\[Delta] }]
