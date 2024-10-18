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


(* ::Subsubsection:: *)
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


(* ::Subsection:: *)
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
(*He02Solver2 Pre WF:*)


(* ::Item:: *)
(*This solver is the YLE solver, but with WhenEvents to catch for specific outcomes - if \[Rho][t] becomes too large (>-\[Delta]t/2)*)


(* ::Item:: *)
(*It additionally appends to McSet and McBounds the value of m^2*)


ParallelEvaluate[Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"]];
Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"]
He02Solver2Pre[t0_, rin_?NumericQ, m_?NumericQ,dim_,wp_,prec_,Ms_]:=Module[{as, d,Vars, Varsoft,ReducedVars,ReducedVarsoft,SetVariables, BSOE,IV,Status,tm, hout,out},
Off[NIntegrate::precw]; Off[NDSolve::precw];
Vars=VarList;
Varsoft=Table[VarList[[i]][t],{i,1,Length[VarList]}] ;
SetVariables=Table[Vars[[i]]-> Varsoft[[i]],{i,1,Length[Vars]}];
d=dim;

ReducedVars=Drop[Vars,{2}];
ReducedVarsoft=Drop[Varsoft,{2}];

BSOE=Join[Table[D[ ReducedVarsoft[[i]],t]==(\[Beta]FYLE2[i,d, Vars, \[Eta]Fys2[dim,Vars,Ms,as],Ms,as] /. SetVariables),{i,1,Length[ReducedVarsoft]}],
{Z'[t]==-(\[Eta]Fys2[d,Vars,Ms,as]/. SetVariables) Z[t]}] /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]};

If[Ms==1, If[LPA==0, BSOE=Drop[BSOE,{uor+1,uor+1+yor}]]]; (* At N=1, flows for Y decople. So, drops these flows as they are unneeded *)

(*If[LPA==1,ReducedVars=ReducedVars, 
ReducedVars=Drop[ReducedVars,{-1-(yor)}];
BSOE=Drop[BSOE /.{Subscript[a, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[a, 2][t],Subscript[c, 0][t]->(1-Subscript[b, 0][t])/\[Rho][t]},{-1-(yor+1)}];];
ReducedVarsoft=Table[ReducedVars[[i]][t],{i,1,Length[ReducedVars]}];*)

tm=0; Status=0;

as=SetPrecision[apmsINTdata[[Ms,2]],150];
IV=Join[InitialConditionsFunc[rin,m, Ms, 150, 150],{Z[tm]==1}];

{time,Broken}=NDSolve[ { BSOE, IV,
(*WhenEvent[\[Rho][t]<-10^-5(*-10^-20*)(*-10*) ,  tm=t; Status=1; "StopIntegration","LocationMethod"->"LinearInterpolation"],*)
WhenEvent[\[Rho][t]>2 (-rin)/ 2  ,  tm=t; Status=2; "StopIntegration","LocationMethod"->"LinearInterpolation"]}, 
Join[ReducedVars,{Z}], {t, tm, t0} ,WorkingPrecision->wp,PrecisionGoal-> prec, AccuracyGoal->prec, 
Method->{"StiffnessSwitching" ,"EquationSimplification"->{Automatic, "TimeConstraint"->\[Infinity]}}, 
MaxSteps-> 10^3]//AbsoluteTiming;

tm=InterpolatingFunctionDomain[First[\[Rho] /. Broken]][[1,1]];
hout=Sqrt[2\[Rho][t]] ((m Exp[-2t]) -2 \[Rho][t] Subscript[u, 2][t])Sqrt[Z[t]]Exp[t (dim+2)/2] /. Broken[[1]] /. t-> tm;
out= Re[Sqrt[(2 \[Rho][t]Exp[t]^(dim-2))/Z[t]]]-Im[Sqrt[(2 \[Rho][t]Exp[t]^(dim-2))/Z[t]]] /. Broken[[1]] /. t->tm;

If[hout \[Element] Reals, AppendTo[McSet, {m,1}],AppendTo[McSet, {m,-1}]];
If[hout \[Element] Reals, McBounds[[3]]=m, McBounds[[2]]= m];

hout=Re[hout]-Im[hout]; (* hout is either real and positive or purely imaginary. If it is purely imaginary, flip the sign. This yields an aove/below zero search *)
(*Echo[{time,tm,hout,out},"{ time, tend, H, \[Phi] } = "];*)
If[v==1,Echo[{time, SetPrecision[m,wp], tm //N,hout //N},"{time, \!\(\*SuperscriptBox[\(m\), \(2\)]\), tend, H } = "]];
hout
  ];


(* ::Section:: *)
(*Crit Data:*)


(*Off[Precision::precsm,FittedModel::precw,NonlinearModelFit::precw,FindRoot::precw, NDSolve::ndsz,NDSolve::mxst ]
ParallelEvaluate[Off[Precision::precsm,FittedModel::precw,NonlinearModelFit::precw,FindRoot::precw, NDSolve::ndsz,NDSolve::mxst ]];*)


(*Off[N::precsm,NumericalMath`LowerPrecision::precsm]*)


(* ::Subsection:: *)
(*Finding mcdata :*)


Print["Just before mcTest"]


si=If[Ns==2, -9, -10]; sf=-8; ss=(sf-si)/((*2*) Ncores-1);
srange=Table[i,{i,si,sf,ss}];
t0=-40;


If[Ns>5, si=-13; sf=-11; ss=(sf-si)/((*2*) Ncores-1);
srange=Table[i,{i,si,sf,ss}];
t0=-50;]


(* ::Subsubsection:: *)
(*Once through:*)


aspms=apmsINTdata[[Ns,2]];
Print[aspms]


Print[" Beginning Step 1"]


wp=60; prec=wp/2;
$MinPrecision=prec;
ParallelEvaluate[$MinPrecision=prec];


FRwp=50; FRprec=FRwp/2;


DistributeDefinitions[wp,prec];
ParallelTable[McSet={ 10^(s) },{s,si,sf,ss}];
ParallelTable[McBounds={ 10^(s), 10^-20, 10^-5 },{s,si,sf,ss}]


{time,mcTestTot}=ParallelTable[{10^s, (m /.FindRoot[He02Solver2Pre[t0, 10^s,m, dim, wp, prec, Ns]==(*If[Ns>3, 10^(2 s), 0]*)10^(2 s)(*0*),{m,McBounds[[2]], McBounds[[3]]},
WorkingPrecision->FRwp, PrecisionGoal->(*prec*)Infinity, AccuracyGoal->FRprec, Method->"Brent"(*, EvaluationMonitor:>Print[m]*)])//Quiet},{s,si,sf,ss}]//AbsoluteTiming;
Print["mc time = ", time/60. , " minutes, or ", time/60./60, " hours " ]


McSetTot=ParallelEvaluate[Length[McSet]];
(*Print[ParallelEvaluate[Print[" For t=", McSet[[1]], " Evaluations = ",Length[McSet]];];];*)
Print[ParallelEvaluate[Length[McSet]]];


McBoundsTot=ParallelEvaluate[McBounds];
McColl={FRprec, time, aspms, McBoundsTot, McSetTot, mcTestTot};


(* ::Subsection:: *)
(*Saving:*)


Save[FileNameJoin[{Directory[], "CritData_Hc/", FileNameOutputCryWF}], McColl]


(* ::Subsection:: *)
(*Gamma computations using flows or the mc data (commented out):*)


(* ::Subsubsection::Closed:: *)
(*Finding \[Gamma] data:*)


(*Print["Just before gamma data"]*)


(*Import[FileNameJoin[{Directory[], "CritData2/", FileNameOutputCryWF}]];*)


(*wp=40; prec=wp/2;
$MinPrecision=prec;
ParallelEvaluate[$MinPrecision=prec];
v=0;*)


(*mcTestTot=SetPrecision[McColl[[-1,-1]],wp]*)


(*t\[Eta]=-5;
mdata=Table[(*m /.*) mcTestTot[[i,2]],{i,1,Length[mcTestTot]}];
mdata=SetPrecision[mdata,wp]*)


(*{time,\[Nu]\[Gamma]Data}=ParallelTable[{si+i ss,Log[E^t0 Sqrt[Subscript[u, 1][t0]+2 \[Rho][t0] Subscript[u, 2][t0]]]/Log[10],
Log[E^(2 t0) Z[t0] (Subscript[u, 1][t0]+2 \[Rho][t0] Subscript[u, 2][t0])]/Log[10],
Log[(3 (E^t0)^(-4-dim) Subscript[u, 2][t0])/(Z[t0]^2 Subscript[u, 1][t0]^4)]/Log[10],
Log[((E^t0)^(-6-2 dim) (90 Subscript[u, 2][t0]^2-15 Subscript[u, 1][t0] Subscript[u, 3][t0]))/(Z[t0]^3 Subscript[u, 1][t0]^7)]/Log[10],
Log[(105 (E^t0)^(-8-3 dim) (72 Subscript[u, 2][t0]^3-24 Subscript[u, 1][t0] Subscript[u, 2][t0] Subscript[u, 3][t0]+Subscript[u, 1][t0]^2 Subscript[u, 4][t0]))/(Z[t0]^4 Subscript[u, 1][t0]^10)]/Log[10],
Log[(945 (E^t0)^(-10-4 dim) (1320 Subscript[u, 2][t0]^4-660 Subscript[u, 1][t0] Subscript[u, 2][t0]^2 Subscript[u, 3][t0]+40 Subscript[u, 1][t0]^2 Subscript[u, 2][t0] Subscript[u, 4][t0]+Subscript[u, 1][t0]^2 (30 Subscript[u, 3][t0]^2-Subscript[u, 1][t0] Subscript[u, 5][t0])))/(Z[t0]^5 Subscript[u, 1][t0]^13)]/Log[10],
(*\[Eta]Fys2[dim,Table[VarList[[i]][t\[Eta]] ,{i,1,Length[VarList]}], Ns, as]*)-Z'[t\[Eta]]/Z[t\[Eta]],
Sqrt[2] E^(1/2 (2+dim) t0) Sqrt[Z[t0]] Sqrt[\[Rho][t0]] Subscript[u, 1][t0]}
/. Subscript[u, 1][t_]->(E^(-2 t) mdata[[i+1]])(*/Z[t]*)-2 \[Rho][t] Subscript[u, 2][t]  
/. YLESolver2[t0, (*rc(1-10^(si+i*ss))*)10^(si+i ss), mdata[[i+1]], dim, wp, prec, Ns][[1]],{i,0,Length[srange]-1} ] //AbsoluteTiming;
Print["nutime=",time]*)


(*hData=Table[\[Nu]\[Gamma]Data[[i,9]],{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Nu]Data=Table[{\[Nu]\[Gamma]Data[[i,1]],\[Nu]\[Gamma]Data[[i,2]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Gamma]Data=Table[{\[Nu]\[Gamma]Data[[i,1]],\[Nu]\[Gamma]Data[[i,3]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Gamma]4Data=Table[{\[Nu]\[Gamma]Data[[i,1]],\[Nu]\[Gamma]Data[[i,4]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Gamma]6Data=Table[{\[Nu]\[Gamma]Data[[i,1]],\[Nu]\[Gamma]Data[[i,5]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Gamma]8Data=Table[{\[Nu]\[Gamma]Data[[i,1]],\[Nu]\[Gamma]Data[[i,6]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Gamma]10Data=Table[{\[Nu]\[Gamma]Data[[i,1]],\[Nu]\[Gamma]Data[[i,7]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Eta]FData=Table[{\[Nu]\[Gamma]Data[[i,8]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
\[Eta]f=\[Eta]FData[[1]];*)


(*Print[(ParallelTable[10^(2 s),{s,si,sf,ss}] -Abs[hData])/ParallelTable[10^(2 s),{s,si,sf,ss}] // Log10[Abs[#]]& ]*)


(*Corr\[Gamma]Data=Table[{\[Nu]\[Gamma]Data[[i,1]],Log10[Abs[\[Nu]\[Gamma]Data[[i,9]]]], \[Nu]\[Gamma]Data[[i,3]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
Corr\[Nu]Data=Table[{\[Nu]\[Gamma]Data[[i,1]],Log10[Abs[\[Nu]\[Gamma]Data[[i,9]]]], \[Nu]\[Gamma]Data[[i,2]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
Corr\[Gamma]4Data=Table[{\[Nu]\[Gamma]Data[[i,1]],Log10[Abs[\[Nu]\[Gamma]Data[[i,9]]]],\[Nu]\[Gamma]Data[[i,4]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
Corr\[Gamma]6Data=Table[{\[Nu]\[Gamma]Data[[i,1]],Log10[Abs[\[Nu]\[Gamma]Data[[i,9]]]],\[Nu]\[Gamma]Data[[i,5]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
Corr\[Gamma]8Data=Table[{\[Nu]\[Gamma]Data[[i,1]],Log10[Abs[\[Nu]\[Gamma]Data[[i,9]]]],\[Nu]\[Gamma]Data[[i,6]]},{i,1,Length[\[Nu]\[Gamma]Data]}];
Corr\[Gamma]10Data=Table[{\[Nu]\[Gamma]Data[[i,1]],Log10[Abs[\[Nu]\[Gamma]Data[[i,9]]]],\[Nu]\[Gamma]Data[[i,7]]},{i,1,Length[\[Nu]\[Gamma]Data]}];*)


(* ::Subsubsubsection::Closed:: *)
(*Gamma:*)


(*Off[N::precsm,NumericalMath`LowerPrecision::precsm]*)


(*Clear[\[Gamma],Cp0];
FitG2=NonlinearModelFit[(*({#[[1]],-Abs[#[[2]]]}&/@ \[Gamma]Data[[1;;-1]])*)\[Gamma]Data[[1;;-1]],\[Gamma]0*x-Log10[Cp0],{\[Gamma]0,{Cp0,10^-4}},x,PrecisionGoal->prec, AccuracyGoal->prec, WorkingPrecision->wp, MaxIterations->10000];
(*Export["gammaplot.jpg",Show[Plot[FitG[x],{x,si-1,sf+1}], ListPlot[\[Gamma]Data] ] ];*)
(*FitG2["ParameterConfidenceIntervalTable"]*)
\[Gamma]=\[Gamma]0/.FitG2["BestFitParameters"];
Cp=Cp0/.FitG2["BestFitParameters"];*)


(*FitG2["ParameterConfidenceIntervalTable"]*)


(*\[Beta]=SetPrecision[\[Gamma]/(\[Delta]-1),wp] /. \[Delta]->5;
\[CapitalDelta]=SetPrecision[\[Beta] \[Delta],wp]/. \[Delta]->5;*)


(*FitG2=NonlinearModelFit[{#[[1]],-Abs[#[[2]]]}&/@ \[Gamma]Data[[1;;-1]],\[Gamma]1*x-Log10[Cp0]-Log10[1 + g 10^(t x)],{{\[Gamma]1,\[Gamma]},{Cp0,Cp},{g,0},{t,1/2}},x,
PrecisionGoal->prec, AccuracyGoal->prec, WorkingPrecision->wp, MaxIterations->10000];
(*Export["gammaplot.jpg",Show[Plot[FitG[x],{x,si-1,sf+1}], ListPlot[\[Gamma]Data] ] ];*)
(*FitG2["ParameterConfidenceIntervalTable"]*)
\[Gamma]=\[Gamma]1/.FitG2["BestFitParameters"];
Cp=Cp0/.FitG2["BestFitParameters"];
\[Theta]\[Gamma]=t/.FitG2["BestFitParameters"];
g2=g/.FitG2["BestFitParameters"];*)


(*FitG2["ParameterConfidenceIntervalTable"]*)


(*FitG2=NonlinearModelFit[Corr\[Gamma]Data[[1;;-1]],\[Gamma]1*x - Log10[Cp0]-Log10[1+f (10^ y 10^(- \[CapitalDelta]0 x))^2+ g2 10^(t x)],{{\[Gamma]1,\[Gamma]},{Cp0,Cp},{g,g2},{t,\[Theta]\[Gamma]},{f,0}, {\[CapitalDelta]0,\[CapitalDelta]}},{x,y}, 
PrecisionGoal->prec, AccuracyGoal->prec, WorkingPrecision->wp, MaxIterations->10^4];
FitG2["ParameterConfidenceIntervalTable"]
\[Gamma]=\[Gamma]1/.FitG2["BestFitParameters"];
Cp=Cp0/.FitG2["BestFitParameters"];
\[Theta]\[Gamma]=t/.FitG2["BestFitParameters"];
g2=g/.FitG2["BestFitParameters"];
f0=f/.FitG2["BestFitParameters"];
f20=f2 /.FitG2["BestFitParameters"];
NumberForm[dt Hc^((1-\[Delta])/(\[Gamma] \[Delta])) (Bc/Cp)^(1/\[Gamma]),5]*)


(* ::Subsubsection::Closed:: *)
(*Using mcTest:*)


(*Fitdata=Log10[{#[[1]],#[[2]]}]&/@mcTestTot;
% //ListPlot*)


(*FitG2=NonlinearModelFit[Fitdata,\[Gamma]1*x-Log10[Cp0],{\[Gamma]1,{Cp0,10^-4}},x,PrecisionGoal->prec, AccuracyGoal->prec, WorkingPrecision->wp, MaxIterations->10000]
(*Export["gammaplot.jpg",Show[Plot[FitG[x],{x,si-1,sf+1}], ListPlot[\[Gamma]Data] ] ];*)
FitG2["ParameterConfidenceIntervalTable"]
\[Gamma]a=\[Gamma]1/.FitG2["BestFitParameters"];
Cpa=Cp0/.FitG2["BestFitParameters"];*)


(*FitG2=NonlinearModelFit[Fitdata,\[Gamma]1*x-Log10[Cp0]-Log10[1 + g 10^(t x)],{{\[Gamma]1,\[Gamma]a},{Cp0,Cpa},{g,0},{t,1/2}},x,
PrecisionGoal->prec, AccuracyGoal->prec, WorkingPrecision->wp, MaxIterations->10000];
(*Export["gammaplot.jpg",Show[Plot[FitG[x],{x,si-1,sf+1}], ListPlot[\[Gamma]Data] ] ];*)
FitG2["ParameterConfidenceIntervalTable"]
\[Gamma]a=\[Gamma]1/.FitG2["BestFitParameters"];
Cpa=Cp0/.FitG2["BestFitParameters"];
\[Theta]\[Gamma]a=t/.FitG2["BestFitParameters"];
g2a=g/.FitG2["BestFitParameters"];*)
