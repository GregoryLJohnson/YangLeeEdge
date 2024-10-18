(* ::Package:: *)

(* ::Subsection:: *)
(*Setup:*)


ClearAll["`*"]


(*ResetDirectory[];
SetDirectory["HPCpack"];*)


(* Arguements 1, 2 ,3 and 4 correspond to math -nopromt- script codename *) 


offset1=1;
offset2=2;


pre1=4+offset1;
pre2=5+offset1;


LPAp=1; LPApp=1;
LPA=0;


Ncores=(*8*)Rationalize[Chop[ToExpression[$CommandLine[[pre1]]]]];


LaunchKernels[Ncores];


first=4+offset2;
second=5+offset2;
third=6+offset2;
fourth=7+offset2;
fifth=8+offset2;
sixth=9+offset2;


Ns=Rationalize[Chop[ToExpression[$CommandLine[[first]]]]];
dim=Rationalize[Chop[ToExpression[$CommandLine[[second]]]]];
uor=Rationalize[Chop[ToExpression[$CommandLine[[third]]]]];
yor=Rationalize[Chop[ToExpression[$CommandLine[[fourth]]]]];
zsor=Rationalize[Chop[ToExpression[$CommandLine[[fifth]]]]];
(*aR=Rationalize[Chop[ToExpression[$CommandLine[[sixth]]]]];*)
(*\[Alpha]=1;*)


(*Print[Ncores]
Print[Ns, dim, uor, yor, zsor]*)


(*Print[Rationalize[Chop[ToExpression[$CommandLine[[first]]]]]]*)


(*Print[$HomeDirectory," and ", Directory[]]*)


(*N\[Phi]s=$CommandLine[[first]];*)
dims=$CommandLine[[second]];
(*\[Alpha]s=$CommandLine[[sixth]];*)


(*ClearAll["`*"]
ResetDirectory[];
SetDirectory["HPCpack"];
Ns=1; dim=3;
uor=6; yor=2; zsor=3;
dims=3;
LPA=0;
(*Ncores=8;*)*)


Off[FindRoot::precw,Precision::precsm]


(* ::Subsection::Closed:: *)
(*Generate File Names:*)


(* ::Input:: *)
(*(* Theory and Truncation String *)*)


(*TheoryString=If[ZsExpansion==0,
StringJoin["N_d_",ToString[dims], "_Uor_",ToString[uor],"_Zor_",ToString[zor],"_Yor_",ToString[yor](*,".wl"*)],
StringJoin["N_d_",ToString[dims],"_Uor_",ToString[uor],"_Zor_",ToString[zor],"_Zsor_",ToString[zsor](*,".wl"*)]];
If[Scalar==1, TheoryString=StringJoin["N_d_",ToString[dims],"_Uor_",ToString[uor],"_Zor_",ToString[zor](*,".wl"*)]];*)


TheoryString=StringJoin["N_",ToString[N\[Phi]s],"_d_",ToString[dims],"_a_",ToString[as], "_Uor_",ToString[uor],"_Yor_",ToString[yor],"_Zsor_",ToString[zsor](*,".wl"*)];


TheoryString=StringJoin[StringReplace[TheoryString,{"."->"p","-"->"m"}],".wl"];


FileNameImportRGeqns=StringJoin[{"RGeqns_",TheoryString}]
FileNameImport\[Phi]RGeqns=StringJoin[{"phi_RGeqns_",TheoryString}]
FileNameImportFPdata=StringJoin[{"FPdata_",TheoryString}]


FileNameOutputCryWF=StringJoin[{"WFCriticalDataYLE_",ToString[Ns],TheoryString}]
FileNameOutputYLEWF=StringJoin[{"WFYLE_",ToString[Ns],TheoryString}]


Print[TheoryString]


(* ::Subsection::Closed:: *)
(*Import FRG:*)


Import[FileNameJoin[{Directory[],"Packages/ONvars_YZ.wl"}]]


Import[FileNameJoin[{Directory[], "RGeqns/", FileNameImportRGeqns}]];
(*Import[FileNameJoin[{Directory[], "RGeqns/", FileNameImport\[Phi]RGeqns}]];*)


Import[FileNameJoin[{Directory[], "RGeqns_R/", FileNameImport\[Phi]RGeqns}]];


Import[FileNameJoin[{Directory[],"FPData/", FileNameImportFPdata}]];
Import[FileNameJoin[{Directory[],"FPData_1/", FileNameImportFPdata}]];


apmsINTdata


apms\[CapitalDelta]Data


HcPMSdata
HcPMSData={#[[1]],#[[2]]} &/@HcPMSdata


FPIntNap5=Cases[FPIntNap5,{_Integer,__}] //NumericalSort;


(* \[Eta]pms vs \[CapitalDelta]pms *)
(* Turn for \[CapitalDelta]pms - remove line for \[Eta]pms *)
(*apmsINTdata=apms\[CapitalDelta]Data*)
(* Turn for Hcpms - remove line for \[Eta]pms *)
apmsINTdata=HcPMSData


(* ::Subsection::Closed:: *)
(*Solvers:*)


SetSystemOptions["NDSolveOptions" -> "DefaultSolveTimeConstraint" -> 100.`]
SetSystemOptions["NDSolveOptions" -> "DefaultScanDiscontinuityTimeConstraint" -> 0.1`]


(* ::Subsubsection::Closed:: *)
(*IC Finder from FP:*)


WFSolver[uo_, zo_, yo_?NumericQ, Ms_, as_, wp_, prec_ ]:=Module[{(*Guess,*)VarInit,VarS}, 
VarInit=Drop[Table[{VarList[[i]],Guess[[i]]},{i,1,Length[VarList]}],{2}];
VarS=FindRoot[Table[(\[Beta]FYLE2[index,dim,VarList, \[Eta]Fys2[dim,VarList,Ms,as],Ms,as]/. Subscript[u, 1]->-2\[Rho] Subscript[u, 2])==0,{index,1,Length[VarList]-1}],
VarInit,WorkingPrecision->wp,PrecisionGoal->prec,AccuracyGoal->prec,MaxIterations->1 10^3];
Guess=VarList /. VarS;
VarS
]


InitialConditionsFunc[\[Delta]t_?NumericQ,m2_?NumericQ, Ms_(*, as_*),wp_,prec_]:=
Module[{EqPrep,Rho0,InitialC,Initial,rhofp,EqPrepZs,EqPrepY,R,PotentialNearYLE,InitialZs, InitialY, VarS, as},

Guess=SetPrecision[VarList /.FPIntNap5[[Ms,3]],
100]; (* Finding fixed point PMS solution for N at given truncation*)
as=SetPrecision[apmsINTdata[[Ms,2]],150];
VarS=WFSolver[uor, yor, zsor, Ms, as, 100, 100];

PotentialNearYLE=Sum[Subscript[u, i](R-\[Rho])^i/(i!),{i,1,uor}]//.Subscript[u, 1]->-2\[Rho] Subscript[u, 2]/.VarS//Expand;
rhofp=\[Rho]/.VarS;
EqPrep=PotentialNearYLE + \[Delta]t R;
EqPrepY=Sum[Subscript[y, i](R-\[Rho])^i/(i!),{i,0,yor}](*/.Subscript[c, 0]->(1-Subscript[b, 0])/\[Rho]*)/.VarS;
EqPrepZs=Sum[Subscript[zs, i](R-\[Rho])^i/(i!),{i,0,zsor}]/.VarS;
Rho0=R/.NSolve[D[EqPrep,R]+2 R D[EqPrep,{R,2}]==m2 && rhofp/5<R<rhofp*5, R, Reals, WorkingPrecision->wp][[1]];
InitialC=Table[Subscript[u, i][0]==D[EqPrep,{R,i}]/.R->Rho0,{i,2,uor}];
InitialY=Table[Subscript[y, i][0]==D[EqPrepY,{R,i}]/.R->Rho0,{i,0,yor}];
InitialZs=Table[Subscript[zs, i][0]==D[EqPrepZs,{R,i}]/.R->Rho0,{i,1,zsor}];
Initial=If[LPA==1,
Flatten[Join[{\[Rho][0]==Rho0},InitialC (*,InitialY,InitialZs*)]],
Flatten[Join[{\[Rho][0]==Rho0},InitialC,InitialY,InitialZs]]];
Initial
]


(* ::Input:: *)
(*(*VarList*)
(*wp=60;prec=60;*)
(*Ns=1; as=Rationalize[apmsINTdata[[Ns,2]],10^-prec];*)
(*Guess=VarList /. FPIntNap5[[Ns+3,3]];*)
(*WFSolver[uor, yor,zsor, Ns, as, wp, prec]//AbsoluteTiming*)
(*Clear[Ns,as, wp, prec]*)*)


(* ::Input:: *)
(*(*wp=60;prec=60;*)
(*Ns=3;*)
(*dt=10^-20; m1=10^-2;*)
(*InitialConditionsFunc[dt,m1,Ns,wp,prec]*)
(*Clear[Ns,as, wp, prec]*)*)


(* ::Subsubsection::Closed:: *)
(*YLE Switching 2 WF:*)


Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"]


YLEswitcher2[t0_, rin_, m_, dim_, wp_, prec_, \[Rho]low_, Ms_]:=Module[{d,as, Vars, Varsoft,Vars\[Phi],Varsoft\[Phi],ReducedVars,ReducedVarsoft,SetVariables,Mappingoft,SetVariables\[Phi],ReducedVars\[Phi],ReducedVarsoft\[Phi], BSOE,YLESOE,YLEphi, Broken, IV,Checker,IVinput,tm,ts, t0p},
Off[NIntegrate::precw]; Off[NDSolve::precw];
Vars=VarList;
Varsoft=Table[VarList[[i]][t],{i,1,Length[VarList]}] ;
SetVariables=Table[Vars[[i]]-> Varsoft[[i]],{i,1,Length[Vars]}];

Vars\[Phi]=VarList\[Phi];
Varsoft\[Phi]=Table[Vars\[Phi][[i]][t],{i,1,Length[Vars\[Phi]]}] ;
SetVariables\[Phi]=Table[Vars\[Phi][[i]]-> Varsoft\[Phi][[i]],{i,1,Length[Vars\[Phi]]}];

d=dim;

ReducedVars=Drop[Vars,{2}];(* Dropping a1 *)
ReducedVarsoft=Drop[Varsoft,{2}];

ReducedVars\[Phi]=Drop[Vars\[Phi],{3}]; (* Dropping a2 *)
ReducedVarsoft\[Phi]=Drop[Varsoft\[Phi],{3}];

BSOE=Join[Table[D[ ReducedVarsoft[[i]],t]==(\[Beta]FYLE2[i,d, Vars, \[Eta]Fys2[dim,Vars,Ms,as],Ms,as] /. SetVariables),{i,1,Length[ReducedVarsoft]}],
{Z'[t]==-(\[Eta]Fys2[d,Vars,Ms,as]/. SetVariables) Z[t]}] /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]};

(*YLESOE=Join[Table[D[ ReducedVarsoft\[Phi][[i]],t]==(\[Beta]FYLE\[Phi]2[i,d, Vars\[Phi], \[Eta]Fys2\[Phi][dim,Vars\[Phi], Ms, as],Ms,as]/. SetVariables\[Phi]),{i,1,Length[ReducedVarsoft\[Phi]]}],
{Z'[t]==-(\[Eta]Fys2\[Phi][d,Vars\[Phi], Ms, as]/.SetVariables\[Phi])Z[t]}] /. {Subscript[u, 2][t] -> m Exp[-2t] };*)

as=SetPrecision[apmsINTdata[[Round[Ms],2]],150];
IV=Join[InitialConditionsFunc[rin,m, Ms, 100, 100],{Z[tm]==1}];

Clear[tm];
tm=0; Checker=1;

While[ tm>t0,

If[ Checker == 0, 
ts=tm;
Mappingoft=(Join[{x,D[u[1/2 x^2],{x,1}]},Table[D[u[1/2 x^2],{x,i}],{i,3,2uor}],If[LPA==1,{},Table[D[y[1/2 x^2],{x,i}],{i,0,2yor}]],
Table[D[zs[1/2 x^2],{x,i}],{i,1,2zsor}]] /. x->Sqrt[2x] /. Truncation /. SetField/. SetVariables) /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]}; (* \[Phi] -> \[Rho] mapping *)
IVinput=Join[Table[ReducedVarsoft[[i]]->IV[[i]],{i,1,Length[IV]-1}],{Z[t]->IV[[-1]]}];

IV=Join[Table[ReducedVarsoft\[Phi][[i]]==(Mappingoft/.IVinput)[[i]] /. t->tm,{i,1,Length[ReducedVarsoft\[Phi]]}],{Z[tm]==IV[[-1]]}];
Print[IV];

YLESOE=Join[Table[D[ ReducedVarsoft\[Phi][[i]],t]==(\[Beta]FYLE\[Phi]2[i,d, Vars\[Phi], \[Eta]Fys2\[Phi][dim,Vars\[Phi], Ms, as],Ms,as]/. SetVariables\[Phi]),{i,1,Length[ReducedVarsoft\[Phi]]}],
{Z'[t]==-(\[Eta]Fys2\[Phi][d,Vars\[Phi], Ms, as]/.SetVariables\[Phi])Z[t]}] /. {Subscript[u, 2][t] -> m Exp[-2t] };
Print["YLESOE Compiled"];

YLEphi=NDSolve[{ YLESOE, IV },
Join[ReducedVars\[Phi],{Z}], {t, tm, t0} , WorkingPrecision->wp , PrecisionGoal->prec, AccuracyGoal-> prec,  
"Method"->{"StiffnessSwitching" ,"EquationSimplification"->{(*"Solve"*)Automatic, "TimeConstraint"->\[Infinity]}}  ,
 StepMonitor :>If[v==1, Print["\[Phi][t]=",\[Phi][t],", t=",t, ", \[Eta]=", \[Eta]Fys2\[Phi][dim,Varsoft\[Phi], Ms, as]/.{Subscript[u, 2][t] -> m Exp[-2t] }]],
 MaxSteps-> 5 10^3];
If[Checker==0, tm=t0];
,
Broken=NDSolve[ { BSOE, IV,WhenEvent[If[LPA==1, \[Rho][t] < \[Rho]low, - (Z'[t]/Z[t]) <\[Rho]low  (*&&   \[Rho][t] < \[Rho]low*)] ,  Checker=0;tm=t;IV=Join[ ReducedVarsoft /. t->tm,{Z[tm]}];Print["\[Eta] < ",\[Rho]low , ", -\!\(\*FractionBox[\(\(Z'\)[t]\), \(Z[t]\)]\)=", -(Z'[t]/Z[t])]; "StopIntegration"]}, 
Join[ReducedVars,{Z}], {t, tm, t0} ,WorkingPrecision->wp,PrecisionGoal-> prec, AccuracyGoal->prec, 
"Method"->{"StiffnessSwitching" ,"EquationSimplification"->{(*"Solve"*)Automatic, "TimeConstraint"->\[Infinity]}}  , 
StepMonitor :>If[v==1,Print["\[Rho][t]=",\[Rho][t],", t=",t, ", \[Eta]=", \[Eta]Fys2[dim,Varsoft, Ms, as]/.{Subscript[u, 1][t] ->m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]}]],
MaxSteps->10^4];
Print[IV,tm];
If[Checker==1, tm=t0];  
  ]];
t0p=InterpolatingFunctionDomain[Last[\[Phi] /. YLEphi]][[1,1]];
(*{{ts,t0p},YLEphi}  *)
{{ts,t0p},Broken,YLEphi}
];


(* ::Subsubsection::Closed:: *)
(*YLE Switching 2 WF Real:*)


YLEswitcher2Real[t0_, rin_, m_, dim_, wp_, prec_, \[Rho]low_, Ms_]:=Module[{d,as, Vars, Varsoft,Vars\[Phi],Varsoft\[Phi],ReducedVars,ReducedVarsoft,
SetVariables,Mappingoft,SetVariables\[Phi],ReducedVars\[Phi],ReducedVarsoft\[Phi], BSOE,YLESOE,YLEphi, Broken, IV,Checker,IVinput,tm,ts, t0p},
Off[NIntegrate::precw]; Off[NDSolve::precw];
Vars=VarList;
Varsoft=Table[VarList[[i]][t],{i,1,Length[VarList]}] ;
SetVariables=Table[Vars[[i]]-> Varsoft[[i]],{i,1,Length[Vars]}];

Vars\[Phi]=VarList\[Phi];
Varsoft\[Phi]=Table[Vars\[Phi][[i]][t],{i,1,Length[Vars\[Phi]]}] ;
SetVariables\[Phi]=Table[Vars\[Phi][[i]]-> Varsoft\[Phi][[i]],{i,1,Length[Vars\[Phi]]}];

d=dim;

ReducedVars=Drop[Vars,{2}];(* Dropping a1 *)
ReducedVarsoft=Drop[Varsoft,{2}];

ReducedVars\[Phi]=Drop[Vars\[Phi],{3}]; (* Dropping a2 *)
ReducedVarsoft\[Phi]=Drop[Varsoft\[Phi],{3}];

Print["Loading BSOE"];

BSOE=Join[Table[D[ ReducedVarsoft[[i]],t]==(\[Beta]FYLE2[i,d, Vars, \[Eta]Fys2[dim,Vars,Ms,as],Ms,as] /. SetVariables),{i,1,Length[ReducedVarsoft]}],
{Z'[t]==-(\[Eta]Fys2[d,Vars,Ms,as]/. SetVariables) Z[t]}] /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]};

Print["Loaded BSOE. Computing IC"];

(*YLESOE=Join[Table[D[ ReducedVarsoft\[Phi][[i]],t]==(\[Beta]FYLE\[Phi]2[i,d, Vars\[Phi], \[Eta]Fys2\[Phi][dim,Vars\[Phi], Ms, as],Ms,as]/. SetVariables\[Phi]),{i,1,Length[ReducedVarsoft\[Phi]]}],
{Z'[t]==-(\[Eta]Fys2\[Phi][d,Vars\[Phi], Ms, as]/.SetVariables\[Phi])Z[t]}] /. {Subscript[u, 2][t] -> m Exp[-2t] };*)

as=SetPrecision[apmsINTdata[[Round[Ms],2]],150];
IV=Join[InitialConditionsFunc[rin,m, Ms, 100, 100],{Z[tm]==1}];

Print["Begining NDSOLVE"];

Clear[tm];
tm=0; Checker=1;

While[ tm>t0,

If[ Checker == 0, 
ts=tm;
Mappingoft=(Join[{x,D[u[1/2 x^2],{x,1}]},Table[D[u[1/2 x^2],{x,i}],{i,3,2uor}],If[LPA==1,{},Table[D[y[1/2 x^2],{x,i}],{i,0,2yor}]],
Table[D[zs[1/2 x^2],{x,i}],{i,1,2zsor}]] /. x->Sqrt[2x] /. Truncation /. SetField/. SetVariables) /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]}; (* \[Phi] -> \[Rho] mapping *)
IVinput=Join[Table[ReducedVarsoft[[i]]->IV[[i]],{i,1,Length[IV]-1}],{Z[t]->IV[[-1]]}];

IV=Join[Table[ReducedVarsoft\[Phi][[i]]==If[((Mappingoft/.IVinput)[[i]] /. t->tm) \[Element] Reals, Re[(Mappingoft/.IVinput)[[i]] /. t->tm], Im[(Mappingoft/.IVinput)[[i]] /. t->tm]]
,{i,1,Length[ReducedVarsoft\[Phi]]}],{Z[tm]==IV[[-1]]}] /. t->tm;
Print[IV];

YLESOE=Join[Table[D[ ReducedVarsoft\[Phi][[i]],t]==(\[Beta]FYLE\[Phi]2C[i,d, Vars\[Phi], \[Eta]Fys2\[Phi]C[dim,Vars\[Phi], Ms, as],Ms,as]/. SetVariables\[Phi])(*//.eRepl*),{i,1,Length[ReducedVarsoft\[Phi]]}],
{Z'[t]==-(\[Eta]Fys2\[Phi]C[d,Vars\[Phi], Ms, as]/.SetVariables\[Phi])Z[t]}] /. {Subscript[u, 2][t] -> m Exp[-2t] };

Print["next"];

YLEphi=NDSolve[{ YLESOE (*//.eRepl*), IV },
Join[ReducedVars\[Phi],{Z}], {t, tm, t0} , WorkingPrecision->wp , PrecisionGoal->prec, AccuracyGoal-> prec,  
"Method"->{"StiffnessSwitching" ,"EquationSimplification"->{(*"Solve"*)Automatic, "TimeConstraint"->\[Infinity]}}  ,
 StepMonitor :>If[v==1, Print["\[Phi][t]=",\[Phi][t],", t=",t, ", \[Eta]=", \[Eta]Fys2\[Phi]C[dim,Varsoft\[Phi], Ms, as]/.{Subscript[u, 2][t] -> m Exp[-2t] }]],
 MaxSteps->  5 10^5];
If[Checker==0, tm=t0];
,
Broken=NDSolve[ { BSOE, IV,WhenEvent[If[LPA==1, \[Rho][t] < \[Rho]low, - (Z'[t]/Z[t]) <\[Rho]low  (*&&   \[Rho][t] < \[Rho]low*)] ,  Checker=0;tm=t;IV=Join[ ReducedVarsoft /. t->tm,{Z[tm]}];Print["\[Eta] < ",\[Rho]low , ", -\!\(\*FractionBox[\(\(Z'\)[t]\), \(Z[t]\)]\)=", -(Z'[t]/Z[t])]; "StopIntegration"]}, 
Join[ReducedVars,{Z}], {t, tm, t0} ,WorkingPrecision->wp,PrecisionGoal-> prec, AccuracyGoal->prec, 
"Method"->{"StiffnessSwitching" ,"EquationSimplification"->{(*"Solve"*)Automatic, "TimeConstraint"->\[Infinity]}}  , 
StepMonitor :>If[v==1,Print["\[Rho][t]=",\[Rho][t],", t=",t, ", \[Eta]=", \[Eta]Fys2[dim,Varsoft, Ms, as]/.{Subscript[u, 1][t] ->m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]}]],
MaxSteps->10^4];
Print[IV,tm];
If[Checker==1, tm=t0];  
  ]];
t0p=InterpolatingFunctionDomain[Last[\[Phi] /. YLEphi]][[1,1]];
(*{{ts,t0p},YLEphi}  *)
{{ts,t0p},Broken,YLEphi}
];


(* ::Subsubsection::Closed:: *)
(*YLE Switching 2 WF Real partial renormed variables:*)


(*eRepl=Join[Table[e^(2n)->I^(2n),{n,1,20}],Table[e^(-2n)->I^(2n),{n,1,20}],
Table[e^(2n+1)-> I^(2n+1) /.{I->e,-I->-e},{n,1,20}],
Table[e^(-2n+1)-> I^(-2n+1) /.{I->e,-I->-e},{n,1,20}]];*)


Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"]


YLEswitcher2RealR[t0_, rin_, m_, dim_, wp_, prec_, \[Rho]low_, Ms_]:=Module[{d,as, Vars, Varsoft,Vars\[Phi],Varsoft\[Phi],ReducedVars,ReducedVarsoft,
SetVariables,Mappingoft,SetVariables\[Phi],ReducedVars\[Phi],ReducedVarsoft\[Phi], BSOE,YLESOE,YLEphi, Broken, IV,Checker,IVinput,tm,ts, t0p},
Off[NIntegrate::precw]; Off[NDSolve::precw];
Vars=VarList;
Varsoft=Table[VarList[[i]][t],{i,1,Length[VarList]}] ;
SetVariables=Table[Vars[[i]]-> Varsoft[[i]],{i,1,Length[Vars]}];

Vars\[Phi]=VarList\[Phi];
Varsoft\[Phi]=Table[Vars\[Phi][[i]][t],{i,1,Length[Vars\[Phi]]}] ;
SetVariables\[Phi]=Table[Vars\[Phi][[i]]-> Varsoft\[Phi][[i]],{i,1,Length[Vars\[Phi]]}];

d=dim;

ReducedVars=Drop[Vars,{2}];(* Dropping a1 *)
ReducedVarsoft=Drop[Varsoft,{2}];

ReducedVars\[Phi]=Drop[Vars\[Phi],{3}]; (* Dropping a2 *)
ReducedVarsoft\[Phi]=Drop[Varsoft\[Phi],{3}];

BSOE=Join[Table[D[ ReducedVarsoft[[i]],t]==(\[Beta]FYLE2[i,d, Vars, \[Eta]Fys2[dim,Vars,Ms,as],Ms,as] /. SetVariables),{i,1,Length[ReducedVarsoft]}],
{Z'[t]==-(\[Eta]Fys2[d,Vars,Ms,as]/. SetVariables) Z[t]}] /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]};

(*YLESOE=Join[Table[D[ ReducedVarsoft\[Phi][[i]],t]==(\[Beta]FYLE\[Phi]2[i,d, Vars\[Phi], \[Eta]Fys2\[Phi][dim,Vars\[Phi], Ms, as],Ms,as]/. SetVariables\[Phi]),{i,1,Length[ReducedVarsoft\[Phi]]}],
{Z'[t]==-(\[Eta]Fys2\[Phi][d,Vars\[Phi], Ms, as]/.SetVariables\[Phi])Z[t]}] /. {Subscript[u, 2][t] -> m Exp[-2t] };*)

as=SetPrecision[apmsINTdata[[Round[Ms],2]],150];
IV=Join[InitialConditionsFunc[rin,m, Ms, 100, 100],{Z[tm]==1}];

Clear[tm];
tm=0; Checker=1;

While[ tm>t0,

If[ Checker == 0, 
ts=tm;
Mappingoft=(Join[{x,D[u[1/2 x^2],{x,1}]},Table[D[u[1/2 x^2],{x,i}],{i,3,2uor}],If[LPA==1,{},Table[D[y[1/2 x^2],{x,i}],{i,0,2yor}]],
Table[D[zs[1/2 x^2],{x,i}],{i,1,2zsor}]] /. x->Sqrt[2x] /. Truncation /. SetField/. SetVariables) /. {Subscript[u, 1][t] -> m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]}; (* \[Phi] -> \[Rho] mapping *)
IVinput=Join[Table[ReducedVarsoft[[i]]->IV[[i]],{i,1,Length[IV]-1}],{Z[t]->IV[[-1]]}];


Mappingoft[[1]]=Mappingoft[[1]] (Exp[tm((dim-2)/2)] Z[t]^(-1/2))^1  ;
Mappingoft[[2]]=Mappingoft[[2]] (Exp[tm((dim+2)/2)] Z[t]^(1/2))^1 ;
For[i=0, i<=2yor ,i++, Mappingoft[[2uor+1+i]]=Mappingoft[[2uor+1+i]]   Z[t]^(1(2+i/2))];



IV=Join[Table[ReducedVarsoft\[Phi][[i]]==If[((Mappingoft/.IVinput)[[i]] /. t->tm) \[Element] Reals, Re[(Mappingoft/.IVinput)[[i]] /. t->tm], Im[(Mappingoft/.IVinput)[[i]] /. t->tm]]
,{i,1,Length[ReducedVarsoft\[Phi]]}],{Z[tm]==IV[[-1]]}] /. t->tm;
Print[IV];

YLESOE=Join[Table[D[ ReducedVarsoft\[Phi][[i]],t]==(\[Beta]FYLE\[Phi]2CB[i,d, Vars\[Phi], \[Eta]Fys2\[Phi]CB[dim,Vars\[Phi], Ms, as],Ms,as]/. SetVariables\[Phi])(*//.eRepl*),{i,1,Length[ReducedVarsoft\[Phi]]}],
{Z'[t]==-(\[Eta]Fys2\[Phi]CB[d,Vars\[Phi], Ms, as]/.SetVariables\[Phi])Z[t]}] /. {Subscript[u, 2][t] -> m Exp[-2t] };

Print["next"];

YLEphi=NDSolve[{ YLESOE (*//.eRepl*), IV },
Join[ReducedVars\[Phi],{Z}], {t, tm, t0} , WorkingPrecision->wp , PrecisionGoal->prec, AccuracyGoal-> prec,  
"Method"->{"StiffnessSwitching" ,"EquationSimplification"->{(*"Solve"*)Automatic, "TimeConstraint"->\[Infinity]}}  ,
 StepMonitor :>If[v==1, Print["\[Phi][t]=",\[Phi][t],", t=",t, ", \[Eta]=", \[Eta]Fys2\[Phi]CB[dim,Varsoft\[Phi], Ms, as]/.{Subscript[u, 2][t] -> m Exp[-2t] }]],
 MaxSteps-> 5 10^5];
If[Checker==0, tm=t0];
,
Broken=NDSolve[ { BSOE, IV,WhenEvent[If[LPA==1, \[Rho][t] < \[Rho]low, - (Z'[t]/Z[t]) <\[Rho]low  (*&&   \[Rho][t] < \[Rho]low*)] ,  Checker=0;tm=t;IV=Join[ ReducedVarsoft /. t->tm,{Z[tm]}];Print["\[Eta] < ",\[Rho]low , ", -\!\(\*FractionBox[\(\(Z'\)[t]\), \(Z[t]\)]\)=", -(Z'[t]/Z[t])]; "StopIntegration"]}, 
Join[ReducedVars,{Z}], {t, tm, t0} ,WorkingPrecision->wp,PrecisionGoal-> prec, AccuracyGoal->prec, 
"Method"->{"StiffnessSwitching" ,"EquationSimplification"->{(*"Solve"*)Automatic, "TimeConstraint"->\[Infinity]}}  , 
StepMonitor :>If[v==1,Print["\[Rho][t]=",\[Rho][t],", t=",t, ", \[Eta]=", \[Eta]Fys2[dim,Varsoft, Ms, as]/.{Subscript[u, 1][t] ->m Exp[-2t] -2 \[Rho][t] Subscript[u, 2][t]}]],
MaxSteps->10^5];
Print[IV,tm];
If[Checker==1, tm=t0];  
  ]];
t0p=InterpolatingFunctionDomain[Last[\[Phi] /. YLEphi]][[1,1]];
(*{{ts,t0p},YLEphi}  *)
{{ts,t0p},Broken,YLEphi}
];


(* ::Subsection:: *)
(*Finding YLE:*)


(* ::Subsubsection:: *)
(*The edge:*)


(*wp=60; prec=wp/2;*)
wp=40; prec=wp/2;
(*wp=20; prec=wp/2;*)
$MinPrecision=prec;


Print["Starting Computation"]


(*(* For Standard Variables *)
t0=-40; 
dt0=10^-10; m0=10^-40;
v=1;
switch=(*-10^-1*)-5 10^-2;
{time,{{ts,t0p},Broken,YLE}}=YLEswitcher2Real[t0, dt0, m0,dim,wp,prec,switch, Ns]//AbsoluteTiming*)


(*(* For Standard Variables *)
Hc=Re[Subscript[u, 1][t0p]Sqrt[Exp[t0p(dim+2)]Z[t0p]]] /. YLE[[1]]
Mc=Re[\[Phi][t0p]Sqrt[(Exp[t0p]^(dim-2))/Z[t0p]]] /. YLE[[1]]
\[Eta]fy=-Z'[t0p]/Z[t0p] /. YLE[[1]] //Re
a3c=Re[Subscript[u, 3][t0p] /. YLE[[1]]]*)


(* For Modified Variables *)
t0=-40; 
dt0=10^-10; m0=10^-40;
v=1;
switch=(*-10^-1*)-5 10^-2;
{time,{{ts,t0p},Broken,YLE}}=YLEswitcher2RealR[t0, dt0, m0,dim,wp,prec,switch, Ns]//AbsoluteTiming


(* For Modified Variables *)
Hc=Subscript[u, 1][t0p] /. YLE[[1]]
Mc=\[Phi][t0p] /. YLE[[1]]
\[Eta]fy=-Z'[t0p]/Z[t0p] /. YLE[[1]] //Re
a3c=Re[Subscript[u, 3][t0p] /. YLE[[1]]]


Print[" YLE time was ", time/60. ," minutes or ",time/60./60, " hours" ]


Print[Hc]
Print[Mc]
Print[\[Eta]fy]
Print[a3c]


(*Plot[-Z'[-t]/Z[-t] /. YLE ,{t,-ts,-t0p}]*)


(*Show[Plot[-Z'[-t]/Z[-t] /. Broken,{t,0,-ts}, PlotRange->All],Plot[-Z'[-t]/Z[-t] /. YLE ,{t,-ts,-t0p},PlotRange->All]]*)


(*Plot[Re[Subscript[u, 1][-t]Sqrt[Exp[-t(dim+2)]Z[-t]]] /. YLE ,{t,-ts,-t0p}]
Plot[Re[\[Phi][-t]Sqrt[(Exp[-t]^(dim-2))/Z[-t]]] /. YLE ,{t,-ts,-t0p}]*)


(* ::Subsubsection:: *)
(*Saving:*)


FileNameOutputYLEWF //Print


Save[FileNameJoin[{Directory[], "YLE_HcR40/", FileNameOutputYLEWF}], {Ns, dim, uor, yor, zsor}]
Save[FileNameJoin[{Directory[], "YLE_HcR40/", FileNameOutputYLEWF}], {Hc, Mc, \[Eta]fy, a3c}]
Save[FileNameJoin[{Directory[], "YLE_HcR40/", FileNameOutputYLEWF}], {t0,dt0, m0, switch}]


\[Phi]range={ts,t0p};
Save[FileNameJoin[{Directory[], "YLE_HcR40/", FileNameOutputYLEWF}], time]
Save[FileNameJoin[{Directory[], "YLE_HcR40/", FileNameOutputYLEWF}], {ts,t0p}]
Save[FileNameJoin[{Directory[], "YLE_HcR40/", FileNameOutputYLEWF}], {Broken,YLE}]


(*FlowColl={{ts,t0p},Broken,YLE};
Save[FileNameJoin[{Directory[], "YLE_deltaR40/", FileNameOutputYLEWF}], FlowColl]*)


(*FRGYLEParameters=*){t0,dt0, m0, switch} //Print
{Hc, Mc, \[Eta]fy, a3c}//Print
{ts,t0p} //Print
(*{time,{{ts,t0p},Broken,YLE}}*)


(*hc=Hc H0;
M0=Mc/(Cp^(-1/\[Gamma]) dt0^(\[Gamma]/(\[Delta]-1)));
zcoR\[Chi]=dt0 Hc^((1-\[Delta])/(\[Gamma] \[Delta])) (Bc/Cp)^(1/\[Gamma])
Fg=Abs[Mc hc^(-1/\[Delta])]
ArgFg=Arg[I Mc (I hc)^(-1/\[Delta])]*)



(*zctilde=dt0 Hc^((1-\[Delta])/(\[Gamma] \[Delta])) (Bc/Cp)^(1/\[Gamma]) ((1+N\[Phi] \[Delta]-\[Delta])/\[Delta])^(1/\[Gamma]);*)
