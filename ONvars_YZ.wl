(* ::Package:: *)

(* ::Subsection:: *)
(*O(N) via (U,Z,Y):*)


(* ::Subsubsection::Closed:: *)
(*Defining Truncation and Variables in \[Rho]:*)


SetField={x-> \[Rho]};
SetDim[dim_]:={d-> dim, vd-> (*(2^(dim+1) \[Pi]^(dim/2) Gamma[dim/2])^-1*)1};

VarListOption[n_,m_,ms_]:=If[m==-1 &&ms==-1 &&LPA==1,
Join[ { \[Rho] }, Table[Subscript[u, i],{i,1,n}] ],
Join[ { \[Rho] }, Table[Subscript[u, i],{i,1,n}],Table[Subscript[y, i],{i,0,m}],Table[Subscript[zs, i],{i,1,ms}] ]];

TruncateOption[n_,m_,ms_]:=If[ ms==-1 && m==-1 &&LPA==1,
Flatten[{Table[D[u[x],{x,i}]-> Subscript[u, i],{i,1,n}], Table[D[u[x],{x,i}]-> 0,{i,n+1,n+8}], 
{zs[x]->1}, Table[D[zs[x],{x,i}]-> 0,{i,1,n+8}],
{y[x]->0}, Table[D[y[x],{x,i}]-> 0,{i,1,n+8}]}],
Flatten[{Table[D[u[x],{x,i}]-> Subscript[u, i],{i,1,n}], Table[D[u[x],{x,i}]-> 0,{i,n+1,n+8}], 
{y[x]->Subscript[y, 0]}, Table[D[y[x],{x,i}]-> Subscript[y, i],{i,1,m}],Table[D[y[x],{x,i}]-> 0,{i,m+1,n+8}],
{zs[x]->1}, Table[D[zs[x],{x,i}]-> Subscript[zs, i],{i,1,ms}],Table[D[zs[x],{x,i}]-> 0,{i,ms+1,n+8}]}]];

VarList=VarListOption[uor,yor,zsor];
Truncation=TruncateOption[uor,yor,zsor];


(* ::Subsubsection::Closed:: *)
(*Defining Truncation and Variables in \[Phi]:*)


(* ::Input:: *)
(*(*\[Phi] variables*)*)


SetField\[Phi]={x-> \[Phi]};
SetDim[dim_]:={d-> dim, vd-> (*(2^(dim+1) \[Pi]^(dim/2) Gamma[dim/2])^-1*)1};

VarListOption\[Phi][n_,m_,ms_]:=If[m==-1 &&ms==-1 &&LPA==1,
Join[ {\[Phi] }, Table[Subscript[u, i],{i,1,n}] ],
Join[ { \[Phi] }, Table[Subscript[u, i],{i,1,n}],Table[Subscript[y, i],{i,0,m}],Table[Subscript[zs, i],{i,1,ms}] ]];

TruncateOption\[Phi][n_,m_,ms_]:=If[ m<0 && ms<0  && LPA==1,
(Flatten[{Table[D[u[x],{x,i}]-> Subscript[u, i],{i,1,n}], Table[D[u[x],{x,i}]-> 0,{i,n+1,n+8}], 
{y[x]->0}, Table[D[y[x],{x,i}]-> 0,{i,1,n+8}],
{zs[x]->1}, Table[D[zs[x],{x,i}]-> 0,{i,1,n+8}]}]),
Flatten[{Table[D[u[x],{x,i}]-> Subscript[u, i],{i,1,n}], Table[D[u[x],{x,i}]-> 0,{i,n+1,n+8}], 
{y[x]->Subscript[y, 0]}, Table[D[y[x],{x,i}]-> Subscript[y, i],{i,1,m}],Table[D[y[x],{x,i}]-> 0,{i,m+1,n+8}],
{zs[x]->1}, Table[D[zs[x],{x,i}]-> Subscript[zs, i],{i,1,ms}],Table[D[zs[x],{x,i}]-> 0,{i,ms+1,n+8}]}]];

VarList\[Phi]=VarListOption\[Phi][2uor,2yor,2zsor];
Truncation\[Phi]=TruncateOption\[Phi][2uor,2yor,2zsor];
