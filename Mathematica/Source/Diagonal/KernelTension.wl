(* ::Package:: *)

BeginPackage["KernelTension`"]

PositionCorrelationTension::usage = "PositionCorrelationTension[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order 0 and frequency k in response to tension modulations."
TangentCorrelationTension::usage = "TangentCorrelationTension[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order 1 and frequency k in response to tension modulations."
CurvatureCorrelationTension::usage = "CurvatureCorrelationTension[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order 2 and frequency k in response to tension modulations."

Begin["`Private`"]


(* ::Text:: *)
(*Import generic kernels that we will here specialize.*)


Needs["KernelGeneric`"->"K`","KernelGeneric.wl"];


(* ::Text:: *)
(*Now we specialize the kernels to make them useful for further calculations. Higher order kernels such as the tangent kernel, the curvature kernel, and all tension kernels, do not yield useful results for generic (\[Lambda],\[Kappa]), so we neglect them.*)


(* Specialize to n == 0. *)
Block[{case, n=0},
case[1] = K`CorrelationTension[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationTension[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Lambda]->\[Lambda]val];
case[4] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Lambda]->\[Lambda]val];
case[5] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Kappa]->\[Kappa]val];
case[6] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Kappa]->\[Kappa]val];

Options[PositionCorrelationTension]={\[Lambda]->0,\[Kappa]->0}; (
PositionCorrelationTension[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
Which[
	\[Kappa]val===0 && \[Lambda]val===0 && k===0, #1,
	\[Kappa]val===0 && \[Lambda]val===0 && k=!=0, #2,
	\[Kappa]val===0 && \[Lambda]val=!=0 && k===0, #3,
	\[Kappa]val===0 && \[Lambda]val=!=0 && k=!=0, #4,
	\[Kappa]val=!=0 && \[Lambda]val===0 && k===0, #5,
	\[Kappa]val=!=0 && \[Lambda]val===0 && k=!=0, #6,
	\[Kappa]val=!=0 && \[Lambda]val=!=0, Throw["The case (\[Lambda]!=0, \[Kappa]!=0) is not defined."]
]];
)&@@Table[case[i],{i,6}];
];


(* Specialize to n == 1. *)
Block[{case, n=1},
case[1] = K`CorrelationTension[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationTension[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Lambda]->\[Lambda]val];
case[4] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Lambda]->\[Lambda]val];

Options[TangentCorrelationTension]={\[Lambda]->0,\[Kappa]->0}; (
TangentCorrelationTension[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
Which[
	\[Kappa]val===0 && \[Lambda]val===0 && k===0, #1,
	\[Kappa]val===0 && \[Lambda]val===0 && k=!=0, #2,
	\[Kappa]val===0 && \[Lambda]val=!=0 && k===0, #3,
	\[Kappa]val===0 && \[Lambda]val=!=0 && k=!=0, #4,
	\[Kappa]val=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@Table[case[i],{i,4}];
];


(* Specialize to n == 2. *)
Block[{case, n=2},
case[1] = K`CorrelationTension[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationTension[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Lambda]->\[Lambda]val];
case[4] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Lambda]->\[Lambda]val];

Options[CurvatureCorrelationTension]={\[Lambda]->0,\[Kappa]->0}; (
CurvatureCorrelationTension[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
Which[
	\[Kappa]val===0 && \[Lambda]val===0 && k===0, #1,
	\[Kappa]val===0 && \[Lambda]val===0 && k=!=0, #2,
	\[Kappa]val===0 && \[Lambda]val=!=0 && k===0, #3,
	\[Kappa]val===0 && \[Lambda]val=!=0 && k=!=0, #4,
	\[Kappa]val=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@Table[case[i],{i,4}];
];


End[]
EndPackage[]
