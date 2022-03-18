(* ::Package:: *)

BeginPackage["KernelActivity`"]

PositionCorrelationActivity::usage = "PositionCorrelationActivity[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order 0 and frequency k in response to activity modulations."
TangentCorrelationActivity::usage = "TangentCorrelationActivity[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order 1 and frequency k in response to activity modulations."
CurvatureCorrelationActivity::usage = "CurvatureCorrelationActivity[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order 2 and frequency k in response to activity modulations."

Begin["`Private`"]


(* ::Text:: *)
(*Import generic kernels that we will here specialize.*)


Needs["KernelGeneric`"->"K`","KernelGeneric.wl"];


(* ::Text:: *)
(*Now we specialize the kernels to make them useful for further calculations. Higher order kernels such as the tangent kernel, the curvature kernel, and all tension kernels, do not yield useful results for generic (\[Lambda],\[Kappa]), so we neglect them.*)


(* Specialize to n == 0. *)
Options[PositionCorrelationActivity]={\[Lambda]->0,\[Kappa]->0};
(PositionCorrelationActivity[\[CapitalDelta]s_,k_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
N[Which[
	\[Kappa]val===0 && \[Lambda]val===0, #1,
	\[Kappa]val===0, #2,
	True, #3
],$MachinePrecision]])&@@{
	K`CorrelationActivity[\[CapitalDelta]s,k,0],
	K`CorrelationActivity[\[CapitalDelta]s,k,0,\[Lambda]->\[Lambda]val],
	K`CorrelationActivity[\[CapitalDelta]s,k,0,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val]
};


(* Specialize to n == 1. *)
Options[TangentCorrelationActivity]={\[Lambda]->0,\[Kappa]->0};
(TangentCorrelationActivity[\[CapitalDelta]s_,k_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
N[Which[
	\[Kappa]val===0 && \[Lambda]val===0, #1,
	\[Kappa]val===0, #2,
	True, Throw["The case (\[Lambda]!=0, \[Kappa]!=0) is not defined."]
],$MachinePrecision]])&@@{
	K`CorrelationActivity[\[CapitalDelta]s,k,1],
	K`CorrelationActivity[\[CapitalDelta]s,k,1,\[Lambda]->\[Lambda]val]
};


(* Specialize to n == 2. *)
Options[CurvatureCorrelationActivity]={\[Lambda]->0,\[Kappa]->0};
(CurvatureCorrelationActivity[\[CapitalDelta]s_,k_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
N[Which[
	\[Kappa]val===0 && \[Lambda]val===0, #1,
	\[Kappa]val===0, #2,
	True, Throw["The case (\[Lambda]!=0, \[Kappa]!=0) is not defined."]
],$MachinePrecision]])&@@{
	K`CorrelationActivity[\[CapitalDelta]s,k,2],
	K`CorrelationActivity[\[CapitalDelta]s,k,2,\[Lambda]->\[Lambda]val]
};


End[]
EndPackage[]
