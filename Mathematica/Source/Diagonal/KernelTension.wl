(* ::Package:: *)

BeginPackage["KernelTension`"]

PositionCorrelationTension::usage = "PositionCorrelationTension[\[CapitalDelta]s_,k_] gives the spectral kernel of order 0 and frequency k in response to tension modulations."
PositionCorrelationTensionC::usage = "PositionCorrelationTensionC[\[CapitalDelta]s_,k_] gives the spectral kernel of order 0 and frequency k in response to tension modulations. Compiled function."
TangentCorrelationTension::usage = "TangentCorrelationTension[\[CapitalDelta]s_,k_] gives the spectral kernel of order 1 and frequency k in response to tension modulations."
TangentCorrelationTensionC::usage = "TangentCorrelationTensionC[\[CapitalDelta]s_,k_] gives the spectral kernel of order 1 and frequency k in response to tension modulations. Compiled function."
CurvatureCorrelationTension::usage = "CurvatureCorrelationTension[\[CapitalDelta]s_,k_] gives the spectral kernel of order 2 and frequency k in response to tension modulations."
CurvatureCorrelationTensionC::usage = "CurvatureCorrelationTensionC[\[CapitalDelta]s_,k_] gives the spectral kernel of order 2 and frequency k in response to tension modulations. Compiled function."

Begin["`Private`"]


(* ::Text:: *)
(*Import generic kernels that we will here specialize.*)


Needs["KernelGeneric`"->"K`","KernelGeneric.wl"];


(* ::Text:: *)
(*Now we specialize the kernels to make them useful for further calculations. Higher order kernels such as the tangent kernel, the curvature kernel, and all tension kernels, do not yield useful results for generic (\[Lambda],\[Kappa]), so we neglect them.*)


cachedir=NotebookDirectory[]<>"Cache";
If[!FileExistsQ[cachedir],CreateDirectory[cachedir]];


compile[expr_]:=Compile[{{\[CapitalDelta]s,_Real},{k,_Real},{lambda,_Real},{kappa,_Real}}, 
	expr, CompilationTarget->"C", RuntimeOptions->{"Speed","EvaluateSymbolically" -> False}];


(* Specialize to n == 0. *)
Block[{case, n=0},
With[{file = cachedir<>"/KernelTension_0.wl"}, If[!FileExistsQ[file],
case[1] = K`CorrelationTension[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationTension[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Lambda]->lambda];
case[4] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Lambda]->lambda];
case[5] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Kappa]->kappa];
case[6] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Kappa]->kappa];

Save[file, case];,
Get[file];
]];

Options[PositionCorrelationTension]={\[Lambda]->0,\[Kappa]->0}; (
PositionCorrelationTension[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1,
	kappa===0 && lambda===0 && k=!=0, #2,
	kappa===0 && lambda=!=0 && k===0, #3,
	kappa===0 && lambda=!=0 && k=!=0, #4,
	kappa=!=0 && lambda===0 && k===0, #5,
	kappa=!=0 && lambda===0 && k=!=0, #6,
	kappa=!=0 && lambda=!=0, Throw["The case (\[Lambda]!=0, \[Kappa]!=0) is not defined."]
]];
)&@@Table[case[i],{i,6}];

Options[PositionCorrelationTensionC]={\[Lambda]->0,\[Kappa]->0}; (
PositionCorrelationTensionC[\[CapitalDelta]s:_?NumericQ,k:_?NumericQ:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda===0 && k=!=0, #2[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k===0, #3[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k=!=0, #4[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0 && lambda===0 && k===0, #5[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0 && lambda===0 && k=!=0, #6[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0 && lambda=!=0, Throw["The case (\[Lambda]!=0, \[Kappa]!=0) is not defined."]
]];
)&@@(compile/@Table[case[i],{i,6}]);
];


(* Specialize to n == 1. *)
Block[{case, n=1},
With[{file = cachedir<>"/KernelTension_1.wl"}, If[!FileExistsQ[file],
case[1] = K`CorrelationTension[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationTension[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Lambda]->lambda];
case[4] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Lambda]->lambda];

Save[file, case];,
Get[file];
]];

Options[TangentCorrelationTension]={\[Lambda]->0,\[Kappa]->0}; (
TangentCorrelationTension[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1,
	kappa===0 && lambda===0 && k=!=0, #2,
	kappa===0 && lambda=!=0 && k===0, #3,
	kappa===0 && lambda=!=0 && k=!=0, #4,
	kappa=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@Table[case[i],{i,4}];

Options[TangentCorrelationTensionC]={\[Lambda]->0,\[Kappa]->0}; (
TangentCorrelationTensionC[\[CapitalDelta]s:_?NumericQ,k:_?NumericQ:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda===0 && k=!=0, #2[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k===0, #3[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k=!=0, #4[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@(compile/@Table[case[i],{i,4}]);
];


(* Specialize to n == 2. *)
Block[{case, n=2},
With[{file = cachedir<>"/KernelTension_2.wl"}, If[!FileExistsQ[file],
case[1] = K`CorrelationTension[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationTension[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationTension[\[CapitalDelta]s,0,n,\[Lambda]->lambda];
case[4] = K`CorrelationTension[\[CapitalDelta]s,k,n,\[Lambda]->lambda];

Save[file, case];,
Get[file];
]];

Options[CurvatureCorrelationTension]={\[Lambda]->0,\[Kappa]->0}; (
CurvatureCorrelationTension[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1,
	kappa===0 && lambda===0 && k=!=0, #2,
	kappa===0 && lambda=!=0 && k===0, #3,
	kappa===0 && lambda=!=0 && k=!=0, #4,
	kappa=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@Table[case[i],{i,4}];

Options[CurvatureCorrelationTensionC]={\[Lambda]->0,\[Kappa]->0}; (
CurvatureCorrelationTensionC[\[CapitalDelta]s:_?NumericQ,k:_?NumericQ:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda===0 && k=!=0, #2[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k===0, #3[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k=!=0, #4[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@(compile/@Table[case[i],{i,4}]);
];


End[]
EndPackage[]
