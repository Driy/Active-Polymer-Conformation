(* ::Package:: *)

BeginPackage["KernelActivity`"]

PositionCorrelationActivity::usage = "PositionCorrelationActivity[\[CapitalDelta]s_,k_] gives the spectral kernel of order 0 and frequency k in response to activity modulations."
PositionCorrelationActivityC::usage = "PositionCorrelationActivityC[\[CapitalDelta]s_,k_] gives the spectral kernel of order 0 and frequency k in response to activity modulations. Compiled function."
TangentCorrelationActivity::usage = "TangentCorrelationActivity[\[CapitalDelta]s_,k_] gives the spectral kernel of order 1 and frequency k in response to activity modulations."
TangentCorrelationActivityC::usage = "TangentCorrelationActivityC[\[CapitalDelta]s_,k_] gives the spectral kernel of order 1 and frequency k in response to activity modulations. Compiled function."
CurvatureCorrelationActivity::usage = "CurvatureCorrelationActivity[\[CapitalDelta]s_,k_] gives the spectral kernel of order 2 and frequency k in response to activity modulations."
CurvatureCorrelationActivityC::usage = "CurvatureCorrelationActivityC[\[CapitalDelta]s_,k_] gives the spectral kernel of order 2 and frequency k in response to activity modulations. Compiled function."

Begin["`Private`"]


(* ::Text:: *)
(*Import generic kernels that we will here specialize.*)


Needs["KernelGeneric`"->"K`","KernelGeneric.wl"];


(* ::Text:: *)
(*Now we specialize the kernels to make them useful for further calculations. Higher order kernels such as the tangent kernel, the curvature kernel, and all tension kernels, do not yield useful results for generic (\[Lambda],\[Kappa]), so we neglect them.*)


cachedir=NotebookDirectory[]<>"Cache";
If[!FileExistsQ[cachedir],CreateDirectory[cachedir]];


compile[expr_]:=Compile[{{\[CapitalDelta]s,_Real},{k,_Integer},{lambda,_Real},{kappa,_Real}}, 
	expr, CompilationTarget->"C", RuntimeOptions->"Speed"];


(* Specialize to n == 0. *)
Block[{case, n=0},
With[{file = cachedir<>"/KernelActivity_0.wl"}, If[!FileExistsQ[file],
case[1] = K`CorrelationActivity[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationActivity[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationActivity[\[CapitalDelta]s,0,n,\[Lambda]->lambda];
case[4] = K`CorrelationActivity[\[CapitalDelta]s,k,n,\[Lambda]->lambda];
case[5] = K`CorrelationActivity[\[CapitalDelta]s,0,n,\[Kappa]->kappa];
case[6] = K`CorrelationActivity[\[CapitalDelta]s,k,n,\[Kappa]->kappa];
case[7] = K`CorrelationActivity[\[CapitalDelta]s,0,n,\[Lambda]->lambda,\[Kappa]->kappa];
case[8] = K`CorrelationActivity[\[CapitalDelta]s,k,n,\[Lambda]->lambda,\[Kappa]->kappa];

Save[file, case];,
Get[file];
]];

Options[PositionCorrelationActivity]={\[Lambda]->0,\[Kappa]->0}; (
PositionCorrelationActivity[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1,
	kappa===0 && lambda===0 && k=!=0, #2,
	kappa===0 && lambda=!=0 && k===0, #3,
	kappa===0 && lambda=!=0 && k=!=0, #4,
	kappa=!=0 && lambda===0 && k===0, #5,
	kappa=!=0 && lambda===0 && k=!=0, #6,
	kappa=!=0 && lambda=!=0 && k===0, #7,
	kappa=!=0 && lambda=!=0 && k=!=0, #8
]];
)&@@Table[case[i],{i,8}];

Options[PositionCorrelationActivityC]={\[Lambda]->0,\[Kappa]->0}; (
PositionCorrelationActivityC[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda===0 && k=!=0, #2[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k===0, #3[\[CapitalDelta]s,k,lambda,kappa],
	kappa===0 && lambda=!=0 && k=!=0, #4[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0 && lambda===0 && k===0, #5[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0 && lambda===0 && k=!=0, #6[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0 && lambda=!=0 && k===0, #7[\[CapitalDelta]s,k,lambda,kappa],
	kappa=!=0 && lambda=!=0 && k=!=0, #8[\[CapitalDelta]s,k,lambda,kappa]
]];
)&@@(compile/@Table[case[i],{i,8}]);
];


(* Specialize to n == 1. *)
Block[{case, n=1},
With[{file = cachedir<>"/KernelActivity_1.wl"}, If[!FileExistsQ[file],
case[1] = K`CorrelationActivity[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationActivity[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationActivity[\[CapitalDelta]s,0,n,\[Lambda]->lambda];
case[4] = K`CorrelationActivity[\[CapitalDelta]s,k,n,\[Lambda]->lambda];

Save[file, case];,
Get[file];
]];

Options[TangentCorrelationActivity]={\[Lambda]->0,\[Kappa]->0}; (
TangentCorrelationActivity[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1,
	kappa===0 && lambda===0 && k=!=0, #2,
	kappa===0 && lambda=!=0 && k===0, #3,
	kappa===0 && lambda=!=0 && k=!=0, #4,
	kappa=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@Table[case[i],{i,4}];

Options[TangentCorrelationActivityC]={\[Lambda]->0,\[Kappa]->0}; (
TangentCorrelationActivityC[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
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
With[{file = cachedir<>"/KernelActivity_2.wl"}, If[!FileExistsQ[file],
case[1] = K`CorrelationActivity[\[CapitalDelta]s,0,n];
case[2] = K`CorrelationActivity[\[CapitalDelta]s,k,n];
case[3] = K`CorrelationActivity[\[CapitalDelta]s,0,n,\[Lambda]->lambda];
case[4] = K`CorrelationActivity[\[CapitalDelta]s,k,n,\[Lambda]->lambda];

Save[file, case];,
Get[file];
]];

Options[CurvatureCorrelationActivity]={\[Lambda]->0,\[Kappa]->0}; (
CurvatureCorrelationActivity[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
Which[
	kappa===0 && lambda===0 && k===0, #1,
	kappa===0 && lambda===0 && k=!=0, #2,
	kappa===0 && lambda=!=0 && k===0, #3,
	kappa===0 && lambda=!=0 && k=!=0, #4,
	kappa=!=0, Throw["The case \[Kappa]!=0 is not defined."]
]];
)&@@Table[case[i],{i,4}];

Options[CurvatureCorrelationActivityC]={\[Lambda]->0,\[Kappa]->0}; (
CurvatureCorrelationActivityC[\[CapitalDelta]s_,k_:0,OptionsPattern[]] := With[{lambda=OptionValue[\[Lambda]],kappa=OptionValue[\[Kappa]]}, 
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
