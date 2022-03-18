(* ::Package:: *)

BeginPackage["Spectral`"]

PositionCorrelationNum::usage = "PositionCorrelationNum[x,y] determines the correlation function between the positions of two points x and y.";
SquaredSeparationNum::usage = "SquaredSeparationNum[x,y] determines the squared separation between the positions of two points x and y.";
TangentCorrelationNum::usage = "TangentCorrelationNum[x,y] determines the correlation function between the tangent vectors at two points x and y.";

Begin["`Private`"]


(* ::Text:: *)
(*Now we specialize the kernels to make them useful for further calculations. The tension kernel does not yield useful results for generic (\[Lambda],\[Kappa]), so we neglect it.*)


Begin["`Kernel`"]

(* Import generic kernels that we will here specialize. *)
Import["Kernel.wl"]

(* Specialize to n == 0. *)
Options[PositionCorrelationActivity]={\[Lambda]->0,\[Kappa]->0};
(PositionCorrelationActivity[\[CapitalDelta]s_,k_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
N[Which[
	\[Kappa]val===0 && \[Lambda]val===0, #1,
	\[Kappa]val===0, #2,
	True, #3
],$MachinePrecision]])&@@{
	Kernel`CorrelationActivity[\[CapitalDelta]s,k,0],
	Kernel`CorrelationActivity[\[CapitalDelta]s,k,0,\[Lambda]->\[Lambda]val],
	Kernel`CorrelationActivity[\[CapitalDelta]s,k,0,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val]
};

(* Specialize to n == 1. *)
Options[TangentCorrelationActivity]={\[Lambda]->0,\[Kappa]->0};
(TangentCorrelationActivity[\[CapitalDelta]s_,k_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]]}, 
N[Which[
	\[Kappa]val===0 && \[Lambda]val===0, #1,
	\[Kappa]val===0, #2,
	True, Throw["The case (\[Lambda]!=0, \[Kappa]!=0) is not defined."]
],$MachinePrecision]])&@@{
	Kernel`CorrelationActivity[\[CapitalDelta]s,k,1],
	Kernel`CorrelationActivity[\[CapitalDelta]s,k,1,\[Lambda]->\[Lambda]val]
};

End[]


(* ::Text:: *)
(*Import Fourier packages since we will use them for numerical Fourier transforms.*)


Needs["FourierSeries`"]; ParallelNeeds["FourierSeries`"];


(* ::Text:: *)
(*Define Functions that use Fourier kernels to map activity and tension profiles to correlation functions.*)


Options[PositionCorrelationNum]={\[Lambda]->0,\[Kappa]->0,"ActivitySpectrum"->None};
PositionCorrelationNum[x_,y_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]],ActFun=OptionValue["ActivitySpectrum"]},
Re@Which[
	ActFun===None, Spectral`Private`Kernel`PositionCorrelationActivity[Abs[x-y],0,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val],
	True, NInverseFourierTransform[
		Spectral`Private`Kernel`PositionCorrelationActivity[Abs[x-y],k,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val]ActFun[k],
		k,(x+y)/2,
		FourierParameters->{1,-1},
		Method->"AdaptiveMonteCarlo",
		MaxRecursion->30, PrecisionGoal->4
	]
]];

Options[SquaredSeparationNum]={\[Lambda]->0,\[Kappa]->0,"ActivitySpectrum"->None};
SquaredSeparationNum[x_,y_,opts:OptionsPattern[]] := With[{fun=Spectral`PositionCorrelationNum},
	fun[x,x,opts]+fun[y,y,opts]-2fun[x,y,opts]
];

Options[TangentCorrelationNum]={\[Lambda]->0,\[Kappa]->0,"ActivitySpectrum"->None};
TangentCorrelationNum[x_,y_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]],ActFun=OptionValue["ActivitySpectrum"]},
Re@Which[
	ActFun===None, Spectral`Private`Kernel`TangentCorrelationActivity[Abs[x-y],0,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val],
	True, NInverseFourierTransform[
		Spectral`Private`Kernel`TangentCorrelationActivity[Abs[x-y],k,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val]ActFun[k],
		k,(x+y)/2,
		FourierParameters->{1,-1},
		Method->"AdaptiveMonteCarlo",
		MaxRecursion->30, PrecisionGoal->4
	]
]];


(* ::Text:: *)
(*Disable NIntegrate warnings for parallel evaluations, to keep the output tidy.*)


Off[NIntegrate::ncvb]; ParallelEvaluate[Off[NIntegrate::ncvb]];
Off[Infinity::indet]; ParallelEvaluate[Off[Infinity::indet]];
Off[Power::infy]; ParallelEvaluate[Off[Power::infy]];


DistributeDefinitions["Spectral`"]


End[]
EndPackage[]
