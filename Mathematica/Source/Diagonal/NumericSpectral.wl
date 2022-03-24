(* ::Package:: *)

BeginPackage["NumericSpectral`"]

PositionCorrelationNum::usage = "PositionCorrelationNum[x,y] determines the correlation function between the positions of two points x and y.";
SquaredSeparationNum::usage = "SquaredSeparationNum[x,y] determines the squared separation between the positions of two points x and y.";
TangentCorrelationNum::usage = "TangentCorrelationNum[x,y] determines the correlation function between the tangent vectors at two points x and y.";

Begin["`Private`"]


(* ::Text:: *)
(*Import specialized kernels.*)


Needs["KernelActivity`"->"K`","KernelActivity.wl"];


(* ::Text:: *)
(*Import Fourier packages since we will use them for numerical Fourier transforms.*)


Needs["FourierSeries`"]; ParallelNeeds["FourierSeries`"];


(* ::Text:: *)
(*Define Functions that use Fourier kernels to map activity and tension profiles to correlation functions.*)


Options[PositionCorrelationNum]={\[Lambda]->0,\[Kappa]->0,"ActivitySpectrum"->None};
PositionCorrelationNum[x_,y_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]],ActFun=OptionValue["ActivitySpectrum"]},
Re@Which[
	ActFun===None, N[K`PositionCorrelationActivity[Abs[x-y],0,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val],$MachinePrecision],
	True, NInverseFourierTransform[
		N[K`PositionCorrelationActivity[Abs[x-y],k,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val]ActFun[k],$MachinePrecision],
		k,(x+y)/2,
		FourierParameters->{1,-1},
		Method->"AdaptiveMonteCarlo",
		MaxRecursion->30, PrecisionGoal->4
	]
]];

Options[SquaredSeparationNum]={\[Lambda]->0,\[Kappa]->0,"ActivitySpectrum"->None};
SquaredSeparationNum[x_,y_,opts:OptionsPattern[]] := With[{fun=NumericSpectral`PositionCorrelationNum},
	fun[x,x,opts]+fun[y,y,opts]-2fun[x,y,opts]
];

Options[TangentCorrelationNum]={\[Lambda]->0,\[Kappa]->0,"ActivitySpectrum"->None};
TangentCorrelationNum[x_,y_,OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]],ActFun=OptionValue["ActivitySpectrum"]},
Re@Which[
	ActFun===None, N[K`TangentCorrelationActivity[Abs[x-y],0,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val],$MachinePrecision],
	True, NInverseFourierTransform[
		N[K`TangentCorrelationActivity[Abs[x-y],k,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val]ActFun[k],$MachinePrecision],
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
