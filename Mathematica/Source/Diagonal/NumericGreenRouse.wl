(* ::Package:: *)

BeginPackage["Green`"]

SquaredSeparationNum::usage = "SquaredSeparationNum[x,y] determines the squared separation between the positions of two points x and y.";
TangentCorrelationNum::usage = "TangentCorrelationNum[x,y] determines the correlation function between the tangent vectors at two points x and y.";

Begin["`Private`"]


(* ::Text:: *)
(*Define Greens function kernels for a Rouse polymer, which are obtained by analytic calculations.*)


Begin["`Kernel`"]

SeparationActivity[s1_,s2_] := 1/(2Pi) (2Log[s1^2+s2^2]-Log[2s1^2]-Log[2s2^2]);
SeparationTension[s1_,s2_] := (1/2) (RealSign[s1]RealSign[s2]-1);
TangentCorrelationActivity[s1_,s2_] := 2/Pi (s1 s2)/(s1^2+s2^2)^2;

End[]


(* ::Text:: *)
(*Define Functions that use Greens function kernels to map activity and tension profiles to mean squared separation and tangent autocorrelation functions.*)


Options[SquaredSeparationNum]={"ActivityProfile"->(1&),"TensionProfile"->(0&)};
SquaredSeparationNum[s1_,s2_,OptionsPattern[]] := With[{ActFun=OptionValue["ActivityProfile"], TenFun=OptionValue["TensionProfile"]},
NIntegrate[
	N[ActFun[x]Green`Private`Kernel`SeparationActivity[s1-x,s2-x] + TenFun[x]Green`Private`Kernel`SeparationTension[s1-x,s2-x],$MachinePrecision],
	{x,-Infinity,Infinity},
	Method->"GlobalAdaptive",
	WorkingPrecision->8,
	AccuracyGoal->8,
	MinRecursion->8,
	Exclusions->{s1,s2}
]];


Options[TangentCorrelationNum]={"ActivityProfile"->(1&)};
TangentCorrelationNum[s1_,s2_,OptionsPattern[]] := With[{ActFun=OptionValue["ActivityProfile"]},
NIntegrate[
	N[ActFun[x]Green`Private`Kernel`TangentCorrelationActivity[s1-x,s2-x],$MachinePrecision],
	{x,-Infinity,Infinity},
	Method->"GlobalAdaptive",
	WorkingPrecision->8,
	AccuracyGoal->8,
	MinRecursion->8,
	Exclusions->{s1,s2}
]];


(* ::Text:: *)
(*Disable NIntegrate warnings for parallel evaluations, to keep the output tidy.*)


Off[NIntegrate::ncvb]; ParallelEvaluate[Off[NIntegrate::ncvb]];
Off[Infinity::indet]; ParallelEvaluate[Off[Infinity::indet]];
Off[Power::infy]; ParallelEvaluate[Off[Power::infy]];


DistributeDefinitions["Green`"]


End[]
EndPackage[]
