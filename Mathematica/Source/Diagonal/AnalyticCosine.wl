(* ::Package:: *)

BeginPackage["AnalyticCosine`"]

PositionCorrelation::usage = "PositionCorrelationNum[x,y] determines the correlation function between the positions of two points x and y.";
SquaredSeparation::usage = "SquaredSeparationNum[x,y] determines the squared separation between the positions of two points x and y.";
TangentCorrelation::usage = "TangentCorrelationNum[x,y] determines the correlation function between the tangent vectors at two points x and y.";
CurvatureCorrelation::usage = "CurvatureCorrelationNum[x,y] determines the correlation function between the curvature vectors at two points x and y."

Begin["`Private`"]


Needs["KernelActivity`"->None,"KernelActivity.wl"]
Needs["KernelTension`"->None,"KernelTension.wl"]


Options[Generator]={\[Lambda]->0,\[Kappa]->0,"ActivityMagnitude"->0,"TensionMagnitude"->0,"WaveMode"->1};
Generator[funActivity_,funTension_,x_,y_,opts:OptionsPattern[]] := With[{\[Lambda]val=OptionValue[\[Lambda]],\[Kappa]val=OptionValue[\[Kappa]],\[CapitalTau]=OptionValue["ActivityMagnitude"],\[CapitalSigma]=OptionValue["TensionMagnitude"],k=OptionValue["WaveMode"]},
	\[CapitalDelta]s=Abs[x-y]; S=x+y;
	funActivity[\[CapitalDelta]s,0,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val] + \[CapitalTau] Cos[k S/2] funActivity[\[CapitalDelta]s,k,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val] - \[CapitalSigma]/2 Cos[k S/2] funTension[\[CapitalDelta]s,k,\[Lambda]->\[Lambda]val,\[Kappa]->\[Kappa]val]
];


Options[PositionCorrelation]={\[Lambda]->0,\[Kappa]->0,"ActivityMagnitude"->0,"TensionMagnitude"->0,"WaveMode"->1};
PositionCorrelation[x_,y_,opts:OptionsPattern[]] := AnalyticCosine`Private`Generator[KernelActivity`PositionCorrelationActivity,KernelTension`PositionCorrelationTension, x,y,opts];


Options[TangentCorrelation]={\[Lambda]->0,\[Kappa]->0,"ActivityMagnitude"->0,"TensionMagnitude"->0,"WaveMode"->1};
TangentCorrelation[x_,y_,opts:OptionsPattern[]] := AnalyticCosine`Private`Generator[KernelActivity`TangentCorrelationActivity,KernelTension`TangentCorrelationTension, x,y,opts];


Options[CurvatureCorrelation]={\[Lambda]->0,\[Kappa]->0,"ActivityMagnitude"->0,"TensionMagnitude"->0,"WaveMode"->1};
CurvatureCorrelation[x_,y_,opts:OptionsPattern[]] := AnalyticCosine`Private`Generator[KernelActivity`CurvatureCorrelationActivity,KernelTension`CurvatureCorrelationTension, x,y,opts];


Options[SquaredSeparation]={\[Lambda]->0,\[Kappa]->0,"ActivityMagnitude"->0,"TensionMagnitude"->0,"WaveMode"->1};
SquaredSeparation[x_,y_,opts:OptionsPattern[]] := With[{fun=AnalyticCosine`PositionCorrelation},
	fun[x,x,opts]+fun[y,y,opts]-2fun[x,y,opts]
];


End[]
EndPackage[]
