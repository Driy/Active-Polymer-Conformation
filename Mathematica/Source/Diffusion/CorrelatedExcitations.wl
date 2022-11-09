(* ::Package:: *)

(* ::Text:: *)
(*Methods for calculating the mean squared distance that a locus travels within some time window. Active processes are here encoded as a region of size \[Lambda] that experiences correlated excitations (correlation coefficient \[Rho]). We track a locus at the center of this region.*)
(**)
(*Time is measured in units of \[Xi]/\[Kappa]. Squared distance is measured in b^2=Subscript[A, 0]/(2\[Kappa]).*)


BeginPackage["DiffusionBox`"]

SquaredDistanceActive::usage = "SquaredDistanceActive[\[Tau],\[Rho],\[Lambda]] indicates the mean squared distance traveled within the time window \[Tau], by the center locus of a polymer segment with length \[Lambda] that experiences correlated excitations with correlation coefficient \[Rho]."
SquaredDistancePassive::usage = "SquaredDistancePassive[\[Tau],\[Rho],\[Lambda]] indicates the mean squared distance traveled within the time window \[Tau], for a passive model that has the same steady-state conformation as a polymer that experiences correlated excitations; see SquaredDistanceActive[\[Tau],\[Rho],\[Lambda]]."
CorrectionActive::usage = "CorrectionActive[\[Tau]] encodes the effect of correlated activity on diffusion."
CorrectionPassive::usage = "CorrectionPassive[\[Tau]] gives, for the passive surrogate model, the correction term to the diffusion dynamics."

Begin["`Private`"]


SquaredDistanceActive[\[Tau]_,\[Rho]_,\[Lambda]_]:=2(Sqrt[\[Tau]/Pi]+\[Rho]*\[Lambda]^2*CorrectionActive[\[Tau]/\[Lambda]^2]);


SquaredDistancePassive[\[Tau]_,\[Rho]_,\[Lambda]_]:=2(Sqrt[\[Tau]/Pi]+\[Rho]*\[Lambda]^2*CorrectionPassive[\[Tau]/\[Lambda]^2]);


CorrectionActive[\[Tau]_]:=2/Pi^2 NIntegrate[
	(1-Exp[-q^2\[Tau]])/(q^2+k^2) (Sin[k/2]Sin[q/2])/(k q),
{q,-Infinity,Infinity},{k,-Infinity,Infinity},
Method->"AdaptiveMonteCarlo"];


CorrectionPassive[\[Tau]_]:=2/Pi^2 NIntegrate[
(q^2 (1-Exp[-k^2\[Tau]])-k^2 (1-Exp[-q^2\[Tau]]))/(q^4-k^4) (Sin[k/2]Sin[q/2])/(k q),
{q,-Infinity,Infinity},{k,-Infinity,Infinity},
Method->"AdaptiveMonteCarlo"];


End[]
EndPackage[]
