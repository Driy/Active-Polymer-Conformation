(* ::Package:: *)

(* ::Text:: *)
(*Methods for calculating the mean squared distance that a locus travels within some time window. Active processes are here encoded as modulations in the magnitude of statistically independent excitations.*)
(**)
(*Time is measured in units of \[Xi]/\[Kappa]. Squared distance is measured in b^2=Subscript[A, 0]/(2\[Kappa]).*)


BeginPackage["DiffusionDiagonal`"]

SquaredDistanceActive::usage = "SquaredDistanceActive[s,\[Tau],\[Epsilon],\[Lambda]] indicates the mean squared distance traveled by locus s within the time window \[Tau], for activity modulations with magnitude \[Epsilon] and wave length \[Lambda]."
SquaredDistancePassive::usage = "SquaredDistancePassive[s,\[Tau],\[Epsilon],\[Lambda]] indicates the mean squared distance traveled by locus s within the time window \[Tau], for a passive model that has the same steady-state conformation as a polymer that experiences activity modulations; see SquaredDistanceActive[s,\[Tau],\[Epsilon],\[Lambda]]."
CorrectionActive::usage = "CorrectionActive[\[Tau]] encodes the effect of activity modulations on diffusion."
CorrectionPassive::usage = "CorrectionPassive[\[Tau]] gives, for the passive surrogate model, the correction term to the diffusion dynamics."

Begin["`Private`"]


SquaredDistanceActive[s_,\[Tau]_,\[Epsilon]_,\[Lambda]_]:=2(Sqrt[\[Tau]/Pi]+\[Epsilon]*\[Lambda]*Cos[s/\[Lambda]]*CorrectionActive[\[Tau]/\[Lambda]^2]);


SquaredDistancePassive[s_,\[Tau]_,\[Epsilon]_,\[Lambda]_]:=2(Sqrt[\[Tau]/Pi]+\[Epsilon]*\[Lambda]*Cos[s/\[Lambda]]*CorrectionPassive[\[Tau]/\[Lambda]^2]);


CorrectionActive[\[Tau]_]:=1/Pi NIntegrate[
	(1-Exp[-q^2\[Tau]])/(q^2+(q+1)^2),
{q,-Infinity,Infinity}]


CorrectionPassive[\[Tau]_]:=1/Pi NIntegrate[
	(q^2 (1-Exp[-(q+1)^2\[Tau]])-(q+1)^2 (1-Exp[-q^2\[Tau]]))/(q^4-(q+1)^4),
{q,-Infinity,Infinity}]


End[]
EndPackage[]
