(* ::Package:: *)

(* ::Text:: *)
(*Define Spectral kernels for a Rouse polymer, which are obtained by analytic calculations. We here define the kernels in generic unevaluated form. This form is not useful per se when evaluated, but it will be reused and specialized in many contexts.*)


BeginPackage["Kernel`"]

CorrelationActivity::usage = "CorrelationActivity[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order n and frequency k in response to activity modulations."
CorrelationTension::usage = "CorrelationTension[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order n and frequency k in response to tension modulations."

Begin["`Private`"]


(* Define constructors *)
Options[J]={\[Lambda]->0,\[Kappa]->0};
J[q_,OptionsPattern[]] := OptionValue[\[Lambda]] + q^2 + OptionValue[\[Kappa]] q^4;


(* This Fourier Transform is fine. *)
CorrelationActivity[\[CapitalDelta]s_,k_,n_Integer:0,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Join[{k,\[CapitalDelta]s},Values@List@opts];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp=InverseFourierTransform[
		(q^2-(k/2)^2)^n/(J[q+k/2,opts]+J[q-k/2,opts]), 
		q, \[CapitalDelta]s,
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp=FullSimplify[tmp,Assumptions->assumptions]/.k->Abs[k];
	Return[tmp];
];


(* Mathematica can only carry out this Fourier Transform in for special (\[Lambda],\[Kappa])... *)
CorrelationTension[\[CapitalDelta]s_,k_,n_Integer:0,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Join[{k,\[CapitalDelta]s},Values@List@opts];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp=InverseFourierTransform[
		(q^2-(k/2)^2)^(n+1)/(J[q+k/2,opts]J[q-k/2,opts]), 
		q, \[CapitalDelta]s,
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp=FullSimplify[tmp,Assumptions->assumptions]/.k->Abs[k];
	Return[tmp];
];


End[]
EndPackage[]
