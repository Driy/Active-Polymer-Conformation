(* ::Package:: *)

(* ::Text:: *)
(*Define Spectral kernels for a Rouse polymer, which are obtained by analytic calculations. We here define the kernels in generic unevaluated form. This form is not useful per se when evaluated, but it will be reused and specialized in many contexts.*)
(**)
(*Note that we here assume that all parameters, including the magnitude of the separation in sequence space \[CapitalDelta]s, are positive. In this way, we exclude any \[Delta]-Distribution contributions and their derivatives.*)


BeginPackage["KernelGeneric`"]

CorrelationActivity::usage = "CorrelationActivity[\[CapitalDelta]s,k,n] gives the spectral kernel of order n and frequency k in response to activity modulations."
CorrelationTension::usage = "CorrelationTension[\[CapitalDelta]s,k,n] gives the spectral kernel of order n and frequency k in response to tension modulations."
StiffnessActivity::usage = "StiffnessActivity[\[CapitalDelta]s,k] calculates a matrix of effective Hookean interactions."

Begin["`Private`"]


(* Define Jacobian *)
(* Measured in terms of characteristic activity over characteristic length squared times inverse friction *)
Options[J]={\[Lambda]->0,\[Kappa]->0};
J[q_,OptionsPattern[]] := (1/2)(OptionValue[\[Lambda]] + q^2 + OptionValue[\[Kappa]] q^4);


(* This Fourier Transform is fine. *)
CorrelationActivity[\[CapitalDelta]s_,k_:0,n_Integer:0,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Join[{\[CapitalDelta]s, k},Values@List@opts];
	vars=Select[vars,Head[#]==Symbol&];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp=InverseFourierTransform[
		(q^2-(k/2)^2)^n/(J[q+k/2,opts]+J[q-k/2,opts]), 
		q, \[CapitalDelta]s,
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp=FullSimplify[tmp,Assumptions->assumptions]/.k->Abs[k];
	Return[tmp];
];


StiffnessActivity[\[CapitalDelta]s_,k_:1,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Join[{\[CapitalDelta]s, k},Values@List@opts];
	vars=Select[vars,Head[#]==Symbol&];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp=InverseFourierTransform[
		(J[q+k/2,opts]J[q-k/2,opts])/(J[q+k/2,opts]+J[q-k/2,opts]), 
		q, \[CapitalDelta]s,
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp=FullSimplify[tmp,Assumptions->assumptions]/.k->Abs[k];
	Return[4 tmp];
];


(* Mathematica can only carry out this Fourier Transform in for special (\[Lambda],\[Kappa])... *)
CorrelationTension[\[CapitalDelta]s_,k_:0,n_Integer:0,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Join[{\[CapitalDelta]s, k},Values@List@opts];
	vars=Select[vars,Head[#]==Symbol&];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp=InverseFourierTransform[
		(q^2-(k/2)^2)^(n+1)/(4 J[q+k/2,opts]J[q-k/2,opts]), 
		q, \[CapitalDelta]s,
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp=FullSimplify[tmp,Assumptions->assumptions]/.k->Abs[k];
	Return[-tmp];
];


End[]
EndPackage[]
