(* ::Package:: *)

(* ::Text:: *)
(*Define Spectral kernels for a Rouse polymer, which are obtained by analytic calculations. We here define the kernels in generic unevaluated form. This form is not useful per se when evaluated, but it will be reused and specialized in many contexts.*)
(**)
(*Note that we here assume that all parameters, including the magnitude of the separation in sequence space \[CapitalDelta]s, are positive. In this way, we exclude any \[Delta]-Distribution contributions and their derivatives.*)


BeginPackage["KernelGeneric`"]

CorrelationActivity::usage = "CorrelationActivity[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order n and frequency k in response to activity modulations."
CorrelationTension::usage = "CorrelationTension[\[CapitalDelta]s_,k_,n] gives the spectral kernel of order n and frequency k in response to tension modulations."
StiffnessActivity::usage = "StiffnessActivity[\[CapitalDelta]s_] calculates a matrix of effective Hookean interactions."

Begin["`Private`"]


(* Define Jacobian *)
Options[J]={\[Lambda]->0,\[Kappa]->0};
J[q_,OptionsPattern[]] := OptionValue[\[Lambda]] + q^2 + OptionValue[\[Kappa]] q^4;


(* This Fourier Transform is fine. *)
CorrelationActivity[s1_,s2_,n_Integer:0,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Values@List@opts;
	vars=Select[vars,Head[#]==Symbol&];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp[1]=InverseFourierTransform[
		(-q k)^(n)/(J[q,opts]+J[k,opts]), 
		{q,k}, {s1,s2},
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp[2]=InverseFourierTransform[
		(-q k)^(n)/(J[q,opts]+J[k,opts]), 
		{k,q}, {s2,s1},
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp[3]=(tmp[1]+tmp[2])/2;
	tmp[3]=FullSimplify[tmp[3],Assumptions->assumptions];
	Return[tmp[3]];
];


StiffnessActivity[s1_,s2_,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Values@List@opts;
	vars=Select[vars,Head[#]==Symbol&];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp[1]=InverseFourierTransform[
		(J[q,opts]J[k,opts])/(J[q,opts]+J[k,opts]), 
		{q,k}, {s1,s2},
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp[2]=InverseFourierTransform[
		(J[q,opts]J[k,opts])/(J[q,opts]+J[k,opts]), 
		{k,q}, {s2,s1},
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp[3]=(tmp[1]+tmp[2])/2;
	tmp[3]=FullSimplify[tmp[3],Assumptions->assumptions];
	Return[tmp[3]];
];


(* Mathematica can only carry out this Fourier Transform in for special (\[Lambda],\[Kappa])... *)
CorrelationTension[s1_,s2_,n_Integer:0,opts:OptionsPattern[]] := Module[{tmp,vars,assumptions},
	vars=Values@List@opts;
	vars=Select[vars,Head[#]==Symbol&];
	assumptions=Join[Map[#\[Element]Reals&,vars],Map[#>0&,vars]];
	tmp[1]=InverseFourierTransform[
		(-q k)^(n+1)/(J[q,opts]J[k,opts]), 
		{q,k}, {s1,s2},
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp[2]=InverseFourierTransform[
		(-q k)^(n+1)/(J[q,opts]J[k,opts]), 
		{k,q}, {s2,s1},
		FourierParameters->{1,-1},
		Assumptions->assumptions];
	tmp[3]=(tmp[1]+tmp[2])/2;
	tmp[3]=FullSimplify[tmp[3],Assumptions->assumptions];
	Return[tmp[3]];
];


End[]
EndPackage[]
