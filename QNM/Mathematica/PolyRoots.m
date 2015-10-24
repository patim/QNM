(* ::Package:: *)

BeginPackage["PolyRoots`"]
\[Omega]up::usage="";


\[Omega]down::usage="";
Mconstruct::usage="";
MySpec::usage="";
Det2::usage="";
DetPlot2::usage="DetPlot2[a1, a2, NN, ind, m, opts:OptionsPattern[]]"<>
" plots the determinant function between a1 and a2 for NN, ind (l) and m.";
Bisec::usage="";
FPM::usage="";
DetZero::usage="";
BraDetZeros::usage="";
DetZeros2::usage="";
DetZerosForN::usage="";
DetZeroOvertone::usage="";


Unprotect[s];
s = -2;
Protect[s];
Protect[RFMethod, a0, fpm, bisec, newton, Eps, Points, Upper];


Begin["Private`"]


\[Omega]up[m_, a_, N_] := (a m - I N Sqrt[1 - a^2])/(
  2 (1 + Sqrt[1 - a^2]));


\[Omega]down[N_] := -((I N)/4);


Mconstruct[m_, a_, N_, A_, up_: True] := 
  Module[{M, \[Sigma]p, \[Sigma]m, \[Zeta], \[Xi], \[Eta], 
    p, \[Gamma], \[Delta], \[Alpha], \[Sigma], g, f, 
    h, \[Eta]sign, \[Omega], q},
   g[k_] := Re[k (k - 4 p + \[Gamma] + \[Delta] - 1) - \[Sigma]];
   f[k_] := -Re[(k + 1) (k + \[Gamma])];
   h[k_] := Re[4 p (k + \[Alpha] - 1)];
   
   \[Omega] = 
    If[up, q = -1 - s + N; \[Omega]up[m, a, N], 
     q = -1 + N; \[Omega]down[N], q = -1 - s + N; \[Omega]up[m, a, N]];
   
   \[Sigma]p = (2 \[Omega] (1 + Sqrt[1 - a^2]) - 
       m a)/(2 Sqrt[
        1 - a^2]); \[Sigma]m = (2 \[Omega] (1 - Sqrt[1 - a^2]) - 
       m a)/(2 Sqrt[1 - a^2]);
   \[Zeta] = I \[Omega]; \[Xi] = (-s - (s + 2 I \[Sigma]p))/2;
   \[Eta]sign = If[up, 1, -1, 1];
   \[Eta] = (-s + \[Eta]sign (s - 2 I \[Sigma]m))/2;
   p = Sqrt[1 - a^2] \[Zeta];
   \[Gamma] = 1 + s + 2 \[Eta];
   \[Delta] = 1 + s + 2 \[Xi];
   \[Alpha] = 
    1 + s + \[Xi] + \[Eta] - 2 \[Zeta] + s (I \[Omega])/\[Zeta];
   \[Sigma] = 
    A + (a \[Omega])^2 - 8 \[Omega]^2 + 
     p (2 \[Alpha] + \[Gamma] - \[Delta]) + (1 + 
        s - (\[Gamma] + \[Delta])/2) (s + (\[Gamma] + \[Delta])/2);
   M = SparseArray[{Band[{2, 1}] -> h[Range[1, q]], 
      Band[{1, 1}] -> g[Range[0, q]], 
      Band[{1, 2}] -> f[Range[0, q - 1]]}];
   M
   ];


(*The first argument trncN is passed as "a reference" to the existing \
variable. It updates its value if there is a need to increase the \
number of terms in the spectral expansion for higher accuracy of the \ 
solution*)
SetAttributes[MySpec, HoldFirst];
MySpec[trncN_, index_, m_, c_, eps_] := Module[{specOut, A, e},
  specOut = KerrQNM`AngularSpectralRootForl[s, index, m, c, trncN];
  e = Abs[specOut[[3, trncN]]];
  While[e > eps,
   (*index variable is sort of like l (mulitpole number)*)
   
   specOut = KerrQNM`AngularSpectralRootForl[s, index, m, c, ++trncN];
   e = Abs[specOut[[3, trncN]]];
   ];
  A = specOut[[1]];
  A
  ]


(*Rather than calculating the determinant of the matrix we extract \
the smallest value from the diagonal after the diagonalization of \
matrix matr into w. We also try to get the valid sign of the determinant 
from the information contained in matrices*)
Det2[matr_] := Module[{u, w, v, wdiag, det, sign},
  {u, w, v} = 
   SingularValueDecomposition[matr, Tolerance -> 0];
  wdiag = Normal[Diagonal[w]];
  sign = Det[SetPrecision[u.v, $MinPrecision]];
  det = sign*wdiag[[Position[Ordering[wdiag], 1][[1]]]][[1]];
  det
  ]


Options[DetPlot2] = 
  Union[{Points -> 250}, {Upper -> True}, Options[ListPlot]];
DetPlot2[a1_, a2_, NN_, ind_, m_, opts : OptionsPattern[]] := 
 Module[{a, detList = {}, prec = $MinPrecision, M, da, det, aroot, c, 
   A, trncN = 4, \[Omega], up},
  up = OptionValue[Upper];
  c[a_] := Module[{\[Omega]},
    \[Omega] = If[up, \[Omega]up[m, a, NN], \[Omega]down[NN]];
    N[a*\[Omega], prec]];
  da = N[(a2 - a1)/OptionValue[Points], prec];
  
  For[a = SetPrecision[a1, $MinPrecision], 
   a <= SetPrecision[a2, $MinPrecision], 
   a = SetPrecision[a + da, $MinPrecision],
   A = MySpec[trncN, ind, m, c[a], 10^-12];
   M = Mconstruct[m, a, NN, A, up];
   det = Det2[M];
   detList = Join[detList, {{a, det}}];
   ];
  (*Print["det: ", Precision[det], " a: ", Precision[a]];*)
  ListPlot[detList, FilterRules[{opts}, Options[ListPlot]]]
  ]


SyntaxInformation[
   Bisec] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SetAttributes[Bisec, HoldAll];

Options[Bisec] = {Eps -> 10^-12};

(*Bisection method of root finding*)

Bisec[f_, intvl_, OptionsPattern[]] := 
 Module[{ff, xroot, xrootNext, dx = 1, x, x1, x2, i},
  x1 = SetPrecision[intvl[[2]], $MinPrecision];
  x2 = SetPrecision[intvl[[3]], $MinPrecision];
  x = SetPrecision[intvl[[1]], $MinPrecision];
  (*ff=SetPrecision[Function[Evaluate[x],f],$MinPrecision];*)
  
  ff = Function[Evaluate[x], f];
  If[ff[x1]*ff[x2] > 0,
	Print["No roots within the range [", x1, ", ", x2, "]"]; 
   Abort[];,Null, Abort[]];
  xroot = x1;
  i = 1;
  While[dx > OptionValue[Eps],
   xrootNext = (x1 + x2)/2;
   dx = Abs[xroot - xrootNext];
   If[dx <= OptionValue[Eps],
    Break[];
    ,
    If[ff[xrootNext]*ff[x1] < 0, x2 = xrootNext, x1 = xrootNext];
    xroot = xrootNext;
    ];
   i++;
   ];
  {xrootNext, dx, ff[xrootNext]}
  ]


SyntaxInformation[
   FPM] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SetAttributes[FPM, HoldAll];

Options[FPM] = {Eps -> 10^-12};

(*False Position Method of root finding*)

FPM[f_, intvl_, OptionsPattern[]] := 
 Module[{ff, xroot, xrootNext, dx = 1, x, x1, x2, i, df, y0, 
		 ff1, ff2, ffNext},
  x1 = SetPrecision[intvl[[2]], $MinPrecision];
  x2 = SetPrecision[intvl[[3]], $MinPrecision];
  x = SetPrecision[intvl[[1]], $MinPrecision];
  (*ff=SetPrecision[Function[Evaluate[x],f],$MinPrecision];*)
  
  ff = Function[Evaluate[x], f];
  If[ff[x1]*ff[x2] > 0,
	Print["No roots within the range [", x1, ", ", x2, "]"]; 
   Abort[];, Null, Abort[]];
  xroot = x1;
  i = 1;
  While[dx > OptionValue[Eps],
   ff1 = ff[x1]; ff2 = ff[x2];
   df = (ff2-ff1)/(x2-x1);
   y0 = (ff1 x2 - ff2 x1)/(x2-x1);
   xrootNext = - y0/df;
   dx = Abs[xroot - xrootNext];
   ffNext = ff[xrootNext];
   If[dx <= OptionValue[Eps],
    Break[];
    ,
    If[ffNext*ff1 < 0, x2 = xrootNext, x1 = xrootNext];
    xroot = xrootNext;
   ];
   i++;
   ];
  {xrootNext, dx, ffNext}
  ]


SyntaxInformation[
   Newton] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SetAttributes[Newton, HoldAll];

Options[Newton] = {Eps -> 10^-12};

(*Newton Method of root finding*)
Newton[f_, intvl_, x0_, OptionsPattern[]]:=
Module[{x, ff, i, maxit=50, dx=1, \[CapitalDelta]x=10^-6, df, xprev, xnext, ffPrev, ffNext},
  
  x = SetPrecision[intvl[[1]], $MinPrecision];
  ff = Function[Evaluate[x], f];
  xprev = SetPrecision[x0, $MinPrecision];
  i=1;
  While[dx > OptionValue[Eps],
	ffNext = ff[xprev+\[CapitalDelta]x];
	ffPrev = ff[xprev];

	df = (ffNext-ffPrev)/\[CapitalDelta]x;
	xnext = xprev-ffPrev/df;
	dx = Abs[xnext-xprev];
	xprev = xnext;
	i++;
	If[i==maxit, 
		Print["Newton: the method has reached the maximum of iterations=", 
		maxit,". Ok, bye!"];Abort[];
	];
  ];
  Print[i];
 {xnext,dx,ffNext}
]


Options[DetZeroOvertone] = {Eps -> 10^-12, Upper -> True, RFMethod -> newton};

(*allows changing the the first element within the function*)
SetAttributes[DetZeroOvertone, HoldFirst]; 
DetZeroOvertone[n_, afin_, ind_, m_, OptionsPattern[]] := 
Module[{fit, i, a, araw1, araw2, aguess1=0, aguess2, bra1, bra2, aroot, e, up, c, aDet, 
		aDetfunc, polyN, trncN=4, prec=$MinPrecision, fitpts, ptsnum=3, polyNi},
  up = OptionValue[Upper];
  e = OptionValue[Eps];
  c[a_,pN_] := Module[{\[Omega]},
    \[Omega] = If[up, \[Omega]up[m, a, pN], \[Omega]down[pN]];
    N[a*\[Omega], prec]
  ];
  aDet[a_,pN_] := Module[{A},
    A = MySpec[trncN, ind, m, c[a,pN], e];
    Det2[Mconstruct[m, a, pN, A, up]]
  ];
  fitpts = Take[n,-ptsnum];
  polyN = n[[-1,3]]+1;(*polynomial N for the next root*)
  fitpts = Table[{fitpts[[i,1]],fitpts[[i,2]]}, {i,1,ptsnum}];
  fit = Fit[fitpts, {1, a}, a];
  araw1 = NSolve[fit==I*\[Omega]up[m, a, polyN], a];
  aguess1 = Max[araw1[[1,1,2]],araw1[[2,1,2]]]; (*extracting the possitve root*)

  For[polyNi = polyN, aguess1 <= afin, polyNi++, 
	araw2 = NSolve[fit==I*\[Omega]up[m, a, polyNi+1], a];
	aguess2 = Max[araw2[[1,1,2]],araw2[[2,1,2]]]; (*extracting the possitve root*)
	
(*brackets are half-way between the roots*)
	bra1 = (aguess1+fitpts[[2,1]])/2;
	bra2 = (aguess1+aguess2)/2;

	aDetfunc[a_]:= aDet[a, polyNi];
	aroot = FPM[aDetfunc[a], {a, bra1, bra2}, Eps->e];
	Print["aroot: ",aroot];
	AppendTo[n, {aroot[[1]],I*\[Omega]up[m, aroot[[1]], polyNi],polyNi}];
	(*swapping variables*)
	aguess1 = aguess2;
	fitpts = {fitpts[[2]], fitpts[[3]], {aroot[[1]], I*\[Omega]up[m, aroot[[1]], polyNi]}};
    fit = Fit[fitpts, {1, a}, a];
  ]
]


(*
  Tolerance Eps is used in MySpec as well as in root finding routines 
*)
Options[DetZero] = {Eps -> 10^-12, Upper -> True, RFMethod -> fpm, 
					a0->0};
DetZero[a1_, a2_, polyN_, ind_, m_, OptionsPattern[]] := 
 Module[{a1hp, a2hp, a, aN, M, prec = $MinPrecision, aroot, c, A, 
   det1, det2, k, aDet, p, up, trncN = 4, e, da, \[Omega]},
  a1hp = SetPrecision[a1, prec];
  a2hp = SetPrecision[a2, prec];
 
  up = OptionValue[Upper];
  e = OptionValue[Eps];
  c[a_] := Module[{\[Omega]},
    \[Omega] = If[up, \[Omega]up[m, a, polyN], \[Omega]down[polyN]];
    N[a*\[Omega], prec]];

  A = MySpec[trncN, ind, m, c[a1hp], e];
  M = Mconstruct[m, a1hp, polyN, A, up];
  det1 = Det2[M];

  A = MySpec[trncN, ind, m, c[a2hp], e];
  M = Mconstruct[m, a2hp, polyN, A, up];
  det2 = Det2[M];

  If[ det1*det2 > 0, 
   Print["No roots within the range [", N[a1hp], ", ", N[a2hp], "]"]; 
   Abort[];, Null, Abort[]];
  aDet[a_] := Module[{A},
    A = MySpec[trncN, ind, m, c[a], e];
    Det2[Mconstruct[m, a, polyN, A, up]]
    ];

  Switch[OptionValue[RFMethod],
	fpm, aroot = FPM[aDet[a], {a, a1hp, a2hp}, Eps -> e];,
	bisec, aroot = Bisec[aDet[a], {a, a1hp, a2hp}, Eps -> e];,
	newton, aroot = Newton[aDet[a], {a, a1hp, a2hp}, OptionValue[a0], Eps->e], 
	_, Print["DetZero: Unknown root finding method"]; Abort[];
  ];

  aN = aroot[[1]];
  A = MySpec[trncN, ind, m, c[aroot[[1]]], e];
  da = a2 - a1;
  \[Omega] = If[up, \[Omega]up[m, aN, polyN], \[Omega]down[polyN]];
  {aN, \[Omega], N[aroot[[3]], prec], A, da}
  ]


(*Bracketing the roots
  Tolerance option Eps is for MySpec
*)
Options[BraDetZeros] = {Eps -> 10^-12, Upper -> True, Points -> 50};
BraDetZeros[a1_, a2_, polyN_, ind_, m_, OptionsPattern[]] := 
 Module[{i, a1hp, a2hp, da, daflag, 
   damin = 10^-8, \[CapitalDelta]a1, \[CapitalDelta]a2, k, aDet, c, 
   trncN = 4, a, A, M, \[Omega], up, 
   prec = $MinPrecision, cntnflag=False},

  up = OptionValue[Upper];
  c[a_] := Module[{\[Omega]},
    \[Omega] = If[up, \[Omega]up[m, a, polyN], \[Omega]down[polyN]];
    N[a*\[Omega], prec]];
  
  aDet[a_] := Module[{A},
    A = MySpec[trncN, ind, m, c[a], OptionValue[Eps]];
    Det2[Mconstruct[m, a, polyN, A, up]]
    ];

  If[Head[Global`brackets]==List && Length[Global`brackets]==1,
	Clear[Global`brackets];
  ];

  (*If the brackets list has enteries*)
  If[Head[Global`brackets]==List && Length[Global`brackets]>1,
	i = Length[Global`brackets];
	da = Global`brackets[[i,2]] - Global`brackets[[i,1]];
	a = Global`brackets[[i,2]];	
	a2hp = SetPrecision[a2, prec];
	cntnflag=True;
  ];

  If[Head[Global`brackets]==Symbol 
	|| (Head[Global`brackets]==List && Length[Global`brackets]==1),
	Global`brackets={{polyN,ind,m}};
	a1hp = SetPrecision[a1, prec];
	a2hp = SetPrecision[a2, prec];
	da = (a2hp - a1hp)/OptionValue[Points];
    i = 2;
	a = a1hp + da;
  ];

  daflag = False;
  For[Null, a < a2hp, a = a + da,
   If[aDet[a]*aDet[a - da] < 0,
     If[cntnflag==True, Null, AppendTo[Global`brackets, {a - da, a}]];
	 cntnflag=False;
     If[da > damin,
      If[i > 3,
        \[CapitalDelta]a2 = a - Global`brackets[[i - 1, 1]];
        \[CapitalDelta]a1 = 
         Global`brackets[[i - 1, 1]] - Global`brackets[[i - 2, 1]];
        k = \[CapitalDelta]a2/\[CapitalDelta]a1;
        If[k < 1,
         a = a + da; da = k da;
         ];
        ];
      , daflag = True;
      ];
     i++;
     ];
   ];
  If[daflag, Print["da reached the minimum: ", da]];
  Global`brackets
  ]


(* Finds all the roots on the interval [a1,a2] by first calling the bracketing routine
   BraDetZeros[] and then refinig it with DetZero (using false position method as a 
   default method) 
*)
Options[DetZeros2] = Union[Options[BraDetZeros], Options[DetZero]];
DetZeros2[a1_, a2_, polyN_, ind_, m_, opts:OptionsPattern[]] := 
 Module[{i, bra1, bra2, detzero, a1hp, a2hp, up},
  a1hp = SetPrecision[a1, $MinPrecision];
  a2hp = SetPrecision[a2, $MinPrecision];
  up = OptionValue[Upper];
  (*Creates Global`brackets*)
(*AppendTo[Global`timeDebug[[1]],AbsoluteTiming[BraDetZeros[a1hp, a2hp, polyN, ind, m, 
    FilterRules[{opts}, Options[BraDetZeros]]]][[1]]];*)
 
 BraDetZeros[a1hp, a2hp, polyN, ind, m, FilterRules[{opts}, Options[BraDetZeros]]]; 
 Print["Bracketing finished for N=", polyN, " bratime=",Global`timeDebug[[1,-1]]];

  If[Head[Global`roots]==List && Length[Global`roots]==1,
	Clear[Global`roots]
  ];

  If[Head[Global`roots]==Symbol, 
	Global`roots={{polyN,ind,m}};i=2;
  ];
  
  If[Head[Global`roots]==List,
	i=Length[Global`roots]+1;
  ];

(*AppendTo[Global`timeDebug[[2]],AbsoluteTiming[*)
  For[Null, i <= Length[Global`brackets], i++,
   bra1 = Global`brackets[[i, 1]];
   bra2 = Global`brackets[[i, 2]];
   detzero = DetZero[bra1, bra2, polyN, ind, m, 
				FilterRules[{opts}, Options[DetZero]]];
   AppendTo[Global`roots, detzero];
   ];
  (*][[1]]];*)
Print["Root refining time=",Global`timeDebug[[2,-1]]];
  AbortProtect[
   Clear[Global`brackets];
   Print["Brackets cleard for N=",polyN];
  ];
  Take[Global`roots, {2,-1}]
  ]


Options[DetZerosForN] = Options[DetZeros2];
DetZerosForN[a1_, a2_, N1_, N2_, ind_, m_, opts : OptionsPattern[]] := 
 Module[{iN, iNroots, \[Omega], points, up, ptsadd = 15, Nstart},
  If[Head[Global`finroots]==Symbol, Global`finroots = {}];
  up = OptionValue[Upper];
  
  Nstart = N1;
  If[Head[Global`finroots]==List && Length[Global`finroots]>0,
	Nstart=Global`finroots[[-1,2,1]]+1;
	If[Nstart<= N2,
		Print["Glboal`finroots has enteries. Starting from N=",Nstart];
		,Print["Glboal`finroots has enteries. All the roots between 
				N1=",N1," and N2=",N2," have been calculated"];
	];
  ];	

  If[Nstart > 4,
   points = OptionValue[Points] + ptsadd (Nstart - 4),
   points = OptionValue[Points];,
   points = OptionValue[Points];
   ];

  For[iN = Nstart, iN <= N2, iN++,
   points = points + ptsadd;
   iNroots = DetZeros2[a1, a2, iN, ind, m, Points->points, FilterRules[{opts}, Options[DetZeros2]]];
   
   AbortProtect[
    AppendTo[Global`finroots, {iNroots, {iN, ind, m}}];
    Print["Finished N=", iN];
    Clear[Global`roots];
    Print["temporary roots cleared for N=", iN];
   ];
   (*If[Mod[iN, 5] == 0,
    Save["AutoSaveRootsN" <> ToString[iN], Global`finroots];
    Print["Saved data for N=", iN];
    ];*)
   ];
  Global`finroots
  ]


End[] (*Private*)


EndPackage[]
