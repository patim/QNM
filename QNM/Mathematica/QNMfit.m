(* ::Package:: *)

(*
QNM multimode ringdown fitting 
Maxim Zalutskiy
Gregory Cook
2013
*)

BeginPackage["QNMfit`"]
(*RD::usage="RD[\[Alpha],\[Beta],lm,modes,pmodes,mass,spin,\[Theta]angle,t,t0] represents ";*)
ReadKerrQNM::usage="ReadKerrQNM[l, m, n] loads quasi-normal modes "<>
"for given l, m and n.";

Clm::usage="Clm[l, m, modes, pmodes, \[Delta], a, t, t0=0] generates "<>
"\!\(\*SubscriptBox[\(C\), \(lm\)]\) of the ringdown
\!\(\*SubscriptBox[\(C\), \(lm\)]\)(t) = \!\(\*SubscriptBox[\(\[CapitalSigma]\), "<>
"\(l' m' n\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(d\), \(l'\)], "<>
"\(mm'\)]\)(\[Theta]){\!\(\*SubscriptBox[\(A\), \(l' m' n\)]\) "<>
"\!\(\*SubscriptBox[\(\[ScriptCapitalA]\), \(ll' m'\)]\)(\!\(\*SubscriptBox[\(c\), "<>
"\(l' m' n\)]\))\!\(\*SuperscriptBox[\(e\), \(\(-\*SubscriptBox[\(i\[Omega]\), "<>
"\(l' m' n\)]\) t + \*SubscriptBox[\(i\[Phi]\), \(l' m' n\)]\)]\) "<>
"\!\(\*SuperscriptBox[\(e\), \(\(-t\)/\*SubscriptBox[\(\[Tau]\), \(l' m' n\)]\)]\)"<>"
	 + (-1\!\(\*SuperscriptBox[\()\), \(l + l'\)]\)\!\(\*SubscriptBox"<>
"[SuperscriptBox[\(A\), \('\)], \(l' \((\(-m'\))\) n\)]\) "<>
"\!\(\*SubscriptBox[SuperscriptBox[\(\[ScriptCapitalA]\), \(*\)], \(ll' \((\(-m'\))\)\)]\)"<>
"(\!\(\*SubscriptBox[\(c\), \(l' \((\(-m'\))\) n\)]\))\!\(\*SuperscriptBox"<>
"[\(e\), "<>"\(\(-\*SubscriptBox[\(i\[Omega]\), \(l' \((\(-m'\))\) n\)]\) t + "<>
"\*SubscriptBox[\(i\[Phi]\), \(l' \((\(-m'\))\) n\)]\)]\) \!\(\*SuperscriptBox"<>
"[\(e\), \(\(-t\)/\*SubscriptBox[\(\[Tau]\), \(l' \((\(-m'\))\) n\)]\)]\)}."<>
"
The first term coresponds to the 'positive m' modes. The second term "<>
"corresponds to 'negative m' primed modes (pmodes)"<>"
modes: {{\!\(\*SubscriptBox[\(ll\), \(1\)]\),m,\!\(\*SubscriptBox[\(n\), "<>
"\(1\)]\),\!\(\*SubscriptBox[\(A\), \(1\)]\),\!\(\*SubscriptBox[\(\[Phi]\), "<>
"\(1\)]\)},{\!\(\*SubscriptBox[\(ll\), \(2\)]\),m,\!\(\*SubscriptBox[\(n\), "<>
"\(2\)]\),\!\(\*SubscriptBox[\(A\), \(2\)]\),\!\(\*SubscriptBox[\(\[Phi]\), "<>
"\(2\)]\)},...}
pmodes:{{\!\(\*SubscriptBox[\(ll\), \(1\)]\),m,\!\(\*SubscriptBox[\(n\), "<>
"\(1\)]\),\!\(\*SubscriptBox[\(Ap\), \(1\)]\),\!\(\*SubscriptBox[\(\[Phi]p\), "<>
"\(1\)]\)},{\!\(\*SubscriptBox[\(ll\), \(2\)]\),m,\!\(\*SubscriptBox[\(n\), "<>
"\(2\)]\),\!\(\*SubscriptBox[\(Ap\), \(2\)]\),\!\(\*SubscriptBox[\(\[Phi]p\), "<>
"\(2\)]\)},...}.
If the values for A's, \[Phi]'s, Ap's or \[Phi]p's are not provided the function "<>
"will generate the symbolic expression.
\[Delta] - mass ratio, a - spin parameter, t - time, t0 - time shift (set to "<>
"zero by default).";

GenData::usage="GenData[lm, modes, pmodes, \[Delta], a, t1, tN, dt] generates "<>
"rindown data by calling \!\(\*SubscriptBox[\(C\), \(lm\)]\)(t) for multiple"<>
" sets of l and m.
lm: {{\!\(\*SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)},"<>
"{\!\(\*SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)},...}
'Positive m' modes
modes:{{\!\(\*SubscriptBox[\(ll\), \(1\)]\),m,\!\(\*SubscriptBox[\(n\), \(1\)"<>
"]\),\!\(\*SubscriptBox[\(A\), \(1\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\)},{\!\("<>
"\*SubscriptBox[\(ll\), \(2\)]\),m,\!\(\*SubscriptBox[\(n\), \(2\)]\),\!\("<>
"\*SubscriptBox[\(A\), \(2\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\)},...}

'Negative m' primed modes
pmodes:{{\!\(\*SubscriptBox[\(ll\), \(1\)]\),m,\!\(\*SubscriptBox[\(n\), \(1\)]\),"<>
"\!\(\*SubscriptBox[\(Ap\), \(1\)]\),\!\(\*SubscriptBox[\(\[Phi]p\), \(1\)]\)},{\!\(\*"<>
"SubscriptBox[\(ll\), \(2\)]\),m,\!\(\*SubscriptBox[\(n\), \(2\)]\),\!\(\*"<>
"SubscriptBox[\(Ap\), \(2\)]\),\!\(\*SubscriptBox[\(\[Phi]p\), \(2\)]\)},...}.

The output of the function:
{{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)}},
{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!"<>
"\(\*SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)}},...}

\[Delta] - mass ratio, a - spin parameter, t1,tN - time frame, dt - time increment.
";

MyFit::usage="MyFit[data, lm, modes, pmodes] fits ringdown waveform to the data.
data:
{{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!"<>
"\(\*SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)}},
{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\("<>
"\*SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)}},...}

lm:
{{\!\(\*SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)},{\!\(\*"<>
"SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)},...}

modes & pmodes:
{{\!\(\*SubscriptBox[\(ll\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\),\!\(\*"<>
"SubscriptBox[\(n\), \(1\)]\), \!\(\*SubscriptBox[\(A\), \(1\)]\),\!\(\*"<>
"SubscriptBox[\(\[Phi]\), \(1\)]\),Fixed},{\!\(\*SubscriptBox[\(ll\), \(2\)]\),\!\("<>
"\*SubscriptBox[\(m\), \(2\)]\),\!\(\*SubscriptBox[\(n\), \(2\)]\), \!\(\*"<>
"SubscriptBox[\(A\), \(2\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\),Fixed},...}

The check is performed whether the indexes from lm list match those in modes "<>
"and pmodes. 
The optional parameters are A, \[Phi] (serve as initial guesses) and Fixed. If the "<>
"flag Fixed is 
provided and set to True, the values are locked and not used as fit parameters.

Options:
OptRange->{\!\(\*SubscriptBox[\(x\), \(start\)]\),\!\(\*SubscriptBox[\(x\), "<>
"\(finish\)]\)} - optimization range (default OptRange->All)
t0->\!\(\*SubscriptBox[\(t\), \(shift\)]\) - time shift in the wave form "<>
"(default t0->0)
Spin->{Value, Fixed} - intial guess for the value of the spin parameter a, "<>
"if Fixed is True, spin is not used as a fit parameter
Mass->{Value, Fixed} - intial guess for the value of mass \[Delta], if Fixed is "<>
"True, mass is not used as a fit parameter";

DataPlot2::usage="DataPlot[data] plots the data as a function of time and "<>
"and the \!\(\*SubscriptBox[\(C\), \(lm\)]\) number";
DataPlot::usage="DataPlot[data] plots the data points.
data:
{{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)"<>
"\!\(\*SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\("<>
"\*SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)}},
{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\("<>
"\*SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)}},...}";

DataFitPlot::usage="DataFitPlot[data, fit, lm] plots both data points and fit.
data:
{{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!"<>
"\(\*SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)}},
{{{\!\(\*SubscriptBox[\(t\), \(1\)]\),Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\("<>
"\*SubscriptBox[\(]\), \(1\)]\),Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*"<>
"SubscriptBox[\(]\), \(1\)]\)},{\!\(\*SubscriptBox[\(t\), \(2\)]\),Re[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\),Im[\!\(\*"<>
"SubscriptBox[\(C\), \(lm\)]\)\!\(\*SubscriptBox[\(]\), \(2\)]\)},...},{\!\(\*"<>
"SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)}},...}

fit:
an object returned by MyFit routine

lm:
{l,m} lm piece of data/fit being plotted";

ModePlot::usage="ModePlot[fit, lm] plots individual modes of the fit for l "<>
"and m.

fit:
an object returned by MyFit routine

lm:
{l,m} lm piece of the fit being plotted";

GetQNMValues::usage="GetQNMValues[fit] outputs QNM with fitted values
Options:
Fixed->True fixes QNM in the output (Fixed->False by default)";

Getaz::usage = "Getaz[fit] returns \!\(\*SubscriptBox[\(a\), \(z\)]\) "<>
"component of the spin";

Getatheta::usage = "Getatheta[fit, (\!\(\*SuperscriptBox[\(a\), \(c\)]\)"<>
"\!\(\*SubscriptBox[\()\), \(x\)]\),(\!\(\*SuperscriptBox[\(a\), \(c\)]\)\!"<>
"\(\*SubscriptBox[\()\), \(y\)]\),(\!\(\*SuperscriptBox[\(a\), \(z\)]\)\!\(\*"<>
"SubscriptBox[\()\), \(z\)]\)] returns a (spin) and \[Theta] (the angle between the "<>
"spin and the z-axis)";
MassSpinConfRegion::usage = "MassSpinConfRegion[fit] Shows the mass and spin "<>
"confidence region graphically";
Protect[Mass,Spin,Theta,Phi,OptRange,Fixed,yrange,xrange,fontsize,imagesize];

Begin["Private`"]


(* Quasi-normal modes loader*)
ReadKerrQNM[l_,m_,n_]:=Module[{nname,mnamea,mname,lname,rawdat,Nelems,NL,
	Lmin,Lmax,i,a,\[Omega],Alm,Allmn},
	nname=If[n<10,"0"<>ToString[n],ToString[n]];
	mnamea=If[Abs[m]<10,"0"<>ToString[Abs[m]],ToString[Abs[m]]];
	mname=If[m<0,"-"<>mnamea,"+"<>mnamea];
	lname = If[l<10,"0"<>ToString[l],ToString[l]];
	Import[path<>"KerrQNM_"<>nname<>".h5",{"HDF5","Datasets",{"/n"<>
		nname<>"/m"<>mname<>"/L"<>lname<>".dat"}}]
]


indx[n_]:=If[n<0,"m",""]<>ToString[Abs[n]];
Varname[v_,l_,m_,n_]:=ToExpression[ToString[v]<>indx[l]<>indx[m]<>indx[n]];


DefineInterpolations[lm_,modes_,neg_]:=Module[{i,k,msize,l,m,ll,mm,nn,
	NeedQNM,Need\[Omega],Need\[ScriptCapitalA],QNMdata,QNMsize,lmin,lmax,Nl,j,lmsize,Re\[Omega],Im\[Omega],spin},
	msize=Length[modes];
	lmsize=Length[lm];
	For[i=1,i<=msize,++i,
		ll=modes[[i,1]];
		mm=neg*modes[[i,2]];
		nn=modes[[i,3]];
(*		NeedQNM=If[Head[F[ll,mm,nn]]==InterpolatingFunction,False,True,True];
Print["NeedQNM=",NeedQNM];
		If[!NeedQNM,
			For[j=1,j<=lmsize,++j,
				l=lm[[j,1]];
				NeedQNM=If[Head[\[ScriptCapitalA][l,ll,mm,nn]]==InterpolatingFunction,False,
												True,True];
				If[NeedQNM,Break[]];
			];
		];*)
		NeedQNM = True; (*always True since F[ll, mm, nn] has been eliminated*)
		If[NeedQNM,
			QNMdata = ReadKerrQNM[ll,mm,nn];
			QNMsize = Length[QNMdata];
			lmin=Max[2,Abs[mm]];
			Nl=(Length[QNMdata[[1]]]-5)/2;
			lmax = lmin+Nl-1;
			Need\[Omega]=If[Head[\[Omega]bar[ll,mm,nn]]==InterpolatingFunction,False,True,True];
Print["Need\[Omega]=", Need\[Omega]];			
			If[Need\[Omega],
				Re\[Omega][i_]:=QNMdata[[i,2]];
				Im\[Omega][i_]:=QNMdata[[i,3]];
				spin[i_]:=QNMdata[[i,1]];

				\[Omega]bar[ll,mm,nn]=Interpolation[Table[{spin[j],Re\[Omega][j]+I*Im\[Omega][j]},{j,
													1,QNMsize}]];

			];
			For[k=1,k<=lmsize,++k,
				l=lm[[k,1]];
				m=lm[[k,2]];
				Need\[ScriptCapitalA]=If[Head[\[ScriptCapitalA][l,ll,mm,nn]]==InterpolatingFunction,False,True,
												True];				
				If[Need\[ScriptCapitalA],
					\[ScriptCapitalA][l,ll,mm,nn]=Interpolation[
						Table[{QNMdata[[j,1]],QNMdata[[j,2(l-lmin)+6]]
								+QNMdata[[j,2(l-lmin)+7]]I},{j,1,QNMsize}]];
				];
			];
		];
	
	];
]


Clm[l_,m_,modes_,pmodes_,\[Delta]_,a_,\[Theta]_,t_,t0_:0]:= 
	Module[{Sum=0,ModeSum},

	ModeSum[modelist_,Aname_,\[Phi]name_,which_]:=
		Module[{i,ll,mm,nn,A,\[Phi]ase,sum=0,size,isize},
		size=Length[modelist];
		For[i=1,i<=size,i++,
			isize=Length[modelist[[i]]];
			ll=modelist[[i,1]];
			mm=modelist[[i,2]];
			nn=modelist[[i,3]];
	
			If[isize>3,
				If[NumericQ[modelist[[i,4]]],
					A=modelist[[i,4]]
				,
					A=Varname[modelist[[i,4]],ll,mm,nn]
				];
				If[NumericQ[modelist[[i,5]]],
					\[Phi]ase=modelist[[i,5]]
				,
					\[Phi]ase=Varname[modelist[[i,5]],ll,mm,nn]
				];
			,
				A=Varname["Private`"<>Aname,ll,mm,nn];
				\[Phi]ase=Varname["Private`"<>\[Phi]name,ll,mm,nn];
			];
			If[which==1,
				sum+=A*WignerD[{ll,-m,-mm},\[Theta]]*\[ScriptCapitalA][l,ll,mm,nn][Re[a]]
						*Exp[-I*\[Omega]bar[ll,mm,nn][Re[a]]*(t-t0)/\[Delta]+I*\[Phi]ase];
			,
				sum+=(-1)^(l+ll)*A*WignerD[{ll,-m,-mm},\[Theta]]
				*Conjugate[\[ScriptCapitalA][l,ll,-mm,nn][Re[a]]
				*Exp[-I*\[Omega]bar[ll,-mm,nn][Re[a]]*(t-t0)/\[Delta]+I*\[Phi]ase]];
			];
		];
		sum
	];
	Sum = ModeSum[modes,"A","\[Phi]",1] + ModeSum[pmodes,"Ap","\[Phi]p",-1];
	Sum
]


GenData[lm_,modes_,pmodes_,\[Delta]_,a_,\[Theta]_,t1_,tN_,dt_,t0_:0]:=
	Module[{C,l,m,mm,t,imodes,ipmodes,lmsize,size,psize,i,j,datalist={}},
	lmsize = Length[lm];
	size=Length[modes];
	psize=Length[pmodes];
	DefineInterpolations[lm,modes,1];
	DefineInterpolations[lm,pmodes,-1];

	For[i=1, i<=lmsize, i++,
		l=lm[[i,1]];
		m=lm[[i,2]];
		(*loop over "positive" modes*)	
		For[j=1;imodes={},j<= size,j++,
			mm=modes[[j,2]];
			If[m!= mm,Continue[];];
			imodes=imodes~Join~{modes[[j]]};
		];

		(*loop over "primed" modes*)
		For[j=1;ipmodes={}, j<=psize, j++,
			mm=pmodes[[j,2]];
			If[m!= -mm,Continue[];];
			ipmodes=ipmodes~Join~{pmodes[[j]]};
		];
		C=Clm[l,m,imodes,ipmodes,\[Delta],a,\[Theta],t,t0];
		(*Print["lm:",l,m," C:",C];*)
		datalist=datalist~Join~{{Table[{t,Re[C],Im[C]},{t,t1,tN,dt}],lm[[i]]}};
	];
	datalist
]


ClmSplit[\[Alpha]_,l_,m_,modes_,pmodes_,\[Delta]_,a_,\[Theta]_,t_,t0_]:=
	KroneckerDelta[\[Alpha],1] Re[Clm[l,m,modes,pmodes,\[Delta],a,\[Theta],t,t0]] + 
	KroneckerDelta[\[Alpha],0] Im[Clm[l,m,modes,pmodes,\[Delta],a,\[Theta],t,t0]]

RD[\[Alpha]_,\[Beta]_,lm_,modes_,pmodes_,mass_,spin_,\[Theta]angle_,t_,t0_]:=
	Module[{size,i,ModeAnalyze},
	ModeAnalyze[mds_]:=Module[{i,newmodes={},modesSize,isize},
		modesSize = Length[mds];
		For[i = 1, i <= modesSize, i++,
			isize = Length[mds[[i]]];
			(*initial guesses provided or they are empty (zero by default)*)
			If[isize==5 || (isize==6 && mds[[i,6]]==False) || isize==3, 
				newmodes=newmodes~Join~{Take[mds[[i]],3]};
			];	
			If[isize==6 && mds[[i,6]]==True, (*values provided*)
				newmodes=newmodes~Join~{Take[mds[[i]],5]};
			];
		];
		newmodes	
	];

	Clear[\[Delta],a,\[Theta]];
	If[(Length[mass]==2 && mass[[2]]==True), (*value provided*)
		\[Delta] = mass[[1]]; 
	];

	If[(Length[spin]==2 && spin[[2]]==True), (*value provided*)
		a = spin[[1]]; 
	];

	If[(Length[\[Theta]angle]==2 && \[Theta]angle[[2]]==True), (*value provided*)
		\[Theta] = \[Theta]angle[[1]]; 
	];

	size=Length[lm];
	(*i goes over lm list *)
	\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(size\)]\(KroneckerDelta[\[Beta], i]*ClmSplit[\[Alpha], lm[\([i, 1]\)], lm[\([i, 2]\)], ModeAnalyze[modes], \n\t\t\tModeAnalyze[pmodes], \[Delta], a, \[Theta], t, t0]\)\)
]


FitParam[modes_,pmodes_,mass_,spin_,\[Theta]angle_]:=
	Module[{A\[Phi]s,Ap\[Phi]ps,\[Delta]a\[Theta]={},ModeAnalyze,modeAnalyzed,\[Delta]a\[Theta]Analyze,fitmodes,fitpmodes},
	ModeAnalyze[mds_,A_,\[Phi]_]:=Module[{i,ll,m,n,size,As={},\[Phi]s={},isize,fitmodes={}},
		size=Length[mds];
		fitmodes = mds;
		For[i=1, i<=size, i++,
			ll=mds[[i,1]];
			m=mds[[i,2]];
			n=mds[[i,3]];
			isize = Length[mds[[i]]];
			If[isize==3,
				As=Append[As,Varname[A,ll,m,n]];
				\[Phi]s=Append[\[Phi]s,Varname[\[Phi],ll,m,n]];
			, (*else isize > or < 3*)
				If[isize==5 || (isize==6 && mds[[i,6]]==False), 
					(*initial guesses provided*)
					As=Append[As,{Varname[A,ll,m,n],mds[[i,4]]}];
					\[Phi]s=Append[\[Phi]s,{Varname[\[Phi],ll,m,n],mds[[i,5]]}];
				];

				If[isize==6 && mds[[i,6]]==True, (*values provided*)
					fitmodes=Delete[fitmodes,i];
					Continue[];
				];
			];
		];
		{As~Join~\[Phi]s,fitmodes}
	];
	modeAnalyzed = ModeAnalyze[modes,"Private`A","Private`\[Phi]"];
	A\[Phi]s = modeAnalyzed[[1]];
	fitmodes=modeAnalyzed[[2]];
	modeAnalyzed = ModeAnalyze[pmodes,"Private`Ap","Private`\[Phi]p"];
	Ap\[Phi]ps = modeAnalyzed[[1]];
	fitpmodes = modeAnalyzed[[2]];
	Clear[\[Delta],a,\[Theta]];
	\[Delta]a\[Theta]Analyze[\[Delta]a\[Theta]_,var_]:=Module[{\[Delta]a\[Theta]Analyzed={}},
		If[(Length[\[Delta]a\[Theta]]==2 && \[Delta]a\[Theta][[2]]==False) || Length[\[Delta]a\[Theta]]==1,
			\[Delta]a\[Theta]Analyzed = {{var,\[Delta]a\[Theta][[1]]}}; 
		];

		If[Length[\[Delta]a\[Theta]]==0,
			\[Delta]a\[Theta]Analyzed={{var,\[Delta]a\[Theta]}};
		];
		\[Delta]a\[Theta]Analyzed	
	];

	\[Delta]a\[Theta] = Join[\[Delta]a\[Theta]Analyze[mass,ToExpression["Private`\[Delta]"]], 
			\[Delta]a\[Theta]Analyze[spin,ToExpression["Private`a"]], \[Delta]a\[Theta]Analyze[\[Theta]angle,
			ToExpression["Private`\[Theta]"]]];
	{\[Delta]a\[Theta]~Join~A\[Phi]s~Join~Ap\[Phi]ps,fitmodes,fitpmodes}
]


ParFrom[fit_,n_]:=fit[[1]]["BestFitParameters"][[n,2]];

MyFit::nodata="no data for mode `1`.";
Options[MyFit]={OptRange->All, t0->0, Mass->{1,False}, Spin->{0,False}, 
				Theta->{0,False}};
MyFit[data_,lm_,modes_,pmodes_,OptionsPattern[]]:=
	Module[{fit,newData={},i,j,datasize,datar,datai,\[Alpha],\[Beta],m,size,psize,lmsize,ll,
			mm,n,fitparameters,fitmodes, fitpmodes,checkedModes,
			checkedPmodes,newmodes={},newpmodes={},tdata,tstart,tend,j1,j2,
			datalm={},ipos,MSTPReturn,MassOut,SpinOut,ThetaOut,index,modesOnset},

	(* building data lm sequence *)
	datalm=Table[data[[i,2]],{i,1,Length[data]}];

	For[i=1,i<=Length[lm],i++,
		If[Length[Position[datalm,lm[[i]]]]==0, 
			Message[MyFit::nodata,lm[[i]]];Return[]];
	
		(*position of data chunck corresponding to lm mode*)
		ipos=Position[datalm,lm[[i]]][[1,1]];
		
		(*bracketing data, in case OptRange is set*)
		datasize=Length[data[[ipos,1]]];
		tstart=OptionValue[OptRange][[1]];
		tend=OptionValue[OptRange][[2]];
		tdata=Table[data[[ipos,1,j,1]],{j,1,datasize}];
		j1=LengthWhile[tdata,#<tstart &]+1;
		j2=LengthWhile[tdata,#<=tend &];
		tdata=.;
		(*newData is an appropriate format for the fitting routine*)
		datar=Table[{1,i,data[[ipos,1,j,1]],data[[ipos,1,j,2]]},{j,j1,j2}];
		datai=Table[{0,i,data[[ipos,1,j,1]],data[[ipos,1,j,3]]},{j,j1,j2}];
		newData=Join[newData,datar,datai];
	];

(******** eliminating wrong elements of the list ********)
	checkedModes=modes;
	size=Length[modes];
	For[j=1,j<= size,j++,
		ll=modes[[j,1]];
		n=modes[[j,3]];
		If[ll<2 || n<0,checkedModes=Delete[checkedModes,j]];
	];
	checkedPmodes=pmodes;
	psize=Length[pmodes];
	For[j=1,j<= psize,j++,
		ll=pmodes[[j,1]];
		mm=pmodes[[j,2]];
		n=pmodes[[j,3]];
		If[ll<2 || n<0,checkedPmodes=Delete[checkedPmodes,j]];
	];
(*****************************************************)

	newmodes=DeleteDuplicates[checkedModes];
	newpmodes=DeleteDuplicates[checkedPmodes];
	fitparameters = FitParam[newmodes,newpmodes,OptionValue[Mass],
					OptionValue[Spin],OptionValue[Theta]];
	fitmodes = fitparameters[[2]];
	fitpmodes = fitparameters[[3]]; 
	DefineInterpolations[lm,newmodes,1];
	DefineInterpolations[lm,newpmodes,-1];
	t1 = AbsoluteTime[];
	fit = Quiet[NonlinearModelFit[newData,{RD[\[Alpha],\[Beta],lm,newmodes,newpmodes,
				OptionValue[Mass],OptionValue[Spin],OptionValue[Theta],
				t,OptionValue[t0]]},fitparameters[[1]],{\[Alpha],\[Beta],t}, ConfidenceLevel -> .9999999]];
	t2 = AbsoluteTime[];
Print["t2-t1: ",t2-t1];
	index=1; 
	(*sets the correspondence between the mass/spin/theta and the indexing 
	of the fit parameters*)
	MSTPReturn[mstp_,num_]:=Module[{output},
		If[Length[mstp]==2 && mstp[[2]]==True,
			output = mstp[[1]];
		,
			output = ParFrom[{fit},num];
			index++;
		];
		output
	];
	MassOut = MSTPReturn[OptionValue[Mass],1];
	SpinOut = MSTPReturn[OptionValue[Spin],index];
	ThetaOut = MSTPReturn[OptionValue[Theta],index];
	modesOnset = index-1;
	
	{fit,lm,newmodes,newpmodes,MassOut,SpinOut,ThetaOut,modesOnset,
	fitmodes,fitpmodes}
]


MakeLabel[name_,l_,m_]:=
	Row[{Subscript[name<>"(",-2],Subscript[C,ToString[l]<>ToString[m]],")"}];

Opts[label_,fontsize_,imsize_]:=
	{Axes->False,Frame->True,FrameLabel->{Text[Style["\!\(\*StyleBox[\"t\",
	\nFontSlant->\"Italic\"]\)/\!\(\*StyleBox[\"M\",\nFontSlant->\"Italic\"]\)",
	FontSize->fontsize]],Text[Style[label,FontSize->fontsize]]},ImageSize->imsize};

Options[DataPlot]={fontsize->12,imagesize->200};
DataPlot[data_,OptionsPattern[]]:=
	Module[{i,j,lmsize,isize,Relabel,Imlabel,l,m,plots={},fosize,imsize},
	lmsize=Length[data];
	fosize=OptionValue[fontsize];
	imsize=OptionValue[imagesize];
	For[i=1,i<= lmsize,i++,
		l=data[[i,2,1]];
		m=data[[i,2,2]];
		Relabel=MakeLabel["Re",l,m];
		Imlabel=MakeLabel["Im",l,m];

		isize=Length[data[[i,1]]];
		plots=plots~Join~{
		(*Re*)
		ListPlot[Table[{data[[i,1,j,1]],data[[i,1,j,2]]},{j,1,isize}],
			Opts[Relabel,fosize,imsize]]
		(*Im*)
		,ListPlot[Table[{data[[i,1,j,1]],data[[i,1,j,3]]},{j,1,isize}],
			Opts[Imlabel,fosize,imsize]]};
	];
	plots
]


Options[DataPlot2]={fontsize->12,imagesize->600};
DataPlot2[data_,OptionsPattern[]]:=Module[{i,t1,tN,j,t,styleScheme,colorScheme,
RDlist,ReRD,ImRD,lmsize,pRelist={},pImlist={},Clmticks,dataSize,PsiPlot},
	dataSize=Length[data];
	For[i=1,i<=dataSize,i++,
		lmsize=Length[data[[i,1]]];
		RDlist=Table[{data[[i,1,j,1]],data[[i,1,j,2]]},{j,1,lmsize}];
		ReRD[i]=Interpolation[RDlist];
		RDlist=Table[{data[[i,1,j,1]],data[[i,1,j,3]]},{j,1,lmsize}];
		ImRD[i]=Interpolation[RDlist];
		pRelist=pRelist~Join~{{i,t,ReRD[i][t]}};
		pImlist=pImlist~Join~{{i,t,ImRD[i][t]}};
	];
	(*The time interval is taken from the last element of data. 
	Should be the same for all elements*)
	t1=data[[Length[data],1,1,1]];
	tN=data[[Length[data],1,lmsize,1]];
	
	(*colorScheme=ColorData[3,"ColorList"];*)
	colorScheme = {RGBColor[0.2,0.2,0.2],RGBColor[76/85,26/255,28/255],
RGBColor[166/255,86/255,8/51],RGBColor[11/51,42/85,184/255],RGBColor[77/255,35/51,74/255],
RGBColor[1,127/255,0],RGBColor[247/255,43/85,191/255],RGBColor[1,1,1/5]};
	styleScheme=Table[{Thickness[0.005], colorScheme[[i]]},{i,1,Length[colorScheme]}];
	Clmticks=Table[i,{i,1,dataSize}];
	PsiPlot[list_,label_]:=ParametricPlot3D[list,{t,t1,tN},BoxRatios->{3,8,1},
		Axes->{True,True,True},
		AxesLabel->{Text[Style["\!\(\*SubscriptBox[\(C\), \(lm\)]\)",
		FontSize->OptionValue[fontsize]]],
		Text[Style["t/M",FontSize->OptionValue[fontsize]]],
		Text[Style[label,OptionValue[imagesize],
		FontSize->OptionValue[fontsize]]]},ImageSize->OptionValue[imagesize],
		PlotStyle->styleScheme, Ticks->{Clmticks,Automatic,None}];
	
	GraphicsColumn[{PsiPlot[pRelist, "Re(\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(4\)]\))"], 
	PsiPlot[pImlist, "Im(\!\(\*SubscriptBox[\(\[CapitalPsi]\), \(4\)]\))"]},
	ImageSize->OptionValue[imagesize]]
]


DataFitPlot::nodata="cannot find a plot for mode `1`";
Needs["PlotLegends`"];
Options[DataFitPlot]={yrange->{-0.01,0.01},xrange->{0,80},
					fontsize->16,imagesize->600};
DataFitPlot[data_,fit_,lm_,OptionsPattern[]]:=
	Module[{lpoptioins,poptions,eroptioins,Reshowoptions,Imshowoptions,dataCnum,
			fitCnum,size,rdata,idata,fsize,isize,datasize,error,rlp,ilp,erp,eip,
			pr,pi,i,fitlm,fitplot,fitfun,p1,p2,errorlegend},
	fitplot = fit[[1]];
	fitlm = fit[[2]];
	datasize=Length[data];
	fsize=OptionValue[fontsize];
	isize=OptionValue[imagesize];
	fitCnum=Position[fitlm,lm][[1,1]];

	For[i=1,i<= datasize,i++,
		If[data[[i,2,1]]==lm[[1]] && data[[i,2,2]]==lm[[2]],
			dataCnum=i;Break[];
		];
	(*Got to the end of the list, but the requested mode is nowhere to be seen*)
	If[i==datasize,Message[DataFitPlot::nodata,lm];Return[]];
	];
	size=Length[data[[dataCnum,1]]];
	rdata=Table[{data[[dataCnum,1,i,1]],data[[dataCnum,1,i,2]]},{i,1,size}];
	idata=Table[{data[[dataCnum,1,i,1]],data[[dataCnum,1,i,3]]},{i,1,size}];
	lpoptioins={PlotStyle->Blue, PlotRange->{OptionValue[xrange],
				OptionValue[yrange]}};
	poptions={PlotStyle->Red,PlotRange->{All,OptionValue[yrange]}};
	eroptioins={PlotStyle->Black,Joined->True, PlotRange->{OptionValue[xrange],
				OptionValue[yrange]}};

	Reshowoptions=Opts[MakeLabel["Re",lm[[1]],lm[[2]]],fsize,isize];
	Imshowoptions=Opts[MakeLabel["Im",lm[[1]],lm[[2]]],fsize,isize];
	
	fitfun[t_]=fitplot[1,fitCnum,t];
	error=Table[{rdata[[i,1]],fitfun[rdata[[i,1]]]-rdata[[i,2]]},{i,1,Length[rdata]}];
    pr=Plot[fitfun[t],{t,OptionValue[xrange][[1]],OptionValue[xrange][[2]]},
			PlotStyle->Red,PlotRange->{All,OptionValue[yrange]}];
    rlp = ListPlot[rdata,lpoptioins];
    erp = ListPlot[error,eroptioins];
	
	fitfun[t_]=fitplot[0,fitCnum,t];
	error=Table[{idata[[i,1]],fitfun[idata[[i,1]]]-idata[[i,2]]},{i,1,Length[idata]}];
	pi=Plot[fitfun[t],{t,OptionValue[xrange][[1]],OptionValue[xrange][[2]]},
		PlotStyle->Red,PlotRange->{All,OptionValue[yrange]}];
	ilp=ListPlot[idata,lpoptioins];
	eip=ListPlot[error,eroptioins]; 
	
	errorlegend=Text[Style["error",FontSize->fsize]];

	p1=ShowLegend[Show[rlp,erp,pr,Reshowoptions],
	{{{Graphics[{Black,Line[{{0,0},{1.5,0}}]}], errorlegend}},
	LegendShadow->None,LegendSize->0.3,LegendTextSpace->2,
	LegendPosition->{0.58,0.42}}];	
	
	p2=ShowLegend[Show[ilp,eip,pi,Imshowoptions],
	{{{Graphics[{Black,Line[{{0,0},{1.5,0}}]}], errorlegend}},
	LegendShadow->None,LegendSize->0.3(*0.2*),LegendTextSpace->2,
	LegendPosition->{0.58,0.42(*0.68,0.48*)}}];
	p1
	p2
]


Options[ModePlot]={yrange->{-0.01,0.01},xrange->{0,80},fontsize-> 18,imagesize->600};
ModePlot[fit_,lm_,OptionsPattern[]]:=
	Module[{plotlistRe={},plotlistIm={},i,j,fpos,size,fitsize,psize,mode,modes={},
			pmodes={},a,\[Delta],\[Theta]angle,A,\[Phi],pRe,pIm,Reoptions,Imoptions,line,legend={},
			legendheight,color=3,ColorScheme,fitmodes,fitpmodes,PlotBuild,PBout,
			legendsize=0},
	a=fit[[5]];
	\[Delta]=fit[[6]];
	\[Theta]angle=fit[[7]];
	fpos=fit[[8]];

	modes=fit[[3]];
	fitmodes=fit[[9]];
	size=Length[modes];
	PlotBuild[fpos_,rawmodes_,fitmodes_]:=
	Module[{i,j,plotlistRe={},plotlistIm={},legend={},line,legmode,mode,A,\[Phi],size,
			fitsize},
		size=Length[rawmodes];
		fitsize=Length[fitmodes];

		For[i=1;j=1, i<=size, i++,
			If[Length[rawmodes[[i]]]==6 && rawmodes[[i,6]]==True,
				A=rawmodes[[i,4]];
				\[Phi]=rawmodes[[i,5]];
			];

			If[(Length[rawmodes[[i]]]==6 && rawmodes[[i,6]]==False) ||
				Length[rawmodes[[i]]]==5 || Length[rawmodes[[i]]]==3,
				A=ParFrom[fit,fpos+j];
				\[Phi]=ParFrom[fit,fpos+fitsize+j];
				j++;		
			];
		
			line=Graphics[{ColorData[color,"ColorList"][[legendsize+i]],
				Line[{{0,0},{2,0}}]}];
			legmode=" ("<>ToString[rawmodes[[i,1]]]<>","<>
					ToString[rawmodes[[i,2]]]<>","<>ToString[rawmodes[[i,3]]]<>")";
			AppendTo[legend,{line,legmode}];
			mode={rawmodes[[i,1]],rawmodes[[i,2]],rawmodes[[i,3]]}~Join~{A,\[Phi]};
			AppendTo[plotlistRe,Re[Clm[lm[[1]],lm[[2]],{mode},{},\[Delta],a,\[Theta]angle,t]]];
			AppendTo[plotlistIm,Im[Clm[lm[[1]],lm[[2]],{mode},{},\[Delta],a,\[Theta]angle,t]]];
		];
		{plotlistRe,plotlistIm,legend}
	];
	
	PBout = PlotBuild[fpos,modes,fitmodes];
	plotlistRe = PBout[[1]];
	plotlistIm = PBout[[2]];
	legend = PBout[[3]];
	legendsize = Length[legend];

	pmodes = fit[[4]];
	psize = Length[pmodes];
	fpos += 2*Length[fitmodes];
	fitpmodes = fit[[10]];

	PBout = PlotBuild[fpos,pmodes,fitpmodes];
	plotlistRe = Join[plotlistRe,PBout[[1]]];
	plotlistIm = Join[plotlistIm,PBout[[2]]];
	legend = Join[legend,PBout[[3]]];

	ColorScheme = If[(size+psize)==1,RGBColor[0,0,0], ColorData[color,"ColorList"]];

	Reoptions=Opts[MakeLabel["Re",lm[[1]],lm[[2]]],OptionValue[fontsize],
			OptionValue[imagesize]];
	Imoptions=Opts[MakeLabel["Im",lm[[1]],lm[[2]]],OptionValue[fontsize],
			OptionValue[imagesize]];
	pRe=Plot[plotlistRe,{t,OptionValue[xrange][[1]],OptionValue[xrange][[2]]},
		PlotRange->{OptionValue[xrange],OptionValue[yrange]},
		PlotRange->OptionValue[yrange],PlotStyle->ColorScheme];
	pIm=Plot[plotlistIm,{t,OptionValue[xrange][[1]],OptionValue[xrange][[2]]},
		PlotRange->{OptionValue[xrange],OptionValue[yrange]},
		PlotRange->OptionValue[yrange],PlotStyle->ColorScheme];

	legendheight=4;
	ShowLegend[Show[pRe,Reoptions],
	{legend,
	LegendShadow->None,LegendSize->0.22,LegendTextSpace->legendheight,
	LegendPosition->{0.67,0.56-0.012*legendheight*Length[legend]}}]

	ShowLegend[Show[pIm,Imoptions],
	{legend,
	LegendShadow->None,LegendSize->0.22,LegendTextSpace->legendheight,
	LegendPosition->{0.67,0.56-0.012*legendheight*Length[legend]}}]
]


Options[GetQNMValues]={Fixed->False};
GetQNMValues[fit_,OptionsPattern[]]:=
	Module[{i,modeSize=0,pmodeSize=0,modes,pmodes,modeStart,Massparam,Spinparam,
			Thetaparam},
	modeSize = Length[fit[[9]]];
	pmodeSize = Length[fit[[10]]];

	modeStart=fit[[8]];
	modes=Table[{fit[[9,i,1]],fit[[9,i,2]],fit[[9,i,3]],
		Abs[ParFrom[fit,modeStart+i]],ParFrom[fit,modeStart+modeSize+i],
		OptionValue[Fixed]},{i,1,modeSize}];

	pmodes=Table[{fit[[10,i,1]],fit[[10,i,2]],fit[[10,i,3]],
			Abs[ParFrom[fit,modeStart+2*modeSize+i]],
			ParFrom[fit,modeStart+2*modeSize+pmodeSize+i],
			OptionValue[Fixed]},{i,1,pmodeSize}];

	Massparam=Mass->{fit[[5]],OptionValue[Fixed]};
	Spinparam=Spin->{fit[[6]],OptionValue[Fixed]};
	Thetaparam=Theta->{fit[[7]],OptionValue[Fixed]};
	Print["********modes*********\n",modes,"\n********pmodes*********\n",pmodes,
		  "\n********parameters*****\n",{Massparam,Spinparam,Thetaparam}];
]


ProjectEllipsoid[ellipsoid_,subdim_]:=
	Module[{center,N,semiaxes,directions,projellipse},
	center=ellipsoid[[1]];
	semiaxes=ellipsoid[[2]];
	N=Length[center];
	directions=ellipsoid[[3]].Transpose[IdentityMatrix[N][[subdim]]];
	projellipse=
	Eigensystem[Transpose[directions].DiagonalMatrix[semiaxes^2].directions];
	{center[[subdim]],Sqrt[projellipse[[1]]],projellipse[[2]]}
]

Needs["MultivariateStatistics`"];

PEllipsoid[ellipsoid_,subdim_]:=Module[{projellipse},
	projellipse=ProjectEllipsoid[ellipsoid,subdim];
	Ellipsoid[projellipse[[1]],projellipse[[2]],projellipse[[3]]]
]

MassSpinConfRegion[fit_]:=Module[{color},
	color=RGBColor[1,0,0];
	Graphics[{color,Opacity[0.5], 
	PEllipsoid[fit[[1]]["ParameterConfidenceRegion"],{2,1}]},
	Axes->True,AxesLabel->{"a","\[Delta]"}]/.Hue[0.67,0.6,0.6]->{}/. Line->Polygon
]


Getaz[fit_]:=Module[{\[Delta],a,ac},
	\[Delta] = fit[[5]];
	a = fit[[6]];
	\[Theta] = fit[[7]];
	ac=\[Delta]^2*a;
	ac*Cos[\[Theta]]
]

Getatheta[\[Delta]_,acx_,acy_,acz_]:=Module[{a,ac,\[Theta]},
	ac = Sqrt[acx^2+acy^2+acz^2];
	a = ac/\[Delta]^2;
	\[Theta] = ArcCos[acz/ac];
	{a,\[Theta]}
]

End[]
EndPackage[]
