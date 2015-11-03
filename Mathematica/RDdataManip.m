(* ::Package:: *)

(*
Data manipulation for QNMfit package input 
Maxim Zalutskiy
Gregory Cook
2013
*)
BeginPackage["RDdataManip`"]
TimeShift::usage="";
NormPsi::usage="";
DataAfter::usage="";
DataBefore::usage="";
ReadNRWaveForm::usage="";
GetData::usage="";
DataCut::usage="";
MergeCLists::usage="";
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
Begin["Private`"]

(* Utility funcitons *)
TimeShift[t_,A_] := Table[{A[[i,1]]-t,A[[i,2]],A[[i,3]]},{i,1,Length[A]}]

NormPsi[data_] := Module[{i,size,isize,norm={}},
	size=Length[data];
	For[i=1,i<=size,i++,
		isize=Length[data[[i,1]]];
		AppendTo[norm,Table[{data[[i,1,j,1]],Sqrt[data[[i,1,j,2]]^2+data[[i,1,j,3]]^2]},{j,1,isize}]];
	];
norm
]

DataAfter[t_,A_,step_:1] := Module[{i,istart},
	istart = Module[{i},
		For[i=1,i<=Length[A],i++,
			If[A[[i,1]] > t,Break[]]
		];
		i];
	Table[{A[[i,1]],A[[i,2]],A[[i,3]]},{i,istart,Length[A],step}]
	]

DataBefore[t_,A_,step_:1]:=Module[{i,iend},
	iend = Module[{i},
		For[i=Length[A],i>=1,i--,
			If[A[[i,1]] < t,Break[]]
		];
		i];
	Table[{A[[i,1]],A[[i,2]],A[[i,3]]},{i,1,iend,step}]
]

ReadNRWaveForm[N_,AnnexDir_,l_,m_]:=Module[{mname,lname,Yname,Gname},
	mname=If[m<0,"-"<>ToString[Abs[m]],ToString[m]];
	lname=ToString[l];
	Yname="Y_l"<>lname<>"_m"<>mname<>".dat";
	Gname=If[N==0,"/OutermostExtraction.dir/",If[N==2,"/Extrapolated_N2.dir/",If[N==3,"/Extrapolated_N3.dir/",If[N==4,"/Extrapolated_N4.dir/"]]]];
	Import[AnnexDir<>"rMPsi4_Asymptotic_GeometricUnits.h5",{"HDF5","Datasets",{Gname<>Yname}}]
]

GetData[path_,lm_,t1_,t2_,step_:1]:=Module[{Y,size,i,l,m,data={}},
	size=Length[lm];
	For[i=1,i<=size,i++,
		l=lm[[i,1]];
		m=lm[[i,2]];
		Y=DataCut[t1,t2,ReadNRWaveForm[0,path,l,m],step];
		If[i==1,
			Print["\!\(\*SubscriptBox[\(t\), \(1\)]\)=",Y[[1,1]],", \!\(\*SubscriptBox[\(t\), \(2\)]\)=",Y[[-1,1]]];
		];
		Y=TimeShift[Y[[1,1]],Y];
		AppendTo[data,{Y,{l,m}}];
	];
	data
]

DataCut[t1_,t2_,A_,step_:1]:=Module[{Atemp},
	If[t1>t2,Print["DataCut: t1 cannot be greater than t2"];Return[-1]];
	Atemp=DataAfter[t1,A,step];
	DataBefore[t2,Atemp,step]
]

(*The intput should be in the following form:
Clist={{{{Subscript[Re[Subscript[C, lm]], Subscript[list, 1]]}, {Subscript[Im[Subscript[C, lm]], Subscript[list, 1]]}, {Subscript[l, 1],Subscript[m, 1]}}, {{{Subscript[Re[Subscript[C, lm]], Subscript[list, 2]]}, {Subscript[Im[Subscript[C, lm]], Subscript[list, 2]]}, {Subscript[l, 2],Subscript[m, 2]}},...}
The output:
newClist = {{{{Subscript[t, 1],Subscript[Re[Subscript[C, lm]], 1],Subscript[Im[Subscript[C, lm]], 1]},{Subscript[t, 2],Subscript[Re[Subscript[C, lm]], 2],Subscript[Im[Subscript[C, lm]], 2]},...},{Subscript[l, 1],Subscript[m, 1]}},
{{{Subscript[t, 1],Subscript[Re[Subscript[C, lm]], 1],Subscript[Im[Subscript[C, lm]], 1]},{Subscript[t, 2],Subscript[Re[Subscript[C, lm]], 2],Subscript[Im[Subscript[C, lm]], 2]},...},{Subscript[l, 1],Subscript[m, 1]}},...}
*)
MergeCLists[Clist_]:=Module[{i,j,Clistsize,Csize,newClist={},lmlist},
	Clistsize=Length[Clist];
	For[i=1,i<= Clistsize,i++,
		Csize=Length[Clist[[i,1]]];
		lmlist=Clist[[i,3]];
		newClist=newClist~Join~{{Table[{Clist[[i,1,j,1]],Clist[[i,1,j,2]],Clist[[i,2,j,2]] },{j,1,Csize}],lmlist}};
	];
	newClist
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


End[]
EndPackage[]
