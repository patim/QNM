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
End[]
EndPackage[]
