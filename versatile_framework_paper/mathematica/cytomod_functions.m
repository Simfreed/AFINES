(* ::Package:: *)

(* ::Subsection:: *)
(*Common Directories*)


outdir="/project/weare-dinner/simonfreedman/cytomod/out";


(* ::Subsection:: *)
(*Reading Dat Files*)


ClearAll[pts,pts3];
pts[dir_,name_]:=
If[
FileExistsQ[dir<>"/txt_stack/"<>name<>".dat"],
ToExpression[Import[dir<>"/txt_stack/"<>name<>".dat","TSV"]],
Print[dir<>"/txt_stack/"<>name<>".dat doesn't exist"];{}
];
pts3[dir_,name_]:=Map[Internal`StringToDouble/@(StringSplit[StringTake[#,{2,-2}],", "])&,Import[dir<>"/txt_stack/"<>name<>".dat","TSV"],{2}];
importCheck[fname_]:=
If[
FileExistsQ[fname],
ToExpression[Import[fname,"TSV"]],
Print[fname<>" doesn't exist"];{}
];
importCheck3[fname_]:=
If[
FileExistsQ[fname],
Map[Internal`StringToDouble/@(StringSplit[StringTake[#,{2,-2}],", "])&,ToExpression[Import[fname,"TSV"]],{2}];,
Print[fname<>" doesn't exist"];{}
];
(*cytosim*)
ptsCS[dir_,name_]:=
If[
FileExistsQ[dir<>"/"<>name<>".dat"],
ToExpression[Import[dir<>"/"<>name<>".dat","TSV"]],
Print[dir<>"/"<>name<>".dat doesn't exist"];{}
];


(* ::Subsection:: *)
(*Drawing Networks*)


ClearAll[link,actin,amotor,pmotor,draw];
link[l_,maxLength_:100,color_:Red,th_:0.007]:=Graphics[
{color,
Thickness[th],
If[Norm[l[[3;;4]]]<maxLength,
Line[{{l[[1]],l[[2]]},{l[[1]]+l[[3]],l[[2]]+l[[4]]}}],
Point[{{l[[1]],l[[2]]},{l[[1]]+l[[3]],l[[2]]+l[[4]]}}]
]
}
];
actin[a_,col_]:=Graphics[
{col,
Disk[{a[[1]],a[[2]]},a[[3]]]
}
];
radactin[a_,col_,rad_]:=Graphics[
{col,
Disk[{a[[1]],a[[2]]},rad]
}
];

amotor[l_,maxLength_:10,col_:Black,th_:0.02]:=Graphics[
{col,
Thickness[th],
If[Norm[l[[3;;4]]]<maxLength,
Line[{{l[[1]],l[[2]]},{l[[1]]+l[[3]],l[[2]]+l[[4]]}}],
Point[{{l[[1]],l[[2]]},{l[[1]]+l[[3]],l[[2]]+l[[4]]}}]
]
}
];
pmotor[l_,col_:Yellow,th_:0.02]:=Graphics[
{col,
Thickness[th],
Line[{{l[[1]],l[[2]]},{l[[1]]+l[[3]],l[[2]]+l[[4]]}}]
}
];

draw[actins_,links_,amotors_,pmotors_,dir_:"",movieFlag_:False,fov_: {50,50},polar_: True,np_: 11,acolor_: Red,movieName_: "movie",maxLinkLength_:-1,barbedEndFactor_: 3]:=Module[{},(*Generate Image Time Series,make into movie*)
mll=If[maxLinkLength==-1,Abs[Min[fov]/2],maxLinkLength];
frames={};
If[actins!={},AppendTo[frames,Map[actin[#,acolor]&,actins,{2}]]];
If[links!={},AppendTo[frames,Map[link[#,mll]&,links,{2}]]];
If[polar,
If[actins!={},AppendTo[frames,Map[actin[#,Blue]&,actins[[All,1;;-1;;np]],{2}]],If[links!={},rad=Mean[Flatten[Map[Norm,links[[All,All,3;;4]],{2}]]]/barbedEndFactor;
AppendTo[frames,Map[radactin[#,Lighter[Blue],rad]&,links[[All,1;;-1;;np-1]],{2}]]]]];
If[amotors!={},AppendTo[frames,Map[amotor,amotors,{2}]]];
If[pmotors!={},AppendTo[frames,Map[pmotor[#,Green]&,pmotors,{2}]]];
frames=Transpose[frames];
minx=0;miny=0;maxx=0;maxy=0;
If[TensorRank[fov]==1,minx=-fov[[1]]/2;maxx=fov[[1]]/2;miny=-fov[[2]]/2;
maxy=fov[[2]]/2,minx=fov[[1,1]];maxx=fov[[1,2]];miny=fov[[2,1]];
maxy=fov[[2,2]]];
frames=Map[Show[#,Frame->True,PlotRange->{{minx,maxx},{miny,maxy}},(*FrameLabel\[Rule]{"x (\[Mu]m)","y (\[Mu]m)","t = "<>ToString[dt*skipFrames*(t-1)]<>" s "},*)FrameTicks->None,Ticks->None,BaseStyle->{FontSize->24,FontColor->Black}]&,frames];
If[movieFlag,Export[dir<>"/data/"<>movieName<>".avi",frames]];
frames];


ClearAll[drawblack];
drawblack[actins_,links_,amotors_,pmotors_,dir_:"",movieFlag_:False,fov_: {50,50},polar_: True,np_: 11,acolor_: Red,movieName_: "movie",maxLinkLength_:-1,barbedEndFactor_: 3,lthick_:0.007,mthick_:0.02,pthick_:0.02]:=Module[{},(*Generate Image Time Series,make into movie*)
mll=If[maxLinkLength==-1,Abs[Min[fov]/2],maxLinkLength];
frames={};
If[actins!={},AppendTo[frames,Map[actin[#,acolor]&,actins,{2}]]];
If[links!={},AppendTo[frames,Map[link[#,mll,acolor,lthick]&,links,{2}]]];
If[polar,
If[actins!={},AppendTo[frames,Map[actin[#,Cyan]&,actins[[All,1;;-1;;np]],{2}]],If[links!={},rad=Mean[Flatten[Map[Norm,links[[All,All,3;;4]],{2}]]]/barbedEndFactor;
AppendTo[frames,Map[radactin[#,Cyan,rad]&,links[[All,1;;-1;;np-1]],{2}]]]]];
If[amotors!={},AppendTo[frames,Map[amotor[#,4,White,mthick]&,amotors,{2}]]];
If[pmotors!={},AppendTo[frames,Map[pmotor[#,Green,pthick]&,pmotors,{2}]]];
frames=Transpose[frames];
minx=0;miny=0;maxx=0;maxy=0;
If[TensorRank[fov]==1,minx=-fov[[1]]/2;maxx=fov[[1]]/2;miny=-fov[[2]]/2;
maxy=fov[[2]]/2,minx=fov[[1,1]];maxx=fov[[1,2]];miny=fov[[2,1]];
maxy=fov[[2,2]]];
frames=Map[Show[#,Frame->True,PlotRange->{{minx,maxx},{miny,maxy}},(*FrameLabel\[Rule]{"x (\[Mu]m)","y (\[Mu]m)","t = "<>ToString[dt*skipFrames*(t-1)]<>" s "},*)FrameTicks->None,Ticks->None,BaseStyle->{FontSize->24,FontColor->Black},Background->Black]&,frames];
If[movieFlag,Export[dir<>"/data/"<>movieName<>".avi",frames]];
frames];
drawStretch[actins_,links_,amotors_,pmotors_,dir_,movieFlag_,fov_:{500,500},dt_:1,acolor_:Red,cell_:{0,0},strainPct_:0,l0_:1,
colorByEnergy_:True,enghilo_:{0,0},scheme_:"BlueGreenYellow",bgcol_:White]:=
Module[
{actinsr=actins, linksr=links,amotorsr=amotors,pmotorsr=pmotors},
(*Generate Image Time Series, make into movie*)
\[Delta]rx=strainPct*fov[[1]]*cell[[2]];
cm=fov*cell;
minx=fov[[1]](cell[[1]]-0.5);
maxx=fov[[1]](cell[[1]]+0.5);
miny=fov[[2]](cell[[2]]-0.5);
maxy=fov[[2]](cell[[2]]+0.5);
frames={};
If[colorByEnergy,
If[enghilo[[1]]==enghilo[[2]],
frcs=(Norm/@Flatten[links[[All,All,3;;4]],1]-l0);
cf=ColorData[{scheme,Quantile[frcs,{0.025,0.975}]}],
cf=ColorData[{scheme,{enghilo[[2]],enghilo[[1]]}}]
],
cf[i_]:="Red"
];
If[actinsr!={},
actinsr[[All,All,1;;2]]=Map[ (#+cm+{\[Delta]rx,0})&,actinsr[[All,All,1;;2]],{2}];AppendTo[frames,Map[actin[#,Blend[{acolor,Red},Part[#,5]/10.]]&,actinsr,{2}]]];
If[linksr!={},
linksr[[All,All,1;;2]]=Map[ (#+cm+{\[Delta]rx,0})&,linksr[[All,All,1;;2]],{2}];
AppendTo[frames,Map[link[#,Abs[0.5fov[[1]]],cf[Abs[Norm[#[[3;;4]]]-l0]]]&,linksr,{2}]]
];
If[amotorsr!={},
amotorsr[[All,All,1;;2]]=Map[(#+cm+{\[Delta]rx,0})&,amotorsr[[All,All,1;;2]],{2}];
AppendTo[frames,Map[amotor,amotorsr,{2}]]];
If[pmotorsr!={},
pmotorsr[[All,All,1;;2]]=Map[(#+cm+{\[Delta]rx,0})&,pmotorsr[[All,All,1;;2]],{2}];AppendTo[frames,Map[pmotor,pmotorsr,{2}]]];
frames=Transpose[frames];
frames=Map[
Show[#,
Frame->True,
PlotRange->{{minx,maxx},{miny,maxy}},
Background->bgcol,
(*FrameLabel\[Rule]{"x (\[Mu]m)","y (\[Mu]m)","t = "<>ToString[dt*skipFrames*(t-1)]<>" s "},*)
FrameTicks->None,
Ticks->None,
BaseStyle->{FontSize->24,FontColor->Black}
]&
,frames
];
If[movieFlag, Export[dir<>"/data/movie_"<>ToString[scheme]<>".avi",frames]];
frames
];


ClearAll[drawbyV];
drawbyV[actins_,links_,amotors_,pmotors_,vs_,vlohi_:{0,0},scheme_:"BlueGreenYellow",dir_:"",movieFlag_:False,fov_:{50,50},polar_: False,np_: 11,acolor_:Red,movieName_: "movie_vcol",maxLinkLength_:-1,barbedEndFactor_: 3]:=
Module[
{},
(*Generate Image Time Series, make into movie*)
cf=ColorData[{scheme,{vlohi[[1]],vlohi[[2]]}}];
If[vlohi[[1]]==vlohi[[2]],cf=ColorData[{scheme,Quantile[Flatten[vs],{0.025,0.975}]}]];
mll=If[maxLinkLength==-1,Abs[Min[fov]/2],maxLinkLength];
frames={};

If[actins!={},AppendTo[frames,Map[actin[#,acolor]&,actins,{2}]]];
If[links!={},AppendTo[frames,MapThread[link[#1,mll,cf[#2]]&,{links,vs},2]]];
If[polar,
If[actins!={},AppendTo[frames,Map[actin[#,Blue]&,actins[[All,1;;-1;;np]],{2}]],If[links!={},rad=Mean[Flatten[Map[Norm,links[[All,All,3;;4]],{2}]]]/barbedEndFactor;
AppendTo[frames,Map[radactin[#,Lighter[Blue],rad]&,links[[All,1;;-1;;np-1]],{2}]]]]];
If[amotors!={},AppendTo[frames,Map[amotor,amotors,{2}]]];
If[pmotors!={},AppendTo[frames,Map[pmotor[#,Green]&,pmotors,{2}]]];
frames=Transpose[frames];
minx=0;miny=0;maxx=0;maxy=0;
If[TensorRank[fov]==1,minx=-fov[[1]]/2;maxx=fov[[1]]/2;miny=-fov[[2]]/2;
maxy=fov[[2]]/2,minx=fov[[1,1]];maxx=fov[[1,2]];miny=fov[[2,1]];
maxy=fov[[2,2]]];
frames=Map[Show[#,Frame->True,PlotRange->{{minx,maxx},{miny,maxy}},FrameTicks->None,Ticks->None,BaseStyle->{FontSize->24,FontColor->Black}]&,frames];
If[movieFlag,Export[dir<>"/data/"<>movieName<>".avi",frames]];
frames
];


(* ::Subsection:: *)
(*Plotting Options*)


Needs["CustomTicks`"] (*Get this package, it's great!*)


SetOptions[{Plot,LogLogPlot,ListPlot,ListLogPlot,ListLogLogPlot,ListLogLinearPlot},BaseStyle->{FontFamily->"Arial",FontSize->22,SingleLetterItalics->False},Frame->True,AspectRatio->0.8,FrameStyle->Black];
SetOptions[{ListPlot,ListLogPlot,ListLogLogPlot,ListLogLinearPlot},Joined->True,PlotMarkers->Automatic,PlotStyle->"BlueGreenYellow"];
SetOptions[{ListVectorPlot,ListVectorDensityPlot},BaseStyle->{FontSize->26,FontFamily->"Arial"},
FrameStyle->Black,VectorStyle->Black];
SetOptions[{ListDensityPlot},BaseStyle->{FontSize->26,FontFamily->"Arial"},
FrameStyle->Black];
ClearAll[logbarleg];
logarithmicscaling[x_,min_,max_]:=Log[x/min]/Log[max/min];
logbarleg[min_,max_,nticks_,cf_:ColorData[{"BlueGreenYellow",{0,1}}],acc_:{4,1},myopts_:{}]:=
BarLegend[
cf,
Ticks->({logarithmicscaling[#,min,max],NumberForm[#,acc]}&/@(min (max/min)^Range[0,1,1/nticks])),
myopts];


(* ::Subsection:: *)
(*Periodic Boundary Conditions*)


rijP[dr_,fov_]:=dr-Round[dr/fov]*fov;


cmP1D[rs_,fov_]:=Module[{},
\[Theta]s=2Pi (rs+fov/2)/fov;
fov ((ArcTan[-Mean[Cos[\[Theta]s]],-Mean[Sin[\[Theta]s]]])/(2Pi))
];
cmP[rs_,fov_]:=MapThread[cmP1D,{rs\[Transpose],fov}];


(* ::Subsection:: *)
(*MSD*)


ClearAll[msdp,msdm,msdc]
msdc:=Compile[{{drs,_Real,2}},
sums=Table[Accumulate[drs[[t;;-1]]],{t,1,Round[Length[drs]]}];
dx2s=sums[[All,All,1]]^2+sums[[All,All,2]]^2;
msd=Table[Mean[dx2s[[1;;-t,t]]],{t,1,Length[dx2s]}];
ClearAll[sums,dx2s];
msd,
{{t,_Integer}}
];
msdm[drs_]:=Module[{},
sums=Table[Accumulate[drs[[t;;]]],{t,Length[drs]}];
dx2s=sums[[All,All,1]]^2+sums[[All,All,2]]^2;
msd=Table[Mean[dx2s[[1;;-t,t]]],{t,1,Length[dx2s]}];
ClearAll[sums,dx2s];
msd
];
(*Note: requires unrolled trajectory; i.e., not periodic boundary conditions?*)
msdD[pos_,\[CapitalDelta]_]:=Module[{},
norms=Table[Norm[pos[[t+\[CapitalDelta]]]-pos[[t]]],{t,1,Length[pos]-\[CapitalDelta]}];
tots=Accumulate[norms];
tfs=Range[1,Length[pos]-\[CapitalDelta]];
tots/tfs
(*{tfs,tots}\[Transpose]*)
];


(* ::Subsection:: *)
(*Radial Distribution Function*)


getSamp[sampSize_:-1,np_:500,nm_:11]:=Module[{},
(*Get random sample of pairs of actins on different filaments*)
nparts=nm*np;
sz=sampSize;
If[sampSize==-1,sz=np (np-1)/2];
rngs=Table[Range[nm(i+1)+1,nparts],{i,0,nparts/nm-1}];
RandomSample[
Flatten[
Table[{i,j},
{i,nparts},
{j,rngs[[Quotient[i-1,nm]+1]]}],1],
sz
]
];

rdf[acts_,fov_,\[Delta]r_,maxr_,sampPairs_]:=Module[{},
(*Calculate their rdf*)
rs=Range[0,maxr,\[Delta]r];
rsmid=(rs[[2;;]]+rs[[;;-2]])0.5;
norms=2Pi \[Delta]r rsmid Length[sampPairs]/(fov[[1]]*fov[[2]]);
bcs=BinCounts[
Map[
Norm[rijP[#,fov]]&,
acts[[sampPairs[[All,1]],1;;2]]-acts[[sampPairs[[All,2]],1;;2]]],
{0,maxr,\[Delta]r}]/norms;
{rsmid,bcs}\[Transpose]
];
rdftraj[acts_,fov_,\[Delta]r_,maxr_,sampPairs_]:=Module[{},
rs=Range[0,maxr,\[Delta]r];
rsmid=(rs[[2;;]]+rs[[;;-2]])0.5;
norms=2Pi \[Delta]r rsmid Length[sampPairs]/(fov[[1]]*fov[[2]]);
rdft=
Table[
BinCounts[
Map[
Norm[rijP[#,fov]]&,
acts[[t,sampPairs[[All,1]],1;;2]]-acts[[t,sampPairs[[All,2]],1;;2]]],
{0,maxr,\[Delta]r}]/norms,
{t,Length[acts]}
];
rdft
];
rdfResamp[acts_,fov_,\[Delta]r_,maxr_,sampSize_:-1,np_:500,nm_:11]:=Module[{},
sampPairs=getSamp[sampSize,np,nm];
(*Calculate their rdf*)
rs=Range[0,maxr,\[Delta]r];
rsmid=(rs[[2;;]]+rs[[;;-2]])0.5;
norms=2Pi \[Delta]r rsmid Length[sampPairs]/(fov[[1]]*fov[[2]]);
bcs=BinCounts[
Map[
Norm[rijP[#,fov]]&,
acts[[sampPairs[[All,1]],1;;2]]-acts[[sampPairs[[All,2]],1;;2]]],
{0,maxr,\[Delta]r}]/norms;
{rsmid,bcs}\[Transpose]
];


(* ::Subsection:: *)
(*Local nematic order*)


nemt[lks_,lscale_,fov_:{50,50}]:=Module[{},
xcm=Round[rijP[lks[[All,All,1]]+0.5lks[[All,All,3]],fov[[1]]],lscale];
ycm=Round[rijP[lks[[All,All,2]]+0.5lks[[All,All,4]],fov[[2]]],lscale];
xydxdy=Table[{xcm[[t]],ycm[[t]],lks[[t,All,3]],lks[[t,All,4]]}\[Transpose],{t,Length[lks]}];
xydxdyg=Table[GatherBy[xydxdy[[t]],#[[1;;2]]&],{t,Length[xydxdy]}];
xydxdyg2=Table[Cases[xydxdyg[[t]],_?((Length[#]>1)&)],{t,1,Length[xydxdyg]}];
directors=Table[Mean[xydxdyg2[[t,gp,All,3;;4]]],{t,1,Length[xydxdyg2]},{gp,1,Length[xydxdyg2[[t]]]}];
costht=Flatten/@Table[
xydxdyg2[[t,gp,All,3;;4]].directors[[t,gp]]/(Norm/@xydxdyg2[[t,gp,All,3;;4]]*Norm[directors[[t,gp]]])
,
{t,Length[directors]},
{gp,Length[directors[[t]]]}
];
s2=2Mean/@(costht^2)-1;
dfield=Table[{xydxdyg2[[t,gp,1,1;;2]],directors[[t,gp]]},{t,Length[directors]},{gp,Length[directors[[t]]]}];
Clear[xcm,ycm,xydxdy,xydxdyg,costht];
{s2,dfield}];


ClearAll[nem];
nem[lks_,lscale_,fov_:{50,50}]:=Module[{},
xcm=Round[rijP[lks[[All,1]]+0.5lks[[All,3]],fov[[1]]],lscale];
ycm=Round[rijP[lks[[All,2]]+0.5lks[[All,4]],fov[[2]]],lscale];
xydxdy={xcm,ycm,lks[[All,3]],lks[[All,4]]}\[Transpose];
xydxdyg=GatherBy[xydxdy,#[[1;;2]]&];
xydxdyg2=Cases[xydxdyg,_?((Length[#]>1)&)];
directors=Table[Mean[xydxdyg2[[gp,All,3;;4]]],{gp,Length[xydxdyg2]}];
costh=Flatten[Table[
xydxdyg2[[gp,All,3;;4]].directors[[gp]]/(Norm/@xydxdyg2[[gp,All,3;;4]]*Norm[directors[[gp]]])
,
{gp,Length[directors]}
]
];
s2=2Mean[costh^2]-1;
dfield=Table[{xydxdyg2[[gp,1,1;;2]],directors[[gp]]},{gp,Length[directors]}];
Clear[xcm,ycm,xydxdy,xydxdyg,costh];
{s2,dfield}];


nem[lks_,lscale_,fov_:{50,50},maxl_:50]:=Module[{},
xcm=Round[rijP[lks[[All,1]]+0.5lks[[All,3]],fov[[1]]],lscale];
ycm=Round[rijP[lks[[All,2]]+0.5lks[[All,4]],fov[[2]]],lscale];
xydxdy={xcm,ycm,lks[[All,3]],lks[[All,4]]}\[Transpose];
xydxdyg=BinLists[xydxdy,{-fov[[1]]/2.0,fov[[1]]/2.0,lscale},{-fov[[2]]/2.0,fov[[2]]/2.0,lscale},{-maxl,maxl,2maxl},{-maxl,maxl,2maxl}][[All,All,1,1]];
xydxdyg2=Cases[Flatten[xydxdyg,1],_?((Length[#]>1)&)];
directors=Table[Mean[xydxdyg2[[gp,All,3;;4]]],{gp,Length[xydxdyg2]}];
costh=Flatten[Table[
xydxdyg2[[gp,All,3;;4]].directors[[gp]]/(Norm/@xydxdyg2[[gp,All,3;;4]]*Norm[directors[[gp]]])
,
{gp,Length[directors]}
]
];
s2=2Mean[costh^2]-1;
dfield=Table[{xydxdyg2[[gp,1,1;;2]],directors[[gp]]},{gp,Length[directors]}];
Clear[xcm,ycm,xydxdy,xydxdyg,costh];
{s2,dfield}];


(* ::Subsection:: *)
(*Filament Strain*)


filstrain[acts_,fov_:{50,50},nparts_:500,arcLength_:10]:=
    1-Map[
        Norm[
        rijP[#,fov]]&,
    acts[[All,Length[acts[[1]]]/nparts;;-1;;Length[acts[[1]]]/nparts,1;;2]]-acts[[All,1;;-1;;Length[acts[[1]]]/nparts,
    1;;2]],
    {2}]/arcLength;


ClearAll[divdens];
divdens[acts_,xyds_,fov_:{50,50},dr_:5,dt_:10,minDiv_:-1000,maxDiv_:1000]:=
Module[{},
bindenss=Table[BinCounts[acts[[t,All,1;;2]],{-fov[[1]]/2,fov[[1]]/2,dr},{-fov[[2]]/2,fov[[2]]/2,dr}]/dr^2,{t,Length[acts]}];
bindenss=MovingAverage[bindenss,dt];
bindenss=0.5(bindenss[[2;;]]+bindenss[[;;-2]]);divlists=Table[BinLists[xyds[[t,All,{1,2,5}]],{-fov[[1]]/2,fov[[1]]/2,dr},{-fov[[2]]/2,fov[[2]]/2,dr},{minDiv,maxDiv,maxDiv-minDiv}],{t,Length[xyds]}];meandivs=Map[Mean,divlists[[All,All,All,1,All,3]],{3}];
maxT=Min[Length[bindenss],Length[meandivs]];
densdivt=dr^2*Total[Flatten[#]]&/@(meandivs[[1;;maxT]]*bindenss[[1;;maxT]]);
ClearAll[bindenss,divlists,meandivs,maxT];
densdivt
];


(* ::Subsection:: *)
(*Density Weighted Divergence*)


divdens2[acts_,divfs_,fov_:{50,50},dr_:5,dt_:10,minDiv_:-1000,maxDiv_:1000]:=
Module[{},

bindenss=Table[BinCounts[acts[[t,All,1;;2]],{-fov[[1]]/2,fov[[1]]/2,dr},{-fov[[2]]/2,fov[[2]]/2,dr}]/dr^2,{t,Length[acts]}];
bindenss=MovingAverage[bindenss,dt];
bindenss=0.5(bindenss[[2;;]]+bindenss[[;;-2]]);

edgex=Range[-fov[[1]]/2,fov[[1]]/2,dr];
edgey=Range[-fov[[2]]/2,fov[[2]]/2,dr];
rngx=0.5(edgex[[2;;]]+edgex[[;;-2]]);
rngy=0.5(edgey[[2;;]]+edgey[[;;-2]]);
bincenters=Table[{x,y},{x,rngx},{y,rngy}];

divsxy=Table[Map[divfs[[t]],bincenters,{2}],{t,Length[divfs]}];
maxT=Min[Length[bindenss],Length[divfs]];

rhodivs=Table[Total[bindenss[[t]]*divsxy[[t]],2],{t,maxT}];
ClearAll[bindenss,divsxy,maxT];
rhodivs

];


(* ::Subsection:: *)
(*Velocity Autocorrelation*)


normdot=Function[{v1,v2},v1.v2/(Norm[v1]Norm[v2])];
ClearAll[vccf];
vccf[acts_,dr_:0.05,fov_:{50,50},dt_:1]:=Module[{},
vs=Map[rijP[#,fov]&,acts[[(dt+1);;,All,1;;2]]-acts[[;;(-dt-1),All,1;;2]],{2}]/dt;
sampPairs=Table[getSamp[],{Length[vs]}];
pairdists=Table[
Norm[rijP[acts[[t,sampPairs[[t,i,1]],1;;2]]-acts[[t,sampPairs[[t,i,2]],1;;2]],fov]],
{t,Length[vs]},
{i,Length[sampPairs[[t]]]}
];
pairdots=Table[normdot[vs[[t,sampPairs[[t,i,1]]]],vs[[t,sampPairs[[t,i,2]]]]],{t,Length[vs]},{i,Length[sampPairs[[t]]]}];
binmeans=Sort[Mean/@GatherBy[{Flatten[pairdists],Flatten[pairdots]}\[Transpose],Round[#[[1]],dr]&]];
(*ClearAll[vs,sampPairs,pairdists,pairdots];*)
binmeans];


(* ::Subsection:: *)
(*Mesh Size*)


mysubdivide[p1_,p2_,n_,fov_]:=Module[{},
dp=p2-p1;
If[dp==rijP[p2-p1,fov],
coords=Subdivide[p1,p2,n],
(*Else, if the problem is in the x coord*)
ps=Sort[{p1,p2}]+{{fov[[1]],0},{0,0}};
If[ps[[2]]-ps[[1]]!=rijP[ps[[2]]-ps[[1]],fov],
(*Else, if the problem is in the y coord*)
ps=Sort[{p1,p2},#1[[2]]<#2[[2]]&]+{{0,fov[[2]]},{0,0}};
If[ps[[2]]-ps[[1]]!=rijP[ps[[2]]-ps[[1]],fov],
(*Else, the problem is in both coords*)
psx=Sort[{p1,p2}];
psy=Sort[{p1,p2},#1[[2]]<#2[[2]]&];
(*If the lower x is also the lower y*)
ps=psx+{fov,{0,0}};
If[psx!=psy,
(*Else:, the lower x is the higher y*)
ps=psx+{fov{1,-1},{0,1}}];
]
];
coords=Map[rijP[#,fov]&,Subdivide[ps[[1]],ps[[2]],n]]
];
coords
];


(* ::Subsection:: *)
(*Persistence Length*)


ClearAll[anglesT, anglesBwT, ccf,th2];
anglesT[rt_]:=ArcTan[rt[[All,3]],rt[[All,4]]];
anglesTCS[rt_]:=ArcTan[rt[[All,1]],rt[[All,2]]];
anglesBwT[angs_]:=Map[(#-2Pi Floor[#/(2Pi)+0.5])&,angs[[2;;]]-angs[[1;;-2]]];
ccf[ths_]:=Module[{},
sums=Table[Accumulate[ths[[i;;]]],{i,Length[ths]}];
Table[Mean[Cos[sums[[1;;Length[sums]-i,i+1]]]],{i,0,Length[sums]-1}]
];
th2[ths_]:=Module[{},
sums=Table[Accumulate[ths[[i;;]]],{i,Length[ths]}];
Table[Mean[sums[[1;;Length[sums]-i,i+1]]^2],{i,0,Length[sums]-1}]
];
plMeanThsq[drs_,dx_,L_,maxL_]:=Module[{},
angs=ArcTan[drs[[All,All,1]],drs[[All,All,2]]];
angsbw=anglesBwT/@angs;
th2s=th2/@angsbw;
muths=Mean[th2s];
dat={Range[dx,L-dx,dx],muths}\[Transpose];
fit=LinearModelFit[dat[[1;;Min[maxl,Length[dat]]]],x,x];
lp=1/D[Normal[fit],x];
ClearAll[angs,angsbw,th2s,muths,dat,fit];
lp
];
plThsq[drs_,maxX_]:=Module[{},
angs=ArcTan[drs[[All,1]],drs[[All,2]]];
angsbwt=anglesBwT[angs];
th2s=th2[angsbwt];
dx=Mean[Norm/@drs];
max=Min[maxX, Length[th2s]+1];
dat={Range[dx,max dx-dx,dx],th2s[[1;;(max-1)]]}\[Transpose];
fit=LinearModelFit[dat,x,x];
lp=1/D[Normal[fit],x];
ClearAll[angs,angsbwt,th2s,dx,dat,fit];
lp
]



(* ::Subsection:: *)
(*Gaussian Radial Basis Function Interpolation From Weights*)


ClearAll[velx,vely,div];
velx[x_,y_,vf_]:=vf[[All,3]].Exp[(-1./vf[[All,5]]^2)*((Total[#^2]&)/@Transpose[{x,y}-{vf[[All,1]],vf[[All,2]]}])];
vely[x_,y_,vf_]:=vf[[All,4]].Exp[(-1./vf[[All,5]]^2)*((Total[#^2]&)/@Transpose[{x,y}-{vf[[All,1]],vf[[All,2]]}])];
div[x_,y_,vf_]:=(-2/vf[[All,5]]^2(vf[[All,3]](x-vf[[All,1]])+vf[[All,4]](y-vf[[All,2]]))).Exp[(-1./vf[[All,5]]^2)*((Total[#^2]&)/@Transpose[{x,y}-{vf[[All,1]],vf[[All,2]]}])];
