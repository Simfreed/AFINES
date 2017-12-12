(* ::Package:: *)

#!/software/mathematica-10.2-x86_ 64/bin/MathematicaScript -script


dir=$ScriptCommandLine[[2]];
Import["~/Code/cytomod/analysis/cytomod_functions.m"];
fov=If[Length[$ScriptCommandLine]>2,ToExpression/@{$ScriptCommandLine[[3]],$ScriptCommandLine[[4]]},{50,50}];
dr=If[Length[$ScriptCommandLine]>4,ToExpression[$ScriptCommandLine[[5]]],1];
dt=If[Length[$ScriptCommandLine]>5,ToExpression[$ScriptCommandLine[[6]]],10];
tf=If[Length[$ScriptCommandLine]>6,ToExpression[$ScriptCommandLine[[7]]],200];
nbins=If[Length[$ScriptCommandLine]>7,ToExpression[$ScriptCommandLine[[8]]],10];
thresh=If[Length[$ScriptCommandLine]>8,ToExpression[$ScriptCommandLine[[9]]],10];


ClearAll[velx,vely,div];
velx[x_,y_,vf_]:=vf[[All,3]].Exp[(-1./vf[[All,5]]^2)*((Total[#^2]&)/@Transpose[{x,y}-{vf[[All,1]],vf[[All,2]]}])];
vely[x_,y_,vf_]:=vf[[All,4]].Exp[(-1./vf[[All,5]]^2)*((Total[#^2]&)/@Transpose[{x,y}-{vf[[All,1]],vf[[All,2]]}])];
div[x_,y_,vf_]:=(-2/vf[[All,5]]^2(vf[[All,3]](x-vf[[All,1]])+vf[[All,4]](y-vf[[All,2]]))).Exp[(-1./vf[[All,5]]^2)*((Total[#^2]&)/@Transpose[{x,y}-{vf[[All,1]],vf[[All,2]]}])];


(*test
dir="/project/weare-dinner/simonfreedman/cytomod/out/vf_r1/contracting/xldens1/seed10000000";*)


(*Spatial Discretization*)
edgex=Range[-fov[[1]]/2,fov[[1]]/2,dr];
edgey=Range[-fov[[2]]/2,fov[[2]]/2,dr];
dx=dr;dy=dr;
xbins=0.5(edgex[[2;;]]+edgex[[;;-2]]);
ybins=0.5(edgey[[2;;]]+edgey[[;;-2]]);
bins=Table[{x,y},{x,xbins},{y,ybins}];


(*divergence at discretization*)
divfile=dir<>"/analysis/div_dr"<>ToString[dr]<>"_nbins"<>ToString[nbins]<>"_thresh"<>ToString[thresh]<>"_dt"<>ToString[dt]<>".dat";
Print["Importing Divergence"];
divts=ToExpression[Import[divfile,"TSV"]];
If[Dimensions[divts]!={tf-dt,Length[xbins],Length[ybins]},
Print["(Re)calculating Divergence at discretization of "<>ToString[dr]<>" because dim = "<>ToString@Dimensions[divts]];
wdir=
StringJoin[ToString/@{dir,"/analysis/vfield_weights_nbins",nbins,"_dt",dt,"_thresh",thresh}];
Print[wdir];
vfs=Table[Import[wdir<>"/t"<>ToString[t]<>".txt","TABLE"],{t,Range[0,tf-dt-1]}];
mydivs=Table[With[{vf=vfs[[t]]},div[#[[1]],#[[2]],vf]&],{t,Length[vfs]}];divts=Table[Map[mydivs[[t]],bins,{2}],{t,Length[mydivs]}];
Export[divfile,divts];
]
Print["calculated divergence"];


(*bincount at same discretization*)
bcfile =dir<>"/analysis/bin_counts_dr"<>ToString[dr]<>".dat";
Print["Importing Bin Counts"];
bcs=ToExpression[Import[bcfile,"TSV"]];
If[Dimensions[bcs]!={tf,Length[xbins],Length[ybins]},
Print["(Re)calculating bin counts at discretization of "<>ToString[dr]<>" because dim = "<>ToString@Dimensions[bcs]];
acts=pts[dir,"actins"];
acts=acts[[1;;Min[tf,Length[acts]]]];
bcs=Table[BinCounts[acts[[t,All,1;;2]],{-fov[[1]]/2,fov[[1]]/2,dr},{-fov[[2]]/2,fov[[2]]/2,dr}],{t,Length[acts]}];
Export[bcfile,bcs];
];
Print["calculated bincounts"];


(* density divergence at this discretization *)
bcsdt=MovingAverage[bcs,dt];
bcsdt=0.5(bcsdt[[2;;]]+bcsdt[[;;-2]]);
rhodivs=Table[Total[bcsdt[[t]]*divts[[t]],2],{t,Length[divts]}];
Export[dir<>"/analysis/dens_div_dr"<>ToString[dr]<>"_nbins"<>ToString[nbins]<>"_thresh"<>ToString[thresh]<>"_dt"<>ToString[dt]<>".dat",rhodivs];
Print["calculated densdiv"];
