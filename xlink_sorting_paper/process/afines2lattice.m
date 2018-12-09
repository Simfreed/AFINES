(* ::Package:: *)

#!/software/mathematica-10.2-x86_ 64/bin/MathematicaScript -script


dir=$ScriptCommandLine[[2]];
Import["~/Code/cytomod/analysis/cytomod_functions.m"];
ti=If[Length[$ScriptCommandLine]>2,ToExpression[$ScriptCommandLine[[3]]],1];
tf=If[Length[$ScriptCommandLine]>3,ToExpression[$ScriptCommandLine[[4]]],10^5];
l0=If[Length[$ScriptCommandLine]>4,ToExpression[$ScriptCommandLine[[5]]],0.037];
lmax=If[Length[$ScriptCommandLine]>5,ToExpression[$ScriptCommandLine[[6]]],25];
xlinkOutputType=If[Length[$ScriptCommandLine]>6,ToExpression[$ScriptCommandLine[[7]]],1];
fidx=If[Length[$ScriptCommandLine]>7,ToExpression[$ScriptCommandLine[[8]]],0];
fov=If[Length[$ScriptCommandLine]>8,ToExpression/@{$ScriptCommandLine[[9]],$ScriptCommandLine[[10]]},{35,20}];




lks=pts4[dir,"links"];
flks=Map[Cases[#,_?((#[[5]]==fidx)&)]&,lks];
Print["imported links"];

conn1s=If[
xlinkOutputType==1,
pts4[dir,"spacers1_bound"],
With[{s1s=pts2[dir,"amotors"]},Table[Cases[s1s[[t,All]],_?((#[[5]]!=-1&&#[[6]]!=-1)&),1,Heads->False],{t,Length[s1s]}]]
];
Print["imported spacers1"];
conn2s=If[
xlinkOutputType==1,
pts4[dir,"spacers2_bound"],
With[{s2s=pts2[dir,"pmotors"]},Table[Cases[s2s[[t,All]],_?((#[[5]]!=-1&&#[[6]]!=-1)&),1,Heads->False],{t,Length[s2s]}]]
];
Print["imported spacers2"];
tf=Min[Join[{tf},Length/@{lks,conn1s,conn2s}]];
ls=Range[0,lmax-l0,l0];


out=OpenWrite[dir<>"/lattice.dat"];



Do[

d1=Table[distFromBarbedEnd[flks[[t]],fidx,conn1s[[t,i]],fov],{i,Length[conn1s[[t]]]}];
d2=Table[distFromBarbedEnd[flks[[t]],fidx,conn2s[[t,i]],fov],{i,Length[conn2s[[t]]]}];

bcs1=BinCounts[d1,{0,lmax,l0}];
bcs2=BinCounts[d2,{0,lmax,l0}];

lattice=Sign[bcs1-bcs2]/.{-1->2};
(*lattice=Table[If[bcs1[[i]]\[Equal]0&&bcs2[[i]]\[Equal]0, 0, If[bcs1[[i]]>bcs2[[i]],1, If[bcs2[[i]]>bcs1[[i]], 2, RandomInteger[]+1]]],{i,Length[bcs1]}];*)
WriteLine[out,StringRiffle[lattice,","]];
Print["t="<>ToString@t],
{t,ti,tf}
];

Close@out;
(*Determine parallelness of filaments*)
f0lks=Map[Cases[#,_?((#[[5]]==0)&)]&,lks];
f1lks=Map[Cases[#,_?((#[[5]]==1)&)]&,lks];
f0e2e=Total[#[[All,3;;4]]]&/@f0lks;
f1e2e=Total[#[[All,3;;4]]]&/@f1lks;
cosths=Table[f0e2e[[t]].f1e2e[[t]]/(Norm[f0e2e[[t]]]Norm[f1e2e[[t]]]),{t,Length[lks]}];
Export[dir<>"/analysis/ang_bw_fils.dat",cosths];
