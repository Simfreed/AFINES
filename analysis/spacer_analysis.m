(* ::Package:: *)

#!/software/mathematica-10.2-x86_ 64/bin/MathematicaScript -script


dir=$ScriptCommandLine[[2]];
Import["~/Code/cytomod/analysis/cytomod_functions.m"];
fov=If[Length[$ScriptCommandLine]>2,ToExpression/@{$ScriptCommandLine[[3]],$ScriptCommandLine[[4]]},{35,20}];
fidx=If[Length[$ScriptCommandLine]>4,ToExpression[$ScriptCommandLine[[5]]],0];
ti=If[Length[$ScriptCommandLine]>5,ToExpression[$ScriptCommandLine[[6]]],1];
tf=If[Length[$ScriptCommandLine]>6,ToExpression[$ScriptCommandLine[[7]]],2000];


s1s=pts2[dir,"amotors"];
s2s=pts2[dir,"pmotors"];
lks=pts2[dir,"links"];

tf=Min[Join[{tf},Length/@{s1s,s2s,lks}]];


conn1s=Table[Cases[s1s[[t,All]],_?((#[[5]]!=-1&&#[[6]]!=-1)&),1,Heads->False],{t,Length[s1s]}];
conn2s=Table[Cases[s2s[[t,All]],_?((#[[5]]!=-1&&#[[6]]!=-1)&),1,Heads->False],{t,Length[s2s]}];


tf=Min[Join[{tf},Length/@{lks,conn1s,conn2s}]];
distfile=dir<>"/analysis/spacer_dists_sorted_f"<>ToString[fidx]<>"_t"<>ToString[ti]<>"-"<>ToString[tf]<>".dat";
doms1file=dir<>"/analysis/spacer1_domains_f"<>ToString[fidx]<>"_t"<>ToString[ti]<>"-"<>ToString[tf]<>".dat";
doms2file=dir<>"/analysis/spacer2_domains_f"<>ToString[fidx]<>"_t"<>ToString[ti]<>"-"<>ToString[tf]<>".dat";


(*dom1s=Table[SequenceCases[dssorted[[t]],{Repeated[{_,1}]}],{t,Length[dssorted]}];
dom2s=Table[SequenceCases[dssorted[[t]],{Repeated[{_,2}]}],{t,Length[dssorted]}];
Export[doms1file, dom1s];
Export[doms2file, dom2s];
*)
out0=OpenWrite[distfile];
out1=OpenWrite[doms1file];
out2=OpenWrite[doms2file];
flks=Map[Cases[#,_?((#[[5]]==fidx)&)]&,lks];

Do[
d1=Table[{distFromBarbedEnd[flks[[t]],fidx,conn1s[[t,i]],fov],1},{i,Length[conn1s[[t]]]}];
d2=Table[{distFromBarbedEnd[flks[[t]],fidx,conn2s[[t,i]],fov],2},{i,Length[conn2s[[t]]]}];
dssorted=Sort@Join[d1,d2];
WriteLine[out0,ToString@dssorted];
WriteLine[out1,ToString@SequenceCases[dssorted,{Repeated[{_,1}]}]];
WriteLine[out2,ToString@SequenceCases[dssorted,{Repeated[{_,2}]}]];
Print["t="<>ToString@t],
{t,ti,tf}
];
Close[out0];
Close[out1];
Close[out2];

(*Determine parallelness of filaments*)
f0lks=Map[Cases[#,_?((#[[5]]==0)&)]&,lks];
f1lks=Map[Cases[#,_?((#[[5]]==1)&)]&,lks];
f0e2e=Total[#[[All,3;;4]]]&/@f0lks;
f1e2e=Total[#[[All,3;;4]]]&/@f1lks;
cosths=Table[f0e2e[[t]].f1e2e[[t]]/(Norm[f0e2e[[t]]]Norm[f1e2e[[t]]]),{t,Length[lks]}];
Export[dir<>"/analysis/ang_bw_fils.dat",cosths];
