(* ::Package:: *)

#!/software/mathematica-10.2-x86_ 64/bin/MathematicaScript -script


dir=$ScriptCommandLine[[2]];
Import["~/Code/cytomod/analysis/cytomod_functions.m"];
ti=If[Length[$ScriptCommandLine]>2,ToExpression[$ScriptCommandLine[[3]]],1];
tf=If[Length[$ScriptCommandLine]>3,ToExpression[$ScriptCommandLine[[4]]],10^5];
nc=If[Length[$ScriptCommandLine]>4,ToExpression[$ScriptCommandLine[[5]]],4];
rc=If[Length[$ScriptCommandLine]>5,ToExpression[$ScriptCommandLine[[6]]],1];
l0=If[Length[$ScriptCommandLine]>6,ToExpression[$ScriptCommandLine[[7]]],0.037];
lmax=If[Length[$ScriptCommandLine]>7,ToExpression[$ScriptCommandLine[[8]]],25];



seqs=Import[dir<>"/lattice.dat","CSV"];
tf=Min[tf,Length[seqs]];
ls=Range[0,lmax-l0,l0];


files=Table[StringJoin[ToString/@{dir,"/analysis/spacer",i,"_domains_t",ti,"-",tf,"_nc",nc,"_rc",rc,".dat"}],{i,0,2}];
outs=OpenWrite/@files;


Do[
lf=Min[Length[seqs[[t]]],Length[ls]];
doms=LatticeDomainMerge[{ls[[1;;lf]],seqs[[t,1;;lf]]}\[Transpose],nc,rc];
dps=Table[Flatten@Position[doms,_?((#[[1,2]]==i)&),{1},Heads->False],{i,0,2}];
Do[WriteLine[outs[[i]],ToString[doms[[dps[[i]],All,1]]]],{i,3}];
(*domsg=Sort[GatherBy[doms,#[[1,2]]&],#1[[1,1,2]]<#2[[1,1,2]]&];*)
(*WriteLine[out1,ToString@doms[[1]]];
WriteLine[out2,ToString@doms[[2]]];*)
Print["t="<>ToString@t],
{t,ti,tf}
];

Close/@outs;
