(* ::Package:: *)

#!/software/mathematica-10.2-x86_ 64/bin/MathematicaScript -script

dir=$ScriptCommandLine[[2]];
fov=If[Length[$ScriptCommandLine]>2,ToExpression/@{$ScriptCommandLine[[3]],$ScriptCommandLine[[4]]},{50,50}];
ti=If[Length[$ScriptCommandLine]>4,ToExpression[$ScriptCommandLine[[5]]],1];
tf=If[Length[$ScriptCommandLine]>5,ToExpression[$ScriptCommandLine[[6]]],-1];
parts=If[Length[$ScriptCommandLine]>6,$ScriptCommandLine[[7]],"amotors_ext"];


Import["/home/simonfreedman/Code/cytomod/analysis/cytomod_functions.m"];
Print["Loaded Functions"];
amots = pts2[dir, parts];
If[Length[Dimensions[amots]] == 1, amots = amots[[1 ;; -2]]];
Print["Imported Motors"];
amots=amots[[ti;;Min[tf,Length[amots]]]]
If[amots[[1,1]]==Null, Print["no motors for msd"]; Quit[]];

drs=amots[[2;;,All,1;;2]]-amots[[;;-2,All,1;;2]];
(*drs=Map[rijP[#,{50,50}]&,drs,{2}];*)


(*Calculate radial distribution of barbed ends*)
msds=Table[msdc[drs[[All,i]]],{i,Length[drs[[1]]]}];
Export[dir<>"/analysis/msd_t"<>ToString[ti]<>"-"<>ToString[tf]<>".dat",N@msds];
Print["Finished calculating msd"];

