(* ::Package:: *)

#!/software/mathematica-10.2-x86_ 64/bin/MathematicaScript -script
(* Given: 
    links.dat file
actins.dat file
field of view
# segments per filament
velocities from divergence calculation*)


dir=$ScriptCommandLine[[2]];
ti=ToExpression[$ScriptCommandLine[[3]]];
tf=ToExpression[$ScriptCommandLine[[4]]];
dth=ToExpression[$ScriptCommandLine[[5]]];
zval=ToExpression[$ScriptCommandLine[[6]]];
dts=Map[ToExpression,$ScriptCommandLine[[7;;-1]]];

Import["/home/simonfreedman/Code/cytomod/analysis/cytomod_functions.m"];
Print["Loaded Functions"];
(*Uses unrolled trajectory, to avoid boundary defects*)
mots=pts2[dir,"amotors_ext"];
Print["imported mots; dimensions: "<>ToString[Dimensions@mots]];
cms=mots[[ti;;tf,All,1;;2]]+0.5*mots[[ti;;tf,All,3;;4]];
Print["calculated center of masses; dimensions: "<>ToString[Dimensions@cms]];
Do[
vs=cms[[(dt+1);;-1;;dt]]-cms[[1;;(-dt-1);;dt]];
Print["calculated velocities; dimensions: "<>ToString[Dimensions@vs]];
ns=Map[Norm,vs,{2}];
Print["calculating angles"];
angs=Flatten@Table[vecangle[vs[[t,i]],vs[[t+1,i]],ns[[t,i]],ns[[t+1,i]],zval],{t,Length[vs]-1},{i,Length[vs[[t]]]}];
bcs=BinCounts[angs/(2Pi),{0,1,dth}];
Export[dir<>"/analysis/vmot_anghist2_t"<>ToString[ti]<>"-"<>ToString[tf]<>"d"<>ToString[dt]<>"dth"<>ToString[dth]<>"_zval"<>ToString[zval]<>".dat",bcs];
,{dt,dts}
];

Print["Finished calculating velocity histogram"];
