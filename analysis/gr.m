(* ::Package:: *)

#!/software/mathematica-10.2-x86_ 64/bin/MathematicaScript -script
(* Given: 
    links.dat file
actins.dat file
field of view
# segments per filament
velocities from divergence calculation*)


dir=$ScriptCommandLine[[2]];
fov=If[Length[$ScriptCommandLine]>2,ToExpression/@{$ScriptCommandLine[[3]],$ScriptCommandLine[[4]]},{50,50}];
dr=If[Length[$ScriptCommandLine]>4,ToExpression[$ScriptCommandLine[[5]]],0.05];
maxr=If[Length[$ScriptCommandLine]>5,ToExpression[$ScriptCommandLine[[6]]],25];
ti=If[Length[$ScriptCommandLine]>6,ToExpression[$ScriptCommandLine[[7]]],1];
tf=If[Length[$ScriptCommandLine]>7,ToExpression[$ScriptCommandLine[[8]]],400];
nsamp=If[Length[$ScriptCommandLine]>8,ToExpression[$ScriptCommandLine[[9]]],10000];

Import["/home/simonfreedman/Code/cytomod/analysis/cytomod_functions.m"];
Print["Loaded Functions"];
acts = pts2[dir, "actins"];
If[Length[Dimensions[acts]] == 1, acts = acts[[1 ;; -2]]];
Print["Imported Actins"];


(*Calculate radial distribution of barbed ends*)
outfile=OpenWrite[dir<>"/analysis/gr_dr"<>ToString[dr]<>"_maxr"<>ToString[maxr]<>"t"<>ToString[ti]<>"-"<>ToString[tf]<>".dat"];
t=ti;
Do[
Print["t = "<>ToString[t]];
WriteString[outfile,rdfResamp[acts[[t]],fov,dr,maxr,nsamp],"\n"];,
{t,ti,tf}];
Close[outfile];
(*Export[dir<>"/analysis/gr_dr"<>ToString[dr]<>"_maxr"<>ToString[maxr]<>"t"<>ToString[ti]<>"-"<>ToString[tf]<>".dat",gr];*)
Print["Finished calculating rdf"];

