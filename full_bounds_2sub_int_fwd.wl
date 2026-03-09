(* ::Package:: *)

(* ::Input::Initialization:: *)
dir=SetDirectory[NotebookDirectory[]];

<<"SDPB.m"
SDPBversion="2.5.1";

runsdpbmac[file_,input_,precision_]:=Module[{judge},RunProcess[{"docker","run","-v",StringJoin[{StringReplace[dir,"\\"->"/"],":/data"}],"wlandry/sdpb:2.5.1","mpirun","--allow-run-as-root","--oversubscribe","-np",nodenum//ToString,"pvm2sdp",precision//ToString,StringJoin[{"/data/",file}],StringJoin[{"/data/",input}]},"StandardOutput",ProcessEnvironment-><|"PATH"->"/usr/local/bin:"<>Environment["PATH"]|>];
judge=RunProcess[{"docker","run","-v",StringJoin[{StringReplace[dir,"\\"->"/"],":/data"}],"wlandry/sdpb:2.5.1","mpirun","--allow-run-as-root","--oversubscribe","-np",nodenum//ToString,"sdpb",StringJoin["--precision=",precision//ToString],StringJoin[{"--procsPerNode=",nodenum//ToString}],"--dualityGapThreshold=1e-6","--maxIterations=5000","-s",StringJoin[{"/data/",input}]},"StandardOutput",ProcessEnvironment-><|"PATH"->"/usr/local/bin:"<>Environment["PATH"]|>]//ToString;
judge];



nodenum=8;

PhypRule=P->Function[{j,x},Hypergeometric2F1[-j,j+d-3,(d-2)/2,(1-x)/2]];


(* ::Input::Initialization:: *)
getDispRelVec[dispR_]:=D[dispR,{{\[Delta],Mp,Mp . W},1}]
(*Args: 
1. list of eft rules, with g being the objective
2. not improved high E dispersion relation dispRelHighE[Mp,W]
3. not improved low E dispersion relation dispLowE[W,R]
4. improved high E symmetric rep dispersion relation, in terms of {t,\[Mu],Mp,Mp.W,P[J,1],P[J,1+(2 t)/\[Mu]],\[Delta],(P^(0,1))[J,1]}
5. improved high E antisymmetric rep dispersion relation, in terms of {t,\[Mu],Mp,Mp.W,P[J,1],P[J,1+(2 t)/\[Mu]],\[Delta],(P^(0,1))[J,1]}
6. improved low E dispersion relation in terms of B[s,t] (vector) and WB[s,t]=W.B[s,t] and (8 G \[Pi] Pt[1])/t
7. fileFcVals file for caching terms for integral calculations
8. file for storing integrals
9. file for storing forward limit equations
10. nMax = max n, power of the smearing polynomial (1-p)p^n in the integrals
11. Jmax=maximum spin for integrals
12. max number of forward limit equations
13. Mp = M plus matrix that maps abcd->adbc
14. W matrix that maps abcd to abdc (should be diagonal with elements +-1)
15. max k of the forward limit equations
16. spacetime dimension 
17. list of mass sampling points, x, with x defined as \[Mu]=1/(1-x)
18. impact parameter step size, \[Delta]b, for finite impact parameter regime
19.  maximum impact parameter, Bmax, for finite impact parameter regime
20. maximum number of subleading terms in the large impact parameter limit (A1 of Sharp Boundaries for the Swampland)
21. filename for saving polynomials (sdpb input)
22. bool flag for loading saved polynomials
23. filename for saving the output of sdpb ( inequality of min, max g, and their optimal functional vector for each eft value)
24. bool flag for including graviton. if False then it calculates forward limit equations for t singlet as well
25. list of excluded reps
26. bool flag for including large impact param limit
*)
calculateBounds[eftVals_,dispRelHighE_,dispLowE_,CimpEvenR_,CimpOddR_,CimpLowE_,fileFcVals_,integralsFileName_,fwdFileName_,nMax_Integer,Jmax_Integer,nFwdMax_Integer,Mp_,W_,kMax_Integer,d_Integer,xs_,\[Delta]b_,Bmax_,mmax_,polsFileName_,savePolsFile_,boundsFileName_,hasGrav_:True,excludeReps_:{},inclLargeImpPar_:True]:= Module[{highEintsSymm,highEintsAnti,finImpParSymm,finImpParAnti,fwdEqnsLowE,fwdEqnsHighE,CimpLowEintegrals,fwdEqnsHighESymm,fwdEqnsHighEAnti,symmReps,antiReps,colorBasis,finImpParImprDispInts,finImpParFwdEqnsHighE},

symmReps=Flatten@Position[Diagonal[W],_?Positive];
antiReps=Flatten@Position[Diagonal[W],_?Negative];


(*get the integrals*)
If[FileExistsQ[integralsFileName],{highEintsSymm,highEintsAnti,finImpParSymm,finImpParAnti}=Import[integralsFileName];
Print["importing integrals from a file"];,
Print["integrals file not found, calculating"];
{highEintsSymm,highEintsAnti,finImpParSymm,finImpParAnti}=
GetHighEintegrals[getDispRelVec[CimpEvenR],getDispRelVec[CimpOddR],nMax,Jmax,fileFcVals,integralsFileName];
];

colorBasis={IdentityMatrix[Length[W]],Mp,Mp . W};

highEintsSymm= highEintsSymm . colorBasis[[;;,symmReps]]; (*J,nint,rep,n_col_b*)
highEintsAnti= highEintsAnti . colorBasis[[;;,antiReps]];
finImpParSymm=finImpParSymm . colorBasis[[;;,symmReps]]; (*nint,rep,n_col_b*)
finImpParAnti=finImpParAnti . colorBasis[[;;,antiReps]];

If[FileExistsQ[fwdFileName],{fwdEqnsLowE,fwdEqnsHighE}=Import[fwdFileName];
Print["importing forward lim eqns from a file"];,
Print["integrals forward lim eqns not found, calculating"];
{fwdEqnsLowE,fwdEqnsHighE}=GetForwardLimEqns[dispRelHighE,dispLowE,Mp,W,kMax,Jmax,d,hasGrav];
Export[fwdFileName,{fwdEqnsLowE,fwdEqnsHighE}];
];

CimpLowEintegrals=GetLowEintegrals[CimpLowE,Mp,W,3,nMax]; (*nint,n_col_b*)(*only k=3 needed for improved disp rel*)
fwdEqnsHighESymm=Table[fwdEqnsHighE[[R]],{R,symmReps},{J,0,Jmax,2}];
fwdEqnsHighEAnti=Table[fwdEqnsHighE[[R]],{R,antiReps},{J,1,Jmax,2}];




Print["getting finite impact param eqns"];
finImpParImprDispInts=Join[finImpParSymm,finImpParAnti,2]; (*nint,rep,n_col_b*)

finImpParFwdEqnsHighE=Table[GetFintImpParLim[fwdEqnsHighE[[i,j]]],{i,Dimensions[fwdEqnsHighE][[1]]},{j,Dimensions[fwdEqnsHighE][[2]]}] (*R,n fwd*);
Print["get bounds"];
getBounds[d,xs,\[Delta]b,Bmax,Mp,W,highEintsSymm,highEintsAnti,CimpLowEintegrals,finImpParImprDispInts,fwdEqnsLowE,fwdEqnsHighESymm,fwdEqnsHighEAnti,finImpParFwdEqnsHighE,If[hasGrav,nMax,0],Jmax,nFwdMax,eftVals,symmReps,antiReps,excludeReps,mmax,polsFileName,savePolsFile,boundsFileName,inclLargeImpPar]
]


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Get the integrals*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input::Initialization:: *)
GetLargeImpParInts[finImpParInts_,dim_,mmax_]:=Module[{lim,A,B,c},lim=Simplify[(Normal[Assuming[b>0,Series[(finImpParInts/.d->dim),{b,Infinity,mmax+5}]]]//Simplify)/(1/b)^(mmax+5)//Cancel//Simplify//ComplexExpand//TrigReduce,b>0]//TrigReduce;
A=lim/. Cos[b]->0/. Sin[b]->0;
B=D[lim/. Cos[b]->q,q]//Simplify;
c=D[lim/. Sin[b]->q,q]//Simplify;
{{A+B,c},{c,A-B}}]


(* ::Input::Initialization:: *)
safeJoin[A_,B_,lvl_]:=Which[A==={}||A===Null,B,B==={}||B===Null,A,True,Join[A,B,lvl]];

safeMult[A_,B_]:=Which[A==={}||A===Null||B==={}||B===Null,{},True,A . B];

getBounds[dim_Integer,xs_,\[Delta]b_,Bmax_,Mp_,W_,highEintsSymm_,highEintsAnti_,CimpLowEintegrals_,finImpParImprDispInts_,fwdEqnsLowE_,fwdEqnsHighESymm_,fwdEqnsHighEAnti_,finImpParFwdEqnsHighE_,nnMax_Integer,Jmax_Integer,nFwdMax_Integer,eftVals_,symmReps_,antisymmReps_,excl_,mmax_,polsFileName_,saveFile_,boundsFileName_,inclLargeImpPar_:True]:=Module[{largeImpParSymm,largeImpParAnti,getPols,lowE,
lowEintsSymmRepMatrices,
lowEintsAntiSymmRepMatrices,
highEintsSymmRepMatrices,
highEintsAntiRepMatrices,
finImpParSymmMatrices,
finImpParAntiMatrices,
lowForwardSymmMatrices,
lowForwardAntiMatrices,
highForwardSymmMatrices,
highForwardAntiMatrices,elimC,nullSp,highESymm,highEAnti,finiteImpParHighE,largeImpPar,nnTake},(*
Args: 
 highEintsSymm[[J,nn,R,n_col_b]]  
highEintsAnti[[J,nn,R,n_col_b]] 
CimpLowEintegrals[[nn,n_col_b]] 
finImpParImprDispInts[[nn, R, n_col]]
fwdEqnsLowE[[nFwd]]
fwdEqnsHighESymm[[R,J,nFwd]] 
fwdEqnsHighEAnti[[R,J,nFwd]] 
finImpParFwdEqnsHighE[[R,nFwd]]
*)
(*nn=number of integrals*)


nnTake=Max[0,nnMax-1];
	
lowE=Join[Take[CimpLowEintegrals,UpTo[nnTake],All]//Flatten,Take[fwdEqnsLowE,UpTo[nFwdMax]]]; (*nn*n_col_b+nFwd*)

highESymm=safeJoin[Flatten[Take[highEintsSymm,All,UpTo[nnTake],All,All],{{3},{1},{2,4}}],Take[fwdEqnsHighESymm,All,All,UpTo[nFwdMax]],3];(*R,J,nn*n_col_b+nFwd*)

highEAnti=safeJoin[Flatten[Take[highEintsAnti,All,UpTo[nnTake],All,All],{{3},{1},{2,4}}],Take[fwdEqnsHighEAnti,All,All,UpTo[nFwdMax]],3];

finiteImpParHighE=safeJoin[Flatten[Take[finImpParImprDispInts,UpTo[nnTake],All,All],{{2},{1,3}}],Take[finImpParFwdEqnsHighE,All,UpTo[nFwdMax]],2];(*R,nn*n_col_b+nFwd*)


largeImpPar=GetLargeImpParInts[finiteImpParHighE,dim,mmax]; (*{2,2,R,nn*n_col_b+nFwd}*)

elimC=Append[Complement[lowE//Variables,Append[Table[Level[eftVals[[1,i]],1][[1]],{i,Length[eftVals[[1]]]}],t]],GGG];

If[Length[elimC]>0,
nullSp=D[lowE,{elimC,1}]//Transpose//NullSpace;
lowE=nullSp . lowE; (*n_indep*)
highESymm=safeMult[highESymm,Transpose[nullSp]]; (*R,J,n_indep*)
highEAnti=safeMult[highEAnti,Transpose[nullSp]];
finiteImpParHighE=safeMult[finiteImpParHighE,Transpose[nullSp]]; (*R,n_indep*)
largeImpPar=safeMult[largeImpPar,Transpose[nullSp]]; (*{2,2,R,n_indep}*)
];
Print["elim"];

getPols[fileName_]:=Module[{Pols1,Pols2,Pols3,Pols,symmRepsUse,antiRepsUse,symmRepsAll,antiRepsAll,symmPos,antiPos},
symmRepsUse=Complement[symmReps,excl];
antiRepsUse=Complement[antisymmReps,excl];

Print["using symmetric reps:",symmRepsUse];
Print["using antisymmetric reps:",antiRepsUse];

symmPos=AssociationThread[symmReps->Range[Length[symmReps]]];

antiPos=AssociationThread[antisymmReps->Range[Length[antisymmReps]]];

symmRepsUse=Complement[symmReps,excl];
antiRepsUse=Complement[antisymmReps,excl];

Pols1=Join[Flatten[Table[With[{i=symmPos[r]},PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x],{{ \[Mu]^2/(1+J^3)highESymm[[i,J/2+1]]/. \[Mu]->1/(1-x)/. x->X/. d->dim}}]],{J,0,Jmax,2},{X,xs},{r,symmRepsUse}],2],Flatten[Table[With[{i=antiPos[r]},PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x],{{\[Mu]^2/(1+J^3)highEAnti[[i,(J+1)/2]]/. \[Mu]->1/(1-x)/. x->X/. d->dim}}]],{J,1,Jmax,2},{X,xs},{r,antiRepsUse}],2]];


Pols1=DeleteCases[Pols1,PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x_],{{ConstantArray[0.,Length[lowE]]}}]/;True];
Pols1=DeleteCases[Pols1,PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x_],{{ConstantArray[0,Length[lowE]]}}]/;True];


(*finImpParSymm[[R, nn]] 
 *)
Pols2=Table[PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x],{{finiteImpParHighE[[r]]/. d->dim}}],{r,Join[symmRepsUse,antiRepsUse]},{b,\[Delta]b,Bmax,\[Delta]b}]//Flatten;


Pols2=DeleteCases[Pols2,PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x_],{{ConstantArray[0.,Length[lowE]]}}]/;True];
Pols2=DeleteCases[Pols2,PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x_],{{ConstantArray[0,Length[lowE]]}}]/;True];



Pols3=Table[PositiveMatrixWithPrefactor[DampedRational[1,{},1/E,x],largeImpPar[[;;,;;,r]]/. b->Bmax+x/.d->dim],{r,Join[symmRepsUse,antiRepsUse]}];
If[inclLargeImpPar,Pols=Join[Pols1,Pols2,Pols3],Pols=Join[Pols1,Pols2]];


Export[fileName,Pols];
Print["file exported ",fileName];
Pols];


(*Join[Pols1,Pols2,Pols3]*)
If[FileExistsQ[polsFileName]&&saveFile,polsFull=Import[polsFileName];
Print["importing pols file"];,polsFull=getPols[polsFileName];
Print["pols file not found, creating"];];
bounds[eftValRule_,Pols_]:=Module[{nvec,bvec},
Print[eftValRule];
nvec=-D[lowE/. eftValRule,g]//Simplify;
bvec=-lowE/. eftValRule/. g->0//Simplify;
Print[nvec];
Print[bvec];
DeleteFile["test.xml"];
DeleteFile["test"];
DeleteDirectory[{"test.ck","test_out"},DeleteContents->True];
WriteBootstrapSDP["test.xml",SDP[bvec,-nvec,N[Pols/.d->dim,50]]];
runsdpbmac["test.xml","test",1024*2];
y1=Insert[ReadList["test_out/y.txt"][[2;;]],x,Position[Abs[nvec],Max[Abs[nvec]]][[1,1]]]/. a_+b_ e->b*10^a;
Print[lowE . y1>0/. Solve[nvec . y1==-1,x][[1]]/. eftValRule//Simplify];
(*Print[lowE.y1>0/. Solve[nvec.y1==-1,x][[1]]//Simplify]*);
DeleteFile["test.xml"];
DeleteFile["test"];
DeleteDirectory[{"test.ck","test_out"},DeleteContents->True];
WriteBootstrapSDP["test.xml",SDP[bvec,nvec,Pols]];
runsdpbmac["test.xml","test",1024*2];
y2=Insert[ReadList["test_out/y.txt"][[2;;]],x,Position[Abs[nvec],Max[Abs[nvec]]][[1,1]]]/. a_+b_ e->b*10^a;
Print[lowE . y2>0/. Solve[nvec . y2==1,x][[1]]/. eftValRule//Simplify];
(*Print[lowE.y2>0/. Solve[nvec.y2==1,x][[1]]//Simplify];*);
{lowE . y1>0/. Solve[nvec . y1==-1,x][[1]]//Simplify,lowE . y2>0/. Solve[nvec . y2==1,x][[1]]//Simplify,y1/. Solve[nvec . y1==-1,x][[1]],y2/. Solve[nvec . y2==1,x][[1]]}];

Table[
(*Print[ToString[StringForm["bounds_u1E_d_``_Jmax_``_dispNumInt_``_analyticInt_G_``_``_reps_higherJ_``_funda.wdx",d,Jmax,dispNumInt, G/.eftVal,StringReplace[ToString[Join[reps1,reps2],InputForm],{"{"|"}"->"",", "->"_"}],higherJ]]];*)
bound=bounds[eftVal,polsFull];
Export[ToString[boundsFileName],bound];

NotebookSave[];bound,{eftVal,eftVals}]]


(* ::Input::Initialization:: *)
GetForwardLimEqns[dispRelHighE_,dispLowE_,Mp_,W_,kMax_,Jmax_,dim_,hasGrav_:True]:=Module[{fwdEqnsHighE,fwdEqnsLowE,$P,xVals,toModInv,getIndepColsFf,evenEqns,oddEqns,indepIdx,symmReps,antiReps,BPtBasis,dispRelPtBasis},
(*returns independent forward limit equations*)
symmReps=Flatten@Position[Diagonal[W],_?Positive];
antiReps=Flatten@Position[Diagonal[W],_?Negative];

BPtBasis=D[getB[Mp,kMax,symmReps],{Array[Pt,Length[Mp]],1}];

dispRelPtBasis=Table[dispLowE[W,R]/.B->Function[{S,T,R},BPtBasis[[R]]/.{s->S,t->T}],{R,Length[Mp]}];

fwdEqnsHighE=Join[Table[SeriesCoefficient[Array[r,Length[Mp]] . Simplify[dispRelHighE[Mp,W]/.PhypRule][[;;,2;;]],{s,0,n1},{t,0,n2}]//Cancel,{n1,2,kMax},{n2,0,kMax-n1}],Table[SeriesCoefficient[Array[r,Length[Mp]] . Simplify[dispRelHighE[Mp,W]/.PhypRule][[;;,1]],{s,0,n1},{t,0,n2}]//Cancel,{n1,If[hasGrav,3,2],kMax},{n2,0,kMax-n1}]]//Flatten;





fwdEqnsLowE=Join[Table[SeriesCoefficient[dispRelPtBasis[[2;;]],{s,0,n1},{t,0,n2}]//Cancel,{n1,2,kMax},{n2,0,kMax-n1}]//Flatten,Table[SeriesCoefficient[dispRelPtBasis[[1]],{s,0,n1},{t,0,n2}]//Cancel,{n1,If[hasGrav,3,2],kMax},{n2,0,kMax-n1}]//Flatten]; (*n eqns*)


(*Assert[Position[fwdEqnsHighE,0]==Position[fwdEqnsLowE,0],"the zeros of low and high energy eqns don't match"];*);

fwdEqnsHighE=D[fwdEqnsHighE,{Array[r,Length[Mp],1]}]//Transpose; (*R,n eqns*)



evenEqns=Flatten[Table[fwdEqnsHighE[[R]],{R,symmReps},{J,0,Jmax,2}],{{1,2},{3}}];(*JR,n eqns*)oddEqns=If[antiReps=={},{}, Flatten[Table[fwdEqnsHighE[[R]],{R,antiReps},{J,1,Jmax}],{{1,2},{3}}]]; (*JR,n eqns*)




$P = 2147483497;
xVals=RandomInteger[{0,$P},100];
toModInv[frac_]:=Mod[Numerator[frac]*ModularInverse[Denominator[frac],$P],$P];

getIndepColsFf[pols_]:=Module[{polsTablex,indepCols,rank,mat,rank1,idx},polsTablex=Map[toModInv,Flatten[Table[pols,{x,xVals}],1]//Cancel,{2}];


idx={FirstPosition[Total[Abs[polsTablex],{1}],_?(#!=0&),Missing["NotFound"]][[1]]};
indepCols={polsTablex[[;;,idx[[1]]]]};

rank=MatrixRank[indepCols,Modulus->$P];
Do[mat=Append[indepCols,polsTablex[[;;,i]]];
rank1=MatrixRank[mat,Modulus->$P];
If[rank<rank1,indepCols=mat;
rank=rank1;
AppendTo[idx,i];];
,{i,2,polsTablex[[1]]//Length}];
idx];

indepIdx=getIndepColsFf[Join[evenEqns,oddEqns]/.d->dim/.\[Mu]->(1+x)];

{fwdEqnsLowE[[indepIdx]],fwdEqnsHighE[[;;,indepIdx]]} (*only independent equations*)


]

PhypRule=P->Function[{j,x},Hypergeometric2F1[-j,j+d-3,(d-2)/2,(1-x)/2]];

GetHighEintegrals[symmRepImprDispHighE_,antisymmRepImprDispHighE_,nMax_Integer,Jmax_Integer,fileFcVals_,exportFileName_]:=Module[{highEintsSymm,highEintsAnti,out,finImpParSymm,finImpParAnti,highE},
(*symmRepImprDispHighE[[neqns]] *)

highEintsSymm=Transpose[ComputeSmearedIntegrals[symmRepImprDispHighE/.PhypRule//Simplify,0,Jmax,nMax,fileFcVals],{3,1,2}];(*neqns,J,nint*)
highEintsAnti=Transpose[ComputeSmearedIntegrals[antisymmRepImprDispHighE/.PhypRule//Simplify,1,Jmax,nMax,fileFcVals],{3,1,2}];(*neqns,J,nint*)
finImpParSymm=GetFintImpParInts[symmRepImprDispHighE,nMax];(*nint,neqns*)
finImpParAnti=GetFintImpParInts[antisymmRepImprDispHighE,nMax];(*nint,neqns*)
out={highEintsSymm,highEintsAnti,finImpParSymm,finImpParAnti};
SetDirectory[NotebookDirectory[]];
Export[exportFileName,out];
out]

GetLowEintegrals[CimpLowE_,Mp_,W_,kMax_,nMax_]:=Module[{BPtBasis,CimpLowE1,highEintsAnti,out,highE,lowE},

BPtBasis=D[getB[Mp,kMax,Flatten@Position[Diagonal[W],_?Positive]],{Array[Pt,Length[Mp]],1}];

CimpLowE1=CimpLowE/.{Pt[1]->D[Pt[1],{Array[Pt,Length[Mp]],1}],B->Function[{S,T},BPtBasis/.{s->S,t->T}//Evaluate],WB->Function[{S,T},W . BPtBasis/.{s->S,t->T}//Evaluate]};

ComputeLowEintegrals[CimpLowE1,nMax]]

ComputeLowEintegrals[lowE_,nMax_Integer]:=Table[Assuming[0<x<1 && nn>1,Integrate[lowE(p^(nn) (1-p))/.t->-p^2,{p,0,1}]
],{nn,2,nMax}]
ComputeSmearedIntegrals[highE_,Jmin_Integer,Jmax_Integer,nMax_Integer,fileFcVals_]:=Module[{intFac,pInt,imprDispRelHighE1,eqns,ans,getIntFac,fcExpandVals,pPowTable,integPn,rescAppR},
SetDirectory[NotebookDirectory[]];
getIntFac[dispR_]:=Module[{maxP,maxM},
maxP=Table[{Exponent[Denominator[term],(t+\[Mu])]},{term,dispR}]//Max;
maxM=Table[{Exponent[Denominator[term],(t-\[Mu])]},{term,dispR}]//Max;
(p^nn) /((-1+p^2 (-1+x))^maxM (1+p^2 (-1+x))^maxP)
];

intFac=getIntFac[highE];
pInt=Assuming[0<x<1 && nn>1,Integrate[intFac,{p,0,1}]];
If[FileExistsQ[fileFcVals],
fcExpandVals=Import[fileFcVals];,
fcExpandVals=ParallelTable[pInt//FunctionExpand//Simplify,{nn,0,500}];
Export[fileFcVals,fcExpandVals];
];


pPowTable[exp_]:=ParallelTable[D[exp,{p,i}]/i!/.p->0,{i,0,Exponent[exp,p]}];
integPn[exp_,Nn_]:=Module[{ta},ta=CoefficientList[exp,p];ta . fcExpandVals[[1+Nn;;Length[ta]+Nn]]];
rescAppR[expr_,intFac_]:=expr(1-p)p^nn /intFac/. \[Mu]->1/(1-x)/.t->-p^2;
eqns=rescAppR[highE,intFac];

ans=ParallelTable[(integPn[eq//FunctionExpand,Nn])//Cancel//Simplify,{eq,eqns},{J,Jmin,Jmax,2},{Nn,2,nMax}];

ans
]

GetFintImpParLim[highE_]:=Module[{lim},lim=((Series[(highE/.J-> b*m/2/.\[Mu]->m^2/.t->-p^2/.{\!\(\*SuperscriptBox[\(P\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[(b m)/2,1]->(b^2 m^2)/(4 (-2+d))}/.P[(b m)/2,1-(2 p^2)/m^2]->Gamma[(d-2)/2]BesselJ[(d-4)/2,b*p]/(b*p/2)^((d-4)/2)/.P[(b m)/2,1]->(Gamma[(d-2)/2]BesselJ[(d-4)/2,b*p]/(b*p/2)^((d-4)/2)/.p->0))/.m->1/minv//Simplify,{minv,0,4}]//Normal)//FullSimplify)/.minv->1/m;
lim=lim*m^4//Cancel; (*absorbing m dependence into spectral density, if there is remaining m dependence something is wrong*)
lim]

GetFintImpParInts[highE_,nMax_Integer]:=Module[{lim},lim=((Series[(highE/.J-> b*m/2/.\[Mu]->m^2/.t->-p^2/.{\!\(\*SuperscriptBox[\(P\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[(b m)/2,1]->(b^2 m^2)/(4 (-2+d))}/.P[(b m)/2,1-(2 p^2)/m^2]->Gamma[(d-2)/2]BesselJ[(d-4)/2,b*p]/(b*p/2)^((d-4)/2)/.P[(b m)/2,1]->(Gamma[(d-2)/2]BesselJ[(d-4)/2,b*p]/(b*p/2)^((d-4)/2)/.p->0))/.m->1/minv//Simplify,{minv,0,4}]//Normal)//FullSimplify)/.minv->1/m;
lim=lim*m^4//Cancel; (*absorbing m dependence into spectral density, if there is remaining m dependence something is wrong*)
On[Assert];
Assert[FreeQ[lim,m]];
ParallelTable[Assuming[0<x<1 && nn>1&&b>0,Integrate[lim(p^(nn) (1-p))/.t->-p^2,{p,0,1}]
],{nn,2,nMax}]
]


(* ::Input::Initialization:: *)
(*1234 s, 1423 t, 1342 u*)

toSTR[expr_,Ps_List]:=Module[{res=0,ci,cr},Do[(*coefficient of the i-th projector/basis element*)ci=Coefficient[Expand[expr],Ps[[i]]];
(*now expand ci in s,t and record powers*)cr=CoefficientRules[Expand[ci],{s[1,2],s[1,4]}];
res+=Total[(#[[2]]*st[Sequence@@#[[1]],i])&/@cr],{i,Length[Ps]}];
res];
stPolsCol[k_,nReps_,iList_]:=Flatten[Table[s[1,iList[[1]]]^(k-q)(s[1,iList[[3]]])^q  P[a[1],a[iList[[3]]],a[iList[[1]]],a[iList[[2]]]][R],{q,0,k},{R,nReps}]]
commonFixedSubspace[A_,B_]:=NullSpace[Join[A-IdentityMatrix[Length[A]],B-IdentityMatrix[Length[B]]]]

stuSymmPolsCol[k_,Mp_,symmReps_]:=(*args: k of the s^k t^(k-q) coeffs, matrices from t,u channels to s channel*)Module[{basis,pols,coeffs},pols=stPolsCol[k,Length[Mp],{2,3,4}]/.s[1,3]->-s[1,2]-s[1,4];
basis=Table[toSTR[x,Array[P[a[1],a[4],a[2],a[3]],Length[Mp]]],{x,pols}];
coeffs=commonFixedSubspace[D[Table[toSTR[x,Array[P[a[1],a[4],a[2],a[3]],Length[Mp]]],{x,stPolsCol[k,Length[Mp],{4,2,3}]/.Table[Array[P[a[1],a[3],a[4],a[2]],Length[Mp]][[i]]->(Mp . Array[P[a[1],a[4],a[2],a[3]],Length[Mp]]//Simplify)[[i]],{i,Length[Mp]}]/.s[1,3]->-s[1,2]-s[1,4]//Simplify}],{basis,1}]//Transpose,D[Table[toSTR[x,Array[P[a[1],a[4],a[2],a[3]],Length[Mp]]],{x,stPolsCol[k,Length[Mp],{3,2,4}]/.P[a[1],a[4],a[3],a[2]][R_]:>If[MemberQ[symmReps,R],P[a[1],a[4],a[2],a[3]][R],-P[a[1],a[4],a[2],a[3]][R]]/.s[1,3]->-s[1,2]-s[1,4]//Simplify}],{basis,1}]//Transpose];

If[Length[coeffs]>0,coeffs . pols/.s[1,2]->s/.s[1,4]->t/.P[a[1],a[4],a[2],a[3]]->Pt//Simplify,{}]]


getB[Mp_,kMax_,symmReps_]:=Module[{pols},Sum[pols=stuSymmPolsCol[k,Mp,symmReps];
If[Length[pols]>0, Table[g[k,i],{i,Length[pols]}] . pols,0],{k,0,kMax}]]
