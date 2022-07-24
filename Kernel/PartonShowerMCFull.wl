(* ::Package:: *)

(* ::Text:: *)
(*This code includes the main utilities for simulating the parton shower*)


BeginPackage["PartonShowerMCFull`"];


PgtoggvacuumLT::usage = "Splitting function g -> gg as function of z LT (~1/z)";


PgtoqqbarvacuumLT::usage = "Splitting function g -> gg as function of z";


PgtoggvacuumNTnopol::usage = "Splitting function g -> gg as function of z with more terms"


PgtoggvacuumNT::usage = "Splitting function g -> gg as function of z with more terms"


IntegralThetaPgtoggvacuum::usage = "Part of Sudakov over theta function g -> gg"


SudakovPgtoggvacuumLT::usage = "Sudakov function g -> gg LT";


SudakovPgtoqqbarvacuumLT::usage = "Sudakov function g -> qqbar LT"


SudakovPgtoggvacuumGavinOrig::usage = "Sudakov function g -> gg as in Gavin's example"


SudakovPgtoggvacuumNT::usage = "Sudakov function g -> gg with more terms";


SudakovPgtoggvacuumNTnopol::usage = "Sudakov function g -> gg with more terms but no pol";


PlotCompareSudakov::usage = "Compare Sudakov factor dependence LT and NT"


ptFromSudakov::usage = "ptFromSudakov function as in Gavin's example"


MakeShowerGavin::usage = "Shower g-> gg using original implementation of Gavin"


SingleSplittingPtOrdered::usage = "Generate single splitting from g-> gg, g->qqbar shower"


PgtoggvacuumLT[z_]:=2*1/z;


PgtoqqbarvacuumLT[z_]:=(z^2+(1-z)^2);


PgtoggvacuumNT[z_]:=2*(1-z)/z +z*(1-z);


PgtoggvacuumNTnopol[z_]:=2*(1-z)/z


IntegralThetaPgtoggvacuum =  Integrate[1/(theta),{theta, pt0/(z*pt1),1},Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]


(* ::Text:: *)
(*The function IntegralQ2Pgtoggvacuum is expected to be equivalent to IntegralThetaPgtoggvacuum after a change of variable theta->Q2. *)


IntegralQ2Pgtoggvacuum =  Integrate[1/(2Q2),{Q2, pt0*pt0/z,z*pt1*pt1},Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]


SudakovPgtoggvacuumGavinOrig= Exp[- 2*\[Alpha]s*CA/Pi*Log[pt1/pt0]^2]


SudakovPgtoggvacuumLT= Exp[-2 \[Alpha]s*CA/Pi*Integrate[PgtoggvacuumLT[z]*IntegralThetaPgtoggvacuum,{z, pt0/pt1,1},Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]]


SudakovPgtoqqbarvacuumLT= Exp[-2 \[Alpha]s*CA/Pi*Integrate[PgtoqqbarvacuumLT[z]*IntegralThetaPgtoggvacuum,{z, pt0/pt1,1},Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]]


SudakovPgtoggvacuumNT=Exp[-2 \[Alpha]s*CA/Pi*Integrate[PgtoggvacuumNT[z]*IntegralThetaPgtoggvacuum,{z, pt0/pt1,1},Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]];


SudakovPgtoggvacuumNTnopol=Exp[-2 \[Alpha]s*CA/Pi*Integrate[PgtoggvacuumNTnopol[z]*IntegralThetaPgtoggvacuum,{z, pt0/pt1,1},Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]];


SudakovPgtoggvacuumNTQ2=Exp[-2 \[Alpha]s*CA/Pi*Integrate[1/(2Q2)* PgtoggvacuumNT[z],{z, pt0/pt1,1},{Q2, pt0*pt0/z,z*pt1*pt1}, Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]];


SudakovPgtoggvacuumNTQ2Medium=Exp[-2 \[Alpha]s*CA/Pi*Integrate[1/(2Q2)* PgtoggvacuumNT[z]+ 1/Q2,{z, pt0/pt1,1},{Q2, pt0*pt0/z,z*pt1*pt1}, Assumptions->{pt0\[Element] Reals,pt1\[Element] Reals , pt1>pt0, pt0>0,z<1, z>0, pt0!=pt1*z}]];


ptFromSudakov[sudakovValue_,CA_,alphas_,pt1_]:= pt1 * Exp[-Sqrt[Log[sudakovValue]/(-2*alphas*CA/Pi)]]


(* ::Text:: *)
(*Here you can find a simple description of how the method SingleSplittingPtOrdered works:*)
(*- the inputs of the function are the Sudakov factors for the two splitting processes considered (g->gg and g->qqbar), the parton momentum (mpt1_),  the parton   *)
(*  identity (parton_),  the pt cut off (cutoffpt_)*)
(*- The function first checks whether the initial parton has pt larger than the cut off (If[mpt1>cutoffpt,..])*)
(*- possibleSplits in case the parton is a gluon creates the following list of lists {{sudakovgtogg,{g,g}},{sudakovgtoqqbar,{q,q}}}*)
(*  where sudakovgtogg and sudakovgtoqqbar are the expression for the two sudakov factors. *)
(*- rnd  is a list of two random numbers between 0. and 1. that are used for the Sudakov extraction of the two processing g->gg and g->qqbar. *)
(*  E.g. {0.2324,0.9144}. *)
(*- It extract the corresponding pt values for the emission for the two processes using the corresponding Sudakov values*)
(*  E.g.  resultsp = {8.03779,1.69913}*)
(*- For the cases where at least one of the two pt values are positive, it finds the index of the larger pt value*)


SingleSplittingPtOrdered[sudakovgtogg_, sudakovgtoqqbar_,mpt1_,parton_:"g",cutoffpt_:1.]:=
  Module[{rnd,possibleSplits,resultsp,nozero,results},
  If[mpt1>cutoffpt,
    possibleSplits=If[parton=="g",{{sudakovgtogg,{"g","g"}}, {sudakovgtoqqbar,{"q","q"}}},0];
    rnd = RandomReal[{0.,1.0},WorkingPrecision->4]&/@possibleSplits;
    resultsp=Table[pt0//.FindRoot[(possibleSplits[[i,1]]//.{pt1->mpt1,CA->3,\[Alpha]s->0.12 })-rnd[[i]],{pt0,cutoffpt,mpt1}],{i,Length[possibleSplits[[;;,1]]]}];
    nozero=Table[Abs[resultsp[[i]]]>10^-5 && Element[resultsp[[i]],Reals],{i,Length[possibleSplits[[;;,1]]]}];
    results = Table[If[nozero[[i]],resultsp[[i]],0],{i,Length[possibleSplits[[;;,1]]]}];
    If[MemberQ[nozero,True],
      splitindex = Ordering[results,-1][[1]];
      Thread[{results[[splitindex]],possibleSplits[[splitindex,2]]}], 
      {{0,parton}}],{{0,parton}}]]


(* ::Text:: *)
(*Here you can find a simple description of how the method Shower works*)
(*1) It first creates two list mylistnsplittings and nquarks, which have both the length of the number of events (or showers) generated. *)
(*     Each slot of the list will contain the number of splittings and of qqbar pairs created in that event. *)
(*2) For each event then, it initialize a list called descendants to a single gluon. The descendant list is made event by event and in each moment *)
(*     represent the entire composition of the shower. In particular, it includes the identity (quark or gluon) and the pT.*)
(*3) A while-loop is then used to perform the shower. The shower continues until all the parton in the shower have pT smaller than cutoffpt.*)
(*    The line MemberQ[Thread[descendants[[;;,1]]<cutoffpt], False] checks if there is at least one parton with pt > cutoffpt. In that case it continues showering and in  *)
(*    particular calling the function "SingleSplittingPtOrdered". Lets see in some detail what the While loop does:*)
(*    - descendants= {{1.5,g}, {0.1, g}}*)
(*    - descendants[[;;,1]] access the first slot for each sublist and gives  {1.5,0.1}*)
(*    - Thread[descendants[[;;,1]]<1.]  gives {False,True}*)
(*    - MemberQ[Thread[descendants[[;;,1]]<1.], False] would give True as an answer, cause the first gluon has 1.5 > 1 GeV*)
(* 4) When the shower is over, it counts the numbers of quarks and gluons and outputs it.*)


Shower[maxevents_,partoninit_:100.,cutoffpt_:1.,sudakovgtogg_:SudakovPgtoggvacuumLT, sudakovgtoqqbar_:SudakovPgtoqqbarvacuumLT,dodebug_:1]:=
Module[{i,j},
  mylistnsplittings = Table[0,{i,maxevents}];
  nquarks=Table[0,{i,maxevents}];
  For[i=0, i<maxevents,i++,
    If[dodebug==1,Print["Generating event number ",i]; If[Mod[i,1]==0,Print["Generating event number ",i]];];
    descendants={{partoninit,"g"}};
    If[dodebug==1, Print["descendants shower initialized to ",descendants];];
    While[MemberQ[Thread[descendants[[;;,1]]<cutoffpt],False],
          descendants=Flatten[Table[SingleSplittingPtOrdered[sudakovgtogg,sudakovgtoqqbar,descendants[[j,1]],descendants[[j,2]]],{j,Length[descendants]}],1];
    ];
    mylistnsplittings[[i+1]]=Length[descendants];
    nquarks[[i+1]]=Count[descendants[[;;,2]],"q"];];
  Print["The shower has produced ",nquarks," qqbar pairs in ",mylistnsplittings, " splittings" ]
]


SingleSplittingPtOrderedFixed[sudakovgtogg_, sudakovgtoqqbar_,mpt1_,parton_:"g",cutoffpt_:1.]:=
Module[{possibleSplits},
  output= {{mpt1,parton}};
  If[mpt1>cutoffpt  &&  parton == "g",(*if condition*)
     possibleSplits={{sudakovgtogg,{"g","g"}}, {sudakovgtoqqbar,{"q","q"}}};
     rnd = RandomReal[{0.,1.0},WorkingPrecision->4]&/@possibleSplits;
     resultsp=Table[pt0//.FindRoot[(possibleSplits[[i,1]]//.{pt1->mpt1,CA->3,\[Alpha]s->0.12 })-rnd[[i]],{pt0,cutoffpt,mpt1}],{i,Length[possibleSplits[[;;,1]]]}];
     nozero=Table[Abs[resultsp[[i]]]>10^-5 && Element[resultsp[[i]],Reals],{i,Length[possibleSplits[[;;,1]]]}];
     results = Table[If[nozero[[i]],resultsp[[i]],0],{i,Length[possibleSplits[[;;,1]]]}];
     processindex = Ordering[results,-1][[1]];
     ptprimsplittee = resultsp[[processindex]];
     typesplittee = possibleSplits[[processindex,2]];
     If[MemberQ[nozero,True],(*if condition*)
     output={{ptprimsplittee,typesplittee[[1]]},{mpt1-ptprimsplittee,typesplittee[[1]]}}](*end of if on MemberQ[nozero,True*)
  ]; (*done if mpt1>cutoffpt  &&  parton == "g" is fullfilled*)
  output
]; (*end of function-module*)


ShowerFixed[partoninit_:100.,cutoffpt_:1.,sudakovgtogg_:SudakovPgtoggvacuumLT, sudakovgtoqqbar_:SudakovPgtoqqbarvacuumLT,dodebug_:1]:=(
  descendants={{partoninit,"g"}};
  iter=0;
  While[MemberQ[Thread[descendants[[;;,1]]<cutoffpt && iter<60],False],
    iter = iter + 1;
    descendants=Flatten[Table[SingleSplittingPtOrderedFixed[sudakovgtogg,sudakovgtoqqbar,descendants[[j,1]],descendants[[j,2]]],{j,Length[descendants]}],1];
  ];
  descendants
);


SingleSplittingPtOrderedFixedUrs[sudakovgtogg_, sudakovgtoqqbar_,fsplitgtogg_,fsplitgtoqqbar_, tscaleinit_:100, parton_:"g",zinit_:1,tscalecutoff_:1.]:=
Module[{possibleSplits,tscalesplitting,splittingfunction},
  output= {{tscaleinit,parton,zinit}};
  If[tscaleinit>tscalecutoff  &&  parton == "g",(*if condition*)
     possibleSplits={{sudakovgtogg,{"g","g"}}, {sudakovgtoqqbar,{"q","q"}}};
     rnd = RandomReal[{0.,1.0},WorkingPrecision->4]&/@possibleSplits;
     resultsp=Table[pt0//.FindRoot[(possibleSplits[[i,1]]//.{pt1->tscaleinit,CA->3,\[Alpha]s->0.12 })-rnd[[i]],{pt0,tscalecutoff,tscaleinit}],{i,Length[possibleSplits[[;;,1]]]}];
     nozero=Table[Abs[resultsp[[i]]]>10^-5 && Element[resultsp[[i]],Reals],{i,Length[possibleSplits[[;;,1]]]}];
     results = Table[If[nozero[[i]],resultsp[[i]],0],{i,Length[possibleSplits[[;;,1]]]}];
     processindex = Ordering[results,-1][[1]];
     tscalesplitting = resultsp[[processindex]];
     typesplittee = possibleSplits[[processindex,2]];
     If [processindex==1, splittingfunction=fsplitgtogtg,splittingfunction=fsplitgtoqqbar];
     If[MemberQ[nozero,True],(*if condition*)
     distrib = ProbabilityDistribution[splittingfunction[z], {z, tscalecutoff/tscaleinit, 1.},Method -> "Normalize"];
     zvalue=RandomVariate[distrib,1][[1]];
     output={{tscalesplitting,typesplittee[[1]],zinit*zvalue},{tscalesplitting,typesplittee[[1]],zinit*(1-zvalue)}}](*end of if on MemberQ[nozero,True*)
  ]; (*done if mpt1>cutoffpt  &&  parton == "g" is fullfilled*)
  output
]; (*end of function-module*)


ShowerFixedUrs[tscaleinit_:100.,tscalecutoff_:1.,sudakovgtogg_:SudakovPgtoggvacuumLT, sudakovgtoqqbar_:SudakovPgtoqqbarvacuumLT,fsplitgtogg_:PgtoggvacuumNT,fsplitgtoqqbar_:PgtoqqbarvacuumLT, dodebug_:1]:=(
  descendants={{tscaleinit,"g",1.}}; 
  iter=0;Print["Here=",descendants];
  While[MemberQ[Thread[descendants[[;;,1]]<tscalecutoff && iter<60],False],
  iter = iter + 1;
  descendants=Flatten[Table[SingleSplittingPtOrderedFixedUrs[sudakovgtogg,sudakovgtoqqbar,fsplitgtogg,fsplitgtoqqbar,descendants[[j,1]],descendants[[j,2]],descendants[[j,3]]],{j,Length[descendants]}],1];];
  descendants
)


(* ::Text:: *)
(*Some instructions on the structure of the shower utility. *)
(*- mylistnsplittings is a list that is initialized to {0} before the shower process starts and contain the complete splitting history*)
(*- nquarks is a list that is initialized to {0} before the shower process starts and contain the number of quarks in the shower*)


(* ::Text:: *)
(*--------------------------------------------------------------------------------------------------------------------------------------------------------*)
(*Plotting utilities*)


PlotCompareSudakovThetaQ2[inputpthigh_,inputCA_,input\[Alpha]s_]:=Module[{functNT,functNTQ2,functNTQ2medium},
functNT= SudakovPgtoggvacuumNT//.{pt1->inputpthigh, CA->inputCA, \[Alpha]s->input\[Alpha]s};
functNTQ2= SudakovPgtoggvacuumNTQ2//.{pt1->inputpthigh, CA->inputCA, \[Alpha]s->input\[Alpha]s};
functNTQ2medium= SudakovPgtoggvacuumNTQ2Medium//.{pt1->inputpthigh, CA->inputCA, \[Alpha]s->input\[Alpha]s};
Plot[{functNT,functNTQ2,functNTQ2medium},{pt0,1,100},Frame->True,FrameStyle->Black,PlotStyle -> {Thickness[0.03], Thickness[0.015],Thickness[0.01]}, PlotLegends->{"P(g to gg)= 2*(1-z)/z +z*(1-z)","P(g to gg)= 2*(1-z)/z +z*(1-z) \!\(\*SuperscriptBox[\(Q\), \(2\)]\)-integrated", "P(g to gg)= 2*(1-z)/z +z*(1-z)+1/\!\(\*SuperscriptBox[\(Q\), \(2\)]\)-integrated \!\(\*SuperscriptBox[\(Q\), \(2\)]\)"}]]


PlotCompareSudakovWithGavin[inputpthigh_,inputCA_,input\[Alpha]s_]:=Module[{functGavin,functLT,functNTnopol,functNT},
functGavin= SudakovPgtoggvacuumGavinOrig//.{pt1->inputpthigh, CA->inputCA, \[Alpha]s->input\[Alpha]s};
functLT= SudakovPgtoggvacuumLT//.{pt1->inputpthigh, CA->inputCA, \[Alpha]s->input\[Alpha]s};
functNTnopol = SudakovPgtoggvacuumNTnopol//.{pt1->inputpthigh, CA->inputCA, \[Alpha]s->input\[Alpha]s};
functNT = SudakovPgtoggvacuumNT//.{pt1->inputpthigh, CA->inputCA, \[Alpha]s->input\[Alpha]s};
Plot[{functGavin,functLT,functNTnopol,functNT},{pt0,1,inputpthigh},Frame->True, AxesLabel->{"p_{T,0}","Sudakov (no splitting prob.)"}, 
                                                       PlotStyle -> {Thickness[0.03], Thickness[0.015],Thickness[0.01],Thickness[0.01]}, FrameStyle->Black, 
                                                       PlotLegends->{"exp{-2*\[Alpha]s*\!\(\*SubscriptBox[\(C\), \(A\)]\)/\[Pi]*Log(\!\(\*SubscriptBox[\(p\), \(T, 1\)]\)/\!\(\*SubscriptBox[\(p\), \(T, 0\)]\)\!\(\*SuperscriptBox[\()\), \(2\)]\)} from Gavin ","P(g to gg)= 2*1/z","P(g to gg)= 2*(1-z)/z", "P(g to gg)= 2*(1-z)/z +z*(1-z)"}]]


(* ::Text:: *)
(*--------------------------------------------------------------------------------------------------------------------------------------------------------*)
(*IMPORTANT: the functions below are not explicitly used for the package, but are provided for reference*)
(*- MakeShower is the original implementation included by Gavin in his tutorial.*)
(*- the basic idea here is that one calculates the pt of the emission (which is in this approach the primary emission)*)
(*  and keep iterating until the emission pt goes below a given pt threshold. *)


MakeShowerSimple[ptCut_,maxevents_,seed_,ptinit_,input\[Alpha]s_,myCA_,Sudakov_]:=(
  mylistnsplittings = {};
  totlistz = {};
  totlistzjet = {};
  SeedRandom[seed];
  RandomReal[];
  randmin=0.;
  randmax=1.;
  Print["Initialization"];
  For[i=0, i<maxevents,i++,
    condition = 1; 
    nsplittings=0;
    evtlistz={};
    inter=0;
    sudakov=1;
    valuezjet=1;
    While[condition == 1 && inter<100, 
      rnd = RandomReal[{randmin,randmax},WorkingPrecision->4];
      ptbeforesplit = pt0//. FindRoot[sudakforcalc-sudakov,{pt0,ptCut,ptinit}];
      sudakov = sudakov * rnd;
      sudakforcalc = Sudakov//. {CA->myCA,pt1->ptinit,\[Alpha]s->input\[Alpha]s};
      pt = pt0//. FindRoot[sudakforcalc-sudakov,{pt0,ptCut,ptinit}];
      If[pt<ptCut, 
         AppendTo[mylistnsplittings,nsplittings];
         AppendTo[totlistzjet,valuezjet];
         AppendTo[totlistz,evtlistz]; 
         condition=0;
      ];
      If[condition==1, 
         nsplittings=nsplittings+1;
         AppendTo[evtlistz,pt/ptbeforesplit];
         valuezjet=valuezjet*pt/ptbeforesplit;
         inter = inter+1;
      ];
    ];
  ];
);


Begin["`Private`"]


End[];


EndPackage[];
