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


SingleSplittingPtOrderedValidated[sudakovgtogg_, sudakovgtoqqbar_,fsplitgtogg_,fsplitgtoqqbar_, tscaleinit_:100, parton_:"g",zinit_:1,tscalecutoff_:1.,dodebug_:1]:=
Module[{possibleSplits,tscalesplitting,splittingfunction,output,zlowcutoff},
  If[dodebug==1,Print["SingleSplittingPtOrderedValidated::Debug, initial scale=",tscaleinit,", parton=",parton," , zinit=",zinit];];
  If[parton == "q",output={{tscalecutoff,parton,zinit}}];
  If[parton == "g" && tscaleinit==tscalecutoff,output={{tscaleinit,parton,zinit}};];
  If[tscaleinit>tscalecutoff  &&  parton == "g",(*if condition*)
     possibleSplits={{sudakovgtogg,{"g","g"},fsplitgtogg,tscalecutoff/tscaleinit}, {sudakovgtoqqbar,{"q","q"},fsplitgtoqqbar,tscalecutoff/tscaleinit}};
     rnd = RandomReal[{0.,1.0},WorkingPrecision->4]&/@possibleSplits;
     resultsp=Table[pt0//.FindRoot[(possibleSplits[[i,1]]//.{pt1->tscaleinit,CA->3,\[Alpha]s->0.12 })-rnd[[i]],{pt0,tscalecutoff,tscaleinit}],{i,Length[possibleSplits[[;;,1]]]}];
     nozero=Table[Abs[resultsp[[i]]]>10^-5 && Element[resultsp[[i]],Reals],{i,Length[possibleSplits[[;;,1]]]}];
     results = Table[If[nozero[[i]],resultsp[[i]],0],{i,Length[possibleSplits[[;;,1]]]}];
     (*processindex = Ordering[results,-1][[1]];*)
     processindex = 1;
     tscalesplitting = resultsp[[processindex]];
     typesplittee = possibleSplits[[processindex,2]];
     splittingfunction=possibleSplits[[processindex,3]];
     funct[z_]:=(1-z)/z+z/(1-z)+z*(1-z);
     If[tscalesplitting>tscalecutoff,
       zlowcutoff=tscalecutoff/tscalesplitting;
       (*Print[zlowcutoff];
       Print[tscalesplitting];
       Print[zlowcutoff];*)
       lowb=100.;
       highb=100.;
       If[zlowcutoff>=0.5,lowb=1-zlowcutoff; highb=zlowcutoff];
       If[zlowcutoff<0.5,lowb=zlowcutoff; highb=1-zlowcutoff];
       If [zlowcutoff<0., Print["ERROR!!! z boundary for extraction is lower than 0"]];
       distrib = ProbabilityDistribution[funct[z], {z, lowb, highb},Method -> "Normalize"];
       zvalue=RandomVariate[distrib,1][[1]];
       If [zvalue<0. || zvalue>1., Print["ERROR!!! z value extracted is not in the correct boundaries [0,1], z=", zvalue]];
       output={{tscalesplitting,typesplittee[[1]],zinit*zvalue},{tscalesplitting,typesplittee[[1]],zinit*(1-zvalue)}};
     ];
     If[tscalesplitting<=tscalecutoff,
        output={{tscalecutoff,"g",zinit}}];
     If[dodebug==1,Print["SingleSplittingPtOrderedValidated::Debug, list of descendants=",descendants];];
  ]; (*done if mpt1>cutoffpt  &&  parton == "g" is fullfilled*)
  output
]; 


ShowerValidated[tscaleinit_:100.,tscalecutoff_:1.,sudakovgtogg_:SudakovPgtoggvacuumLT, sudakovgtoqqbar_:SudakovPgtoqqbarvacuumLT,fsplitgtogg_:PgtoggvacuumNT,fsplitgtoqqbar_:PgtoqqbarvacuumLT, dodebug_:1,maxiteration_:200]:=(
  descendants={{tscaleinit,"g",1.}}; 
  iter=0;
  If[dodebug==1,Print["Starting from=",descendants];];
  While[MemberQ[Thread[descendants[[;;,1]]<=tscalecutoff && iter<maxiteration],False],
  If[dodebug==1,Print["Descendant list at iteration="iter," is= ",descendants];];
  If[iter==maxiteration-1, Print["Check for an infinite loop or a very high-multiplicity event!!!"]];
  iter = iter + 1;
  descendants=Flatten[Table[SingleSplittingPtOrderedValidated[sudakovgtogg,sudakovgtoqqbar,fsplitgtogg,fsplitgtoqqbar,descendants[[j,1]],descendants[[j,2]],descendants[[j,3]],tscalecutoff,dodebug],{j,Length[descendants]}],1];];
  nquarks=Count[descendants[[;;,2]],"q"];
  If [Mod[nquarks,2]==1, Print["This event has a odd number of quarks= ", nquarks];];
  If[dodebug==1,Print["Final list="iter," is= ",descendants];];
  descendants
)


RunShowerMulti[maxevents_:100.,tscaleinit_:100.,tscalecutoff_:1.,sudakovgtogg_:SudakovPgtoggvacuumLT, sudakovgtoqqbar_:SudakovPgtoqqbarvacuumLT,fsplitgtogg_:PgtoggvacuumNT,fsplitgtoqqbar_:PgtoqqbarvacuumLT, dodebug_:1,maxiteration_:200,path_]:=(
  descendtot = {};
  zvaluestot = {};
  nmultiplicitytot = {};
  nmultiplicityquarkstot = {};
  nquarkstot = {};
  For[i = 0, i < maxevents, i++,
    Clear[descendentevent, zvaluesevent]; 
    descendentevent = ShowerValidated[tscaleinit, tscalecutoff, SudakovPgtoggvacuumLT, SudakovPgtoqqbarvacuumLT, PgtoggvacuumLT, PgtoqqbarvacuumLT, 0, 200];
    zvaluesevent = Table[descendentevent[[jendex]][[3]], {jendex, Length[descendentevent]}];
    nquarks = Count[descendentevent[[;; , 2]], "q"];
    AppendTo[descendtot, descendentevent];
    AppendTo[zvaluestot, zvaluesevent];
    AppendTo[nmultiplicitytot, Length[descendentevent]];
    AppendTo[nmultiplicityquarkstot, nquarks];
  ];  
  PlotShowerQuantities[nmultiplicitytot, nmultiplicityquarkstot, zvaluestot,tscaleinit];
  Export[StringJoin[path,StringTemplate["histomultinit_maxevents`1`_tscalecutoff`2`_scale`3`GeVc.pdf"][maxevents,tscalecutoff,tscaleinit]], histomult]
  Export[StringJoin[path,StringTemplate["histonquarks_maxevents`1`_tscalecutoff`2`_scale`3`GeVc.pdf"][maxevents,tscalecutoff,tscaleinit]], histonquarks]
  Export[StringJoin[path,StringTemplate["histolog1overz_maxevents`1`_tscalecutoff`2`_scale`3`GeVc.pdf"][maxevents,tscalecutoff,tscaleinit]], histolog1overz]
  histomult
  histonquarks
  histolog1overz
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


PlotShowerQuantities[nmultiplicitytot_,nmultiplicityquarkstot_,zvaluestot_,ptinitial_]:=(
histomult= Histogram[nmultiplicitytot,{-0.5,80.5,1.},"Probability", Frame -> True, AxesLabel->{HoldForm["Parton multiplicity"],
          HoldForm[Entries]}, PlotLabel->StringTemplate["Distribution of parton multiplicity for initial gluon with t=`1` GeV"][ptinitial],
          LabelStyle->{FontFamily->"Helvetica", 12, GrayLevel[0]},ImageSize->Large];
histonquarks= Histogram[nmultiplicityquarkstot,{-0.5,10.5,1.},"Probability", AxesLabel->{HoldForm["Quark multiplicity"],
          HoldForm[Entries]}, PlotLabel->StringTemplate["Distribution of quark multiplicity for initial gluon of with t=`1` GeV"][ptinitial],
          LabelStyle->{FontFamily->"Helvetica", 12, GrayLevel[0]},ImageSize->Large];
histolog1overz = Histogram[Log[1/Flatten[zvaluestot]],{0.,30.,0.5},"Probability", AxesLabel->{HoldForm[Log[1/x]],HoldForm[Entries]}, 
                           PlotLabel->StringTemplate["Distribution of Log[\!\(\*FractionBox[\(1\), \(x\)]\)] for initial gluon with t=`1` GeV"][ptinitial],
                           LabelStyle->{FontFamily->"Helvetica", 12, GrayLevel[0]},ImageSize->Large];
)


(* ::Text:: *)
(*Utilities*)


SilentErrors[]:=(
Off[RandomReal::precw];
Off[Greater::nord];
Off[LessEqual::nord];
Off[FindRoot::cvmit];
Off[FindRoot::frmp];
Off[FindRoot::sszero];
Off[FindRoot::lstol];
Off[Part::partd];
Off[General::stop];
Off[General::szero];
Off[General::ovfl];)


SetStandardSeeds[]:=(
SeedRandom[1234];
RandomReal[];
)


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
