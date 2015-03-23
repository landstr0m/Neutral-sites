(* ::Package:: *)

(* ::Section:: *)
(*Introduction*)


BeginPackage["ProteomicsAnalysis`", {"HierarchicalClustering`"}];


(* ::Subsection:: *)
(*Help messages*)


TrypsinDigestion::usage = "TrypsinDigestion[P] or TrypsinDigestion[P, m] gives \
a list of the peptides resulting from the tryptic digestion of P.  The parameter m \
is the number of missed cleavages. It is a \"Listable\" function too.";


PrecursorMass::usage = "";


ProteomicsDendrogram::usage = "";


ProteomicsHeatmap::usage = "";


Begin["Private`"];


(* ::Section:: *)
(*Code*)


(* ::Subsection:: *)
(*Precursor mass calculations*)


PrecursorMass[p_String, q_:1] := (Total[(Characters[p] /. $AminoAcidMasses)] +
	  18 + q - StringLength[p] + 46 * StringCount[p, "C"])/q


$AminoAcidMasses = {
	"A"->72.04498`,   "C"->104.01706`,
	"D"->116.03481`,  "E"->130.05046`,
	"F"->148.07628`,  "G"->58.02933`,
	"H"->138.06678`,  "I"->114.09193`,
	"K"->129.10283`,  "L"->114.09193`,
	"M"->132.04836`,  "N"->115.0508`,
	"P"->98.06063`,   "Q"->129.06645`,
	"R"->157.10898`,  "S"->88.0399`,
	"T"->102.05555`,  "V"->100.07628`,
	"W"->187.08718`,  "Y"->164.0712`,
	_String->Indeterminate};


(* ::Subsection:: *)
(*Trypsin Digestion*)


TrypsinDigestion[s_String, n_Integer:0] :=
  With[{p = Join[{0}, DeleteCases[StringPosition[s, "R" | "K" ~~
          Except["P"]][[All, 1]], 1], {StringLength[s]}]}, 
     Flatten[Table[StringTake[s, p[[{i, i + m}]] + {1, 0}], {m, n + 1}, 
       {i, Length[p] - m}], 1]]


TrypsinDigestion[sl_List, n_Integer:0] := Map[TrypsinDigestion[#, n] &, sl]


PDFext[fn_] := fn <> If[! StringMatchQ[fn, __ ~~ ".pdf" ~~ EndOfString], ".pdf", ""]


exportButton[g_] :=
	Row[{"Export:", Button["PDF", Export[
		SystemDialogInput["FileSave", WindowTitle -> "Export PDF"] // PDFext,
		g, "PDF"], ImageSize -> Automatic]}]


Options[ProteomicsDendrogram] = {
	PlotLabel -> Automatic
};


Options[ProteomicsHeatmap] = {
	"Labels" -> Automatic
};


Options[heatmap] = {
	"Labels" -> Automatic
};


color[x_, {m_, M_}] := If[x < 0, RGBColor[0, x/m, 0], RGBColor[x/M, 0, 0]]


heatmap[d_, All, opts:OptionsPattern[]] := Module[{min, max},
	min = Min[d];
	max = Max[d];
	Graphics[
		Table[Text[OptionValue["Labels"][[j]], {j + 1/2, 0}],
			{j, Length[OptionValue["Labels"]]}]
		~Join~
		Table[{color[d[[i, j]], {min,max}],
			Rectangle[{j, -i}, {j + 1, -i - 1}]},
			{i, Length[d]}, {j, Length[d[[1]]]}],
		AspectRatio -> Length[d]/40, ImageSize -> 300]]


heatmap[d_, Row, opts:OptionsPattern[]]:= Module[{min, max},
	Graphics[
		Table[Text[OptionValue["Labels"][[j]], {j + 1/2, 0}],
			{j, Length[OptionValue["Labels"]]}]
		~Join~
		Table[
			min = Min[d[[i]]];
			max=Max[d[[i]]];
			{color[d[[i, j]], {min, max}], Rectangle[{j, -i}, {j + 1, -i - 1}]},
			{i, Length[d]}, {j, Length[d[[1]]]}],
		AspectRatio -> Length[d]/40, ImageSize -> 300]]


heatmap[d_, Column, opts:OptionsPattern[]] := Module[{min, max},
	min = Map[Min, d // Transpose];
	max = Map[Max, d // Transpose];
	Graphics[
		Table[Text[OptionValue["Labels"][[j]], {j + 1/2, 0}],
			{j, Length[OptionValue["Labels"]]}]
		~Join~
		Table[{color[d[[i, j]], {min[[j]], max[[j]]}],
			Rectangle[{j, -i}, {j + 1, -i - 1}]},
			{i, Length[d]}, {j, Length[d[[1]]]}],
		AspectRatio -> Length[d]/40, ImageSize -> 300]]


ProteomicsDendrogram[data_List, OptionsPattern[]] := Module[{g, nlist, pl},
	nlist = Union[Quantile[data[[2 ;; ,1]], Range[0.1,0.9,0.1]] // Round];
	pl = OptionValue[PlotLabel];
	Manipulate[
		With[{selected = Select[data[[2 ;; ]], First[#] > n && (And @@ (NumericQ /@ Rest[#])) &]},
			g = Show[
				DendrogramPlot[
					Transpose[selected[[All, 2 ;; ]]],
					Linkage -> m,
					LeafLabels -> data[[1, 2 ;; ]],
					PlotLabel -> If[pl === Automatic,
						"Number of proteins: " <> ToString[Length[selected]], pl]],
				Axes -> {False, True},
				AxesOrigin -> {0, 0},
				PlotRange -> All,
				PlotRangePadding -> {{0, 0}, {0, Scaled[0.1]}}]],
		{n, nlist, ControlType -> SetterBar},
		{m, {"Ward", "Single", "Average", "Complete", "WeightedAverage", "Centroid", "Median"}},
		Item[exportButton[g]],
		SaveDefinitions -> True,
		TrackedSymbols :> True]]


ProteomicsHeatmap[data_List, OptionsPattern[]] := Module[
	{xdata, nlist, g, l},

	nlist = Union[Quantile[data[[2 ;; ,1]], Range[0.1,0.9,0.1]] // Round];
	xdata = Cases[data[[2 ;; ]], {_, __Real}];
	l = OptionValue["Labels"];
	If[l === Automatic, l = data[[1, 2 ;; ]],
		If[l === None, l = {}]];

	Manipulate[
		Module[{selected, xy, st, ind, rules},
			selected = Select[xdata, First[#] > n &];
			xy = Map[# - Mean[#] &, selected[[All, 2 ;; ]]];
			rules = Thread[xy -> (pos /@ Range[Length[xy]])];
			ind = Cases[Agglomerate[xy, Linkage -> "Ward"] /. rules, pos[_], +\[Infinity]]
					/. pos[x_] :> x;
			st = FindShortestTour[xy[[ind]]] // Last;
			g  = heatmap[xy[[ind]][[st]], t, "Labels" -> l]],

		{{n, nlist[[-1]], "n"}, nlist, ControlType -> SetterBar},
		{{t, Column}, {All, Column, Row}},
		Item[exportButton[g]],
		TrackedSymbols :> True,
		SynchronousUpdating -> False,
		SaveDefinitions -> True]
]


End[];


EndPackage[];
