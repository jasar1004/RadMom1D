#!MC 1410
$!ReadDataSet  '"./M1/test1D.dat" '
  ReadDataOption = New
  ResetStyle = Yes
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"x" "kp" "E" "Fx" "Sr"'
$!ReadDataSet  '"./P1/test1D.dat" '
  ReadDataOption = Append
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"x" "kp" "E" "Fx" "Sr"'
$!ReadDataSet  '"./P3/test1D.dat" '
  ReadDataOption = Append
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"x" "kp" "E" "Fx" "Sr"'

$!AlterData
  Equation = '{E (kW m^-2 sr^-1)} = {E}/1000'

$!ActiveLineMaps += [2]
$!ActiveLineMaps += [3]
$!LineMap [1]  Assign{YAxisVar = 6}
$!LineMap [1]  Assign{Zone = 1}
$!LineMap [2]  Assign{YAxisVar = 6}
$!LineMap [2]  Assign{Zone = 2}
$!LineMap [3]  Assign{YAxisVar = 6}
$!LineMap [3]  Assign{Zone = 3}

$!LineMap [1]  Name = 'M1'
$!LineMap [2]  Name = 'P1'
$!LineMap [3]  Name = 'P3'
$!View Fit

$!Pick AddAtPosition
  X = 1.01060358891
  Y = 0.615415986949
  ConsiderStyle = Yes
$!FrameControl ActivateByNumber
  Frame = 1
$!Pick AddAtPosition
  X = 1.01060358891
  Y = 0.615415986949
  CollectingObjectsMode = HomogeneousAdd
  ConsiderStyle = Yes
$!FrameLayout ShowBorder = No
$!Pick AddAtPosition
  X = 7.91435562806
  Y = 3.10807504078
  ConsiderStyle = Yes
$!Pick AddAtPosition
  X = 7.9013050571
  Y = 3.10807504078
  CollectingObjectsMode = HomogeneousAdd
  ConsiderStyle = Yes

$!LineMap [1]  Lines{LinePattern = Solid}
$!LineMap [2]  Lines{LinePattern = LongDash}
$!LineMap [3]  Lines{LinePattern =  Dashed}

$!LinePlotLayers ShowSymbols = Yes

$!LineMap [1]  Symbols{SymbolShape{GeomShape = Square}}
$!LineMap [2]  Symbols{SymbolShape{GeomShape = Diamond}}
$!LineMap [3]  Symbols{SymbolShape{GeomShape = Grad}}
$!LineMap [1-3]  Symbols{Size = 1.2}

$!LineMap [1-3]  Lines{LineThickness = 0.8}

$!LineMap [1-3]  Symbols{Skipping = 10}

$!FrameLayout ShowBorder = No

$!GlobalLinePlot Legend{Show = Yes}
$!GlobalLinePlot Legend{Box{BoxType = None}}
$!GlobalLinePlot Legend{TextShape{Height = 2}}
$!XYLineAxis XDetail 1 {Title{TitleMode = UseText}}
$!XYLineAxis XDetail 1 {Title{Offset = 5}}
$!XYLineAxis XDetail 1 {Title{Text = 'x (m)'}}
$!XYLineAxis YDetail 1 {Title{TitleMode = UseText}}
$!XYLineAxis YDetail 1 {Title{Offset = 9}}
$!XYLineAxis YDetail 1 {Title{Text = 'E (kW.m-2)'}}

$!XYLineAxis GridArea{DrawBorder = Yes}
$!XYLineAxis GridArea{LineThickness = 0.20}
$!XYLineAxis XDetail 1 {Ticks{LineThickness = 0.20}}
$!XYLineAxis YDetail 1 {Ticks{LineThickness = 0.20}}
$!XYLineAxis YDetail 1 {AxisLine{LineThickness = 0.20}}
$!XYLineAxis XDetail 1 {AxisLine{LineThickness = 0.20}}
$!XYLineAxis XDetail 1 {Ticks{ShowOnGridBorderMax = Yes}}
$!XYLineAxis YDetail 1 {Ticks{ShowOnGridBorderMax = Yes}}
$!XYLineAxis XDetail 1 {Title{TextShape{IsBold = No}}}
$!XYLineAxis XDetail 1 {Title{TextShape{IsItalic = Yes}}}
$!XYLineAxis YDetail 1 {Title{TextShape{IsBold = No}}}
$!XYLineAxis YDetail 1 {Title{TextShape{IsItalic = Yes}}}

$!XYLineAxis XDetail 1 {TickLabel{TextShape{Height = 4}}}
$!XYLineAxis YDetail 1 {TickLabel{TextShape{Height = 4}}}
$!XYLineAxis YDetail 1 {Title{TextShape{Height = 4.6}}}
$!XYLineAxis XDetail 1 {Title{TextShape{Height = 4.6}}}
$!GlobalLinePlot Legend{TextShape{Height = 4}}

$!RedrawAll
