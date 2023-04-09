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

$!ActiveLineMaps += [2]
$!ActiveLineMaps += [3]
$!LineMap [1]  Assign{YAxisVar = 3}
$!LineMap [1]  Assign{Zone = 1}
$!LineMap [2]  Assign{YAxisVar = 3}
$!LineMap [2]  Assign{Zone = 2}
$!LineMap [3]  Assign{YAxisVar = 3}
$!LineMap [3]  Assign{Zone = 3}

$!LineMap [1]  Name = 'M1'
$!LineMap [2]  Name = 'P1'
$!LineMap [3]  Name = 'P3'
$!View Fit
$!GlobalLinePlot Legend{Show = Yes}
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
$!GlobalLinePlot Legend{Box{BoxType = None}}
$!View Fit
