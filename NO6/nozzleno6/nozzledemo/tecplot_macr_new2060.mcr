#!MC 1410
$!VarSet |MFBD| = 'C:\Program Files\Tecplot\Tecplot 360 EX 2015 R2'
$!READDATASET  '"D:\study\lerning\gitlearning\homework\NO6\nozzleno6\nozzledemo\nozzledemo\2060\2.dat" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  VARNAMELIST = '"x" "y" "u" "v" "p" "den" "mach" "entropy" "vmu"'
$!GLOBALRGB REDCHANNELVAR = 3
$!GLOBALRGB GREENCHANNELVAR = 3
$!GLOBALRGB BLUECHANNELVAR = 3
$!SETCONTOURVAR 
  VAR = 5
  CONTOURGROUP = 1
  LEVELINITMODE = RESETTONICE
$!SETCONTOURVAR 
  VAR = 6
  CONTOURGROUP = 2
  LEVELINITMODE = RESETTONICE
$!SETCONTOURVAR 
  VAR = 7
  CONTOURGROUP = 3
  LEVELINITMODE = RESETTONICE
$!SETCONTOURVAR 
  VAR = 8
  CONTOURGROUP = 4
  LEVELINITMODE = RESETTONICE
$!SETCONTOURVAR 
  VAR = 9
  CONTOURGROUP = 5
  LEVELINITMODE = RESETTONICE
$!SETCONTOURVAR 
  VAR = 5
  CONTOURGROUP = 6
  LEVELINITMODE = RESETTONICE
$!SETCONTOURVAR 
  VAR = 5
  CONTOURGROUP = 7
  LEVELINITMODE = RESETTONICE
$!SETCONTOURVAR 
  VAR = 5
  CONTOURGROUP = 8
  LEVELINITMODE = RESETTONICE
$!FIELDLAYERS SHOWCONTOUR = YES
$!REDRAWALL 
$!PICK ADDATPOSITION
  X = 8.5993006993
  Y = 2.72272727273
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!GLOBALCONTOUR 1  LEGEND{ISVERTICAL = NO}
$!REDRAWALL 
$!PICK ADDATPOSITION
  X = 4.33636363636
  Y = 2.21923076923
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
$!REDRAWALL 
$!PICK ADDATPOSITION
  X = 5.71258741259
  Y = 2.2527972028
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!GLOBALCONTOUR 1  LEGEND{AUTORESIZE = YES}
$!REDRAWALL 
$!PICK SHIFT
  X = -3.02097902098
  Y = 0.637762237762
$!PICK ADDATPOSITION
  X = 2.16573426573
  Y = 1.53671328671
  CONSIDERSTYLE = YES
$!PICK ADDATPOSITION
  X = 2.16573426573
  Y = 1.53671328671
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!TWODAXIS YDETAIL{SHOWAXIS = NO}
$!TWODAXIS XDETAIL{SHOWAXIS = NO}
$!REDRAWALL 
$!PICK ADDATPOSITION
  X = 4.64965034965
  Y = 2.99125874126
  CONSIDERSTYLE = YES
$!PICK SHIFT
  X = -0.593006993007
  Y = 0.156643356643
$!PICK ADDATPOSITION
  X = 4.85104895105
  Y = 3.09195804196
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!SETCONTOURVAR 
  VAR = 3
  CONTOURGROUP = 1
  LEVELINITMODE = RESETTONICE
$!REDRAWALL 
$!PICK SHIFT
  X = -0.145454545455
  Y = 0.0223776223776
$!PICK ADDATPOSITION
  X = 4.71678321678
  Y = 1.92832167832
  CONSIDERSTYLE = YES
$!PICK ADDATPOSITION
  X = 4.28041958042
  Y = 2.81223776224
  CONSIDERSTYLE = YES
$!PICK SHIFT
  X = -0.0671328671329
  Y = 0.0559440559441
$!PICK ADDATPOSITION
  X = 4.12377622378
  Y = 1.96188811189
  CONSIDERSTYLE = YES
$!EXPORTSETUP IMAGEWIDTH = 804
$!EXPORTSETUP EXPORTFNAME = 'D:\study\lerning\gitlearning\homework\NO6\nozzleno6\nozzledemo\nozzledemo\2060\u.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!PICK ADDATPOSITION
  X = 4.14615384615
  Y = 3.0472027972
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!SETCONTOURVAR 
  VAR = 4
  CONTOURGROUP = 1
  LEVELINITMODE = RESETTONICE
$!REDRAWALL 
$!EXPORTSETUP EXPORTFNAME = 'D:\study\lerning\gitlearning\homework\NO6\nozzleno6\nozzledemo\nozzledemo\2060\v.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!PICK ADDATPOSITION
  X = 3.2958041958
  Y = 3.05839160839
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!SETCONTOURVAR 
  VAR = 5
  CONTOURGROUP = 1
  LEVELINITMODE = RESETTONICE
$!REDRAWALL 
$!EXPORTSETUP EXPORTFNAME = 'D:\study\lerning\gitlearning\homework\NO6\nozzleno6\nozzledemo\nozzledemo\2060\p.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!PICK ADDATPOSITION
  X = 2.85944055944
  Y = 3.12552447552
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!SETCONTOURVAR 
  VAR = 6
  CONTOURGROUP = 1
  LEVELINITMODE = RESETTONICE
$!REDRAWALL 
$!PICK SHIFT
  X = -0.324475524476
  Y = 0.134265734266
$!PICK ADDATPOSITION
  X = 2.84825174825
  Y = 2.40944055944
  CONSIDERSTYLE = YES
$!EXPORTSETUP EXPORTFNAME = 'D:\study\lerning\gitlearning\homework\NO6\nozzleno6\nozzledemo\nozzledemo\2060\den.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!PICK ADDATPOSITION
  X = 3.92237762238
  Y = 3.23741258741
  CONSIDERSTYLE = YES
$!PICK ADDATPOSITION
  X = 3.92237762238
  Y = 3.23741258741
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!SETCONTOURVAR 
  VAR = 7
  CONTOURGROUP = 1
  LEVELINITMODE = RESETTONICE
$!REDRAWALL 
$!EXPORTSETUP EXPORTFNAME = 'D:\study\lerning\gitlearning\homework\NO6\nozzleno6\nozzledemo\nozzledemo\2060\mach.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!PICK ADDATPOSITION
  X = 3.58671328671
  Y = 3.19265734266
  COLLECTINGOBJECTSMODE = ALWAYSADD
  CONSIDERSTYLE = YES
$!SETCONTOURVAR 
  VAR = 8
  CONTOURGROUP = 1
  LEVELINITMODE = RESETTONICE
$!REDRAWALL 
$!RemoveVar |MFBD|
