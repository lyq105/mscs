######################################################################
# Automatically generated by qmake (2.01a) ?? 9? 13 17:45:07 2015
######################################################################

TEMPLATE = app
TARGET = 
DEPENDPATH += . forms resources src src/gui
INCLUDEPATH += . src/gui src/pmcell src/elastic
unix:INCLUDEPATH +=/usr/include/vtk-5.8/
unix:LIBS = -lQVTK -lvtkRendering -lvtkCommon -lvtkWidgets -lvtkHybrid -lvtkIO -lvtkFiltering -lvtkGraphics
unix:QMAKE_CXXFLAGS += -Wno-deprecated -Wno-unused-parameter

win32:INCLUDEPATH += "C:/Program Files/VTK_vs2008/include/vtk-5.10/"
win32:LIBS +=  "C:/Program Files/VTK_vs2008/lib/vtk-5.10/QVTK.lib" \
        "C:/Program Files/VTK_vs2008/lib/vtk-5.10/vtkRendering.lib" \
        "C:/Program Files/VTK_vs2008/lib/vtk-5.10/vtkWidgets.lib" \
        "C:/Program Files/VTK_vs2008/lib/vtk-5.10/vtkHybrid.lib" \
        "C:/Program Files/VTK_vs2008/lib/vtk-5.10/vtkIO.lib" \
        "C:/Program Files/VTK_vs2008/lib/vtk-5.10/vtkFiltering.lib" \
        "C:/Program Files/VTK_vs2008/lib/vtk-5.10/vtkCommon.lib" \
        "C:/Program Files/VTK_vs2008/lib/vtk-5.10/vtkGraphics.lib"

# Input
HEADERS += src/gui/aboutdialog.h \
           src/gui/analysistype.h \
           src/gui/cellmodel.h \
           src/gui/mainwindow.h \
           src/gui/materiel.h \
           src/gui/meshpara.h \
           src/gui/setcelldialog.h \
           src/gui/sotsinterface.h \
           src/gui/tspara.h \
           src/gui/mylogger.h \
    src/gui/vtkviewer.h \
    src/elastic/AnalyticalSolution.h \
    src/elastic/analyzer.h \
    src/elastic/BeamSecShape.h \
    src/elastic/CanXBeamConZMomAnaSol.h \
    src/elastic/CanYBeamConXMomAnaSol.h \
    src/elastic/CanZBeamConYMomAnaSol.h \
    src/elastic/CMCell.h \
    src/elastic/EllipseGen.h \
    src/elastic/EllipseMade.h \
    src/elastic/EllipseSurface.h \
    src/elastic/Fem.h \
    src/elastic/FixRodXTenAnaSol.h \
    src/elastic/FixRodXTwistAnaSol.h \
    src/elastic/FixRodYTenAnaSol.h \
    src/elastic/FixRodYTwistAnaSol.h \
    src/elastic/FixRodZTenAnaSol.h \
    src/elastic/FixRodZTwistAnaSol.h \
    src/elastic/Gauss.h \
    src/elastic/Geometry.h \
    src/elastic/Gloloaded_vector.h \
    src/elastic/GloStiffMatrix.h \
    src/elastic/Hns.h \
    src/elastic/HomoPara.h \
    src/elastic/HomoSolver.h \
    src/elastic/Matbase.h \
    src/elastic/MathMatrix.h \
    src/elastic/MatPro.h \
    src/elastic/Mesher.h \
    src/elastic/Mscm.h \
    src/elastic/PCMCell.h \
    src/elastic/SolveEqu.h \
    src/elastic/Vector2D.h \
    src/elastic/Vector3D.h \
FORMS += forms/aboutdialog.ui \
         forms/analysistype.ui \
         forms/cellmodel.ui \
         forms/mainwindow.ui \
         forms/materiel.ui \
         forms/meshpara.ui \
         forms/setcelldialog.ui
SOURCES += src/mscs_main.cpp \
           src/gui/aboutdialog.cpp \
           src/gui/analysistype.cpp \
           src/gui/cellmodel.cpp \
           src/gui/mainwindow.cpp \
           src/gui/materiel.cpp \
           src/gui/meshpara.cpp \
           src/gui/setcelldialog.cpp \
           src/gui/sotsinterface.cpp \
           src/gui/tspara.cpp \
    src/gui/vtkviewer.cpp \
    src/elastic/AnalyticalSolution.cpp \
    src/elastic/analyzer.cpp \
    src/elastic/BeamSecShape.cpp \
    src/elastic/CanXBeamConZMomAnaSol.cpp \
    src/elastic/CanYBeamConXMomAnaSol.cpp \
    src/elastic/CanZBeamConYMomAnaSol.cpp \
    src/elastic/CMCell.cpp \
    src/elastic/EllipseGen.cpp \
    src/elastic/EllipseMade.cpp \
    src/elastic/EllipseSurface.cpp \
    src/elastic/Fem.cpp \
    src/elastic/FixRodXTenAnaSol.cpp \
    src/elastic/FixRodXTwistAnaSol.cpp \
    src/elastic/FixRodYTenAnaSol.cpp \
    src/elastic/FixRodYTwistAnaSol.cpp \
    src/elastic/FixRodZTenAnaSol.cpp \
    src/elastic/FixRodZTwistAnaSol.cpp \
    src/elastic/Gauss.cpp \
    src/elastic/Geometry.cpp \
    src/elastic/Gloloaded_vector.cpp \
    src/elastic/GloStiffMatrix.cpp \
    src/elastic/Hns.cpp \
    src/elastic/HomoPara.cpp \
    src/elastic/HomoSolver.cpp \
    src/elastic/Matbase.cpp \
    src/elastic/MathMatrix.cpp \
    src/elastic/MatPro.cpp \
    src/elastic/Mesher.cpp \
    src/elastic/Mscm.cpp \
    src/elastic/PCMCell.cpp \
    src/elastic/SolveEqu.cpp \
    src/elastic/Vector2D.cpp \
    src/elastic/Vector3D.cpp
RESOURCES += resources/Mscs.qrc
