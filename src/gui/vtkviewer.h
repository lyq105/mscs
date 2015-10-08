#ifndef VTKVIEWER_H
#define VTKVIEWER_H


class QVTKWidget;
class VTKActor;

class VTKviewer
{
public:
    VTKviewer();
    void render();

private:
    QVTKWidget* qvtk;
};

#endif // VTKVIEWER_H
