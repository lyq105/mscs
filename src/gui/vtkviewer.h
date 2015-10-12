#ifndef VTKVIEWER_H
#define VTKVIEWER_H


class QVTKWidget;
class VTKActor;

class VTKviewer
{
public:
    VTKviewer();
    ~VTKviewer();

    void render();
    QVTKWidget* getVTKWidget(){return qvtk;}

private:
    QVTKWidget* qvtk;
};

#endif // VTKVIEWER_H
