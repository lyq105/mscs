#ifndef PMMESHER_H
#define PMMESHER_H
#include "Mesher.h"

class PMMesher : public Mesher
{
public:
    PMMesher(double ms, int c, int ln, double tr, CMCell* cell);

    ///  生成颗粒增强单胞网格
    int generate_mesh();
private:
    double mesh_size;
    int coat;
    int layer_num;
    double thick_ratio;
};

#endif // PMMESHER_H
