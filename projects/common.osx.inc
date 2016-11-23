UTILSRC = ../../src/util/SIMPLE_PARSER.cpp \
          ../../src/util/TIMER.cpp \
          ../../src/util/TIMING_BREAKDOWN.cpp \
          ../../src/util/glvu.cpp

GEOMETRY = ../../src/geometry/TRIANGLE.cpp \
           ../../src/geometry/OBJ.cpp \
           ../../src/geometry/TRIANGLE_BVH.cpp \
           ../../src/dtgrid/SPARSE_SDF.cpp \
           ../../src/deformCD/aabb.cpp \
           ../../src/deformCD/aap.cpp \

INTEGRATOR = ../../src/integrator/DEFOR_SOLID_3D.cpp \
             ../../src/integrator/COUPLED_INTEGRATOR_3D.cpp

MATERIAL = ../../src/material/COROTATION.cpp

CC = /usr/bin/clang++

GLFLAGS = -L/usr/lib/ -framework GLUT -framework OpenGL -framework Cocoa -framework Accelerate
PNGFLAGS = -L/opt/local/lib/ -lpng
ZLIBFLAGS = -L/opt/local/lib/ -lz

CPPFLAGS = -m64 -O2 -c -Wno-deprecated-declarations -Wno-deprecated-register -std=c++0x \
            -I/opt/local/include \
            -I/usr/local/include \
            -I../../include/glvu \
            -I../../include \
            -I../../include/dtgrid

LDFLAGS_COMMON = ${GLFLAGS} ${ZLIBFLAGS} ${PNGFLAGS}
