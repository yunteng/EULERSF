include ../common.osx.openmp.inc

EXECUTABLE = ../../bin/rubberBall
regression: $(EXECUTABLE)
ALL: $(EXECUTABLE)

SOURCEC    = main.cpp \
             $(UTILSRC) \
             $(GEOMETRY) \
             $(INTEGRATOR) \
             $(MATERIAL)

OBJSC = $(SOURCEC:.cpp=.o)

.cpp.o:
	$(CC) $(CPPFLAGS) ${IPOPTINCL} $< -o $@
	
regression: ALL
ALL: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJSC)
	${CC} $(OBJSC) ${LDFLAGS_COMMON} ${IPOPTLIBS} -o ${EXECUTABLE}

clean:
	rm -rf *.o ${EXECUTABLE}

cleanall:
	rm -rf $(OBJSC) ${EXECUTABLE}