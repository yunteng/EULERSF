include ../common.osx.inc

EXECUTABLE = ../../bin/rasterizeSolid
regression: $(EXECUTABLE)
ALL: $(EXECUTABLE)

SOURCEC    = main.cpp \
             $(UTILSRC) \
             $(GEOMETRY)

OBJSC = $(SOURCEC:.cpp=.o)

.cpp.o:
	$(CC) $(CPPFLAGS) $< -o $@
	
regression: ALL
ALL: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJSC)
	${CC} $(OBJSC) ${LDFLAGS_COMMON} -o ${EXECUTABLE}

clean:
	rm -rf *.o ${EXECUTABLE}

cleanall:
	rm -rf $(OBJSC) ${EXECUTABLE}