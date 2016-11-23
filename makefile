SHELL := /bin/bash -e

BINS = rasterizeSolid rubberBall

SRCS = deformCD dtgrid geometry integrator material util

all:
	mkdir -p bin
	-for d in $(BINS); do (echo -e "\n==== Entering $$d ====\n";cd ./projects/$$d; make -j12;cd ../..); done

clean:
	-for d in $(BINS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./projects/$$d; make clean;cd ../..); done 
  -for d in $(SRCS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./src/$$d; rm -f *.o; cd ../..); done

objclean:
	-for d in $(SRCS); do (echo -e "\n==== Cleaning $$d ====\n";cd ./src/$$d; rm -f *.o; cd ../..); done

ctags:
	ctags -R *
