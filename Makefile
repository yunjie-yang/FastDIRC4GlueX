CFLAGS_BASE = -O3 -mtune=native -Wno-comment -L ./lib/ `root-config --cflags` `root-config --glibs` -lMinuit
CFLAGS_OPT = -g
CFLAGS_BASE += $(CFLAGS_OPT)
INCLUDE = -I./src/

CFLAGS = $(CFLAGS_BASE)

LIBLOC = ./lib

OUT = fastdirc_exe

#DRIVER = fastdirc_minimal
DRIVER = fastdirc_tree_minimal
#DRIVER = fastdirc_tree_full_KinBins
#DRIVER = fastdirc_tree_observables

OBJFILES = dirc_base_sim.o
OBJFILES += dirc_threesegbox_sim.o
OBJFILES += dirc_rect_digitizer.o
OBJFILES += dirc_spread_gaussian.o
OBJFILES += GlueXUserOptions.o

OBJFILES += dirc_spread_gaussian_new.o
OBJFILES += dirc_spread_gaussian_exp.o

DRCLIB = ./src/DrcHit_cc.so
DRCLIB += ./src/DrcEvent_cc.so

OBJLOC = $(patsubst %,$(LIBLOC)/%,$(OBJFILES))
LIBFILES = $(LIBLOC)
vpath %.o ./lib/
vpath %.cpp ./src/

%.o : %.cpp
	g++ -Wall $(CFLAGS) $(INCLUDE) -g -o $@ -c $<
	mv $@ $(LIBLOC)

.PHONY : all
all: $(DRIVER).cpp $(OBJFILES)
	g++ -Wall $(DRIVER).cpp $(OBJLOC) $(DRCLIB) $(CFLAGS) $(INCLUDE) -o $(OUT)

.PHONY : libs
libs: $(OBJFILES)
	echo libraries built

.PHONY : clean
clean:
	rm -f lib/*.o
	rm -f $(OUT)

.PHONY : cleanall
cleanall:
	rm -f lib/*.o
	rm -f *.gcda
	rm -f $(OUT)
