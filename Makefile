CFLAGS_BASE = -O3 -mtune=native -Wno-comment -rpath $(ROOTSYS)/lib -L ./lib/ `root-config --cflags` `root-config --glibs` -lMinuit
CFLAGS_OPT = -g
CFLAGS_BASE += $(CFLAGS_OPT)
INCLUDE = -I./src/

CFLAGS = $(CFLAGS_BASE)

LIBLOC = ./lib

OUT = fastdirc_exe

#DRIVER = fastdirc_tree

#DRIVER = fastdirc_tree_tune_geometry
#DRIVER = fastdirc_tree_tune_supports

#DRIVER = fastdirc_tree_track_selection
#DRIVER = fastdirc_tree_KinBins

#DRIVER = fastdirc_minimal
#DRIVER = fastdirc_tree_minimal

#DRIVER = fastdirc_tree_objective_functions_study

DRIVER = fastdirc_tree_full_KinBins
#DRIVER = fastdirc_tree_full_simple

#DRIVER = fastdirc_tree_observables

#DRIVER = fastdirc_tree_output_tree

OBJFILES = dirc_base_sim.o
OBJFILES += dirc_threesegbox_sim.o
OBJFILES += dirc_rect_digitizer.o
OBJFILES += dirc_spread_gaussian.o
OBJFILES += GlueXUserOptions.o

OBJFILES += dirc_spread_gaussian_new.o
OBJFILES += dirc_spread_gaussian_exp.o

#DRCLIB = /w/halld-scifs17exp/home/yunjiey/build_master/halld_recon/Linux_CentOS7-x86_64-gcc4.8.5/plugins/pid_dirc.so
#DRCLIB = /home/yunjiey/Documents/gluex_install/gluex_top/halld_recon/Linux_Ubuntu16.04-x86_64-gcc5.4.0/plugins/pid_dirc.so
DRCLIB = ./src/DrcHit_cc.so
DRCLIB += ./src/DrcEvent_cc.so

#OBJFILES += dirc_optical_sim.o
#OBJFILES += dirc_babar_sim.o
#OBJFILES += dirc_lut_enum.o
#OBJFILES += dirc_lut.o
#OBJFILES += dirc_gluex_lut_enum.o
#OBJFILES += dirc_babar_digitizer.o
#OBJFILES += dirc_probability_spread.o
#OBJFILES += dirc_spread_relative.o
#OBJFILES += dirc_spread_radius.o
#OBJFILES += dirc_spread_linear_soft.o
#OBJFILES += dirc_probability_separation.o
#OBJFILES += dirc_progressive_separation.o

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
	
