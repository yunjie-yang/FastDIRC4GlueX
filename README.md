This repository hosts the development code for adapting the generic FastDIRC code (see https://github.com/jmhardin/FastDIRC) for the GlueX DIRC detector.

# Intro and structure
* `src` directory: it contains the source code performing the core functions of FastDIRC simulation and reconstruction; the compiled library files will be put into the `lib` directory
* `geometry_files` directory: it contains two csv files specifying the GlueX DIRC detector geometry in the internal FastDIRC parameterization. `FastDIRC_HDDS_Nominal.csv` is Yunjie Yang's attempt to translate the nominal HDDS geometry into the FastDIRC parameterization, and can be used as the default.
* `macros` directory: it contains a few macros for post-processing the FastDIRC outputs: creating DLL plots, draw ROC curves, evaluate separation power etc.
* `*.cpp`: these are the so-called _drivers_, which performs analysis of the DIRC trees by calling the core functions of FastDIRC; in principle, one driver file could perform multiple analysis tasks, but it is probably easier to have (and make the code more readable)
* `*.in`: these are configuration files to control the behavior of the driver, such as controlling the number of tracks to process, specifying the offsets and rotations etc.
* `Makefile`: this is the Makefile, which also specifies which driver to compile via the `DRIVER` variable (not including the .cpp part), and where to put the driver output file via the `OUT` variable.


# Compile FastDIRC to process GlueX DIRC data trees
There exists a plugin called `dirc_tree` which selects pion and kaon events from rho to pi+pi- and phi to K+K- reactions and create trees. Each entry of the tree is a `DrcEvent` which includes the track information (momentum, position etc.), some combo level information (invariant mass, missing mass squared, number of drift chamber hits etc.) and the PMT hits (specified by `DrcHit` class). Therefore, in order to read such trees, one needs to build the libraries for `DrcEvent` and `DrcHit` classes, and link them during the compilation of FastDIRC.

1. Create the libraries for `DrcEvent` and `DrcHit` classes. First, make sure the `DrcHit.*` and `DrcEvent.*` in the `src` directory are compatible with the trees that you will be processing. Next, go to `src`, and do:
```
$ root -l
$ .L DrcHit.cc+
$ .L DrcEvent.cc+
```
This creates `*.d`, `*.pcm` and `*.so` objects.

2. Compile a driver. Specify the desired `DRIVER` (_e.g._, fastdirc_tree_minimal.cpp) and `OUT`, then simply do:
```
$ make (-j4)
```
This creates an executable specified by `OUT`, _e.g._, `fastdirc_exe`.

# Run FastDIRC and analyze output
1. Run the executable with a configuration file:
```
$ ./fastdirc_exe config_tree_minimal.in
```
This creates a root file (_e.g._,output_hists_example.root) with a few histograms in it.

2. Analyze the output:
```
$ cd macros
$ python get_res_and_plot.py
```
This creates a csv file in the top directory and some plots in the output_plots directory. The `get_res_and_plot.py` has some configurable labels and runs two root macros: `graphicHistos.C` which computes the equivalent per-track resolution and `draw_DLLs_arg.C` which draws the DLL plots.

# Core FastDIRC functions
