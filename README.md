This repository hosts the development code for adapting the generic FastDIRC code (see https://github.com/jmhardin/FastDIRC) for the GlueX DIRC detector.

## Compile FastDIRC
To compile FastDIRC, and do:
```
$ make
```
This will compile the driver specified by `DRIVER` variable in the `Makefile` and the binary executable `fastdirc_exe` will be put into the location specified by `OUT`.

### Running FastDIRC with GlueX DIRC data trees
There exists a plugin called `dirc_tree` which selects pion and kaon events from rho to pi+pi- and phi to K+K- reactions and create trees. Each entry of the tree is a `DrcEvent` which includes the track information (momentum, position etc.), some combo level information (invariant mass, missing mass squared, number of drift chamber hits etc.) and the PMT hits (specified by `DrcHit` class). Therefore, in order to read such trees, one needs to build the libraries for `DrcEvent` and `DrcHit` classes. The procedure to do that is:

- copy `DrcHit.*` and `DrcEvent.*` from halld_recon/src/plugins/Analysis/dirc_tree to FastDIRC_dev/src
- go to src
- do:
```
$ root -l
$ .L DrcHit.cc+
$ .L DrcEvent.cc+
```
This creates `*.d`, `*.pcm` and `*.so` objects. Then one can include the two `.so` files in in the `Makefile` of FastDIRC when compiling the driver. One also needs to include the header files in the driver. See `fastdirc_tree_minimal.cpp` driver for an example.
