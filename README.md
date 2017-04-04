# MultiCumulants
General purpose code for computing cumulants

#### Checking out this code for development
```bash
git clone git@github.com:jdbrice/MultiCumulants.git
cd MultiCumulants
git submodule init
git submodule update
```

#### Development application
```bash
make dev
./bin/dev.app
```

log file is bin/debug.log - will show all logging output during development running

#### Toy MonteCarlo application
```bash
make bin/toymc.app
./bin/toymc.app
```
to create the ROOT shared library for interactive use:
```bash
make rootlib
```

then the ROOT shared library is in lib/ToyMC.so, you can load it into ROOT with:
```c++
gSystem->Load( "lib/ToyMC.so");
```




#### Contribution Guidelines

