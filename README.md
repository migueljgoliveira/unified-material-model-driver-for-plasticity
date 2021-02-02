# UMMDp - Unified Material Model Driver for Plasticity

## Pre-requisites

* Fortran compiler

## Usage

### Preparation of Program Source Files

Concatenate the UMMDp source files into one single file with the plug-in file first. Simply use the batch files (.sh/.bat) or run each command separately.

#### Unix/Linux

````
$ compile.sh
````

  or

````
$ cp source/plug_ummdp_abaqus.f source/tmp.f
$ cat source/ummdp*.f >> source/tmp.f
$ mv source/tmp.f compiled/UMMDp.f
````

#### Windows

````
$ compile.bat
````

  or

````
> copy "source\plug_ummdp_abaqus.f" "source\tmp.f"
> type "source\ummdp*.f" >> "source\tmp.f"
> move "source\tmp.f" "compiled\UMMDp.f"
````
