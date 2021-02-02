# UMMDp - Unified Material Model Driver for Plasticity

## :wrench: Pre-requisites

* Fortran compiler

## :rocket: Usage

<details><summary><b>Preparation of program source files</b></summary>

Concatenate the UMMDp source files into one single file with the plug-in file first. Simply use the batch files (.sh/.bat) or run each command separately.

##### Unix/Linux

```sh
$ compile.sh
```

  or

```sh
$ cp source/plug_ummdp_abaqus.f source/tmp.f
$ cat source/ummdp*.f >> source/tmp.f
$ mv source/tmp.f compiled/ummdp.f
```

##### Windows

```cmd
> compile.bat
```

  or

```cmd
> copy "source\plug_ummdp_abaqus.f" "source\tmp.f"
> type "source\ummdp*.f" >> "source\tmp.f"
> move "source\tmp.f" "compiled\ummdp.f"
```
</details>

<details><summary><b>Preparation of the input file</b></summary>

This section describes the keywords in Abaqus input data file for use in the UMMDp.

1. Definition of the principal axis for the material anisotropy (for more information, please refer to Abaqus's manual)
    ````
    *ORIENTATION, NAME=ORI-1
    1., 0., 0., 0., 1., 0.
    3, 0.
    ````

2. Definition of the material model (more details are provided later)
    ````
    *MATERIAL, NAME=UMMDp
    *USER MATERIAL, CONSTANTS=27
    0, 0, 1000.0, 0.3, 2, -0.069, 0.936, 0.079,
    1.003, 0.524, 1.363, 0.954, 1.023, 1.069, 0.981, 0.476,
    0.575, 0.866, 1.145, -0.079, 1.404, 1.051, 1.147, 8.0,
    0, 1.0, 0
    ````

3. Define the number of internal state variables (SDV)

    Set the number of state variables equal to 1+NTENS, where NTENS is the number of
    components of the tensor variables. NTENS=3 for plane stress or a shell element,
    and NTENS=6 for a solid element. The 1st state variable is reserved for the equivalent plastic
    strain, and NTENS is reserved for the plastic strain components. The following ex-
    ample corresponds to a solid element without kinematic hardening:
    ````
    *DEPVAR
    7,
    ````
    In the case of kinematic hardening, the number of internal state variables corresponds
    to the equivalent plastic strain, plastic strain components and components of each
    partial back-stress tensor.

4. Define the user output variables (UVARM)

    UMMDp can output three user output variables:

    - UVARM(1): current equivalent stress (the value calculated by substituting the stress com-
    ponents for the yield function)

    - UVARM(2): current yield stress (the value calculated by substituting the equivalent plastic
    strain for the function of the isotropic hardening curve)

    - UVARM(3:8): current components of the total back-stress tensor
    ````
    *USER OUTPUT VARIABLES
    8,
    ````

5. Define output variables for post processing

    This keyword controls the output variables (e.g. equivalent plastic strain and equiv-
    alent stress) for post processing.
    ````
    *OUTPUT, FIELD
    *ELEMENT OUTPUT
    SDV, UVARM
    ````

</details>

<details><summary><b>Execution of the program</b></summary>
  
To execute the program there are two options: 1. link the user subroutine in source code
or 2. link the user subroutine previously compiled:

1. To execute the program with the user subroutine in source code, execute the command:
    ```
    $> abaqus job=jobname user=ummdp.f
    ```

2. To execute the program with the user subroutine previously compiled, execute the commands:
    ````
    > abaqus job=jobname user=ummdp.obj
    ````
    ````
    $ abaqus job=jobname user=ummdp.o
    ````

To compile the file ummdp.obj/o use
    ````
    $> abaqus make library=ummdp.f
    ````

</details>

## :computer: Setup

<details><summary><b>Debug and print</b></summary>

The first input parameter corresponds to the definition of debug and print mode, defined by the variable nvbs0. It is a mandatory parameter and the options are:
  * 0 - Error messages only
  * 1 - Summary of multistage return mapping
  * 2 - Detail of multistage return mapping and summary of Newton-Raphson
  * 3 - Detail of Newton-Raphson
  * 4 - Input/Output
  * 5 - All status for debug and print
  
</details>

<details><summary><b>Elastic properties</b></summary>

  * prela(1) - ID for elastic properties
  * prela(2~) - Data depends on ID

Only isotropic Hooke elastic properties can be defined. There are 2 ways to set them:

* Young's Modulus and Poisson's Ratio
  * ID = 0
  * prela(1) = 0
  * prela(2) = 200.0E+3 (Young's modulus)
  * prela(3) = 0.3 (Poisson's ratio)

* Bulk Modulus and Modulus of Rigidity
  * ID = 1
  * prela(1) = 0
  * prela(2) = 166666.7 (Bulk modulus)
  * prela(3) = 76923.08 (Modulus of rigidity)

</details>

<details><summary><b>Yield criterion</b></summary>

</details>

<details><summary><b>Isotropic hardening</b></summary>

</details>

<details><summary><b>Kinematic hardening</b></summary>

</details>
