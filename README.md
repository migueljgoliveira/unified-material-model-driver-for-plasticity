# UMMDp - Unified Material Model Driver for Plasticity

## :wrench: Pre-requisites

* Fortran compiler

## :rocket: Setup

### Source Files

Concatenate the UMMDp source files into one single file with the plug-in file first. Simply use the batch files (.sh/.bat) or run each command separately.

##### Unix/Linux

```
$ compile.sh
```

  or

```
$ cp source/plug_ummdp_abaqus.f source/tmp.f
$ cat source/ummdp*.f >> source/tmp.f
$ mv source/tmp.f compiled/ummdp.f
```

##### Windows

```
> compile.bat
```

  or

```
> copy "source\plug_ummdp_abaqus.f" "source\tmp.f"
> type "source\ummdp*.f" >> "source\tmp.f"
> move "source\tmp.f" "compiled\ummdp.f"
```

### Abaqus Input File

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

### Program Execution
  
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

## :computer: Setup

### Debug & Print

The first input parameter corresponds to the definition of debug and print mode, defined by the variable nvbs0. It is a mandatory parameter and the options are:

  * 0 - Error messages only
  * 1 - Summary of multistage return mapping
  * 2 - Detail of multistage return mapping and summary of Newton-Raphson
  * 3 - Detail of Newton-Raphson
  * 4 - Input/Output
  * 5 - All status for debug and print
  
### Elastic Properties

  * prela(1) - ID for elastic properties
  * prela(2~) - Data depends on ID

Only isotropic Hooke elastic properties can be defined. There are 2 ways to set them:

* Young's Modulus and Poisson's Ratio

  ```
  ID = 0
  prela(1) = 0
  prela(2) = 200.0E+3 (Young's modulus)
  prela(3) = 0.3 (Poisson's ratio)
  ```
  
* Bulk Modulus and Modulus of Rigidity

  ```
  ID = 1
  prela(1) = 0
  prela(2) = 166666.7 (Bulk modulus)
  prela(3) = 76923.08 (Modulus of rigidity)
  ```

### Yield Criterion

* von Mises (1913) :heavy_check_mark:
* Hill48 (1948) :heavy_check_mark:
* Yld2004-18p (2005) :heavy_check_mark: <code><sup>12</sup></code>
* CPB (2006) :heavy_check_mark:
* Karafillis-Boyce (1993) :grey_question:
* Hu (2005) :grey_question:
* Yoshida 6th Polynomial (2011) :grey_question:
* Gotoh Biquadratic (1978) :grey_question:
* Yld2000-2d (2003) :heavy_check_mark:
* Vegter (2006) :grey_question:
* BBC2005 (2005) :grey_question:
* Yld89 (1989) :grey_question:
* BBC2008 (2008) :grey_question:
* Hill 1990 (1990) :grey_question:
  
### Isotropic Hardening

  * Perfectly Plastic :heavy_check_mark:
  * Linear Hardening :heavy_check_mark:
  * Swift :heavy_check_mark:
  * Ludwick :heavy_check_mark:
  * Voce :heavy_check_mark:
  * Voce + Linear :heavy_check_mark:
  * Voce + Swift :heavy_check_mark:

### Kinematic Hardening

  * No Kinematic Hardening :heavy_check_mark:
  * Prager (1949) :heavy_check_mark:
  * Ziegler (1959) :heavy_check_mark:
  * Armstrong-Frederick (1966) :heavy_check_mark:
  * Chaboche (1979) :heavy_check_mark:
  * Chaboche (1979) - Ziegler Type :heavy_check_mark:
  * Yoshida-Uemori


## :books: References
<sup>1</sup> R. von Mises. 1913. Mechanik der festen Korper im plastisch deformablen Zustand. Gottin. Nachr. Math. Phys., 1: 582-592.

<sup>2</sup> R. Hill. 1948. A theory of the yielding and plastic flow of anisotropic metals. Proc. Roy. Soc. London, 193:281-297.

<sup>3</sup> F. Barlat, H. Aretz, J.W. Yoon, M.E. Karabin, J.C. Brem, R.E. Dick. 2005. Linear transformation-based anisotropic yield functions. International Journal of Plasticity 21:1009-1039.

<sup>4</sup> O. Cazacu, B. Plunkett, F. Barlat. 2006. Orthotropic yield criterion for hexagonal close packed metals. International Journal of Plasticity 22:1171-1194.

<sup>5</sup> A.P. Karafillis, M.C. Boyce. 1993. A general anisotropic yield criterion using bounds and a transformation weighting tensor. Journal of the Mechanics of Physics and Solids 41:1859-1886.

<sup>6</sup> W. Hu. 2005. An orthotropic yield criterion in a 3-D general stress state. International Journal of Plasticity 21:1771-1796.

<sup>7</sup> F. Yoshida, H. Hamasaki, T. Uemori. 2013. A user-friendly 3D yield function to describe anisotropy of steel sheets. International Journal of Plasticity 45:119-139.

<sup>8</sup> M. Gotoh. 1977. A theory of plastic anisotropy based on a yield function of fourth order (plane stress state) - I. International Journal of Mechanical Sciences 19-9:505-512.

<sup>9</sup> F. Barlat, J.C. Brem, J.W. Yoon, K. Chung, R.E. Dick, D.J. Lege, F. Pourboghrat, S.H. Choi, E. Chu. 2003. Plane stress yield function for aluminium alloy sheets-part 1: theory. International Journal of Plasticity 19:1297-1319.

<sup>10</sup> H. Vegter, A.H. van den Boogaard. 2006. A plane stress yield function for anisotropic sheet material by interpolation of biaxial stress states. International Journal of Plasticity 22:557-580.

<sup>11</sup> D. Banabic, D.S. Aretz, H. Comsa, L. Paraianu. 2005. An improved analytical description of orthotropy in metallic sheets. International Journal of Plasticity 21:493-512.

<sup>12</sup> F. Barlat, J. Lian. 1989. Plastic behavior and stretchability of sheet metals. Part I: a yield function for orthotropic sheets under plane stress conditions. International Journal of Plasticity. 5:51-66.

<sup>13</sup> D.S. Comsa, D. Banabic. 2008. Plane-stress yield criterion for highly-anisotropic sheet metals. Proceedings of NUMISHEET 2008.

<sup>14</sup> R. Hill. 1990. Constitutive modelling of orthotropic plasticity in sheet metals. Journal of the Mechanics and Physics of Solids.
