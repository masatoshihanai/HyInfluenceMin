# _**Hypergraph-based Influence Minimization**_ : Influence Minimization Solver for Hypergraph 

## Reference

- M. Hanai, et.al. _"Influence Minimization of Virus Infection via Social Gathering Activities: Towards Effective Spread Control of COVID-19"_ (2020)

## 1, Quick Start

```bash
$ git clone git@github.com:masatoshihanai/HyInfluenceMin.git
$ mkdir build; cd build
$ cmake .. ; make
```

## Outline of the Page

1. [__Graph data format__](#data)
2. [__How to compile__](#compile)
3. [__How to run minimization solver__](#run)
5. [__Related tools__](#related)
6. [__Acknowledgement / Contact__](#ack)

---

### 1. Graph data format {#data}

HyMinSolver supports the check-in format provided in [SNAP](https://snap.stanford.edu/data/loc-Gowalla.html).

```bash
#[user] [check-in time]         [latitude]      [longitude]     [location id]
196514  2010-07-24T13:45:06Z    53.3648119      -2.2723465833   145064
196514  2010-07-24T13:44:58Z    53.360511233    -2.276369017    1275991
196514  2010-07-24T13:44:46Z    53.3653895945   -2.2754087046   376497
196514  2010-07-24T13:44:38Z    53.3663709833   -2.2700764333   98503
196514  2010-07-24T13:44:26Z    53.3674087524   -2.2783813477   1043431
196514  2010-07-24T13:44:08Z    53.3675663377   -2.278631763    881734
196514  2010-07-24T13:43:18Z    53.3679640626   -2.2792943689   207763
196514  2010-07-24T13:41:10Z    53.364905       -2.270824       1042822	
```

### 2. How to compile {#compile}

First, prepare the standard development tools for C/C++

__Requirements__: `C/C++ compiler`  `Make` `CMake` `Git`

- __C++ compiler__: [GCC](https://gcc.gnu.org/install/), [IntelC++](https://software.intel.com/en-us/c-compilers), [Clang/LLVM](https://clang.llvm.org/index.html), or etc.
- [Make](https://www.gnu.org/software/make/), [CMake](https://cmake.org/), [Git](https://git-scm.com/)

Next, get the code from github

```bash
### get the code from git repository
$ git clone git@github.com:masatoshihanai/HyInfluenceMin.git
```

Configure and compile

```bash
### make new build directory
$ cd HyInfluenceMin
$ mkdir build
$ cd build
```

```bash
### configure and compile
$ cmake ..
$ make
```

We have successed to compile and run the program in these environments:

##### Tested Environments

- Compiler: `GCC 7.4`
- OS:  `Ubuntu 18.04`

### 3. How to run HyMinSolver {#run}

```bash
$ ./HyMin <checkin-data> <# of restriction activities>
```

Run with `-h` for the detail information.
```bash
$ ./HyMin -h
```
```bash
 HyMin - Hypergraph-based Influence Minimization Soler
   Usage: $ ./HyMin [-h] [-v] <checkin-file> <# of restriction activities>      
      -h: Help
      -v: Verbose. Output details
       
    Example to use:
      $ ./HyMin -v ../data/Brightkite.checkin
```

### TODOs
- 