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

HyMinSoler supports the check-in format provided by [SNAP](https://snap.stanford.edu/data/index.html).

```bash
# userID location
0    1
1    0
1    2
1    3
2    1
3    2
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

