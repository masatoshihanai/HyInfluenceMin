/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef HYMINSOLVER_WTIME_HPP
#define HYMINSOLVER_WTIME_HPP

#ifdef HAS_OMP
#include <omp.h>
double getTime() {return omp_get_wtime();};
#else
#include <time.h>
inline double getTime() {return (double) clock() / CLOCKS_PER_SEC;}
#endif

#endif //HYMINSOLVER_WTIME_HPP
