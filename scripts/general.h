//////
// Intro

// Some bit hacks and permutations

// (c) 2011..2014 by Jasper L. Neumann
// www.sirrida.de / programming.sirrida.de
// E-Mail: info@sirrida.de

// Granted to the public domain
// First version: 2013-02-15
// Last change: 2013-04-12

// Include file to supply general types and functions.

// At program start you should call
//   init_general
// for some general checks and for the creation of the mutex.
// Also, the random generator (random_int) is seeded.
//
// if (!init_general()) {
//   return 1;
//   }

#ifndef unit__general
#define unit__general

#include <stdint.h>

#if 0
  // Allow for Posix thread synchronization primitives
  #define USE_MUTEX
#endif

#define CONSOLE_PROGRAM
  // Has the program a console output facility?
  // Allow init_general() to dump error messages onto the console

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// Tell the compiler of an unused parameter
#ifdef __GNUC__
  #define UNUSED(x) x __attribute__((unused))
#else
  #define UNUSED(x) x
#endif

#ifndef mycall
  // Define your favourite function options here
  // #define mycall __attribute__((fastcall))
  // #define mycall __attribute__((regparm(3)))
  #define mycall
#endif

#if 0

  // Types via stdint.h
  #include <stdint.h>

  // Fixed size integer types
  typedef uint8_t t_8u;
  typedef uint16_t t_16u;
  typedef uint32_t t_32u;
  typedef uint64_t t_64u;
  // typedef __uint128_t t_128u;  // not a standard type

  // General integer type; >= 16 bit
  typedef int_fast16_t t_int;
  // General unsigned integer type; >= 16 bit
  typedef uint_fast16_t t_uint;
  // General integer type; >= 32 bit
  typedef int_fast32_t t_longint;

  typedef char t_char;

  // The boolean type should be predefined as bool (C++) or _Bool (C99, stdbool.h)
  // Since there is no usable standard way, we do it on our own. Sigh!
  typedef unsigned char t_bool;  // Poor replacement for boolean type
  #define false (0!=0)
  #define true (0==0)

#else

  // Types self-made
  // Fixed size integer types
  typedef signed char t_8s;
  typedef unsigned char t_8u;
  typedef signed short int t_16s;
  typedef unsigned short int t_16u;
  typedef signed int t_32s;
  typedef unsigned int t_32u;
  typedef signed long long int t_64s;
  typedef unsigned long long int t_64u;
  // typedef __int128_t t_128s;  // not a standard type
  // typedef __uint128_t t_128u;  // not a standard type

  // General integer type; >= 16 bit
  typedef int t_int;
  // General unsigned integer type; >= 16 bit
  typedef unsigned int t_uint;
  // General integer type; >= 32 bit
  typedef int t_longint;

  typedef char t_char;

  // The boolean type; might be predefined
  typedef unsigned char t_bool;  // Poor replacement for boolean type
  #ifndef false
    #define false ((t_bool)(0!=0))
  #endif
  #ifndef true
    #define true ((t_bool)(0==0))
  #endif

#endif

#ifdef USE_MUTEX
  #include <pthread.h>
  extern pthread_mutex_t global_mutex;
  // Access to a global general purpose mutex
  #define GLOBAL_MUTEX_WAIT   pthread_mutex_lock(&global_mutex);
  #define GLOBAL_MUTEX_SIGNAL pthread_mutex_unlock(&global_mutex);
  // Memory fence / barrier
  #define SYNC_SYNCHRONIZE       __sync_synchronize();
  #define SYNC_SYNCHRONIZE_LOAD  __sync_synchronize();
  #define SYNC_SYNCHRONIZE_STORE __sync_synchronize();
#else
  #define GLOBAL_MUTEX_WAIT
  #define GLOBAL_MUTEX_SIGNAL
  #define SYNC_SYNCHRONIZE
  #define SYNC_SYNCHRONIZE_LOAD
  #define SYNC_SYNCHRONIZE_STORE
#endif

mycall t_int random_int(t_int x);
mycall t_bool init_general();

#endif

// eof.
