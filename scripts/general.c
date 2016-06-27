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

#include "general.h"

#ifdef USE_MUTEX
  pthread_mutex_t global_mutex;
#endif

mycall t_int random_int(t_int x) {
// 0..x-1, x quite small, typically <=32767
// Replace with your own generator if you want

  return rand() % x;
  }

mycall t_bool init_general() {
// General initializations and checks.
// Results false if some error is detected, also an error message is printed.

  t_bool res;

  res = true;
  if (sizeof(t_8s)*8 != 8) {
#ifdef CONSOLE_PROGRAM
    printf("t_8s defective\n");
#endif
    res = false;
    }
  if (sizeof(t_8u)*8 != 8) {
#ifdef CONSOLE_PROGRAM
    printf("t_8u defective\n");
#endif
    res = false;
    }
  if (sizeof(t_16s)*8 != 16) {
#ifdef CONSOLE_PROGRAM
    printf("t_16s defective\n");
#endif
    res = false;
    }
  if (sizeof(t_16u)*8 != 16) {
#ifdef CONSOLE_PROGRAM
    printf("t_16u defective\n");
#endif
    res = false;
    }
  if (sizeof(t_32s)*8 != 32) {
#ifdef CONSOLE_PROGRAM
    printf("t_32s defective\n");
#endif
    res = false;
    }
  if (sizeof(t_32u)*8 != 32) {
#ifdef CONSOLE_PROGRAM
    printf("t_32u defective\n");
#endif
    res = false;
    }
  if (sizeof(t_64s)*8 != 64) {
#ifdef CONSOLE_PROGRAM
    printf("t_64s defective\n");
#endif
    res = false;
    }
  if (sizeof(t_64u)*8 != 64) {
#ifdef CONSOLE_PROGRAM
    printf("t_64u defective\n");
#endif
    res = false;
    }
  if (sizeof(t_int)*8 < 16) {
#ifdef CONSOLE_PROGRAM
    printf("t_int defective\n");
#endif
    res = false;
    }
  if (sizeof(t_uint)*8 < 16) {
#ifdef CONSOLE_PROGRAM
    printf("t_uint defective\n");
#endif
    res = false;
    }
  if (sizeof(t_longint)*8 < 32) {
#ifdef CONSOLE_PROGRAM
    printf("t_longint defective\n");
#endif
    res = false;
    }

  srand((unsigned int)(time(NULL)));

#ifdef USE_MUTEX
  pthread_mutex_init(&global_mutex, 0);
#endif

  return res;
  }

// eof.
