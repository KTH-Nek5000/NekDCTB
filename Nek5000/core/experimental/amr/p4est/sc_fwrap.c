
/*
 * sc_fwrap.c
 * Fortran interface for sc library.
 *
 *  Created on: Feb 20, 2016
 *      Author: Adam Peplinski
 */
#ifdef MPI
#define SC_ENABLE_MPI
#endif

#ifdef MPIIO
#define SC_ENABLE_MPIIO
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "../../../experimental/amr/nekamr.h"

#include <sc.h>

#include "name.h"

#include "../../../experimental/amr/p4est/sc_fwrap.h"

/* initialize/finalize */
void fsc_check(int * mpiret)
{
	SC_CHECK_MPI (*mpiret);
}

void fsc_init(MPI_Fint * fmpicomm,
	      int * catch_signals, int * print_backtrace,
	      int * log_threshold)
{
  MPI_Comm mpicomm;
  mpicomm = MPI_Comm_f2c(*fmpicomm);
  sc_init (mpicomm, *catch_signals, *print_backtrace,
  	   NULL, *log_threshold);
}

void fsc_finalize()
{
  sc_finalize ();
}

/* package registration for logging */
void fsc_pkg_reg(int * package_id, int * log_threshold, char * name, int len_n)
{
  char *full=" ";
  *package_id =sc_package_register (NULL, *log_threshold,
		  name, full);
}

void fsc_pkg_is_reg(int * package_id, int * package_reg)
{
  *package_reg = sc_package_is_registered (*package_id);
}

void fsc_pkg_unreg(int * package_id)
{
  sc_package_unregister (*package_id);
}

void fsc_pkg_print(int * log_priority)
{
  sc_package_print_summary (*log_priority);
}

/* logging */
void fsc_log (int * package, int * category, int * priority,
	      char *msg, int len_msg)
{
  char str[len_msg];
  strcpy(str, msg);
  strcat(str, "\n");
  sc_log (__FILE__, __LINE__, *package, *category, *priority, str);
}

void fsc_set_log(int * log_threshold)
{
  sc_set_log_defaults (stdout, NULL, *log_threshold);
}

/* aborting simulation */
void fsc_abort (char *msg, int len_msg)
{
  char str[len_msg];
  strcpy(str, msg);
  strcat(str, "\n");
  SC_ABORT (str);
}

void fsc_check_abort (int * q, char *msg, int len_msg)
{
  char str[len_msg];
  strcpy(str, msg);
  strcat(str, "\n");
  SC_CHECK_ABORT ( *q, str);
}

/* simulation parameters */
