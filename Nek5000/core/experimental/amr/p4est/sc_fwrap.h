/** @file sc_fwrap.h
 * @brief Fortran77 interface of sc library
 *
 * @author Adam Peplinski
 * @date Feb 20, 2016
 *
 * @ingroup nekamr
 */
#ifndef NEKP4EST_SC_FWRAP_H_
#define NEKP4EST_SC_FWRAP_H_

#if !defined(NEKAMR_H_)
#error "sc_fwrap.h requires nekamr.h"
#endif

#if !defined(SC_H)
#error "sc_fwrap.h requires sc.h"
#endif

#if !defined(NAME_H)
#error "sc_fwrap.h requires name.h"
#endif

/* FORTRAN interface
* #define fsc_ FORTRAN_NAME(fsc_,FSC_)
*/
#define fsc_check     FORTRAN_NAME(fsc_check,FSC_CHECK)
#define fsc_init      FORTRAN_NAME(fsc_init,FSC_INIT)
#define fsc_finalize  FORTRAN_NAME(fsc_finalize,FSC_FINALIZE)
#define fsc_pkg_reg    FORTRAN_NAME(fsc_pkg_reg,FSC_PKG_REG)
#define fsc_pkg_is_reg FORTRAN_NAME(fsc_pkg_is_reg,FSC_PKG_IS_REG)
#define fsc_pkg_unreg  FORTRAN_NAME(fsc_pkg_unreg,FSC_PKG_UNREG)
#define fsc_pkg_print  FORTRAN_NAME(fsc_pkg_print,FSC_PKG_PRINT)
#define fsc_set_log FORTRAN_NAME(fsc_set_log,FSC_SET_LOG)
#define fsc_log     FORTRAN_NAME(fsc_log,FSC_LOG)
#define fsc_abort FORTRAN_NAME(fsc_abort,FSC_ABORT)
#define fsc_check_abort FORTRAN_NAME(fsc_check_abort,FSC_CHECK_ABORT)

/**
 *
 * @param mpiret
 */
void fsc_check(int * mpiret);

/**
 *
 * @param fmpicomm
 * @param catch_signals
 * @param print_backtrace
 * @param log_threshold
 */
void fsc_init(MPI_Fint * fmpicomm,int * catch_signals,
		int * print_backtrace, int * log_threshold);

/**
 *
 */
void fsc_finalize();

/**
 *
 * @param package_id
 * @param log_threshold
 * @param name
 */
void fsc_pkg_reg(int * package_id, int * log_threshold, char * name, int len_n);

/**
 *
 * @param package_id
 * @param package_reg
 */
void fsc_pkg_is_reg(int * package_id, int * package_reg);

/**
 *
 * @param package_id
 */
void fsc_pkg_unreg(int * package_id);

/**
 *
 * @param log_priority
 */
void fsc_pkg_print(int * log_priority);

/**
 *
 * @param package
 * @param category
 * @param priority
 * @param msg
 * @param len_msg
 */
void fsc_log (int * package, int * category, int * priority,
	      char *msg, int len_msg);

/** Controls the default SC log behavior.
 *
 * @param [in] log_threshold   Set default SC log threshold
 */
void fsc_set_log(int * log_threshold);

/** The main sc log function
 *
 * @param msg
 * @param len_msg     Fortran77-C required; character length
 */
void fsc_abort (char *msg, int len_msg);

/**
 *
 * @param q
 * @param msg
 * @param len_msg
 */
void fsc_check_abort (int * q, char *msg, int len_msg);

#endif /* NEKP4EST_SC_FWRAP_H_ */
