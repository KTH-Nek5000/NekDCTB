/** @file nekamr.h
 *  @brief Definitions for C-Fortran77 interface in combination with p4est
 *
 * @author Adam Peplinski
 * @date Feb 20, 2016
 *
 * @ingroup nekamr
 */
/** @defgroup nekamr Non-conforming Nek5000 with p4est
 *
 * AMR version of nek5000
 */

#ifndef NEKAMR_H_
#define NEKAMR_H_

/* check preprocessing flags */
#if !defined(N_DIM)
#error "Undefined N_DIM; check makenek.inc"
#endif

#if N_DIM < 2 || N_DIM > 3
#error "Wrong N_DIM value; check makenek.inc"
#endif

/* Numbers designating the type and action for refinement  */
#define AMR_RM_NONE       0     /* No refinement action */
#define AMR_RM_H_REF      1     /* Refine by splitting element */
#define AMR_RM_H_CRS     -1     /* Coarsen by merging element */
#define AMR_RM_P_REF      2     /* Refine by rising polynomial order in element */
#define AMR_RM_P_CRS     -2     /* Coarsen by lowering polynomial order in element */

/* For mesh partitioning */
/* #define AMR_PRT_GRPH */
#undef AMR_PRT_GRPH

/* For ParMETIS */
/* do we use physical coordinates of nodes to speed up partitioning */
/* #define AMR_PRT_CRD */
#undef AMR_PRT_CRD

/* do we use weights for partitioning */
/* #define AMR_PRT_WGT */
#undef AMR_PRT_WGT

#endif /* NEKAMR_H_ */
