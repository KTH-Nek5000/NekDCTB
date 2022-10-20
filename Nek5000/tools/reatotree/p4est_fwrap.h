/** @file p4est_fwrap.h
 *  @brief Fortran77 interface of p4est library
 *
 * @author Adam Peplinski
 * @date Feb 20, 2016
 *
 * @ingroup nekp4est
 */
#ifndef NEKP4EST_P4EST_FWRAP_H_
#define NEKP4EST_P4EST_FWRAP_H_

#if !defined(NEKP4EST_H_)
#error "p4est_fwrap.h requires nekp4est.h"
#endif

#if !defined(P4EST_H) && !defined(P8EST_H)
#error "p4est_fwrap.h requires p4est.h or p8est.h"
#endif

#if !defined(NAME_H)
#error "p4est_fwrap.h requires name.h"
#endif

/* FORTRAN interface
 * #define fp4est_ FORTRAN_NAME(fp4est_,FP4EST_)
 */
/* wrappers */
/* initialize */
#define fp4est_init          FORTRAN_NAME(fp4est_init,FP4EST_INIT)
/* Connectivity */
#define fp4est_cnn_new       FORTRAN_NAME(fp4est_cnn_new,FP4EST_CNN_NEW)
#define fp4est_cnn_del       FORTRAN_NAME(fp4est_cnn_del,FP4EST_CNN_DEL)
#define fp4est_cnn_attr      FORTRAN_NAME(fp4est_cnn_attr,FP4EST_CNN_ATTR)
#define fp4est_cnn_valid     FORTRAN_NAME(fp4est_cnn_valid,FP4EST_CNN_VALID)
#define fp4est_cnn_complete  FORTRAN_NAME(fp4est_cnn_complete,FP4EST_CNN_COMPLETE)
#define fp4est_cnn_save      FORTRAN_NAME(fp4est_cnn_save,FP4EST_CNN_SAVE)
#define fp4est_cnn_load      FORTRAN_NAME(fp4est_cnn_load,FP4EST_CNN_LOAD)
/* tree_ management */
#define fp4est_tree_new      FORTRAN_NAME(fp4est_tree_new,FP4EST_TREE_NEW)
#define fp4est_tree_del      FORTRAN_NAME(fp4est_tree_del,FP4EST_TREE_DEL)
#define fp4est_tree_valid    FORTRAN_NAME(fp4est_tree_valid,FP4EST_TREE_VALID)
#define fp4est_tree_save     FORTRAN_NAME(fp4est_tree_save,FP4EST_TREE_SAVE)
#define fp4est_tree_load     FORTRAN_NAME(fp4est_tree_load,FP4EST_TREE_LOAD)
/* nekp4est internal load balance */
#define fp4est_part          FORTRAN_NAME(fp4est_part,FP4EST_PART)
/* I/O (VTK) */
#define fp4est_vtk_write     FORTRAN_NAME(fp4est_vtk_write,FP4EST_VTK_WRITE)
#define fp4est_vtk_iscalar   FORTRAN_NAME(fp4est_vtk_iscalar,FP4EST_VTK_ISCALAR)
/* testing */
#define fp4est_test_print    FORTRAN_NAME(fp4est_test_print,FP4EST_TEST_PRINT)
/* Initialize element data */
#define nekp4est_init_msh_dat FORTRAN_NAME(nekp4est_init_msh_dat,NEKP4EST_INIT_MSH_DAT)

/** Data type for user variables; required by p4est */
typedef struct user_data_s {
  int imsh; /**< velocity (0) and temperature (1) mesh indicator */
  int igrp; /**< element group */
  int crv[P4EST_FACES]; /**< curvature data; concerns external faces; requred by face projection */
  int bc[P4EST_FACES]; /**< boundary condition data; concerns external faces */
  int ref_mark; /**< integer to store refinement mark; definition in nekp4est.h */
  // to keep track of changes of nek5000 global element numbering and element distribution
  // element
  int el_gln; /**< element global numbering; nek5000 side */
  int el_ln; /**< element local numbering; nek5000 side */
  int el_nid; /**< mpi rank owning element; nek5000 side */
  // parent
  int parent_gln; /**< parent global numbering; nek5000 side */
  int parent_ln; /**< parent local numbering; nek5000 side */
  int parent_nid; /**< mpi rank owning parent; nek5000 side */
  // children
  int children_gln[P4EST_CHILDREN]; /**< children global numbering; nek5000 side */
  int children_ln[P4EST_CHILDREN]; /**< children local numbering; nek5000 side */
  int children_nid[P4EST_CHILDREN]; /**< mpi rank owning children; nek5000 side */
} user_data_t;

/* Wrappers */
/** Initialize p4est package setting log verbosity
 *
 * @param [in] log_threshold  Log level
 */
void fp4est_init(int * log_threshold)
;

/** Initialize elements connectivity
 *
 * @param num_vertices
 * @param num_trees
 * @param num_edges
 * @param num_corners
 * @param vertices
 * @param tree_to_vertex
 * @param tree_to_tree
 * @param tree_to_face
 * @param tree_to_edge
 * @param ett_offset
 * @param edge_to_tree
 * @param edge_to_edge
 * @param tree_to_corner
 * @param ctt_offset
 * @param corner_to_tree
 * @param corner_to_corner
 */
#ifdef P4_TO_P8
void fp4est_cnn_new(int * num_vertices, int * num_trees, int * num_edges,
		int * num_corners, double *vertices, int * tree_to_vertex,
		int * tree_to_tree, int8_t * tree_to_face, int * tree_to_edge,
		int * ett_offset, int * edge_to_tree, int8_t * edge_to_edge,
		int * tree_to_corner, int * ctt_offset, int * corner_to_tree,
		int8_t * corner_to_corner)
;
#else
void fp4est_cnn_new(int * num_vertices, int * num_trees,
		int * num_corners,
		double *vertices,
		int * tree_to_vertex, int * tree_to_tree,
		int8_t * tree_to_face,
		int * tree_to_corner, int * ctt_offset,
		int * corner_to_tree, int8_t * corner_to_corner)
;
#endif

/** Destroy mesh connectivity */
void fp4est_cnn_del()
;

/** Allocate or free the attribute fields in a connectivity
 *
 * @param enable_tree_attr
 */
void fp4est_cnn_attr(int * enable_tree_attr)
;

/** Check connectivity consistency
 *
 * @param [out] is_valid   non zero for correct connectivity
 */
void fp4est_cnn_valid(int * is_valid)
;

/** Internally connect a connectivity based on tree_to_vertex information.
 * Only internal connectivity; nor mesh periodicity
 */
void fp4est_cnn_complete();

/** Save connectivity in a file
 *
 * @param filename
 * @param len_f
 */
void fp4est_cnn_save(char *filename, int len_f)
;

/** Load a connectivity structure from a file
 *
 * @param filename
 * @param len_f
 */
void fp4est_cnn_load(char * filename, int len_f)
;

/** Generate forest
 *
 * @param fmpicomm
 * @param min_level
 */
void fp4est_tree_new(MPI_Fint * fmpicomm, int * min_level)
;

/** Destroy tree */
void fp4est_tree_del()
;

/** Check tree consistency
 *
 * @param [out] is_valid     non zero for correct tree
 */
void fp4est_tree_valid(int * is_valid)
;

/** Save tree to the file
 *
 * @param [in] save_data      if non zero save user data
 * @param filename
 * @param len_f
 */
void fp4est_tree_save(int *save_data, char * filename, int len_f)
;

/** Load tree from a file
 *
 * @param fmpicomm
 * @param [in] load_data           if non zero read user data
 * @param filename
 * @param len_f
 */
void fp4est_tree_load(
		MPI_Fint * fmpicomm, int *load_data, char * filename, int len_f)
;

/** Forest partitioning for p4est
 *
 * @param partforcoarsen   partitioning strategy:
 * 0 - equal element count
 * 1 - octants families prepared for coarsening
 */
void fp4est_part(int * partforcoarsen)
;

/** Write tree structure to VTK file
 *
 * @param filename
 * @param len_f
 */
void fp4est_vtk_write(char * filename, int len_f)
;

/** Write integer scalar field to VTK file
 *
 * @param iscalar
 * @param num
 * @param filename
 * @param len_f
 */
void fp4est_vtk_iscalar(int * iscalar,int *num,char * filename, int len_f)
;

/** Print tree information for testing
 */
void fp4est_test_print()
;

#endif /* NEKP4EST_P4EST_FWRAP_H_ */
