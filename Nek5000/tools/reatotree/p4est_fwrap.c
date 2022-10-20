/*
 * p4est_fwrap.c
 * Fortran interface for p4est library.
 *
 *  Created on: Feb 21, 2016
 *      Author: Adam Peplinski
 */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "nekp4est.h"

#if N_DIM == 2
#undef P4_TO_P8
#else
#include <p4est_to_p8est.h>
#endif

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_iterate.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_iterate.h>
#endif

#include "name.h"

#include "p4est_fwrap.h"

/* Global variables
 *  notice I use tree_nek->user_pointer to store ghost quadrants data
 */
static p4est_connectivity_t *connect_nek = NULL; /**< Nek5000 connectivity structure */
static p4est_t *tree_nek = NULL; /**< Nek5000 tree structure */

/** @brief Initialize element data for 0 level tree
 *
 * @details This is dummy routine required by fp4est_tree_new.
 * I do not initialize tree data here, as there are special routines dedicated
 * to transfer data between Nek5000 and p4est, however it could be used if
 * one would like to generate connectivity and tree directly in Nek5000
 *
 * @param p4est
 * @param which_tree
 * @param quadrant
 */
void init_msh_dat(p4est_t * p4est, p4est_topidx_t which_tree,
		  p4est_quadrant_t * quadrant) {
  user_data_t *data = (user_data_t *) quadrant->p.user_data;
  int iwt, il;

  extern void nekp4est_init_msh_dat(int * iwt, int * imsh, int * igrp,
			       int (*)[P4EST_FACES], int (*)[P4EST_FACES]);

  iwt = (int) which_tree;
  nekp4est_init_msh_dat(&iwt, &data->imsh,&data->igrp,&data->bc,&data->crv);
  
  // reset refinement mar, nek5000 elemnt distribution information and refinement history
  data->ref_mark = 0;
  data->el_gln = -1;
  data->el_ln = -1;
  data->el_nid = -1;
  data->parent_gln = -1;
  data->parent_ln = -1;
  data->parent_nid = -1;
  for(il=0;il<P4EST_CHILDREN;++il){
    data->children_gln[il] = -1;
    data->children_ln[il] = -1;
    data->children_nid[il] = -1;
  }
}

/* Wrappers */
/* initialize */
void fp4est_init(int * log_threshold) {
	p4est_init(NULL, *log_threshold);
}

/* Connectivity */
#ifdef P4_TO_P8
void fp4est_cnn_new(int * num_vertices, int * num_trees, int * num_edges,
		int * num_corners, double *vertices, int * tree_to_vertex,
		int * tree_to_tree, int8_t * tree_to_face, int * tree_to_edge,
		int * ett_offset, int * edge_to_tree, int8_t * edge_to_edge,
		int * tree_to_corner, int * ctt_offset, int * corner_to_tree,
		int8_t * corner_to_corner)
#else
void fp4est_cnn_new(int * num_vertices, int * num_trees,
		int * num_corners,
		double *vertices,
		int * tree_to_vertex, int * tree_to_tree,
		int8_t * tree_to_face,
		int * tree_to_corner, int * ctt_offset,
		int * corner_to_tree, int8_t * corner_to_corner)
#endif
{
#ifdef P4_TO_P8
  connect_nek = p4est_connectivity_new_copy(*num_vertices, *num_trees, *num_edges, *num_corners,
                                            vertices, tree_to_vertex, tree_to_tree,
                                            tree_to_face, tree_to_edge, ett_offset, edge_to_tree,
                                            edge_to_edge, tree_to_corner, ctt_offset, corner_to_tree,
                                            corner_to_corner);
#else
  connect_nek = p4est_connectivity_new_copy(*num_vertices, *num_trees, *num_corners,
                                            vertices, tree_to_vertex, tree_to_tree, tree_to_face,
                                            tree_to_corner, ctt_offset, corner_to_tree, corner_to_corner);
#endif
}

void fp4est_cnn_del() {
	if (connect_nek) p4est_connectivity_destroy(connect_nek);
	connect_nek =  NULL;
}

void fp4est_cnn_attr(int * enable_tree_attr) {
	p4est_connectivity_set_attr(connect_nek, *enable_tree_attr);
}

void fp4est_cnn_valid(int * is_valid) {
	*is_valid = p4est_connectivity_is_valid(connect_nek);
}

void fp4est_cnn_complete() {
	p4est_connectivity_complete(connect_nek);
}

void fp4est_cnn_save(char *filename, int len_f) {
	p4est_connectivity_save(filename, connect_nek);
}

void fp4est_cnn_load(char * filename, int len_f) {
	connect_nek = p4est_connectivity_load(filename, NULL);
}

/* tree_ management */
void fp4est_tree_new(MPI_Fint * fmpicomm, int * min_level) {
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	tree_nek = p4est_new_ext(mpicomm, connect_nek, 0, *min_level, 1,
			sizeof(user_data_t), init_msh_dat, NULL);
}

void fp4est_tree_del() {
	if (tree_nek) p4est_destroy(tree_nek);
	tree_nek = NULL;
}

void fp4est_tree_valid(int * is_valid) {
	*is_valid = p4est_is_valid(tree_nek);
}

void fp4est_tree_save(int *save_data, char * filename, int len_f) {
	p4est_save_ext(filename, tree_nek, *save_data,0);
}

void fp4est_tree_load(MPI_Fint * fmpicomm, int *load_data, char * filename,
		int len_f) {
	MPI_Comm mpicomm;
	mpicomm = MPI_Comm_f2c(*fmpicomm);
	tree_nek = p4est_load_ext(filename, mpicomm, sizeof(user_data_t), *load_data,1,0,
	NULL, &connect_nek);
	tree_nek->user_pointer = NULL;
}

/* p4est internal load balance */
void fp4est_part(int * partforcoarsen) {
	p4est_partition(tree_nek, *partforcoarsen, NULL);
}

/* I/O (VTK) */
void fp4est_vtk_write(char * filename, int len_f) {
	p4est_vtk_write_file(tree_nek, NULL, filename);
}

void fp4est_vtk_iscalar(int * iscalar, int *num, char * filename, int len_f) {
	int il, jl, err;
	double rscalar[*num * P4EST_CHILDREN];

	for (il = 0; il < *num; ++il) {
		for (jl = 0; jl < P4EST_CHILDREN; jl++) {
			rscalar[il * P4EST_CHILDREN + jl] = (double) iscalar[il];
		}
	}

	err = p4est_vtk_write_header(tree_nek, NULL, 1.0, 0, 0, 0, 0, "scalar",
	NULL, filename);
	err = p4est_vtk_write_point_scalar(tree_nek, NULL, filename, "scalar",
			rscalar);
	err = p4est_vtk_write_footer(tree_nek, filename);
}


/* Iterate over element volumes to print element data */
void iter_print(p4est_iter_volume_info_t * info, void *user_data) {
  user_data_t *data = (user_data_t *) info->quad->p.user_data;

  // which quad (local and global element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  int iwlt, iwg;
  int il;// loop index

  // get quad number
  tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
  // local quad number
  iwl = info->quadid + tree->quadrants_offset;
  iwlt = (int) iwl;
  // global quad number
  iwg = (int) info->p4est->global_first_quadrant[info->p4est->mpirank] + iwlt;
  printf("Elem %i %i %i \n",iwg,data->imsh,data->igrp);
  for(il=0;il<P4EST_FACES;++il){
     printf("   Face %i %i %i\n",il,data->crv[il],data->bc[il]);
  }
}

/* print tree information for testing */
void fp4est_test_print(){
  p4est_ghost_t *ghost_nek = NULL;
  ghost_nek = p4est_ghost_new(tree_nek, P4EST_CONNECT_FULL);
#ifdef P4_TO_P8
  p4est_iterate(tree_nek,ghost_nek,NULL,iter_print,NULL, NULL, NULL);
#else
  p4est_iterate(tree_nek,ghost_nek,NULL,iter_print,NULL, NULL);
#endif
  p4est_ghost_destroy(ghost_nek);
}

