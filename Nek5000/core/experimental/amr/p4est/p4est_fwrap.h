/** @file p4est_fwrap.h
 *  @brief Fortran77 interface of p4est library
 *
 * @author Adam Peplinski
 * @date Feb 20, 2016
 *
 * @ingroup nekamr
 */
#ifndef NEKP4EST_P4EST_FWRAP_H_
#define NEKP4EST_P4EST_FWRAP_H_

#if !defined(NEKAMR_H_)
#error "p4est_fwrap.h requires nekamr.h"
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
/* tree and grid info */
#define fp4est_ghost_new     FORTRAN_NAME(fp4est_ghost_new,FP4EST_GHOST_NEW)
#define fp4est_ghost_del     FORTRAN_NAME(fp4est_ghost_del,FP4EST_GHOST_DEL)
#define fp4est_mesh_new      FORTRAN_NAME(fp4est_mesh_new,FP4EST_MESH_NEW)
#define fp4est_mesh_del      FORTRAN_NAME(fp4est_mesh_del,FP4EST_MESH_DEL)
#define fp4est_nodes_new     FORTRAN_NAME(fp4est_nodes_new,FP4EST_NODES_NEW)
#define fp4est_nodes_del     FORTRAN_NAME(fp4est_nodes_del,FP4EST_NODES_DEL)
#define fp4est_lnodes_new    FORTRAN_NAME(fp4est_lnodes_new,FP4EST_LNODES_NEW)
#define fp4est_lnodes_del    FORTRAN_NAME(fp4est_lnodes_del,FP4EST_LNODES_DEL)
/* p4est internal load balance */
#define fp4est_part          FORTRAN_NAME(fp4est_part,FP4EST_PART)
/* refinement, coarsening, balance */
#define fp4est_refine        FORTRAN_NAME(fp4est_refine,FP4EST_REFINE)
#define fp4est_coarsen       FORTRAN_NAME(fp4est_coarsen,FP4EST_COARSEN)
#define fp4est_balance       FORTRAN_NAME(fp4est_balance,FP4EST_BALANCE)
/* I/O (VTK) */
#define fp4est_vtk_write     FORTRAN_NAME(fp4est_vtk_write,FP4EST_VTK_WRITE)
#define fp4est_vtk_iscalar   FORTRAN_NAME(fp4est_vtk_iscalar,FP4EST_VTK_ISCALAR)
/* routines required by Nek5000 for data exchange */
#define fp4est_msh_get_size  FORTRAN_NAME(fp4est_msh_get_size,FP4EST_MSH_GET_SIZE)
#define fp4est_msh_get_dat   FORTRAN_NAME(fp4est_msh_get_dat,FP4EST_MSH_GET_DAT)
#define fp4est_msh_get_hst   FORTRAN_NAME(fp4est_msh_get_hst,FP4EST_MSH_GET_HST)
#define fp4est_msh_get_node  FORTRAN_NAME(fp4est_msh_get_node,FP4EST_MSH_GET_NODE)
#define fp4est_msh_get_tplg  FORTRAN_NAME(fp4est_msh_get_tplg,FP4EST_MSH_GET_TPLG)
#define fp4est_msh_get_lnode FORTRAN_NAME(fp4est_msh_get_lnode,FP4EST_MSH_GET_LNODE)
#define fp4est_msh_get_algn  FORTRAN_NAME(fp4est_msh_get_algn,FP4EST_MSH_GET_ALGN)
#define fp4est_msh_get_graph FORTRAN_NAME(fp4est_msh_get_graph,FP4EST_MSH_GET_GRAPH)
#define fp4est_refm_put      FORTRAN_NAME(fp4est_refm_put,FP4EST_REFM_PUT)
#define fp4est_egmap_put     FORTRAN_NAME(fp4est_egmap_put,FP4EST_EGMAP_PUT)
#define fp4est_crd_cnt_get   FORTRAN_NAME(fp4est_crd_cnt_get,FP4EST_CRD_CNT_GET)
#define fp4est_crd_vrt_get   FORTRAN_NAME(fp4est_crd_vrt_get,FP4EST_CRD_VRT_GET)
#define fp4est_bc_check      FORTRAN_NAME(fp4est_bc_check,FP4EST_BC_CHECK)
#define fp4est_tree_copy     FORTRAN_NAME(fp4est_tree_copy,FP4EST_TREE_COPY)
#define fp4est_tree_check    FORTRAN_NAME(fp4est_tree_check,FP4EST_TREE_CHECK)
#define fp4est_family_get    FORTRAN_NAME(fp4est_family_get,FP4EST_FAMILY_GET)

/** Data type for user variables; required by p4est */
typedef struct user_data_s {
	int imsh; /**< velocity (0) and temperature (1) mesh indicator */
	int igrp; /**< element group */
	int crv[P4EST_FACES]; /**< curvature data; concerns external faces; requred by face projection */
	int bc[P4EST_FACES]; /**< boundary condition data; 0 internal, -1 periodic, otherwise external */
	int ref_mark; /**< integer to store refinement mark; definition in nekamr.h */
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
/** Initialise p4est package setting log verbosity
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

/** Build ghost layer */
void fp4est_ghost_new()
;

/** Destroy ghost layer */
void fp4est_ghost_del()
;

/** Generate mesh information */
void fp4est_mesh_new()
;

/** Destroy mesh information */
void fp4est_mesh_del()
;

/** Generate new node information */
void fp4est_nodes_new()
;

/** Destroy node information */
void fp4est_nodes_del()
;

/** Generate new global nodes (GLL points) numbering
 *
 * @param degree   polynomial degree
 */
void fp4est_lnodes_new(int * degree)
;

/** Destroy global node numbering */
void fp4est_lnodes_del()
;

/** Forest partitioning for p4est
 *
 * @param partforcoarsen   partitioning strategy:
 * 0 - equal element count
 * 1 - octants families prepared for coarsening
 */
void fp4est_part(int * partforcoarsen)
;

/** Perform tree refinement
 *
 * @param max_level    max refinement level
 */
void fp4est_refine(int *max_level)
;

/** Perform tree coarsening
 */
void fp4est_coarsen()
;

/** Perform 2:1 tree balancing
 */
void fp4est_balance()
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

/** Get mesh size information to Nek5000
 *
 * @param [out] nelgt   global element number
 * @param [out] nelgit  element offset (number of elements on lower nid's)
 * @param [out] nelt    number of T-type elements
 * @param [out] nelv    number of V-type elements
 * @param [out] maxl    current max level
 */
void fp4est_msh_get_size(
		int * nelgt,int * nelgit, int * nelt, int * nelv, int * maxl)
;

/** Get mesh data to Nek5000
 *
 * @param[in]  nelv  global number of V-type elements
 * @param[in]  lelt  array dimension for bc, cbc
 * @param[out] igrp  element group
 * @param[out] level element level
 * @param[out] crv   face projection flag
 * @param[out] bc    boundary surface flag
 */
void fp4est_msh_get_dat(int * nelv, int *lelt, int * igrp,
		int * level, int * crv, int* bc)
;

/** Get refinement/coarsening data to Nek5000
 *
 * @param map_nr     local number of unchanged elements
 * @param rfn_nr     local number of refined elements
 * @param crs_nr     local number of coarsened elements
 * @param elgl_map   element mapping info for unchanged elements
 * @param elgl_rfn   element mapping info for refined elements
 * @param elgl_crs   element mapping info for coarsened elements
 */
void fp4est_msh_get_hst(int * map_nr, int * rfn_nr, int * crs_nr, int *elgl_map,
		int * elgl_rfn, int * elgl_crs)
;

/** Get global information on hanging vertices to Nek5000
 *
 * @param [out] lnelt   number of local elements
 * @param [out] node    array with hanging vertices information
 *                      (0-independent, 1-face hanging, 2- edge hanging)
 */
void fp4est_msh_get_node(int * lnelt, int * node)
;

/** Get global vertex, face and edge numbering to Nek5000 for mesh topology
 *
 * @param lnelt      local number of elements
 * @param lnoden     local number of nodes (vert., fac., edg.)
 * @param gnoden     local number of owned nodes
 * @param lnodes     global node numbering list
 * @param hang_elm   is any face/edge hanging list
 * @param hang_fsc   hanging face list
 * @param hang_edg   hanging edge list
 */
#ifdef P4_TO_P8
void fp4est_msh_get_tplg(
		int * lnelt, int * lnoden, int * gnoden,p4est_gloidx_t * lnodes,
		int * hang_elm, int * hang_fsc, int * hang_edg);
#else
void fp4est_msh_get_tplg(
		int * lnelt, int * lnoden, int * gnoden,p4est_gloidx_t * lnodes,
		int * hang_elm, int * hang_fsc);
#endif

/** Get global vertex, face and edge numbering to Nek5000
 *
 * @param lnelt      local number of elements
 * @param degree     polynomial order
 * @param lnodes     global node numbering list
 */
void fp4est_msh_get_lnode(
		int * lnelt, int * degree, p4est_gloidx_t * lnodes);


/** Get face orientation
 *
 * @param fcs_algn   face alignment
 * @param lnelt      local number of elements
 */
void fp4est_msh_get_algn(int * fcs_algn, int * lnelt);

/** Get adjacency graph for partitioning
 *
 * @param node_num      node number in the graph (to create vtxdist in ParMetis notation)
 * @param graph         graph (adjncy in ParMetis notation)
 * @param graph_offset  graph_offset (xadj in ParMetis notation)
 * @param vrt_wgt       vertex weights (vwgt in ParMetis notation)
 * @param edg_wgt       edge weights (vwgt in ParMetis notation)
 * @param npt           number of points on the edge (polynomial order + 1)
 */
void fp4est_msh_get_graph(int * node_num, int * graph,
		int * graph_offset, int * vrt_wgt, int * edg_wgt, int * npt)
;

/** Fill ref_mark in p4est block
 *
 * @param ref_mark   refinement mark array
 */
void fp4est_refm_put(int * ref_mark)
;

/** Fill element global mapping in p4est block
 *
 * @param el_gnum   element global mapping array
 * @param el_lnum   element global mapping array
 * @param el_nid   element global mapping array
 */
void fp4est_egmap_put(int * el_gnum,int * el_lnum,int * el_nid)
;


/** Get approximate physical coordinates of the p4est block centre
 *
 * @param coord   physical coordinates
 */
void fp4est_crd_cnt_get(double * coord)
;

/** Get approximate physical coordinates of the p4est block vertices
 *
 * @param coord   physical coordinates
 */
void fp4est_crd_vrt_get(double * coord)
;

/** Check boundary conditions for V- and T-type mesh
 */
void fp4est_bc_check()
;

/** Make tree copy for later comparison
 *
 * @param quad_data   do we test quadrant data
 */
void fp4est_tree_copy(int *quad_data)
;

/** Check if tree was modified
 *
 * @param check       tree modification marker
 * @param quad_data   do we test quadrant data
 */
void fp4est_tree_check(int * check, int *quad_data)
;

/** Provide information about element families
 *
 * @param family array containing global element numbers in the family
 * @param nelf   number of entrances in family array
 */
void fp4est_family_get(int * family, int * nelf)
;

#endif /* NEKP4EST_P4EST_FWRAP_H_ */
