#ifndef INC_MOLSURF_H
#define INC_MOLSURF_H
#ifdef __cplusplus
extern "C" {
#endif
/*! \file molsurf.h
    \brief This is an adaptation of Paul Beroza's "molsurf" program;
 */
// ---------- DEFINES ----------------------------------------------------------
#define REAL_T double
#define MAXAT_EDGE 30
#define NAME_SIZE 8
// size limitations: these are multiplied by the number of atoms
#define NUM_NEIGHBOR 100
#define NUM_PROBE 100
#define NUM_TORUS 10
#define NUM_CIRCLE 30
#define NUM_FACE 20
#define NUM_CYCLE 20
#define NUM_VERTEX 20
#define NUM_EDGE 20
#define NUM_CUSP 20

#define MAXTOR_EDGE       20
#define MAX_FACE_EDGE     20
#define MAX_FACE_CYCLES    4
#define MAXTOR_PROBE      12

// ---------- DATA structures --------------------------------------------------
typedef REAL_T POINT[3];
typedef struct atom {
        POINT pos;
        REAL_T q, rad;
        char anam[NAME_SIZE];
        char rnam[NAME_SIZE];
        int anum;
        int rnum;
        int buried;                   // 1 = buried
        int neighbor_start, n_neighbors;
        int upper_start, n_upper;     // points to neighbor_torus struct
        int torus_start;              // points to torus struct
        // note that the toruslist[] only contains pairs of atoms 
        // in ascending order if you want all tori that contain an 
        // atom you have to look through the whole toruslist 
        int ntorus; 
        int n_convex_edges;           // convex edges associated with atom
        int convex_edges[MAXAT_EDGE]; // convex edges associated with atom 
        int n_cycles;                 // cycles of edges associated with atom
        int cycle_start;              // points to start of cyclelist for the atom
        REAL_T area;                  // accessible surf area associated with the atom
} ATOM;

typedef struct neighbor_torus {
        int iatom;
        int nprobes;    /* -1 = buried, 0 = free, +1 = partially free */
} NEIGHBOR_TORUS;

typedef struct neighbor{
        int iatom;
} NEIGHBOR;

// ---------- torus and probe data structures 
typedef struct {
        POINT center, uv;       /* torus center and axis unit vector */
        REAL_T rad;
        int a1, a2;
        int concave_edges[MAXTOR_EDGE], n_concave_edges;
        int convex_edges[MAXTOR_EDGE], n_convex_edges;
        int circle1, circle2;
        int low;                /* >1 = low, 0 = normal */
} TORUS;

typedef struct probe {
    POINT pos;
        int a1, a2, a3;        /* atoms associated with probe */
        int c1, c2, c3;         /* circles associated with probe */
        REAL_T height;          /* height of probe above base plane */
        int low;                /* >1 = low, 0 = normal */
} PROBE;
/************* edge data structures ************/
typedef struct vertex {
        POINT pos;
        int iatom;      /* atom associated with vertex */
    int iprobe;   /* probe associated with vertex */
        REAL_T beta;   /* concave triangle angle */
} VERTEX;

typedef struct edge {
    int vert1, vert2;
        int circle;
        int alive;      /* 0 = dead (not part of surface) */
} EDGE;

typedef struct circle {
    int torus;
        int atom_or_probe_num;          /* concave = probe  convex = atom */
        REAL_T rad;
        POINT center;
        POINT axis;
} CIRCLE;
/************* face data structures ************/

typedef struct concave_face {  /* each concave face has 3 edges, */
    int e1, e2, e3;            /* edges                              */
        int probe;             /* probe associated with concave face */
        int alive;              /* 1 = active  0 = dead */
        REAL_T area;
} CONCAVE_FACE;

typedef struct saddle_face {   /* edges oriented clockwise when view from "above" torus */
                               /* starting with concave edge */
    int e1_concave, e2_convex, e3_concave, e4_convex;
        int torus;              /* torus associated with saddle face */
        REAL_T area;
        int alive;
} SADDLE_FACE;

typedef struct cycle {
    int nedges;
        int edge[MAX_FACE_EDGE];
        int atom;
        REAL_T area;
} CYCLE;

typedef struct convex_face {
      int n_cycles;           /* 0 or more cycles border a convex face */
      int cycle[MAX_FACE_CYCLES];
      int atom;       /* atom associated with convex face */
      REAL_T area;
} CONVEX_FACE;
/************** Cusp trimming data structures *********************/

typedef struct low_torus {      /* torus that intersects itself */
      int itorus;               /* index in  TORUS[] array */
      int vert1, vert2;/* 2 vertices where the torus intersects itself */
      int nfaces;
      int face[MAXTOR_PROBE];  /* broken_concave_faces associated with low_torus */
      int ncones;
      int cone[MAXTOR_PROBE];  /* cone faces associated with low torus*/
} LOW_TORUS;

typedef struct cone_face {      /* cones are remaining part of saddles from
                                   a self-intersecting torus */
        int e1_convex;          /* only one convex edge */
        int e2_concave, e3_concave;     /* two edges joining convex edge to cusp point */
        int itorus;             /* index in TORUS array         */
        int cusp_vertex;        /* vertex that is cone apex */
        REAL_T area;
} CONE_FACE;

typedef struct broken_concave_face {
        int itorus[3];
        int probe;
        int n_cycles;
        int concave_cycle[MAX_FACE_CYCLES];
        int alive;
        REAL_T area;
} BROKEN_CONCAVE_FACE;
/* concave face resulting from CONCAVE_FACE intersections  */
/* it's possible for face to be split up */
/* edges that form a cycle */

typedef struct concave_cycle {
    int nedges;
        int edge[MAX_FACE_EDGE];
        int edge_direction[MAX_FACE_EDGE];      /* 1 right hand rule; -1 left hand */
        int cusp_edge[MAX_FACE_EDGE];           /* cusp edge index, -1 = ordinary edge */
        int iprobe;
        int iface;
        int intersects_self;                    /* 1 = 2 cusps intersect 0 = no intersecting cusps */
        REAL_T area;
} CONCAVE_CYCLE;

typedef struct cusp_edge {
        int cycle1;             /* first  cycle */
        int cycle2;             /* second cycle */
        int edge;                               /* concave edge index */
        int probe1, probe2;                     /* probes that form cusp edge */
        int alive;                              /* 1 = alive; 0 = dead */
        int concentric_pair;                    /* 1 = intersects w/ a conc. cusp -> */
                                                /* can't intersect w/ anyone else    */
}CUSP_EDGE;
/* new_cusp is for storing non_axial cusps, which may or may not
   become real cusps, depending on the non_axial cusps that it intersects */

typedef struct cusp_pair {
        POINT circle_center, circle_axis, vert1, vert2;
        REAL_T circle_rad;
        int cycle1, cycle2, cycle3;
        int cusp1, cusp2;               /*  cusps that intersect at new_cusp */
        int group;                      /* group to which new_cusp belogs */
} CUSP_PAIR;

// ---------- FUNCTIONS --------------------------------------------------------
void Set_molsurf_debug(int);
void memory_usage (); 
REAL_T molsurf(REAL_T, ATOM*, int,NEIGHBOR_TORUS *, NEIGHBOR *,
               TORUS *, PROBE *, CONCAVE_FACE *,
               SADDLE_FACE *, CONVEX_FACE *,
               CONE_FACE *, BROKEN_CONCAVE_FACE *,
               CONCAVE_CYCLE *, VERTEX *,
               EDGE *, EDGE *,
               CIRCLE *, CIRCLE *,
               CYCLE *, LOW_TORUS *, CUSP_EDGE *,
               CUSP_PAIR *);

#ifdef __cplusplus
}
#endif
#endif
