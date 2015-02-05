#ifndef INC_ACTION_MOLSURF_H
#define INC_ACTION_MOLSURF_H
#include "Action.h"
#include "molsurf.h"
// Class: Action_Molsurf
/// Wrapper for the molsurf routine in molsurf.c
/** \author Original C code by Paul Beroza.
  * \author C++ adaptation by Dan Roe.
  * This is the cpptraj wrapper for the molsurf routine originally written
  * by Paul Beroza. This calculates the Connolly surface area of the
  * molecule. See "M.L. Connolly, Analytical molecular surface calculation,
  * J. Appl. Cryst., 16, p. 548-558, 1983."
  */
class Action_Molsurf: public Action {
  public:
    Action_Molsurf();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Molsurf(); }
    static void Help();
    ~Action_Molsurf();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* sasa_;
    AtomMask Mask1_;
    ATOM* atom_;
    double probe_rad_;
    double rad_offset_;
    typedef std::vector<AtomMask> Marray;
    Marray SubMasks_;
    typedef std::vector<DataSet*> DSarray;
    DSarray SubData_;
    std::vector<int> mask1idx_;
    // Molsurf internal data structs
  /* neighbor arrays:  these are big so amount of data stored must be small
   * upper_neighbors is of the NEIGHBOR_TORUS type, which contains 2
   * small ints, one for the index of the neighbor the other for the
   * torus.  This last index is key.  The torus associated with
   * the neighbor is unique (the index of the upper_neighbor atom is always 
   * larger than the atom that points to it).  If the index for this torus
   * is -1, then no torus has been instantiated yet, so we don't have to
   * allocated the memory for it.  It the index is -2 the torus is buried.
   * Any other index corresponds to an index in the torus array that has
   * already been instantiated.
   */
    NEIGHBOR_TORUS *upper_neighbors; ///< contains atoms and torus indices
    NEIGHBOR *neighbors;             ///< contains atom indices for all neighbors
    TORUS *toruslist;
    PROBE *probelist;

    CONCAVE_FACE *concave_face;
    SADDLE_FACE *saddle_face;
    CONVEX_FACE *convex_face;
    CONE_FACE *cone_face;
    BROKEN_CONCAVE_FACE *broken_concave_face;
    CONCAVE_CYCLE *concave_cycle;

    VERTEX *vertexlist;
    EDGE *concave_edge_list;
    EDGE *convex_edge_list;
    CIRCLE *convex_circle_list;
    CIRCLE *concave_circle_list;

    CYCLE *cyclelist;
    LOW_TORUS *low_torus;
    CUSP_EDGE *cusp_edge;
    CUSP_PAIR *cusp_pair;
    // -----------------------------
    int AllocateMemory();
    void ClearMemory();
};
#endif
