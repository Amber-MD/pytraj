#ifndef INC_ACTION_VOLMAP_H
#define INC_ACTION_VOLMAP_H
#include "Action.h"
#include "DataSet_GridFlt.h"
class Action_Volmap : public Action {
  public:
    Action_Volmap();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Volmap(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    static void RawHelp();
    /// gets the LJ radius for a given atom from a topology
    double GetRadius_(Topology const&, int);

    int ensembleNum_;
    /// grid resolutions
    double dx_, dy_, dz_;
    /// minimum values in the x-, y-, and z-dimensions
    double xmin_, ymin_, zmin_;
    /// number of frames we analyzed so we can average at the end
    int Nframes_;
    /// If true, set up the grid on first frame based on centermask.
    bool setupGridOnMask_;
    /// mask to center the grid on
    AtomMask centermask_;
    /// mask of atoms to grid
    AtomMask densitymask_;
    /// the grid we are using
    DataSet_GridFlt* grid_;
    /// file name with the peak locations as Carbons in XYZ file format
    std::string peakfilename_;
    /// The value below which to ignore all peaks
    double peakcut_;
    /// the atomic radii of each atom in the gridded selection
    std::vector<float> halfradii_;
    /// the clearance between the edges of our grid and centermask_
    double buffer_;
    /// the scaling factor to divide all radii by
    double radscale_;
    static const double sqrt_8_pi_cubed;
    static const double one_over_6;
};
#endif
