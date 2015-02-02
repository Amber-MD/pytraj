#ifndef INC_ACTION_DIFFUSION_H
#define INC_ACTION_DIFFUSION_H
#include "Action.h"
class Action_Diffusion : public Action {
  public:
    Action_Diffusion();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Diffusion(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    Frame initial_;
    std::vector<double> previousx_;
    std::vector<double> previousy_;
    std::vector<double> previousz_;
    bool printIndividual_;
    double time_;
    bool hasBox_;
    int debug_;
    std::vector<double> distancex_;
    std::vector<double> distancey_;
    std::vector<double> distancez_;
    std::vector<double> distance_;
    std::vector<double> deltax_;
    std::vector<double> deltay_;
    std::vector<double> deltaz_;
    AtomMask mask_;
    CpptrajFile outputx_;
    CpptrajFile outputy_;
    CpptrajFile outputz_;
    CpptrajFile outputr_;
    CpptrajFile outputa_;
    Vec3 boxcenter_; ///< Hold center of box each frame
};
#endif
