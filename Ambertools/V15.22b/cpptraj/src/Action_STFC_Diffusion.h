#ifndef INC_ACTION_STFC_DIFFUSION_H
#define INC_ACTION_STFC_DIFFUSION_H
#include "Action.h"
/** \author Hannes H. Loeffler
  * \author C++ adaptation by Daniel R. Roe
  */
class Action_STFC_Diffusion : public Action {
  public:
    Action_STFC_Diffusion();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_STFC_Diffusion(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    void calculateMSD(const double*,int,int,Vec3 const&);

    int ensembleNum_;
    bool printDistances_; // iarg1
    enum CalcType { DEFAULT = 0, COM, DIST };
    CalcType calcType_; // iarg2
    enum DirectionType { DX = 0, DY, DZ, DXY, DXZ, DYZ, DXYZ };
    DirectionType direction_; // iarg3
    AtomMask mask_;
    AtomMask mask2_;
    CpptrajFile output_;
    CpptrajFile outputnw_;
    std::string outputAverDist_;
    std::string outputNumWat_;
    double time_;
    double lowerCutoff_;
    double upperCutoff_;
    bool hasBox_;
    int n_atom_;

    typedef std::vector<double> Darray;
    Darray initialxyz_;
    Darray distancexyz_;
    Darray distance_;
    Darray deltaxyz_;
    Darray previousxyz_;

    Darray dSum1_;
    Darray dSum2_;
    std::vector<int> nInside_;
    int elapsedFrames_;
};
#endif    
