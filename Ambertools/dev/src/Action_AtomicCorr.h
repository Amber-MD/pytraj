#ifndef INC_ACTION_ATOMICCORR_H
#define INC_ACTION_ATOMICCORR_H
#include "Action.h"
class Action_AtomicCorr : public Action {
  public:
    Action_AtomicCorr();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_AtomicCorr(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    class AtomVector {
      public:
        AtomVector() : idx_(0) {}
        AtomVector(std::string const& sIn, int idxIn) : lbl_(sIn), idx_(idxIn) {}
        int operator-(const AtomVector& rhs) { return idx_ - rhs.idx_; }
        bool empty()               { return vec_.empty();    }
        size_t size()              { return vec_.size();     }
        void push_back(float fval) { vec_.push_back( fval ); }
        std::string const& Label() { return lbl_;            }
        Vec3 VXYZ(int idx)         { return Vec3((double)vec_[idx  ], 
                                                 (double)vec_[idx+1], 
                                                 (double)vec_[idx+2]); }
      private:
        std::vector<float> vec_;
        std::string lbl_;
        int idx_;
    };

    enum AcorrModeType { ATOM = 0, RES };
    static const char* ModeString[];
    AcorrModeType acorr_mode_;
    double cut_;
    int min_;
    int debug_;
    DataSet* dset_;
    DataFile* outfile_;

    typedef std::vector< AtomVector > ACvector;
    ACvector atom_vectors_;
    AtomMask mask_;
    std::vector<AtomMask> resmasks_;
    Frame refframe_;
};
#endif
