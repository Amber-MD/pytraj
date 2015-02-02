#ifndef INC_ACTION_HBOND_H
#define INC_ACTION_HBOND_H
#include <vector>
#include <map>
#include <set>
#include "Action.h"
#include "ImagedAction.h"
#include "DataSet_integer.h"
// Class: Action_Hbond
/// Action to calculate the Hbonds present in each frame.
class Action_Hbond : public Action {
  public:
    Action_Hbond();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Hbond(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    struct HbondType {
      int A;        ///< Acceptor atom#
      int H;        ///< Hydrogen atom#
      int D;        ///< Donor atom#
      int Frames;   ///< # frames this hbond has been present
      double dist;  ///< Used to calc avg distance of this hbond
      double angle; ///< Used to calc avg angle of this hbond
      DataSet_integer* data_; ///< If series, keep track of frames hbond is present.
    };

    ImagedAction Image_;
    Matrix_3x3 ucell_, recip_;
    int debug_;
    int ensembleNum_;
    int Nframes_;
    std::string avgout_;
    std::string solvout_;
    std::string bridgeout_;
    bool useAtomNum_;
    typedef std::map<int,HbondType> HBmapType;
    HBmapType HbondMap_;   ///< Track all solute-solute hbonds found
    HBmapType SolventMap_; ///< Track all solute-solvent hbonds found
    typedef std::map< std::set<int>, int > BridgeType;
    BridgeType BridgeMap_; ///< Track all combos of residues bridged by solvent.
    typedef std::vector<int> HBlistType;
    HBlistType Donor_;                 ///< Array of hbond donor atoms (D0, H0, D1, H1, ...)
    HBlistType Acceptor_;              ///< Array of hbond acceptor atoms (A0, A1, ...)
    HBlistType SolventDonor_;
    HBlistType SolventAcceptor_;
    AtomMask Mask_;
    AtomMask DonorMask_;
    AtomMask DonorHmask_;
    AtomMask AcceptorMask_;
    AtomMask SolventDonorMask_;
    AtomMask SolventAcceptorMask_;
    bool hasDonorMask_;
    bool hasDonorHmask_;
    bool hasAcceptorMask_;
    bool hasSolventDonor_;
    bool hasSolventAcceptor_;
    bool calcSolvent_;
    bool noIntramol_;
    double acut_;
    double dcut2_;
    Topology* CurrentParm_;

    bool series_;
    std::string hbsetname_;
    DataSet* NumHbonds_;
    DataSet* NumSolvent_;
    DataSet* NumBridge_;
    DataSet* BridgeID_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    /// Return true if the first hbond has more frames than the second.
    /** If both have the same # of frames, pick something arbitrary just to
      * be sure that we have a well-defined ordering (otherwise we could get
      * spurious test failures) -- Order equivalent frames based on atom number
      * of Acceptor.
      */
    struct hbond_cmp {
      inline bool operator()(HbondType const& first, HbondType const& second) const {
        if (first.Frames == second.Frames)
          return (first.dist < second.dist);
        else
          return (first.Frames > second.Frames);
      }
    };
    /// Return true if first bridge has more frames than second.
    struct bridge_cmp {
      inline bool operator()(std::pair< std::set<int>, int> const& first, 
                             std::pair< std::set<int>, int> const& second) const
      {
        if (first.second > second.second)
          return true;
        else if (first.second < second.second)
          return false;
        else
          return (first.second < second.second);
      }
    };

    void SearchAcceptor(HBlistType&,AtomMask&,bool);
    void SearchDonor(HBlistType&,AtomMask&,bool,bool);
    inline int AtomsAreHbonded(Frame const&, int, int, int, int, int,bool);
    inline void HbondTypeCalcAvg(HbondType&);
};
#endif
