#ifndef INC_ACTION_NATIVECONTACTS_H
#define INC_ACTION_NATIVECONTACTS_H
#include <map>
#include "Action.h"
#include "ImagedAction.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_integer.h"
/// Calculate the number of native/non-native contacts based on distance
/** Intended to combine and replace contacts, mindist, and maxdist actions.
  */
class Action_NativeContacts : public Action {
  public:
    Action_NativeContacts();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NativeContacts(); }
    static void Help();
  private:
    typedef std::vector<int> Iarray;
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    Iarray SetupContactIndices(AtomMask const&, Topology const&);
    int SetupContactLists(Topology const&, Frame const&);
    int DetermineNativeContacts(Topology const&, Frame const&);
    inline bool ValidContact(int, int, Topology const&) const;

    double distance_;     ///< Cutoff distance
    float pdbcut_;        ///< Only print pdb atoms with bfac > pdbcut.
    int debug_;           ///< Action debug level.
    int ensembleNum_;
    int matrix_min_;      ///< Used for map output
    int resoffset_;       ///< When byResidue, ignore residues spaced this far apart
    unsigned int nframes_;///< Number of frames, for normalizing map
    bool first_;          ///< If true use first frame as reference
    bool byResidue_;      ///< If true calculate distances by residue
    bool includeSolvent_; ///< If true include solvent residues
    bool series_;         ///< If true save time series of native contacts.
    bool usepdbcut_;      ///< If true only print pdb atoms with bfac > pdbcut.
    ImagedAction image_;  ///< Hold imaging-related info/routines.
    AtomMask Mask1_;      ///< First mask in which to search
    AtomMask Mask2_;      ///< Second mask in which to search
    Iarray contactIdx1_;  ///< Hold atom/residue indices for Mask1 (for map)
    Iarray contactIdx2_;  ///< Hold atom/residue indices for Mask2 (for map)
    std::string cfile_;   ///< File to write native contact list to.
    std::string pfile_;   ///< File to write contact PDB to.
    std::string rfile_;   ///< File to write total fraction frames for res pairs.
    DataSet* numnative_;  ///< Hold # of native contacts
    DataSet* nonnative_;  ///< Hold # of non-native contacts
    DataSet* mindist_;    ///< Hold minimum observed distance among contacts
    DataSet* maxdist_;    ///< Hold maximum observed distance among contacts
    DataSet_MatrixDbl* nativeMap_; ///< Hold native contacts map
    DataSet_MatrixDbl* nonnatMap_; ///< Hold non-native contacts map
    Topology* CurrentParm_;
    Frame refFrame_;      ///< For printing out contact PDB.
    const Topology* refParm_;   ///< For printing out contact PDB.
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    Matrix_3x3 ucell_, recip_;
    /// Define contact, atom pair.
    class contactType;
    /// Define contact pair.
    typedef std::pair<int,int> Cpair;
    /// Define map pair for inserting into map.
    typedef std::pair<Cpair, contactType> Mpair;
    /// Define list of contacts.
    typedef std::map<Cpair, contactType> contactListType;
    contactListType nativeContacts_; ///< List of native contacts.
    /// Hold residue total contact frames and total # contacts.
    class resContact {
    public:
    resContact() : nframes_(0), ncontacts_(0) {}
    resContact(int nf) : nframes_(nf), ncontacts_(1) {}
    void Increment(int nf) { nframes_ += nf; ++ncontacts_; }
    int Nframes() const { return nframes_; }
    int Ncontacts() const { return ncontacts_; }
    bool operator<(resContact const& rhs) const {
      if (nframes_ == rhs.nframes_)
        return (ncontacts_ > rhs.ncontacts_);
      else
        return (nframes_ > rhs.nframes_);
    }
    private:
    int nframes_, ncontacts_;
    };
    /// For holding residue pair and total fraction contact.
    typedef std::pair<Cpair, resContact> Rpair;
    /// For sorting residue contact pairs.
    struct res_cmp {
      inline bool operator()(Rpair const& first, Rpair const& second) const {
        return (first.second < second.second);
      }
    };
};
// ----- PRIVATE CLASS DEFINITIONS ---------------------------------------------
class Action_NativeContacts::contactType {
  public:
    contactType() : dist_(0.0), dist2_(0.0), data_(0), nframes_(0), res1_(-1), res2_(-1) {}
    contactType(std::string const& id, int r1, int r2) : dist_(0.0), dist2_(0.0), data_(0),
                                         id_(id), nframes_(0), res1_(r1), res2_(r2) {}
    const char* id() const { return id_.c_str(); }
    int Nframes()    const { return nframes_;    }
    int Res1()       const { return res1_;       }
    int Res2()       const { return res2_;       }
    double Avg()     const { return dist_;       }
    double Stdev()   const { return dist2_;      }
    DataSet_integer& Data() { return *data_;     }
    void Increment(int fnum, double d, double d2) {
      nframes_++;
      dist_ += d;
      dist2_ += d2;
      if (data_!=0) data_->AddVal(fnum, 1);
    }
    void Finalize();
    bool operator<(contactType const& rhs) const {
      if (nframes_ == rhs.nframes_)
        return (dist_ < rhs.dist_);
      else
        return (nframes_ > rhs.nframes_);
    }
    void SetData(DataSet* ds) { data_ = (DataSet_integer*)ds; }
  private:
    double dist_;    ///< (For avg) contact distance when present.
    double dist2_;   ///< (For stdev) contact distance^2 when present.
    DataSet_integer* data_;  ///< If series, keep track of frames contact is present.
    std::string id_; ///< Contact ID.
    int nframes_;    ///< Number of frames contact is present.
    int res1_;
    int res2_;
};
#endif
