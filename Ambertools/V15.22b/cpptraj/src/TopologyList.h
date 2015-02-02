#ifndef INC_TOPOLOGYLIST_H
#define INC_TOPOLOGYLIST_H
#include "Topology.h"
#include "ArgList.h"
// Class: TopologyList
/// Holds a list of Topology classes.
/** Can either add new topology by filename, or add existing topology by 
  * address. Can search for topology in list by index, full/base filename,
  * or tag.
  */
class TopologyList {
  public:
    static const char* ParmArgs;
    TopologyList();
    ~TopologyList();
    void Clear();
    void SetDebug(int);
    Topology* GetParm(int) const;
    Topology* GetParmByIndex(ArgList&) const;
    Topology* GetParm(ArgList&) const;
    int AddParmFile(std::string const&);
    int AddParmFile(std::string const&,ArgList&);
    void AddParm(Topology* pIn) { pIn->SetPindex(TopList_.size()); TopList_.push_back( pIn ); }
    void List() const;
  private:
    std::vector<Topology*> TopList_;
    int debug_;
};
#endif
