#ifndef INC_REFERENCEFRAME_H
#define INC_REFERENCEFRAME_H
#include "Topology.h"
#include "ArgList.h"
/// Hold single frame coordinates to be used as a reference.
/** The frame and topology are stored as pointers instead of classes to
  * save memory; this way multiple actions can use the same reference 
  * structure without each having to have a different copy. Because
  * of this, memory is not freed in ReferenceFrame destructor to avoid
  * potential double-frees when ReferenceFrame is used in e.g. vectors.
  * Freeing must be accomplished with the ClearRef function.
  */
class ReferenceFrame {
  public:
    ReferenceFrame() : frame_(0), parm_(0), num_(0) {}
    ReferenceFrame(int) : frame_(0), parm_(0), num_(-1) {}
    ~ReferenceFrame();
    Frame const& Coord()        const { return *frame_;      }
    Topology const& Parm()      const { return *parm_;       }
    bool error()                const { return num_ == -1;   }
    bool empty()                const { return frame_ == 0;  }
    FileName const& FrameName() const { return name_;        }
    std::string const& Tag()    const { return tag_;         }
    int LoadRef(std::string const&, Topology*, int);
    int LoadRef(std::string const&, ArgList&, Topology*, std::string const&, int);
    int StripRef( AtomMask const& );
    void RefInfo() const;
    void ClearRef();
  private:
    Frame* frame_;      ///< Reference coords, allocated.
    Topology* parm_;    ///< Pointer to assiociated parm in TopologyList.
    FileName name_;     ///< Ref structure filename.
    std::string tag_;   ///< Ref structure optional tag.
    int num_;           ///< Frame number.
};
#endif
