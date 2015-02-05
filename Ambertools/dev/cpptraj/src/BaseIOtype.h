#ifndef INC_BASEIOTYPE_H
#define INC_BASEIOTYPE_H
/// This abstract base class is inherited by any IO object using the FileTypes framework.
class BaseIOtype {
  public:
    typedef BaseIOtype* (*AllocatorType)();
    typedef void (*HelpType)();
};
#endif
