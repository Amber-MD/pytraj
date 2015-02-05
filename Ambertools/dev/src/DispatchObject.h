#ifndef INC_DISPATCHOBJECT_H
#define INC_DISPATCHOBJECT_H
/// Abstract base class that all dispatchable objects will inherit.
class DispatchObject {
  public:
    /// Function pointer for allocating a DispatchObject. 
    typedef DispatchObject* (*DispatchAllocatorType)();
};
#endif
