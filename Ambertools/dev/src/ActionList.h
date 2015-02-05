#ifndef INC_ACTIONLIST_H
#define INC_ACTIONLIST_H
#include "Action.h"
// Class: ActionList
/// Hold actions that will be performed every frame.
/** This class is responsible for holding all actions that will be performed
  * during the course of trajectory processing.
  */
class ActionList {
  public:
    ActionList();
    ~ActionList();
    /// Clear the list
    void Clear();
    /// Set the debug level for actions.
    void SetDebug(int);
    void SetSilent(bool b) { actionsAreSilent_ = b; }
    int Debug() const { return debug_; }
    /// Add given action to the action list and initialize.
    int AddAction(DispatchObject::DispatchAllocatorType, ArgList&,
                  TopologyList*,DataSetList*,DataFileList*);
    /// Set up actions for the given parm.
    int SetupActions(Topology **);
    /// Perform actions on the given frame.
    bool DoActions(Frame **, int);
    /// Call print for each action.
    void Print();
    /// List all actions in the action list.
    void List() const;
    bool Empty() const { return actionlist_.empty(); }
    // The functions below help set up actions when ensemble processing.
    /// \return the number of actions in the list.
    int Naction() const { return actionlist_.size(); }
    /// \return command string for corresponding action.
    std::string const& CmdString(int i) const { return actioncmd_[i];              }
    /// \return new action corresponding to existing action.
    DispatchObject::DispatchAllocatorType
      ActionAlloc(int i)                const { return actionAlloc_[i]; }
  private:
    /// Action initialization and setup status.
    enum ActionStatusType { NO_INIT=0, INIT, SETUP, INACTIVE };
    typedef std::vector<Action*> Aarray;
    /// List of actions
    Aarray actionlist_;
    /// List of action commands
    std::vector<std::string> actioncmd_;
    /// List of action allocators (for ensemble).
    std::vector<DispatchObject::DispatchAllocatorType> actionAlloc_;
    /// List of action statuses
    std::vector<ActionStatusType> actionstatus_;
    /// Default debug level for actions
    int debug_;
    /// If true suppress all init/setup output from actions.
    bool actionsAreSilent_;
};
#endif
