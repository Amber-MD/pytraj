#ifndef INC_ARGLIST_H
#define INC_ARGLIST_H
#include <vector>
#include <string>
// Class: ArgList
/// Hold a list of string arguments and keeps track of their usage.
/** Can be set from an input line using SetList(), with arguments separated 
  * by a specified delimiter, or arguments can be added one-by-one with AddArg.
  * Arguments can be accessed with the various getX routines,
  * where X is specific to certain types, e.g. getNextDouble returns
  * the next double, getNextMask returns an atom mask expression (i.e.
  * it has :, @, % characters etc). All of the getX routines (along with
  * the hasKey routine) mark the argument they access as used, so that
  * subsequent calls with these functions will not return the same
  * argument over and over. 
  */
class ArgList {
  public:
    // Constructors
    ArgList() {}
    ArgList(const char*);
    ArgList(std::string const&);
    ArgList(std::string const&, const char*);
    // Copy/Assignment
    ArgList(const ArgList&);
    ArgList& operator=(const ArgList &);
    /// \return the argument at the given position
    std::string const& operator[](int) const;
    /// \return Internal argument list as vector of strings
    std::vector<std::string> const& List() const { return arglist_; }
    // Iterators
    typedef std::vector<std::string>::const_iterator const_iterator;
    const_iterator begin() const { return arglist_.begin();     }
    const_iterator end()   const { return arglist_.end();       }
    /// \return the number of arguments
    int Nargs()            const { return (int)arglist_.size(); }
    /// \return true if no arguments in list.
    bool empty()           const { return arglist_.empty();     }
    /// \return the argument string
    const char *ArgLine()  const { return argline_.c_str();     }
    std::string const& ArgString() const { return argline_;     }
    /// Clear list
    void ClearList();
    /// Set up argument list from string and given separators
    int SetList(std::string const&, const char *);
    /// \return an argument list of remaining unmarked args.
    ArgList RemainingArgs();
    /// Add argument to the list
    void AddArg(std::string const&);
    /// Mark given argument
    void MarkArg(int);
    /// Print a warning if not all arguments are marked
    bool CheckForMoreArgs() const;
    /// Print the argument list
    void PrintList() const;
    /// Print detailed info for arg list
    void PrintDebug() const;
    /// Remove the first argument
    void RemoveFirstArg();
    /// \return the first argument
    const char *Command() const;
    /// \return true if the first argument matches key
    bool CommandIs(const char*) const;
    /// \return the next unmarked string
    std::string const& GetStringNext();
    /// \return the next unmarked mask
    std::string const& GetMaskNext();
    /// \return the next unmarked tag
    std::string const& getNextTag();
    /// \return true if arg at position is valid integer.
    bool ValidInteger(int) const;
    /// \return true if arg at position is valid double.
    bool ValidDouble(int) const;
    /// \return the next unmarked integer
    int getNextInteger(int);
    /// \return the next unmarked double
    double getNextDouble(double);
    /// \return the string following the given key
    std::string const& GetStringKey(const char *);
    /// \return the integer following the given key 
    int getKeyInt(const char *, int);
    /// \return the double following the given key
    double getKeyDouble(const char*, double);
    /// \return true if the key is present in the list
    bool hasKey(const char*);
    /// \return true if they key is in the list but do not mark.
    bool Contains(const char*) const;
  private:
    /// Empty string to return when args not found
    static const std::string emptystring;
    /// The original argument string (complete list)
    std::string argline_;
    /// List of arguments
    std::vector<std::string> arglist_;
    /// Mark which arguments have been used
    std::vector<bool> marked_;
};
#endif
