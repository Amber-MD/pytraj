#ifndef INC_NAMETYPE_H
#define INC_NAMETYPE_H
#include <cstddef> // size_t
#include <string>
class NameType {
  public:
    NameType();
    NameType(const NameType&);
    NameType(const char*);
    NameType(std::string const&);
    NameType& operator=(const NameType&);

    void ToBuffer(char*) const;
    bool Match(NameType const&) const;
    bool operator==(const NameType&) const;
    bool operator==(const char*) const;
    bool operator!=(const NameType&) const;
    bool operator!=(const char*) const;
    const char* operator*() const { return c_array_; }
    char operator[](int) const;
    std::string Truncated() const;
    void ReplaceAsterisk();

  private:
    const size_t NameSize_;
    char c_array_[6];

    void FormatName();
};  
#endif
