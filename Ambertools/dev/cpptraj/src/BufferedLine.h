#ifndef INC_BUFFEREDLINE_H
#define INC_BUFFEREDLINE_H
#include <vector>
#include "CpptrajFile.h"
/// Used to buffer text files that will be read line-by-line
class BufferedLine : private CpptrajFile {
  public:
    BufferedLine();
    ~BufferedLine();
    /// Get the next line in the buffer.
    const char* Line();
    /// Convert current line into tokens
    int TokenizeLine(const char*);
    /// \return next token, null-delimited.
    const char* NextToken();
    /// \return specified token, not null-delimited.
    inline const char* Token(int);
    /// Open file for reading, set up buffer.
    int OpenFileRead( std::string const& fname ) {
      if ( OpenRead( fname ) ) return 1;
      return ResetBuffer();
    }
    int LineNumber()          const { return nline_;          }
    const char* Buffer()      const { return buffer_;         }
    // Pointer to current buffer position.
    const char* CurrentLine() const { return bufferPosition_; }
    inline std::string GetLine();
    // Members of CpptrajFile that should be public
    using CpptrajFile::Filename;
    using CpptrajFile::CloseFile;
  private:
    int ResetBuffer();
    static const size_t DEFAULT_BUFFERSIZE = 16384;

    char* buffer_;         ///< Character buffer
    char* bufferPosition_; ///< Position in buffer/start of current line.
    /// Array of pointers to beginning and ends of tokens in current line. 
    std::vector<char*> tokens_;
    size_t tokenidx_;      ///< Current position in tokens array
    char saveChar_;        ///< Saved last char of current token
    char* lineEnd_;        ///< End of current line in buffer
    char endChar_;         ///< Character that was at *lineend
    char* endBuffer_;      ///< End position of buffer
    size_t nline_;         ///< Current line number.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
const char* BufferedLine::Token(int idx) {
  if (idx < 0 || idx >= (int)tokens_.size()) return 0;
  return tokens_[idx];
}

std::string BufferedLine::GetLine() {
  const char* ptr = Line();
  if (ptr == 0) return std::string();
  return std::string(ptr);
}
#endif
