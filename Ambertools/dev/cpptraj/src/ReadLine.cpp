#include <cstdio> // required for readline
#include <cstdlib> // free
#include <cstring>
#define READLINE_LIBRARY
#include <readline.h>
#include <history.h>
#include "ReadLine.h"
#include "Command.h"

// duplicate_string()
static char* duplicate_string(const char* s) {
  char* r = (char*)malloc(strlen(s) + 1);
  strcpy(r, s);
  return r;
}

// command_generator()
/** Generator function for command completion.  STATE lets us know whether
  * to start from scratch; without any state (i.e. STATE == 0), then we
  * start at the top of the list.
  */
static char* command_generator(const char* text, int state) {
  static int list_index, len;
  const char *name;

  // If this is a new word to complete, initialize now. This includes
  // saving the length of TEXT for efficiency, and initializing the index
  // variable to 0.
  if (!state) {
    list_index = 0;
    len = strlen(text);
  }

  // Return the next name which partially matches from the command list.
  while ( (name = Command::CmdToken(list_index).Cmd) != 0 )
  {
    list_index++;
    if (strncmp(name, text, len) == 0)
      return (duplicate_string(name));
  }

  // If no names matched, then return NULL.
  return 0;
}

// cpptraj_completion()
static char** cpptraj_completion(const char* text, int start, int end) {
  char** matches = 0;
  // If this word is at the start of the line, assume it is a command.
  if (start == 0 || (strncmp(rl_line_buffer, "help ", 5)==0))
    matches = rl_completion_matches(text, command_generator);
  return matches;
}

// CONSTRUCTOR
ReadLine::ReadLine() {
  // Tell the completer that we want a crack first.
  rl_attempted_completion_function = cpptraj_completion;
}

// -----------------------------------------------------------------------------
/** Get next input line with readline. Lines terminated with a backslash
  * will be concatenated. Comments will be ignored.
  */
int ReadLine::GetInput() {
  input_.clear();
  char* line = readline("> ");
  if (line == 0) return 1; // EOF
  input_ += line;
  // Terminal backslash requests a continuation of the line
  size_t end = strlen( line );
  while (end > 1 && line[end - 1] == '\\') {
    // Remove that backlash
    size_t bs_pos = input_.find_last_of('\\');
    input_.erase( bs_pos, 1 );
    free( line );
    line = readline("");
    if (line == 0) break;
    input_ += line;
    end = strlen( line );
  }
  // Remove leading whitespace.
  std::string::iterator beg = input_.begin();
  while ( beg != input_.end() && isspace(*beg) )
    beg = input_.erase(beg);
  // Find '#' not preceded by blackslash; indicates comment.
  // Remove it and all after.
  end = input_.find_first_of('#');
  if (end != std::string::npos) {
    if (end == 0 || (end > 0 && input_[end-1] != '\\'))
      input_.erase( input_.begin() + end, input_.end() );
  }
  // Add line to history
  if (!input_.empty()) AddHistory(input_.c_str());
  if (line != 0) free( line );
  return 0;
}

void ReadLine::AddHistory(const char* line) {
  if (line != 0) add_history( line );
}

bool ReadLine::YesNoPrompt(const char* prompt) {
  char* line = readline(prompt);
  if (line == 0 || strlen( line ) < 1) return false;
  if (line[0] == 'y' || line[0] == 'Y') return true;
  return false;
}
