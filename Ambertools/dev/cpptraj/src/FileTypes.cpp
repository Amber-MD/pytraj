#include <set>
#include "FileTypes.h"
#include "CpptrajStdio.h"
// FileTypes::GetFormatFromArg()
FileTypes::FileFormatType FileTypes::GetFormatFromArg(KeyPtr begin, ArgList& argIn,
                                                      FileFormatType def)
{
  for (KeyPtr token = begin; token->Key != 0; ++token)
    if (argIn.hasKey( token->Key )) return token->Type;
  return def;
}
// FileTypes::GetFormatFromString()
FileTypes::FileFormatType FileTypes::GetFormatFromString(KeyPtr begin, std::string const& fmt,
                                                         FileFormatType def)
{
  for (KeyPtr token = begin; token->Key != 0; ++token)
    if ( fmt.compare( token->Key )==0 ) return token->Type;
  return def;
}
// FileTypes::GetExtensionForType()
std::string FileTypes::GetExtensionForType(KeyPtr begin, FileFormatType typeIn) {
  for (KeyPtr token = begin; token->Extension != 0; ++token)
    if ( token->Type == typeIn )
      return std::string( token->Extension );
  return std::string();
}
// FileTypes::GetTypeFromExtension()
FileTypes::FileFormatType FileTypes::GetTypeFromExtension(KeyPtr begin, std::string const& extIn,
                                                          FileFormatType def)
{
  for (KeyPtr token = begin; token->Extension != 0; ++token)
    if ( extIn.compare( token->Extension ) == 0 ) return token->Type;
  return def;
}
// -----------------------------------------------------------------------------
// FileTypes::FormatDescription()
const char* FileTypes::FormatDescription(AllocPtr allocArray, FileFormatType typeIn) {
  return allocArray[ typeIn ].Description;
}
// FileTypes::AllocIO()
/** \param silent If true do not print error message if allocator does not exist.
  */
BaseIOtype* FileTypes::AllocIO(AllocPtr allocArray, FileFormatType typeIn, bool silent) {
  if (allocArray[typeIn].Alloc == 0) {
    if (!silent)
      mprinterr("Error: CPPTRAJ was compiled without support for %s files.\n",
                allocArray[typeIn].Description);
    return 0;
  }
  return allocArray[typeIn].Alloc();
}
// FileTypes::FormatKeywords()
std::string FileTypes::FormatKeywords(KeyPtr begin, FileFormatType ftype) {
  std::string keywords;
  std::set<std::string> Keys;
  for (KeyPtr token = begin; token->Key != 0; ++token)
    if (token->Type == ftype) Keys.insert( std::string(token->Key) );
  if (!Keys.empty()) {
    keywords.assign("Keywords:");
    for (std::set<std::string>::const_iterator key = Keys.begin();
                                               key != Keys.end(); ++key)
      keywords.append(" " + *key);
  }
  return keywords;
}
// FileTypes::FormatExtensions()
std::string FileTypes::FormatExtensions(KeyPtr begin, FileFormatType ftype) {
  std::string extensions;
  std::set<std::string> Exts;
  for (KeyPtr token = begin; token->Extension != 0; ++token)
    if (token->Type == ftype) Exts.insert( std::string(token->Extension) );
  if (!Exts.empty()) {
    extensions.assign("Extensions:");
    for (std::set<std::string>::const_iterator ext = Exts.begin();
                                               ext != Exts.end(); ++ext)
      extensions.append(" '" + *ext + "'");
  }
  return extensions;
}
// FileTypes::ReadOptions()
void FileTypes::ReadOptions(KeyPtr begin, AllocPtr allocArray, FileFormatType UNK) {
  for (int i = 0; i < UNK; i++) {
    std::string fmtExtensions = FormatExtensions(begin, i);
    if (allocArray[i].ReadHelp || !fmtExtensions.empty()) {
      mprintf("    Options for %s:", allocArray[i].Description);
      if (!fmtExtensions.empty()) mprintf(" %s", fmtExtensions.c_str());
      mprintf("\n");
      if (allocArray[i].ReadHelp != 0) allocArray[i].ReadHelp();
    }
  }
}
// FileTypes::WriteOptions()
void FileTypes::WriteOptions(KeyPtr begin, AllocPtr allocArray, FileFormatType UNK) {
  for (int i = 0; i < UNK; i++) {
    std::string fmtExtensions = FormatExtensions(begin, i);
    std::string fmtKeywords =  FormatKeywords(begin, i);
    if (allocArray[i].WriteHelp || !fmtExtensions.empty() || !fmtKeywords.empty()) {
      mprintf("    Options for %s:", allocArray[i].Description);
      if (!fmtKeywords.empty()) mprintf(" %s,", fmtKeywords.c_str());
      if (!fmtExtensions.empty()) mprintf(" %s", fmtExtensions.c_str()); 
      mprintf("\n");
      if (allocArray[i].WriteHelp != 0) allocArray[i].WriteHelp();
    }
  }
}
