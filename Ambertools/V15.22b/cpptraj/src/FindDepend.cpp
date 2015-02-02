/*
 * FindDepend 
 * Dan Roe 2010
 * Given a file with #include directives, print a list of dependencies.
 */
#include <cstdio>
#include <cctype>
#include <cstring>
#include <list>
#include <string>

using namespace std;

#define BUFFERSIZE 1024

// Return true if first goes before second
bool compareNames( string first, string second) {
  if ( first.compare(second) < 0)
    return true;
  else
    return false;
}

/*
 * Given a file, put all headers into a list.
 */
list<string> *PeakHeader(char *filename, int indent, bool includeStdHeaders) {
  FILE *infile;
  char buffer[BUFFERSIZE];
  char headername[BUFFERSIZE];
  char *ptr, *includePosition;
  size_t pos;
  list<string> *HeaderList;
  list<string> *SecondList;
  list<string> SpliceList;
  list<string>::iterator it;
  string temp;

  // Safety valve
  if (indent > 100) {
    fprintf(stdout,"Recursion > 100! Bailing out!\n");
    return NULL;
  }

  infile = fopen(filename,"r");
  if (infile==NULL) {
    fprintf(stderr,"Could not open %s\n",filename);
    return NULL;
  }
  HeaderList = new list<string>();

  // Go through file and grab includes
  while ( fgets(buffer, BUFFERSIZE, infile)!=NULL ) {
    ptr=buffer;
    // skip leading whitespace - necessary??
    while ( isspace(*ptr) ) {
      ptr++;
      if (*ptr=='\n' || *ptr=='\0') break;
    }
    // First non-whitespace charcter must be '#'
    if (*ptr != '#') continue;
    // Ignore commented lines
//    if (strncmp(ptr,"//",2)==0) continue;
    // Check if this is an include line - if not, skip
    includePosition = strstr(ptr,"include");
    if ( includePosition==NULL ) continue;
//    if ( includePosition==NULL || strchr(ptr,'#')==NULL ) continue;
    //if ( strncmp(ptr,"#include",8)!=0 ) continue; 
    // Record standard libs but dont process them
    if ( strchr(ptr, '<')==NULL ) {
      // Get header name - assume it is second string - dont include first "
      sscanf(includePosition, "%*s \"%s", headername);
      // Get rid of last "
      pos=strlen(headername);
      if (headername[pos-1]=='"') headername[pos-1]='\0';
    } else {
      if (includeStdHeaders)
        sscanf(includePosition,"%*s %s", headername);
      else
        continue;
    }
    // Check for system headers that might be in quotes
    if (!includeStdHeaders) {
      if (strcmp(headername,"mpi.h")==0) continue;
      if (strcmp(headername,"zlib.h")==0) continue;
      if (strcmp(headername,"bzlib.h")==0) continue;
      if (strcmp(headername,"netcdf.h")==0) continue;
      if (strcmp(headername,"omp.h")==0) continue;
      if (strncmp(headername,"readline",8)==0) continue;
    }
    temp.assign(headername);
    HeaderList->push_back(temp);
  }
  fclose(infile);

  // DEBUG
  //fprintf(stdout,"%i:%s ",indent,filename);
  //for (it=HeaderList->begin(); it!=HeaderList->end(); it++) 
  //  fprintf(stdout,"[%s]",(*it).c_str());
  //fprintf(stdout,"\n");

  // Go through each header in the list, skipping standard headers
  // Use a copy of HeaderList so the iterator wont go on forever
  SpliceList = *HeaderList;
  for (it=SpliceList.begin(); it!=SpliceList.end(); it++) {
    if ( (*it).find('<')==string::npos ) {
      SecondList = PeakHeader((char*)(*it).c_str(),indent+1,includeStdHeaders);
      if (SecondList!=NULL) {
        HeaderList->splice( HeaderList->end(), *SecondList );
        delete SecondList;
      }
    }
  }
  // Assign this filename, only if indent>0
  if (indent>0) {
    temp.assign(filename);
    HeaderList->push_front(temp);
  }

  return HeaderList;
}

// Print a list of dependencies for the given file
void GetDependencies(char *filename) {
  list<string> *HeaderList;
  list<string>::iterator it;
  size_t pos;
  char tempname[BUFFERSIZE];


  HeaderList = PeakHeader(filename,0,false);
  if (HeaderList==NULL) return;
  HeaderList->sort(compareNames);
  // Remove duplicates
  HeaderList->unique();

  // Show the initial filename
  strcpy(tempname, filename);
  // Replace cpp with o
  pos = strlen(filename);
  if ( filename[pos-1] == 'p' &&
       filename[pos-2] == 'p' &&
       filename[pos-3] == 'c' &&
       filename[pos-4] == '.'    ) {
    tempname[pos-3]='o';
    tempname[pos-2]=' ';
    tempname[pos-1]=':';
    tempname[pos]=' ';
    tempname[pos+1]='\0';
    strcat(tempname,filename);
  }
  // Replace c with o
  else if ( filename[pos-1] == 'c' &&
            filename[pos-2] == '.' ) { 
    tempname[pos-1]='o';
    tempname[pos]=' ';
    tempname[pos+1]=':';
    tempname[pos+2]=' ';
    tempname[pos+3]='\0';
    strcat(tempname,filename);
  }
  // Skip f files for now
  else if ( filename[pos-1] == 'f' &&
            filename[pos-2] == '.' ) {
    return;
  }
  fprintf(stdout,"%s",tempname);

  // Print the headers
  for (it=HeaderList->begin(); it!=HeaderList->end(); it++)
    fprintf(stdout," %s",(*it).c_str());
  fprintf(stdout,"\n");
  delete HeaderList;
}

// M A I N
int main(int argc, char **argv) {
  int i;

  if (argc<2) return 0;

  for (i=1; i<argc; i++)
    GetDependencies(argv[i]);

  return 0;
}  
