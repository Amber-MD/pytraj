#include "Residue.h"
#include <cctype> // tolower

char Residue::ConvertResName(std::string const& r) {
  if (r.compare(0,3,"ALA")==0) return 'A';
  if (r.compare(0,3,"ARG")==0) return 'R';
  if (r.compare(0,3,"ASN")==0) return 'N';
  if (r.compare(0,3,"ASP")==0) return 'D';
  if (r.compare(0,3,"ASH")==0) return 'D'; // Protonated ASP
  if (r.compare(0,3,"CYS")==0) return 'C';
  if (r.compare(0,3,"CYM")==0) return 'C'; // Deprotonated CYS
  if (r.compare(0,3,"CYX")==0) return 'C';
  if (r.compare(0,3,"GLN")==0) return 'Q';
  if (r.compare(0,3,"GLU")==0) return 'E';
  if (r.compare(0,3,"GLH")==0) return 'E'; // Protonated GLU
  if (r.compare(0,3,"GLY")==0) return 'G';
  if (r.compare(0,3,"HIS")==0) return 'H';
  if (r.compare(0,3,"HIE")==0) return 'H'; // NE-protonated (HIS)
  if (r.compare(0,3,"HID")==0) return 'H'; // ND-protonated
  if (r.compare(0,3,"HIP")==0) return 'H'; // NE/ND protonated
  if (r.compare(0,3,"ILE")==0) return 'I';
  if (r.compare(0,3,"LEU")==0) return 'L';
  if (r.compare(0,3,"LYS")==0) return 'K';
  if (r.compare(0,3,"LYN")==0) return 'K'; // Deprotonated (neutral) LYS 
  if (r.compare(0,3,"MET")==0) return 'M';
  if (r.compare(0,3,"PHE")==0) return 'F';
  if (r.compare(0,3,"PRO")==0) return 'P';
  if (r.compare(0,3,"SER")==0) return 'S';
  if (r.compare(0,3,"THR")==0) return 'T';
  if (r.compare(0,3,"TRP")==0) return 'W';
  if (r.compare(0,3,"TYR")==0) return 'R';
  if (r.compare(0,3,"VAL")==0) return 'V';
  // Make lower case letter when unrecognized.
  if (!r.empty()) return tolower(r[0]);
  return ' ';
}
