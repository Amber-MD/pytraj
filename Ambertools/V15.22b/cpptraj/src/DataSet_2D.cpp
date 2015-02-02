#include "DataSet_2D.h"

const DataSet_2D::MatrixToken DataSet_2D::TokenArray[] = {
  { "UNDEFINED",                       "UNKNOWN"   },
  { "distance matrix",                 "DIST"      },
  { "covariance matrix",               "COVAR"     },
  { "mass-weighted covariance matrix", "MWCOVAR"   },
  { "correlation matrix",              "CORREL"    },
  { "distance covariance matrix",      "DISTCOVAR" },
  { "IDEA matrix",                     "IDEA"      },
  { "IRED matrix",                     "IRED"      },
  { "dihedral covariance matrix",      "DIHCOVAR"  },
  { 0,                                 0           }
};
