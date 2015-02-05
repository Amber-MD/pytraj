#ifndef INC_ACTION_MATRIX_H
#define INC_ACTION_MATRIX_H
#include "Action.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_Vector.h"
#include "Array1D.h"
#include "ActionFrameCounter.h"
/// Calculate various types of matrices
class Action_Matrix : public Action, ActionFrameCounter {
  public:
    Action_Matrix();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Matrix(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    typedef DataSet_MatrixDbl::Darray Darray;     ///< Mass/vector array type.
    typedef Darray::iterator          v_iterator; ///< Iterator over vector.
    typedef Darray::const_iterator    M_iterator; ///< Iterator over mass.

    DataSet_MatrixDbl* Mat_;
    DataFile* outfile_;
    AtomMask mask1_;
    AtomMask mask2_;
    std::string filename_;
    enum OutputType { BYATOM=0, BYRESIDUE, BYMASK };
    OutputType outtype_;
    int snap_;
    int debug_;
    int ensembleNum_;
    // IRED only
    int order_;                                ///< Legendre order
    std::vector<DataSet_Vector*> IredVectors_; ///< IRED vectors
    // DIHCOVAR only
    Array1D DihedralSets_;
    // MWcovar only
    Darray mass1_;
    Darray mass2_;

    Darray vect2_; ///< Hold diagonal elements squared.
#   ifdef _OPENMP
    /// For OPENMP only, save coord indices (X-Y) for speed 
    std::vector<int> crd_indices_;
#   endif
    bool useMask2_;
    bool useMass_;
    Topology* CurrentParm_; // For ByResidue output

    Darray FillMassArray(Topology const&, AtomMask const&) const;
    void CalcIredMatrix(int);
    void CalcDistanceMatrix(Frame const&);
    inline void StoreVec(v_iterator&, v_iterator&, const double*) const;
    void CalcCovarianceMatrix(Frame const&);
    void StoreXY(v_iterator&, v_iterator&, const double*) const;
    void CalcDihedralCovariance( int );
    void CalcIdeaMatrix(Frame const&);
    void CalcCorrelationMatrix(Frame const&);
    void CalcDistanceCovarianceMatrix(Frame const&);
    void Vect2MinusVect();
    void FinishCovariance(size_t);
    inline void DotProdAndNorm(DataSet_MatrixDbl::iterator&, v_iterator&,
                               v_iterator&, v_iterator&, v_iterator&) const;
    void FinishCorrelation();
    void FinishDistanceCovariance();
};
#endif
