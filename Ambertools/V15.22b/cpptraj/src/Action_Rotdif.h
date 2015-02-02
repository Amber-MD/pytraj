#ifndef INC_ACTION_ROTDIF_H
#define INC_ACTION_ROTDIF_H
#include "Action.h"
#include "Random.h"
#include "DataSet_Vector.h"
// Class: Action_Rotdif
/// Estimate rotational diffusion tensors from MD simulations
/** To estimate rotational diffusion tensors from MD simulations along the
  * lines described by Wong & Case, (Evaluating rotational diffusion from
  * protein MD simulations, J. Phys. Chem. B 112:6013, 2008) using the 
  * following procedure:
  * - Create random vectors
  * - Calculate rotation matrices via RMS fitting to a reference
  * - Apply the rotation matrices to the random vectors
  * - computes a best-fit rotational diffusion tensor, plus various statistics 
  *   on this. There are two optimizations: one giving the optimal diffusion
  *   tensor in the small-anisotropy limit, and the second optimizing the 
  *   diffusion tensor with the full expression.
  * \author Original code: Vance Wong, George Giambasu
  * \author Adapted by: DRR
  */
class Action_Rotdif: public Action {
  public:
    Action_Rotdif();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Rotdif(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    int debug_;
    int rseed_;          ///< Random seed
    int nvecs_;          ///< Number of random vectors to generate
    double tfac_;        ///< time step
    double ti_;          ///< initial time (ns)
    double tf_;          ///< final time (ns), should be less than ncorr * tfac
    int NmeshPoints_;    ///< # cubic spline mesh points 
    int itmax_;          ///< Max number of iterations
    double delmin_;      ///< Convergence criterion
    double d0_;          ///< Initial guess for iso diffusion tensor
    int olegendre_;      ///< order of Legendre polynomial in the correlation function
    int ncorr_;          ///< Max length to compute time correlation fns (# frames)
    double delqfrac_;    ///< how to scale simplexes
    double amoeba_ftol_; ///< Simplex min tolerance
    int amoeba_itmax_;   ///< Simplex min iterations
    int amoeba_nsearch_; ///< Number of simplex min searches
    bool do_gridsearch_; ///< If true perform grid search after simplex min.
    bool useMass_;       ///< Use mass in rms fitting.
    bool usefft_;

    // Workspace for LAPACK functions
    Matrix_3x3 D_tensor_;
    Vec3 D_XYZ_;

    std::string randvecOut_;
    std::string randvecIn_;
    std::string rmOut_;
    std::string deffOut_;
    std::string corrOut_;

    Frame SelectedRef_;
    AtomMask TargetMask_;
    Frame SelectedTgt_;
    CpptrajFile outfile_;

    // Variables used by the random number generator
    Random_Number RNgen_;

    std::vector<Matrix_3x3> Rmatrices_; ///< Store rotation matrices
    DataSet_Vector random_vectors_;     ///< Hold nvecs random vectors
    std::vector<double> D_eff_;         ///< Hold calculated effective D values for each vector
//    std::vector<double> sumc2_;      

    DataSet_Vector RandomVectors();
    int direct_compute_corr(DataSet_Vector const&, int, std::vector<double>&);
    int fft_compute_corr(DataSet_Vector const&, int, std::vector<double>&);
    double calcEffectiveDiffusionConst(double );

    static void PrintMatrix(CpptrajFile&, const char*, Matrix_3x3 const&);
    static void PrintVector(CpptrajFile&, const char*, Vec3 const&);
    static void PrintVec6(CpptrajFile&, const char*, std::vector<double> const&);
    void PrintTau( std::vector<double> const& );
    int Tensor_Fit(std::vector<double>&);
    int DetermineDeffs();
    void PrintDeffs(std::string const&) const;

    int DetermineDeffsAlt();
};
#endif  
