#ifndef SOAX_SOLVER_H_
#define SOAX_SOLVER_H_

#include "itkFEMLinearSystemWrapperItpack.h"
#include "global.h"

namespace soax {

/*
 * A bank of linear system solver with various orders. It is used for
 * evolution of snakes with different sizes. Each solver is randomly
 * accessible by the order of the system. This solver bank can be
 * dynamically expanded to accomodate the linear system with the
 * largest order.
 */
class SolverBank {
 public:
  SolverBank();
  ~SolverBank();

  void Reset();
  /*
   * Solve the linear system with order of vectors.size().
   */
  void SolveSystem(const VectorContainer &vectors, unsigned dim, bool open);

  double GetSolution(unsigned order, unsigned index, bool open);

  double alpha() const {return alpha_;}
  void set_alpha(double a) {alpha_ = a;}

  double beta() const {return beta_;}
  void set_beta(double b) {beta_ = b;}

  double gamma() const {return gamma_;}
  void set_gamma(double g) {gamma_ = g;}

 private:
  typedef itk::fem::LinearSystemWrapperItpack SolverType;
  typedef std::vector<SolverType *> SolverContainer;

  /*
   * Delete all solvers in a SolverContainer and release the memory.
   */
  void ClearSolvers(SolverContainer &solvers);

  void ExpandSolverContainer(SolverContainer &solvers, unsigned position);

  void InitializeSolver(SolverType *solver, unsigned order, bool open);

  void FillMatrixOpen(SolverType *solver, unsigned order);

  void FillMatrixClosed(SolverType *solver, unsigned order);

  void DestroySolutionVectors(SolverContainer &solvers);

  /*
   * Solvers for open snakes.
   */
  SolverContainer open_solvers_;
  /*
   * Solvers for closed snakes.
   */
  SolverContainer closed_solvers_;

  /*
   * Weight for first order continuity of snakes.
   */
  double alpha_;

  /*
   * Weight for second order continuity of snakes.
   */
  double beta_;

  /*
   * Step size of iteration.
   */
  double gamma_;

  DISALLOW_COPY_AND_ASSIGN(SolverBank);
};

} // namespace soax

#endif // SOAX_SOLVER_H_
