#ifndef PHOTON_SOLVER_H
#define PHOTON_SOLVER_H

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <mpi.h>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>


#include <Amesos2.hpp>
#include <Amesos2_Solver.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosSolverManager.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>

#include <vector>

#include "mesh_data.h"
#include "constants.h"
#include "species.h"
#include "field_solver.h"
#include "run_parameters.h"

class PhotonSolver
{
	typedef Tpetra::Vector<>::scalar_type scalar_type;
	typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
	typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
	typedef Tpetra::Vector<>::node_type node_type;
	typedef Tpetra::Map<> map_type;
	typedef Tpetra::CrsMatrix<> crs_matrix_type;
	typedef Tpetra::MultiVector<> mv_type;
	typedef Belos::LinearProblem<scalar_type, mv_type, Tpetra::Operator<> > problem_type;
	typedef Belos::SolverManager<scalar_type, mv_type, Tpetra::Operator<> > belos_manager;
	typedef Belos::PseudoBlockGmresSolMgr<scalar_type, mv_type, Tpetra::Operator<> > belos_block_gm_res_manager;
	typedef Ifpack2::Preconditioner<> prec_type;
	typedef Tpetra::RowMatrix<> row_matrix_type;

	public:
		// default constructor
		PhotonSolver();
		PhotonSolver(const Teuchos::RCP<const Teuchos::Comm<int> >&, RunParameters *);
		void InitializeSolver();
		void ConfigureSolver();
		void updateBC(scalar_type);
		void Solve();
		void PrintMatrixA();
		void PrintMatrixB();
		void PrintMatrixX();
		void UpdateIonizationRate(std::vector<Species *> &, FieldSolver *);
		int GetYSize() { return my_y_size; }

		scalar_type GetPhotonRateAt(local_ordinal_type index) { return (data[index]+data2[index]+data3[index]); }
		scalar_type GetPhotonRateAt(local_ordinal_type x, local_ordinal_type y) { return (data[x + y*rp->x_size]+data2[x+y*rp->x_size]+data3[x+y*rp->x_size]); }

		int getNumberIterationsForSolve() { return belos_solver->getNumIters(); }	


	private:
		Teuchos::RCP<crs_matrix_type> A;
		Teuchos::RCP<mv_type> X;
		Teuchos::RCP<mv_type> B;
		Teuchos::RCP<crs_matrix_type> A2;
		Teuchos::RCP<mv_type> X2;
		Teuchos::RCP<mv_type> B2;
		Teuchos::RCP<crs_matrix_type> A3;
		Teuchos::RCP<mv_type> X3;
		Teuchos::RCP<mv_type> B3;
		Teuchos::RCP<const map_type> global_map;

		Tpetra::global_size_t num_global_entries;
		Teuchos::RCP<const Teuchos::Comm<int> > comm;
		scalar_type one;
		scalar_type zero;
		scalar_type negative_one;
		scalar_type negative_two;
		scalar_type negative_four;

		// Belos Solver
		Teuchos::RCP<Teuchos::ParameterList> pl;
		Teuchos::RCP<Teuchos::ParameterList> pl2;
		Teuchos::RCP<Teuchos::ParameterList> pl3;
		Teuchos::RCP<problem_type> problem;
		Teuchos::RCP<problem_type> problem2;
		Teuchos::RCP<problem_type> problem3;
		Teuchos::RCP<belos_manager> belos_solver;
		Teuchos::RCP<belos_manager> belos_solver2;
		Teuchos::RCP<belos_manager> belos_solver3;
		Teuchos::RCP<Ifpack2::Factory> prec_factory;
		Teuchos::RCP<prec_type> belos_preconditioner;
		Teuchos::RCP<prec_type> belos_preconditioner2;
		Teuchos::RCP<prec_type> belos_preconditioner3;
		Teuchos::ArrayRCP<const scalar_type> data;
		Teuchos::ArrayRCP<const scalar_type> data2;
		Teuchos::ArrayRCP<const scalar_type> data3;

		RunParameters *rp;

		int my_rank;
		int num_procs;
		int bc_value;
		int my_y_size;
		int my_num_elements;
};


#endif

