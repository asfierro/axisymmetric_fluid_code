#ifndef FIELD_SOLVER_H
#define FIELD_SOLVER_H

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
#include <Xpetra_CrsMatrix.hpp>

#include <vector>

#include "mesh_data.h"
#include "constants.h"
#include "species.h"
#include "run_parameters.h"

class FieldSolver
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
		FieldSolver();
		FieldSolver(const Teuchos::RCP<const Teuchos::Comm<int> >&, scalar_type, int, RunParameters *);
		void InitializeSolver();
		void PerformInitialSolve();
		void ConfigureInitialSolve();
		void UpdateBC(scalar_type);
		void Solve();
		void PrintMatrixA();
		void PrintMatrixB();
		void PrintMatrixX();
		void ComputeElectricField();
		void UpdateCharge(std::vector<Species *> &);
		void UpdateCharge();
		double GetChargeDensity(int x, int y) { return charge_density[y][x]; }
		int GetYSize() { return my_y_size; }

		scalar_type GetPotentialAt(local_ordinal_type index) { return data[index]; }
		scalar_type GetPotentialAt(local_ordinal_type x, local_ordinal_type y) { return data[x + y*rp->x_size]; }
		scalar_type GetElectricFieldXAt(local_ordinal_type x, local_ordinal_type y);
		scalar_type GetElectricFieldYAt(local_ordinal_type x, local_ordinal_type y);

		//void SetPotentialAt(local_ordinal_type x, local_ordinal_type y, scalar_type val) { data[x + y*XSIZE] = val; }
		void SetPotentialAt(local_ordinal_type x, local_ordinal_type y, scalar_type val);
		void SetElectricFieldYAt(local_ordinal_type x, local_ordinal_type y, scalar_type val) { electric_field_y[y][x] = val; }
		void SetElectricFieldXAt(local_ordinal_type x, local_ordinal_type y, scalar_type val) { electric_field_x[y][x] = val; }

		void GetPreviousandNextElectricFields();
		int getNumberIterationsForSolve() { return belos_solver->getNumIters(); }	
		void ComputeChargeDensity(std::vector<Species *> &, int);


	private:
		Teuchos::RCP<crs_matrix_type> A;
		Teuchos::RCP<mv_type> X;
		Teuchos::RCP<mv_type> B;
		Teuchos::RCP<const map_type> global_map;

		Tpetra::global_size_t num_global_entries;
		Teuchos::RCP<const Teuchos::Comm<int> > comm;
		scalar_type one;
		scalar_type zero;
		scalar_type negative_one;
		scalar_type negative_two;
		scalar_type negative_four;

		// Amesos Solver
		Teuchos::RCP<Amesos2::Solver<crs_matrix_type, mv_type> > amesos_solver;

		// Belos Solver
		Teuchos::RCP<Teuchos::ParameterList> pl;
		Teuchos::RCP<problem_type> problem;
		Teuchos::RCP<belos_manager> belos_solver;
		Teuchos::RCP<Ifpack2::Factory> prec_factory;
		Teuchos::RCP<prec_type> belos_preconditioner;
		Teuchos::ArrayRCP<const scalar_type> data;

		RunParameters *rp;

		double **electric_field_x;
		double **electric_field_y;
		double **charge_density;

		double *electric_field_next;
		double *electric_field_y_plus_2;
		double *electric_field_previous;

		double *electric_field_x_next;
		double *electric_field_x_plus_2;
		double *electric_field_x_previous;

		double *buffer_row_next;
		double *buffer_row_previous;

		int solver_select;

		int my_rank;
		int num_procs;
		int bc_value;
		int my_y_size;
		int my_num_elements;
};


#endif

