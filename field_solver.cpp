#include "field_solver.h"
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <string>
#include <BelosTpetraOperator.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include <BelosXpetraAdapterOperator.hpp>
#include <Xpetra_Matrix.hpp>
#include <BelosOperatorT.hpp>

FieldSolver::FieldSolver(const Teuchos::RCP<const Teuchos::Comm<int> > &n_comm, scalar_type voltage, int solver_type, RunParameters *new_rp)
{
	// Values that will go in the matrix
	one = static_cast<scalar_type> (1.0);
	negative_two = static_cast<scalar_type> (-2.0);
	negative_one = static_cast<scalar_type> (-1.0);
	zero = static_cast<scalar_type> (0.0);
	negative_four = static_cast<scalar_type> (-4.0);
	bc_value = voltage;

	comm = n_comm;
	my_rank = comm->getRank();
	num_procs = comm->getSize();
	solver_select = solver_type;
	rp = new_rp;
}

void FieldSolver::InitializeSolver()
{
	const int myRank = comm->getRank();
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::ArrayView;
	using Teuchos::tuple;

	my_y_size = rp->y_size / num_procs;
	if(my_rank == (num_procs-1))
	{
		int remainder_size = rp->y_size % num_procs;
		my_y_size = rp->y_size / num_procs + remainder_size;
	}

	my_num_elements = my_y_size * rp->x_size;
	//std::cout << "my rank = " << my_rank << "\t" << my_y_size << "\t" << std::endl;
	global_map = rcp(new map_type (rp->y_size*rp->x_size, my_y_size*rp->x_size, 0, comm));

	if(my_rank == 0)
	{
		std::cout << "total number of processors: " << num_procs << std::endl;
		std::cout << "YSIZE " << rp->y_size << std::endl;
		std::cout << "XSIZE " << rp->x_size << std::endl;
		}

	A = rcp(new crs_matrix_type (global_map,5));

	//my_y_size = global_map->getNodeNumElements();
	//std::cout << "my rank = " << my_rank << "\t" << global_map->getNodeNumElements() << "\t" << std::endl;

	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (my_num_elements); local_row++)
	{
		const global_ordinal_type global_row = global_map->getGlobalElement(local_row);
		//boundary condition insert
		if( (global_row % rp->x_size == 0 || ((global_row+1) % rp->x_size) == 0))
		{
			A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row), tuple<scalar_type>(one));
		}
		else
		{
			if( (global_row + rp->x_size) > (rp->x_size*rp->y_size) )
			{
				A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row),
																					tuple<scalar_type>(negative_one, one));
			}
			else if( (global_row - rp->x_size) < 0)
			{
				A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row, global_row+rp->x_size),
																					tuple<scalar_type>(negative_one, one));	
			}
			else
			{
				//A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-XSIZE, global_row-1, global_row, global_row+1, global_row+XSIZE),
																					//tuple<scalar_type> (one, one, negative_four, one, one));
				
				global_ordinal_type x_pos = global_row % rp->x_size;
				global_ordinal_type y_pos = global_row / rp->x_size;

				scalar_type t1 = static_cast<scalar_type> (1.0-SPACE_STEP/(2.0*((double)y_pos)*SPACE_STEP));
				scalar_type t2 = static_cast<scalar_type> (1.0+SPACE_STEP/(2.0*((double)y_pos)*SPACE_STEP));
				A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row-1, global_row, global_row+1, global_row+rp->x_size),
																					tuple<scalar_type> (t1, one, negative_four, one, t2));
			}
		}
	}

	A->fillComplete();
	X = rcp(new mv_type (A->getDomainMap(),1));
	B = rcp(new mv_type (A->getRangeMap(),1));

	/*
	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (B->getMap()->getNodeNumElements()); local_row++)
	{
		const global_ordinal_type global_row = global_map->getGlobalElement(local_row);

		if( (global_row % XSIZE) == 0)
		{
			B->replaceGlobalValue(global_row, 0, bc_value);
		}
	}*/

	// initialize electric field
	electric_field_x = new double*[my_y_size];
	electric_field_y = new double*[my_y_size];
	charge_density = new double*[my_y_size];
	for(int i = 0; i < (my_y_size); i++)
	{
		electric_field_x[i] = new double[rp->x_size];
		electric_field_y[i] = new double[rp->x_size];
		charge_density[i] = new double[rp->x_size];
	}

	for(int j = 0; j < my_y_size; j++)
	{
		for(int i = 0; i < rp->x_size; i++)
		{
			electric_field_x[j][i] = 0.0;
			electric_field_y[j][i] = 0.0;
			charge_density[j][i] = 0.0;
		}
	}

	electric_field_next = new double[rp->x_size];
	electric_field_previous = new double[rp->x_size];
	electric_field_x_next = new double[rp->x_size];
	electric_field_x_previous = new double[rp->x_size];
	electric_field_x_plus_2 = new double[rp->x_size];
	electric_field_y_plus_2 = new double[rp->x_size];
	buffer_row_next = new double[rp->x_size];
	buffer_row_previous = new double[rp->x_size];
	data = X->get1dView();
}

void FieldSolver::SetPotentialAt(local_ordinal_type x, local_ordinal_type y, scalar_type val)
{
	local_ordinal_type local_row = x + rp->x_size*y;
	const global_ordinal_type global_row = global_map->getGlobalElement(local_row);
	X->replaceGlobalValue(global_row, 0, val);
}

void FieldSolver::UpdateBC(scalar_type new_v)
{
	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (B->getMap()->getLocalNumElements()); local_row++)
	{
		//const global_ordinal_type global_row = global_map->getGlobalElement(local_row);

		if( (local_row % rp->x_size) == 0)
		{
			B->replaceLocalValue(local_row, 0, new_v);
		}
		else if ( ((local_row+1) % rp->x_size) == 0)
		{
			B->replaceLocalValue(local_row, 0, 0.0);
		}
	}
}

void FieldSolver::ConfigureInitialSolve()
{
	switch(solver_select)
	{
		case AMESOS:
			amesos_solver = Amesos2::create<crs_matrix_type, mv_type>("KLU2", A, X, B);
			break;
		case BELOS:
			Teuchos::ParameterList pl1;
			// ** ILUT
			//pl1.set("fact: iluk level-of-fill", 10);
			//pl1.set("fact: drop tolerance", 1e-6);
			//belos_preconditioner = Ifpack2::Factory::create<row_matrix_type>("ILUT", A);
			// ** ILUT
			//
			// ** RELAXATION
			//pl1.set("relaxation: type", "Gauss-Seidel");
			//pl1.set("relaxation: sweeps", 100);
			//pl1.set("relaxation: damping factor", 1.2);
			//belos_preconditioner = Ifpack2::Factory::create<row_matrix_type>("RELAXATION", A);
			// ** RELAXATION
			//
			// **ADDITIVE SCHWARZ
			/*
			Teuchos::ParameterList inner;
			inner.set("fact: iluk level-of-fill", 10);
			inner.set("fact: drop tolerance", 1e-6);
			pl1.set("inner preconditioner name", "ILUT");
			pl1.set("inner preconditioner parameters", inner);
			pl1.set("schwarz: overlap level", 1);
			pl1.set("schwarz: combine mode", "ADD");
			belos_preconditioner = Ifpack2::Factory::create<row_matrix_type>("SCHWARZ", A);
			// **ADDITIVE SCHWARZ
			//
			belos_preconditioner->setParameters(pl1);
			belos_preconditioner->initialize();
			belos_preconditioner->compute();
			pl = Teuchos::rcp(new Teuchos::ParameterList);
			pl->set("Maximum Iterations", 6000);
			pl->set("Maximum Restarts", 30);
			pl->set("Convergence Tolerance", 1e-8);
			pl->set("Implicit Residual Scaling", "Norm of RHS");
			pl->set("Explicit Residual Scaling", "Norm of RHS");
			pl->set("Num Blocks", 200);
			//pl->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
			//pl->set("Output Frequency", 500);
			problem = Teuchos::rcp(new problem_type(A, X, B));
			problem->setLeftPrec(belos_preconditioner);
			problem->setProblem();
			belos_solver = Teuchos::rcp( new belos_block_gm_res_manager(problem, pl) );*/


			// ** MUELU
			//Teuchos::Ptr<Teuchos::ParameterList> pl2;
			//Teuchos::updateParametersFromXmlFile("muelu/muelu.xml",pl2);
			//Teuchos::RCP<Xpetra::Matrix> At = Teuchos::rcp(new Xpetra::CrsMatrixWrap(A));
			//Teuchos::RCP<Teuchos::ParameterXMLFileReader> reader = Teuchos::rcp(Teuchos::ParameterXMLFileReader("muelue/muelu.xml"));
			//Teuchos::RCP<Belos::OperatorT<mv_type> > belosOp = Teuchos::rcp(new Belos::XpetraOp<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(A));
			//
			/*	
			Teuchos::RCP<mv_type> nullspace = Teuchos::rcp(new mv_type (A->getRowMap(),1));
			nullspace->putScalar(one);
			Teuchos::RCP<MueLu::HierarchyManager> mueLuFactory = Teuchos::rcp(new MueLu::ParameterListInterpreter("muelu/muelu.xml",comm));
			//std::string param_file = "muelu/muelu.xml";
			Teuchos::RCP<MueLu::Hierarchy> H;
			H = mueLuFactory->CreateHierarchy();
			H->GetLevel(0)->Set("A", A);
			H->GetLevel(0)->Set("Nullspace", nullspace);
			mueLuFactor->SetHierarchy(*H);
			typedef Belos::OperatorT<mv_type> OP;
			H->IsPreconditioner(true);
			Teuchos::RCP<OP> belos_preconditioner = Teuchos::rcp(new Belos::MueLuOp<scalar_type, local_ordinal_type, global_ordinal_type>(H));
			*/
			//
			
			std::cout << "preconditioner setup begin" << std::endl;
			Teuchos::RCP<MueLu::TpetraOperator<> > preconditioner = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<Tpetra::Operator<> >) A, "muelu/muelu.xml");
			pl = Teuchos::rcp(new Teuchos::ParameterList);
			pl->set("Maximum Iterations", 6000);
			pl->set("Maximum Restarts", 30);
			pl->set("Convergence Tolerance", 1e-6);
			pl->set("Implicit Residual Scaling", "Norm of RHS");
			pl->set("Explicit Residual Scaling", "Norm of RHS");
			//pl->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
			pl->set("Num Blocks", 200);
			problem = Teuchos::rcp(new problem_type(A, X, B));
			problem->setLeftPrec(preconditioner);
			problem->setProblem();
			belos_solver = Teuchos::rcp( new belos_block_gm_res_manager(problem, pl) );
			// ** MUELU
			std::cout << "preconditioner setup complete" << std::endl;
			break;
	}
}

void FieldSolver::PerformInitialSolve()
{
	switch(solver_select)
	{
		case AMESOS:
			amesos_solver->symbolicFactorization().numericFactorization().solve();
			break;
		case BELOS:
			belos_solver->solve();
			if(my_rank == 0)
				std::cout << "num iterations (INITIAL SOLVE) = " << belos_solver->getNumIters() << std::endl;
	}
	ComputeElectricField();
}

void FieldSolver::Solve()
{
	//std::cout << "solving" << std::endl;
	switch(solver_select)
	{
		case AMESOS:
			amesos_solver->symbolicFactorization().numericFactorization().solve();
			break;
		case BELOS:
			belos_solver->reset(Belos::ResetType::Problem);

			Belos::ReturnType ret = belos_solver->solve();

			if( ret != Belos::Converged )
				std::cout << "SOLVER NOT CONVERGED" << std::endl;
			break;
	}
	ComputeElectricField();
}


void FieldSolver::UpdateCharge()
{

	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (B->getMap()->getLocalNumElements()); local_row++)
	{
		if( (local_row % rp->x_size) != 0 && ((local_row+1) % rp->x_size) != 0)
		{
			local_ordinal_type x_pos = local_row % rp->x_size;
			local_ordinal_type y_pos = local_row / rp->x_size;

			double cd = charge_density[y_pos][x_pos] * SPACE_STEP * SPACE_STEP / EPSILON_0;
			B->replaceLocalValue(local_row, 0, -cd);
		}
	}
}

void FieldSolver::UpdateCharge(std::vector<Species *> &species_list)
{
	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (B->getMap()->getLocalNumElements()); local_row++)
	{
		if( (local_row % rp->x_size) != 0 && ((local_row+1) % rp->x_size) != 0)
		{
			local_ordinal_type x_pos = local_row % rp->x_size;
			local_ordinal_type y_pos = local_row / rp->x_size;

			double cd = 0.0;
			for(int i = 0; i < species_list.size(); i++)
			{
				//if(species_list[i]->GetSpeciesType() == ELECTRON || species_list[i]->GetSpeciesType() == ION || species_list[i]->GetSpeciesType() == NEGATIVE_ION)
				if(species_list[i]->isCharged())
					cd += species_list[i]->GetDensity(x_pos, y_pos+1) * species_list[i]->GetSpeciesCharge();
			}
			cd = cd * SPACE_STEP * SPACE_STEP / EPSILON_0;
			B->replaceLocalValue(local_row, 0, -cd);
		}
	}
}

void FieldSolver::ComputeElectricField()
{	
 	// send data to next processor - even processors to odd
	if((my_rank % 2) == 0)
	{
		MPI_Send(&data[rp->x_size*(my_y_size-1)], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&buffer_row_previous[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// send data to next processor - odd processors to even
	if((my_rank % 2) == 1 && my_rank != (num_procs-1))
	{
		MPI_Send(&data[rp->x_size*(my_y_size-1)], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 && my_rank != (num_procs-1))	
	{
		MPI_Recv(&buffer_row_previous[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// send data to previous processor - even processors to odd
	if((my_rank % 2 ) == 0 && my_rank != 0)
	{
		MPI_Send(&data[0], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 && my_rank != (num_procs-1))
	{
		MPI_Recv(&buffer_row_next[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// send data to previous processor - odd processors to even
	if((my_rank % 2) == 1)
	{
		MPI_Send(&data[0], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&buffer_row_next[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	/*
	// general field equation
	for(int j = 0; j < my_y_size-1; j++)
	{
		for(int i = 0; i < (XSIZE-1); i++)
		{
			electric_field_x[j][i] = -1.0 * (data[i+1+j*XSIZE] - data[i+j*XSIZE])/SPACE_STEP;
			electric_field_y[j][i] = -1.0 * (data[i+(j+1)*XSIZE] - data[i+j*XSIZE])/SPACE_STEP;
		}
		electric_field_x[j][XSIZE-1] = -1.0 * (data[(XSIZE-1)+j*XSIZE] - data[(XSIZE-2)+j*XSIZE])/SPACE_STEP;
	}*/
	/*	
	// correct last row y field
	if(my_rank == (num_procs - 1))
	{
		for(int i = 0; i < XSIZE; i++)
			buffer_row_next[i] = data[i+(my_y_size-1)*XSIZE];		// this is empty on the last processor, so fill it
	}
	for(int i = 0; i < (XSIZE); i++)
	{
		electric_field_y[my_y_size-1][i] = -1.0 * (buffer_row_next[i] - data[i+(my_y_size-1)*XSIZE])/SPACE_STEP;
	}

	//  correct last row x field
	for(int i = 0; i < (XSIZE-1); i++)
	{
		electric_field_x[my_y_size-1][i] = -1.0 * (data[i+1+(my_y_size-1)*XSIZE] - data[i+(my_y_size-1)*XSIZE])/SPACE_STEP;
	}
	electric_field_x[my_y_size-1][XSIZE-1] = -1.0 * (data[(XSIZE-1)+(my_y_size-1)*XSIZE] - data[(XSIZE-2)+(my_y_size-1)*XSIZE])/SPACE_STEP;

	//if(my_rank == 0)
	//{
		//for(int i = 0; i < (XSIZE); i++)
			//electric_field_y[0][i] = 0.0;
	//}*/
	// central difference scheme x_field
	for(int j = 0; j < my_y_size; j++)
	{
		for(int i = 1; i < (rp->x_size-1); i++)
		{
			electric_field_x[j][i] = -1.0 * (data[i+1+j*rp->x_size] - data[i-1+j*rp->x_size])/(2.0*SPACE_STEP);
		}
	}
	// central difference scheme y_field
	for(int j = 1; j < my_y_size-1; j++)
	{
		for(int i = 0; i < (rp->x_size); i++)
		{
			electric_field_y[j][i] = -1.0 * (data[i+(j+1)*rp->x_size] - data[i+(j-1)*rp->x_size])/(2.0*SPACE_STEP);
		}
	}
	if(my_rank != (num_procs-1))
	{
		for(int i = 0; i < rp->x_size; i++)
			electric_field_y[my_y_size-1][i] = -1.0 * (buffer_row_next[i] - data[i+(my_y_size-2)*rp->x_size])/(2.0*SPACE_STEP);
	}
	else
	{
		for(int i = 0; i < rp->x_size; i++)
			electric_field_y[my_y_size-1][i] = 0.0;
	}
	// first row y field on all procs
	if(my_rank != 0)
	{
		for(int i = 0; i < rp->x_size; i++)
			electric_field_y[0][i] = -1.0 * (data[i+rp->x_size] - buffer_row_previous[i])/(2.0*SPACE_STEP);
	}
	else
	{
		for(int i = 0; i < rp->x_size; i++)
			electric_field_y[0][i] = 0.0;
	}

	// fix first/last row x field on all procs
	// fix first/last row y field on all procs
	for(int j = 0; j < my_y_size; j++)
	{
		electric_field_x[j][0] = electric_field_x[j][1];	
		electric_field_x[j][rp->x_size-1] = electric_field_x[j][rp->x_size-2];	

		electric_field_y[j][0] = 0.0;
		electric_field_y[j][rp->x_size-1] = 0.0;
	}
	GetPreviousandNextElectricFields();
}

FieldSolver::scalar_type FieldSolver::GetElectricFieldYAt(local_ordinal_type x, local_ordinal_type y)
{
	if(my_rank == (num_procs-1))
	{
		if(y == -1)
			return electric_field_previous[x];
		else if(y == my_y_size || y == (my_y_size+1))
			return electric_field_y[my_y_size-1][x];
		else
			return electric_field_y[y][x];
	}
	else if(my_rank == 0)
	{
		if(y == -1)
			return electric_field_y[1][x];
		else if(y == my_y_size)
			return electric_field_next[x];
		else if(y == (my_y_size+1))
			return electric_field_y_plus_2[x];
		else
			return electric_field_y[y][x];
	}
	else
	{
		if(y == -1)
			return electric_field_previous[x];
		else if(y == my_y_size)
			return electric_field_next[x];
		else if(y == (my_y_size+1))
			return electric_field_y_plus_2[x];
		else
			return electric_field_y[y][x];

	}
}

FieldSolver::scalar_type FieldSolver::GetElectricFieldXAt(local_ordinal_type x, local_ordinal_type y)
{
	if(my_rank == (num_procs-1))
	{
		if(y == -1)
			return electric_field_x_previous[x];
		else if(y == my_y_size || y == (my_y_size+1))
			return electric_field_x[my_y_size-1][x];
		else
			return electric_field_x[y][x];
	}
	else if(my_rank == 0)
	{
		if(y == -1)
			return electric_field_x[1][x];
		else if(y == my_y_size)
			return electric_field_x_next[x];
		else if(y == (my_y_size+1))
			return electric_field_x_plus_2[x];
		else
			return electric_field_x[y][x];
	}
	else
	{
		if(y == -1)
			return electric_field_x_previous[x];
		else if(y == my_y_size)
			return electric_field_x_next[x];
		else if(y == (my_y_size+1))
			return electric_field_x_plus_2[x];
		else
			return electric_field_x[y][x];
	}
}

void FieldSolver::GetPreviousandNextElectricFields()
{
 	// send data to next processor - even processors to odd
	if((my_rank % 2) == 0)
	{
		MPI_Send(electric_field_y[my_y_size-1], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_x[my_y_size-1], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&electric_field_previous[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_x_previous[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// send data to next processor - odd processors to even
	if((my_rank % 2) == 1 && my_rank != (num_procs-1))
	{
		MPI_Send(electric_field_y[my_y_size-1], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_x[my_y_size-1], rp->x_size, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 && my_rank != (num_procs-1))	
	{
		MPI_Recv(&electric_field_previous[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_x_previous[0], rp->x_size, MPI_DOUBLE, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// send data to previous processor - even processors to odd
	if((my_rank % 2 ) == 0 && my_rank != 0)
	{
		MPI_Send(electric_field_y[0], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_x[0], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_y[1], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_x[1], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else if(my_rank != 0 && my_rank != (num_procs-1))
	{
		MPI_Recv(&electric_field_next[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_x_next[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_y_plus_2[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_x_plus_2[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}

	// send data to previous processor - odd processors to even
	if((my_rank % 2) == 1)
	{
		MPI_Send(electric_field_y[0], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_x[0], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_y[1], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(electric_field_x[1], rp->x_size, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&electric_field_next[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_x_next[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_y_plus_2[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&electric_field_x_plus_2[0], rp->x_size, MPI_DOUBLE, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
}

void FieldSolver::ComputeChargeDensity(std::vector<Species *> &species_list, int beg_y)
{
	for(int i = 0; i < rp->x_size; i++)
	{
		for(int j = 0; j < my_y_size; j++)
			charge_density[j][i] = 0.0;
	}
	int y_begin = 1;
	int y_end = my_y_size + 1;

	if(my_rank == 0)
	{
		y_begin = 2;
	}
	else if(my_rank == (num_procs-1))
	{
		y_end = my_y_size;
	}

	for(int j = y_begin; j < y_end; j++)
	{
		double w1 = 0.25;
		double w2 = 0.25;

		for(int i = 1; i < (rp->x_size-1); i++)
		{
			for(int q = 0; q < species_list.size(); q++)
			{
				if(species_list[q]->isCharged())
				{
					charge_density[j-1][i] += w1 * (species_list[q]->GetDensity(i-1,j-1)) * species_list[q]->GetSpeciesCharge();
					charge_density[j-1][i] += w1 * (species_list[q]->GetDensity(i,j-1)) * species_list[q]->GetSpeciesCharge();
					charge_density[j-1][i] += w2 * (species_list[q]->GetDensity(i-1,j)) * species_list[q]->GetSpeciesCharge();
					charge_density[j-1][i] += w2 * (species_list[q]->GetDensity(i,j)) * species_list[q]->GetSpeciesCharge();
				}
			}
		}
	}

	if(my_rank == (num_procs-1))
	{
		double w2 = 0.5;
		for(int i = 1; i < (rp->x_size-1); i++)
		{
			for(int q = 0; q < species_list.size(); q++)
			{
				if(species_list[q]->isCharged())
				{
					charge_density[my_y_size-1][i] += 1.0 * w2 * (species_list[q]->GetDensity(i-1,my_y_size-1)) * species_list[q]->GetSpeciesCharge();
					charge_density[my_y_size-1][i] += 1.0 * w2 * (species_list[q]->GetDensity(i,my_y_size-1)) * species_list[q]->GetSpeciesCharge();
					charge_density[my_y_size-1][i] = 0.0;
				}
			}
		}
	}
	// axis of symmetry
	if(my_rank == 0)
	{
		double w2 = 0.5;
		for(int i = 1; i < (rp->x_size-1); i++)
		{
			for(int q = 0; q < species_list.size(); q++)
			{
				if(species_list[q]->isCharged())
				{
					charge_density[0][i] += 1.0 * w2 * (species_list[q]->GetDensity(i-1,1)) * species_list[q]->GetSpeciesCharge();
					charge_density[0][i] += 1.0 * w2 * (species_list[q]->GetDensity(i,1)) * species_list[q]->GetSpeciesCharge();
				}
			}
		}
	}

	// correct first x row and last x row
	for(int j = 1; j < y_end; j++)
	{
		for(int q = 0; q < species_list.size(); q++)
		{
			if(species_list[q]->isCharged())
			{
				charge_density[j-1][0] = 0.0;
				charge_density[j-1][rp->x_size-1] = 0.0;
			}
		}
	}
}


void FieldSolver::PrintMatrixA()
{
	Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
	A->describe(*fos, Teuchos::VERB_EXTREME);
	*fos << std::endl;
}

void FieldSolver::PrintMatrixX()
{
	Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
	X->describe(*fos, Teuchos::VERB_EXTREME);
	*fos << std::endl;
}

void FieldSolver::PrintMatrixB()
{
	Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
	B->describe(*fos, Teuchos::VERB_EXTREME);
	*fos << std::endl;
}


