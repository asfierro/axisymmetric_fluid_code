#include "photon_solver.h"
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <MueLu_Hierarchy.hpp>
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
#include "reaction.h"

PhotonSolver::PhotonSolver(const Teuchos::RCP<const Teuchos::Comm<int> > &n_comm, RunParameters *new_rp)
{
	// Values that will go in the matrix
	one = static_cast<scalar_type> (1.0);
	negative_two = static_cast<scalar_type> (-2.0);
	negative_one = static_cast<scalar_type> (-1.0);
	zero = static_cast<scalar_type> (0.0);
	negative_four = static_cast<scalar_type> (-4.0);

	comm = n_comm;
	my_rank = comm->getRank();
	num_procs = comm->getSize();
	rp = new_rp;
}

void PhotonSolver::InitializeSolver()
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
	A2 = rcp(new crs_matrix_type (global_map,5));
	A3 = rcp(new crs_matrix_type (global_map,5));

	//my_y_size = global_map->getNodeNumElements();
	//std::cout << "my rank = " << my_rank << "\t" << global_map->getNodeNumElements() << "\t" << std::endl;

	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (my_num_elements); local_row++)
	{
		const global_ordinal_type global_row = global_map->getGlobalElement(local_row);
		//boundary condition insert
		if( (global_row % rp->x_size) == 0 || ((global_row+1) % rp->x_size) == 0)
		{
			A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row), tuple<scalar_type>(one));
			A2->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row), tuple<scalar_type>(one));
			A3->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row), tuple<scalar_type>(one));
		}
		else
		{
			if( (global_row + rp->x_size) > (rp->x_size*rp->x_size) )
			{
			/*
				A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row),
																					tuple<scalar_type>(negative_one, one));
				A2->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row),
																					tuple<scalar_type>(negative_one, one));
				A3->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row),
																					tuple<scalar_type>(negative_one, one));*/
				A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row), tuple<scalar_type>(one));
				A2->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row), tuple<scalar_type>(one));
				A3->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row), tuple<scalar_type>(one));
			}
			else if( (global_row - rp->x_size) < 0)
			{
				A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row, global_row+rp->x_size),
																					tuple<scalar_type>(negative_one, one));	
				A2->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row, global_row+rp->x_size),
																					tuple<scalar_type>(negative_one, one));	
				A3->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row, global_row+rp->x_size),
																					tuple<scalar_type>(negative_one, one));	
			}
			else
			{
				//A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-XSIZE, global_row-1, global_row, global_row+1, global_row+XSIZE),
																					//tuple<scalar_type> (one, one, negative_four, one, one));
				
				global_ordinal_type x_pos = global_row % rp->x_size;
				global_ordinal_type y_pos = global_row / rp->x_size;

				scalar_type t1 = static_cast<scalar_type> (1.0-SPACE_STEP/(2.0*((float)y_pos)*SPACE_STEP));
				scalar_type t2 = static_cast<scalar_type> (1.0+SPACE_STEP/(2.0*((float)y_pos)*SPACE_STEP));
				scalar_type t3 = static_cast<scalar_type> (-SPACE_STEP*SPACE_STEP*688070.25 - 4.0);
				A->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row-1, global_row, global_row+1, global_row+rp->x_size),
																					tuple<scalar_type> (t1, one, t3, one, t2));

				scalar_type t4 = static_cast<scalar_type> (-SPACE_STEP*SPACE_STEP*4796100.0 - 4.0);
				A2->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row-1, global_row, global_row+1, global_row+rp->x_size),
																					tuple<scalar_type> (t1, one, t4, one, t2));

				scalar_type t5 = static_cast<scalar_type> (-SPACE_STEP*SPACE_STEP*178222500.0 - 4.0);
				A3->insertGlobalValues(global_row, tuple<global_ordinal_type> (global_row-rp->x_size, global_row-1, global_row, global_row+1, global_row+rp->x_size),
																					tuple<scalar_type> (t1, one, t5, one, t2));
			}
		}
	}

	A->fillComplete();
	A2->fillComplete();
	A3->fillComplete();
	X = rcp(new mv_type (A->getDomainMap(),1));
	B = rcp(new mv_type (A->getRangeMap(),1));
	X2 = rcp(new mv_type (A2->getDomainMap(),1));
	B2 = rcp(new mv_type (A2->getRangeMap(),1));
	X3 = rcp(new mv_type (A3->getDomainMap(),1));
	B3 = rcp(new mv_type (A3->getRangeMap(),1));

	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (B->getMap()->getLocalNumElements()); local_row++)
	{
		const global_ordinal_type global_row = global_map->getGlobalElement(local_row);

		if( (global_row % rp->x_size) == 0)
		{
			B->replaceGlobalValue(global_row, 0, 0.0);
			B2->replaceGlobalValue(global_row, 0, 0.0);
			B3->replaceGlobalValue(global_row, 0, 0.0);
		}
	}

	data = X->get1dView();
	data2 = X2->get1dView();
	data3 = X3->get1dView();
}

void PhotonSolver::ConfigureSolver()
{
			Teuchos::ParameterList pl1;
			// ** ILUT
			//pl1.set("fact: iluk level-of-fill", 10);
			//pl1.set("fact: drop tolerance", 1e-2);
			//belos_preconditioner = Ifpack2::Factory::create<row_matrix_type>("ILUT", A);
			// ** ILUT
			//
			// ** RELAXATION
			//pl1.set("relaxation: type", "Gauss-Seidel");
			//pl1.set("relaxation: sweeps", 200);
			//pl1.set("relaxation: damping factor", 1.2);
			//belos_preconditioner = Ifpack2::Factory::create<row_matrix_type>("RELAXATION", A);
			// ** RELAXATION
			//
			// **ADDITIVE SCHWARZ
			/*
			Teuchos::ParameterList inner;
			inner.set("fact: iluk level-of-fill", 10);
			pl1.set("inner preconditioner name", "ILUT");
			pl1.set("inner preconditioner parameters", inner);
			pl1.set("schwarz: overlap level", 1);
			pl1.set("schwarz: combine mode", "ADD");
			belos_preconditioner = Ifpack2::Factory::create<row_matrix_type>("SCHWARZ", A);
			belos_preconditioner2 = Ifpack2::Factory::create<row_matrix_type>("SCHWARZ", A2);
			// **ADDITIVE SCHWARZ
			belos_preconditioner->setParameters(pl1);
			belos_preconditioner->initialize();
			belos_preconditioner->compute();

			// **ADDITIVE SCHWARZ
			belos_preconditioner2->setParameters(pl1);
			belos_preconditioner2->initialize();
			belos_preconditioner2->compute();*/


			Teuchos::RCP<MueLu::TpetraOperator<> > preconditioner1 = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<Tpetra::Operator<> >) A, "muelu/muelu.xml");
			pl = Teuchos::rcp(new Teuchos::ParameterList);
			pl->set("Maximum Iterations", 6000);
			pl->set("Maximum Restarts", 20);
			pl->set("Convergence Tolerance", 1e-10);
			pl->set("Implicit Residual Scaling", "Norm of RHS");
			pl->set("Explicit Residual Scaling", "Norm of RHS");
			pl->set("Num Blocks", 300);
			//pl->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
			//pl->set("Output Frequency", 500);
			problem = Teuchos::rcp(new problem_type(A, X, B));
			problem->setLeftPrec(preconditioner1);
			problem->setProblem();
			belos_solver = Teuchos::rcp( new belos_block_gm_res_manager(problem, pl) );

			Teuchos::RCP<MueLu::TpetraOperator<> > preconditioner2 = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<Tpetra::Operator<> >) A2, "muelu/muelu.xml");
			pl2 = Teuchos::rcp(new Teuchos::ParameterList);
			pl2->set("Maximum Iterations", 6000);
			pl2->set("Maximum Restarts", 20);
			pl2->set("Convergence Tolerance", 1e-10);
			pl2->set("Implicit Residual Scaling", "Norm of RHS");
			pl2->set("Explicit Residual Scaling", "Norm of RHS");
			pl2->set("Num Blocks", 300);
			problem2 = Teuchos::rcp(new problem_type(A2, X2, B2));
			problem2->setLeftPrec(preconditioner2);
			problem2->setProblem();
			belos_solver2 = Teuchos::rcp( new belos_block_gm_res_manager(problem2, pl2) );

			pl3 = Teuchos::rcp(new Teuchos::ParameterList);
			pl3->set("Maximum Iterations", 6000);
			pl3->set("Maximum Restarts", 20);
			pl3->set("Convergence Tolerance", 1e-10);
			pl3->set("Implicit Residual Scaling", "Norm of RHS");
			pl3->set("Explicit Residual Scaling", "Norm of RHS");
			pl3->set("Num Blocks", 300);
			Teuchos::RCP<MueLu::TpetraOperator<> > preconditioner3 = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<Tpetra::Operator<> >) A3, "muelu/muelu.xml");
			problem3 = Teuchos::rcp(new problem_type(A3, X3, B3));
			problem3->setLeftPrec(preconditioner3);
			problem3->setProblem();
			belos_solver3 = Teuchos::rcp( new belos_block_gm_res_manager(problem3, pl3) );
}

void PhotonSolver::Solve()
{
	//if(my_rank == 0)
		//std::cout << "****************** FIRST MATRIX SOLVE ***************** " << std::endl;
	belos_solver->reset(Belos::ResetType::Problem);
	belos_solver->solve();

	//if(my_rank == 0)
		//std::cout << "****************** SECOND MATRIX SOLVE ***************** " << std::endl;
	belos_solver2->reset(Belos::ResetType::Problem);
	belos_solver2->solve();

	belos_solver3->reset(Belos::ResetType::Problem);
	belos_solver3->solve();
}

void PhotonSolver::UpdateIonizationRate(std::vector<Species *> &species_list, FieldSolver *p)
{
	/*
	TransportData *tp;
	int electron_index = -1;
	for(int q = 0; q < species_list.size(); q++)
	{
		if(species_list[q]->GetSpeciesType() == ELECTRON)
		{
			electron_index = q;
			break;
		}
	}

	std::vector<Reaction *> g = species_list[electron_index]->GetGainReactions();
	tp = species_list[electron_index]->GetTransportData();
	Reaction *ionization_reaction = NULL;
	for(int q = 0; q < g.size(); q++)
	{
		if(g[q]->GetReactionName() == "n2_ionization")
		{
			ionization_reaction = g[q];
			break;
		}
	}

	if(ionization_reaction == NULL)
	{
		std::cout << "cannot find ionization reaction for photon solver model" << std::endl;
		return;
	}*/
	
	for(local_ordinal_type local_row = 0; local_row < static_cast<local_ordinal_type> (B->getMap()->getLocalNumElements()); local_row++)
	{
		const global_ordinal_type global_row = global_map->getGlobalElement(local_row);

		if( (local_row % rp->x_size) == 0)
		{
			B->replaceLocalValue(local_row, 0, 0.0);
			B2->replaceLocalValue(local_row, 0, 0.0);
			B3->replaceLocalValue(local_row, 0, 0.0);
		}
		else if ( ((local_row+1) % rp->x_size) == 0)
		{
			B->replaceLocalValue(local_row, 0, 0.0);
			B2->replaceLocalValue(local_row, 0, 0.0);
			B3->replaceLocalValue(local_row, 0, 0.0);
		}
		else if ( (global_row + rp->x_size) > (rp->x_size*rp->x_size) )
		{
			B->replaceLocalValue(local_row, 0, 0.0);
			B2->replaceLocalValue(local_row, 0, 0.0);
			B3->replaceLocalValue(local_row, 0, 0.0);
		}
		else
		{
			local_ordinal_type x_pos = local_row % rp->x_size;
			local_ordinal_type y_pos = local_row / rp->x_size;
			/*
			double field_x_middle = (p->GetElectricFieldXAt(x_pos,y_pos) + p->GetElectricFieldXAt(x_pos+1,y_pos) +
																p->GetElectricFieldXAt(x_pos,y_pos+1) + p->GetElectricFieldXAt(x_pos+1,y_pos+1)) / 4.0;

			double field_y_middle = (p->GetElectricFieldYAt(x_pos,y_pos) + p->GetElectricFieldYAt(x_pos+1,y_pos) +
																p->GetElectricFieldYAt(x_pos,y_pos+1) + p->GetElectricFieldYAt(x_pos+1,y_pos+1)) / 4.0;

			double mag_e = sqrt(pow(field_x_middle,2) + pow(field_y_middle,2));
			double E_n = mag_e / N_B / 1e-21;*/

			double z_r = x_pos * SPACE_STEP;
			double r_r = (double)((int)global_row / rp->x_size) * SPACE_STEP;
			double z_0 = rp->x_size / 2.0 * SPACE_STEP;
			double SD = 1e-8;

			double intensity = 3.5e28 * exp(-pow(z_r - z_0,2) / SD - pow(r_r,2) / SD);

			// ---- Using A1 - these already have 0.06 * pq / (pq+p) baked into the numbers out front
			//double val = -126.50194 * ionization_reaction->GetReactionRate(E_n) * species_list[electron_index]->GetDensity(x_pos, y_pos+1) 
						//* N_B * SPACE_STEP * SPACE_STEP ;
			double val = -44685.0 * intensity * SPACE_STEP * SPACE_STEP;
			if(val > 0.0)
				val = 0.0;
			B->replaceLocalValue(local_row, 0, val);

			// ---- Using A2
			//val = -3248.53938 * ionization_reaction->GetReactionRate(E_n) * species_list[electron_index]->GetDensity(x_pos, y_pos+1)
						//* N_B * SPACE_STEP * SPACE_STEP;
			val = -1147500.0 * intensity * SPACE_STEP * SPACE_STEP;
			if(val > 0.0)
				val = 0.0;
			B2->replaceLocalValue(local_row, 0, val);

			// ---- Using A3
			//val = -311222.812192 * ionization_reaction->GetReactionRate(E_n) * species_list[electron_index]->GetDensity(x_pos, y_pos+1)
						//* N_B * SPACE_STEP * SPACE_STEP;
			val = -109935000.0 * intensity * SPACE_STEP * SPACE_STEP;
			if(val > 0.0)
				val = 0.0;
			B3->replaceLocalValue(local_row, 0, val);


			// for use in alpha model
			//val = -92016.0 * ionization_reaction->GetReactionRate(E_n) * species_list[electron_index]->GetDensity(x_pos, y_pos+1) 
						//* (tp->GetMobility(E_n)/N_B) * mag_e * SPACE_STEP * SPACE_STEP;
		}
	}
}

void PhotonSolver::PrintMatrixA()
{
	Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
	A->describe(*fos, Teuchos::VERB_EXTREME);
	*fos << std::endl;
}

void PhotonSolver::PrintMatrixX()
{
	Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
	X->describe(*fos, Teuchos::VERB_EXTREME);
	*fos << std::endl;
}

void PhotonSolver::PrintMatrixB()
{
	Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
	B->describe(*fos, Teuchos::VERB_EXTREME);
	*fos << std::endl;
}


