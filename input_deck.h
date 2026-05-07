#ifndef INPUT_DECK_H
#define INPUT_DECK_H

#include <string>
#include <fstream>
#include <vector>
#include "species.h"
#include "reaction.h"
#include "mesh_data.h"
#include "run_parameters.h"

class InputDeck
{
	private:
		std::string filename;
		std::string GetToken(std::string, std::string);
		std::vector<std::string> input_lines;
		std::vector<int> input_lines_read;

		bool OpenFile(std::ifstream &);
		void CloseFile(std::ifstream &);
		void CheckNonReadLines(mesh_data *);

		void ReadSpeciesInteractions(mesh_data *, std::vector<Species *> &);
		void ReadSpeciesInput(mesh_data *, std::vector<Species *> &, RunParameters *);
		int ReadXSize(mesh_data *);
		int ReadYSize(mesh_data *);
		double ReadVoltage(mesh_data *);
		void ReadVoltageMode(mesh_data *, RunParameters *);
		int ReadOutputStride(mesh_data *);
		int ReadStatusStride(mesh_data *);
		int ReadRestart(mesh_data *);
		double ReadFinalTime(mesh_data *);
		void ReadSpeciesInitialDensity(mesh_data *, RunParameters *, std::vector<Species *> &);
		std::string ReadOutputDirectory(mesh_data *);
		std::string ReadPhotoionizationMethod(mesh_data *);

	public:
		InputDeck() { }
		InputDeck(std::string f) { filename = f; }

		void ReadInputDeck(RunParameters *, mesh_data *, std::vector<Species *> &);
		void Finish(mesh_data *, std::vector<Species *> &, RunParameters *);
};

#endif
