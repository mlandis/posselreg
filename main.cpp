#include <iostream>
#include <string>
#include "Alignment.h"
#include "iomanager.h"
#include "Settings.h"
#include "MbRandom.h"
#include "Mcmc.h"
#include "Model.h"
#include "Util.h"



int main (int argc, char *argv[]) {

	// get user settings
	Settings userSettings(argc, argv);
	
	// instantiate the file managers
	IoManager inputFileMngr( userSettings.getAlignmentFileName() );
	IoManager treeFileMngr( userSettings.getTreeFileName() );
	IoManager outputFileMngr( userSettings.getOutputFileName() );
	
	// read in the alignment
	Alignment alignment( &inputFileMngr );
	alignment.translate();

	// read in the tree
	std::string treeStr = "";
	if (userSettings.getTreeFileName() != "")
		treeStr = Util::getLineFromFile(treeFileMngr.getFilePathName(), 1);
	
	// instantiate the random variable object
	MbRandom myRandom( 1 );
	
	// set up the phylogenetic model
	Model myModel( &alignment, &myRandom, &userSettings, treeStr );
	
	// run the MCMC analysis
	Mcmc myMcmc( &alignment, &userSettings, &myRandom, &myModel, &outputFileMngr );
	
    return 0;
}
