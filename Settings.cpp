#include <cstdlib> 
#include <iostream>
#include "Settings.h"



Settings::Settings(int argc, char *argv[]) {

	// set parameters to some default state
	alignmentFileName  = "";
	treeFileName       = "";
	outputFileName     = "";
	numOmegaCats       = 3;
	brlenLambda        = 10.0;
	numCycles          = 1000000;
	printFrequency     = 100;
	sampleFrequency    = 100;
	summarizeFrequency = 10000;
	
#	if 1
	/* set up fake command-line argument string */
	char *cmdString[30];
	cmdString[ 0] = (char *)"mmm";
	cmdString[ 1] = (char *)"-aln_path";
	cmdString[ 2] = (char *)"/Users/mlandis/Documents/code/posselreg/gap.p1.in";
	cmdString[ 3] = (char *)"-tree_path";
	cmdString[ 4] = (char *)"/Users/mlandis/Documents/code/posselreg/gap.p1.m3.tree";
	cmdString[ 5] = (char *)"-output_name";
	cmdString[ 6] = (char *)"/Users/mlandis/Documents/code/posselreg/gap.p1.out";
	cmdString[ 7] = (char *)"-brlen_lambda";
	cmdString[ 8] = (char *)"10.0";
	cmdString[ 9] = (char *)"-chain_length";
	cmdString[10] = (char *)"1000000";
	cmdString[11] = (char *)"-print_freq";
	cmdString[12] = (char *)"1";
	cmdString[13] = (char *)"-sample_freq";
	cmdString[14] = (char *)"50";
	cmdString[15] = (char *)"-sum_freq";
	cmdString[16] = (char *)"1000";
	
	argc = 17;
	argv = cmdString;
#	endif


	if (argc > 1)
		{
		if (argc % 2 == 0)
			{
			std::cout << "Usage:" << std::endl;
			std::cout << "   -aln_path <PATH>               : Path to file containing the alignment" << std::endl;
			std::cout << "   -tree_path <PATH>              : Path to file containing the newick tree file" << std::endl;
			std::cout << "   -output_name <NAME>            : Output file name" << std::endl;
			std::cout << "   -brlen_lambda <NUMBER>         : Exponential parameter for branch length prior" << std::endl;
			std::cout << "   -chain_length <NUMBER>         : Number of MCMC cycles" << std::endl;
			std::cout << "   -print_freq <NUMBER>           : Frequency with which output is printed to screen" << std::endl;
			std::cout << "   -sample_freq <NUMBER>          : Frequency with which chain is sampled to a file" << std::endl;
			std::cout << "   -sum_freq <NUMBER>             : Frequency with which chain is summarized" << std::endl;
			std::cout << std::endl;
			std::cout << "Example:" << std::endl;
			std::cout << "   mmm -aln_path <input file> -output_name <output file>" << std::endl;
			exit(0);
			}
			
		/* read the command-line arguments */
		std::string status = "none";
		for (int i=1; i<argc; i++)
			{
			std::string cmd = argv[i];
			//cout << cmd << endl;
			if (status == "none")
				{
				/* read the parameter specifier */
				if ( cmd == "-aln_path" )
					status = "aln_path";
				else if ( cmd == "-tree_path" )
					status = "tree_path";
				else if ( cmd == "-output_name" )
					status = "output_name";
				else if ( cmd == "-brlen_lambda" )
					status = "brlen_lambda";
				else if ( cmd == "-chain_length" )
					status = "chain_length";
				else if ( cmd == "-print_freq" )
					status = "print_freq";
				else if ( cmd == "-sample_freq" )
					status = "sample_freq";
				else if ( cmd == "-sum_freq" )
					status = "sum_freq";
				else
					{
					std::cerr << "Could not interpret option \"" << cmd << "\"." << std::endl;
					exit(1);
					}
				}
			else
				{
				/* read the parameter */
				if ( status == "aln_path" )
					{
					alignmentFileName = argv[i];
					}
				else if ( status == "tree_path" )
					{
					treeFileName = argv[i];
					}
				else if ( status == "output_name" )
					{
					outputFileName = argv[i];
					}
				else if ( status == "brlen_lambda" )
					{
					brlenLambda = atof(argv[i]);
					}
				else if ( status == "chain_length" )
					{
					numCycles = atoi(argv[i]);
					}
				else if ( status == "print_freq" )
					{
					printFrequency = atoi(argv[i]);
					}
				else if ( status == "sample_freq" )
					{
					sampleFrequency = atoi(argv[i]);
					}
				else if ( status == "sum_freq" )
					{
					summarizeFrequency = atoi(argv[i]);
					}
				else
					{
					std::cerr << "Unknown status reading command line information" << std::endl;
					exit(1);
					}
				status = "none";
				}
			}
		}
	else
		{
		std::cerr << "ERROR: You must provide command line arguments" << std::endl;
		}
}
