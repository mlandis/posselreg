#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>

class Settings {

	public:
                            Settings(int argc, char *argv[]);
			  std::string   getAlignmentFileName(void) { return alignmentFileName; }
			  std::string   getTreeFileName(void) { return treeFileName; } 
			  std::string   getOutputFileName(void) { return outputFileName; }
			          int   getNumOmegaCats(void) { return numOmegaCats; }
				   double   getBrlenLambda(void) { return brlenLambda; }
				      int   getChainLength(void) { return numCycles; }
					  int   getPrintFrequency(void) { return printFrequency; }
					  int   getSampleFrequency(void) { return sampleFrequency; }
					  int   getSummarizeFrequency(void) { return summarizeFrequency; }
			         void   setAlignmentFileName(std::string s) { alignmentFileName = s; }
					 void   setTreeFileName(std::string s) { treeFileName = s; }
					 void   setOutputFileName(std::string s) { outputFileName = s; }
					 void   setBrlenLambda(double x) { brlenLambda = x; }

	private:
			  std::string   alignmentFileName;
			  std::string   treeFileName;
			  std::string   outputFileName;
			          int   numOmegaCats;
				   double   brlenLambda;
				      int   numCycles;
					  int   printFrequency;
					  int   sampleFrequency;
					  int   summarizeFrequency;
};


#endif