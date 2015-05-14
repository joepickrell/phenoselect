/*
 * PhenoSelect.cpp
 *
 *  Created on: Dec 17, 2012
 *      Author: pickrell
 */


#include "PhenoBF.h"

using namespace std;


void printopts(){
        cout << "\nPhenoSelect v. 0.1\n";
        cout << "-c [file name] counts at control SNPs\n";
        cout << "-i [file name] counts at phenotype SNPs (risk,nonrisk alleles)\n";
        cout << "-o [file name] output file\n";
        cout << "-pops [string] [string] [string] the three population names\n";
        cout << "-w [integer] which population under selection (1, 2, 3)\n";
        cout << "-null no selection\n";
        cout << "-nit [integer] number of MCMC iterations\n";
        cout << "-sp [float] sd of normal prior on s\n";
        cout << "\n";
}

string controlfile;
string phenofile;
string outfile;
int whichpop = 1;
int nit = 500000;
double sp = 0.2;
bool null = false;

int main(int argc, char *argv[]){

    CCmdLine cmdline;
    PhyloPop_params p;
    if (cmdline.SplitLine(argc, argv) < 1){
        printopts();
        exit(1);
    }
    if (cmdline.HasSwitch("-i")) phenofile = cmdline.GetArgument("-i", 0).c_str();
    else{
        printopts();
        exit(1);
    }
    if (cmdline.HasSwitch("-pops")) {
    	if (cmdline.GetArgumentCount("-pops") < 3 or cmdline.GetArgumentCount("-pops") > 3){
    		cerr<< "ERROR: need three arguments to -pops\n";
    		exit(1);
    	}
    	p.pops[0] = cmdline.GetArgument("-pops", 0);
    	p.pops[1] = cmdline.GetArgument("-pops", 1);
    	p.pops[2] = cmdline.GetArgument("-pops", 2);
    }
      else{
    	  cerr << "ERROR: need -pops\n";
          printopts();
          exit(1);
      }
    if (cmdline.HasSwitch("-c")) controlfile = cmdline.GetArgument("-c", 0).c_str();
     else{
         printopts();
         exit(1);
     }

    if (cmdline.HasSwitch("-o")) outfile = cmdline.GetArgument("-o", 0).c_str();
      else{
          printopts();
          exit(1);
      }
    if (cmdline.HasSwitch("-w")) whichpop = atoi(cmdline.GetArgument("-w", 0).c_str());
    if (cmdline.HasSwitch("-nit")) nit = atoi(cmdline.GetArgument("-nit", 0).c_str());
    if (cmdline.HasSwitch("-null")) null  = true;
    if (cmdline.HasSwitch("-sp")) sp = atof(cmdline.GetArgument("-sp", 0).c_str());
    PhenoBF data(controlfile, phenofile, &p);
    data.s_prior = sp;
    data.whichpop = whichpop;
    if (null) data.models = false;
    data.run_MCMC(nit);
    data.print_stored(outfile);
    return 0;
}

