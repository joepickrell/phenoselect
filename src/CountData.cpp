/*
 * CountData.cpp
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */
#include "CountData.h"



CountData::CountData(string infile, PhyloPop_params* p){
	params = p;
	read_counts(infile);
	cout << "npop:"<< npop<< " nsnp:"<<nsnp<< "\n";
}



double CountData::get_freq(int whichpop, int whichsnp){
	double toreturn;
	pair<int, int> cc = get_counts(whichpop, whichsnp);
	toreturn = (double) cc.first/  (  (double) cc.second + (double) cc.first );
 	//toreturn = gsl_matrix_get(alfreqs, whichsnp, whichpop);
	return toreturn;
}


pair<int, int> CountData::get_counts(int whichpop, int whichsnp){
	pair<int, int> toreturn;
	toreturn = allele_counts.at(whichsnp).at(whichpop);
	return toreturn;
}


void CountData::read_counts(string infile){
	//cout << "here?\n"; cout.flush();
    allele_counts.clear();
    head2id.clear();
    id2head.clear();
    npop = 0;
    nsnp = 0;
    string ext = infile.substr(infile.size()-3, 3);
    if (ext != ".gz"){
    	std::cerr << infile << " is not gzipped (only .gz files accepted)\n";
    	exit(1);
    }
	igzstream in(infile.c_str()); //only gzipped files
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    string st, buf;

    intStat = stat(infile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << infile << "\n";
            exit(1);
    }

    /*
     * header contains population names
     */
    getline(in, st);
    stringstream ss(st);
    line.clear();
    while (ss>> buf){
    	line.push_back(buf);
     }

    /*
     * make map from header, number populations according to order
     */
    int start = 0;
    for(int i = start; i < line.size(); i++) {
    	//cout << line.size() << "\n"; cout.flush();
    	head2id.insert(make_pair(line[i], i-start));
    	id2head.insert(make_pair(i-start, line[i]));
    	npop ++;
    }
    int headsize = line.size();

    /*
     * read counts, store in allele_counts
     */
    while(getline(in, st)){
            buf.clear();
            stringstream ss(st);
            line.clear();
            while (ss>> buf){
                    line.push_back(buf);
            }
            vector<pair<int, int> > topush;
            if (line.size() != headsize){
            	cerr << "ERROR: Line "<< nsnp <<" has "<< line.size() << " entries. Header has "<< headsize <<"\n";
            	exit(1);
            }
            for(vector<string>::iterator it = params->pops.begin(); it != params->pops.end(); it++){
            	if (head2id.find(*it) == head2id.end()){
            		cerr << "ERROR: can't find population "<< *it << "\n";
            		exit(1);
            	}
            	int i = head2id.find(*it)->second;
            	//cout << *it << " "<< i <<"\n"; cout.flush();
                typedef boost::tokenizer<boost::char_separator<char> >
                tokenizer;
                boost::char_separator<char> sep(",");
                tokenizer tokens(line[i], sep);
                vector<int> tmpcounts;
                for (tokenizer::iterator tok_iter = tokens.begin();  tok_iter != tokens.end(); ++tok_iter){
                        int tmp = atoi(tok_iter->c_str());
                        tmpcounts.push_back(tmp);
                }
                if (tmpcounts.size() != 2){
                	std::cerr << "ERROR: "<< line[i] << " does not have two alleles (expecting SNP data)\n";
                	exit(1);
                }
                topush.push_back(make_pair(tmpcounts[0], tmpcounts[1]));
            }
            allele_counts.push_back(topush);
            nsnp++;
    }
}


string CountData::get_pop_in_index(int index){
	string toreturn;
	map<string , int>::iterator it = head2id.begin();
	while (it != head2id.end()){
		if (it->second == index) return it->first;
		it++;
	}
	if (it == head2id.end()) {
		cerr << "Trying to get index "<< index << " in CountData, none found\n";
		exit(1);
	}
	return toreturn;
}

