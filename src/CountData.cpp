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
	alfreqs = gsl_matrix_alloc(nsnp, npop);
	set_alfreqs();
}



double CountData::get_freq(int whichpop, int whichsnp){
	double toreturn;
	toreturn = gsl_matrix_get(alfreqs, whichsnp, whichpop);
	return toreturn;
}

string CountData::get_pops(){
	string toreturn = "(";
	map<string , int>::iterator it = pop2id.begin();
	map<string , int>::iterator it2 = pop2id.end();
	it2--;
	toreturn+= it->first +":0.1";
	it++;
	while (it != it2){
		toreturn+=",("+it->first+":0.1";
		it++;
	}
	toreturn+=","+it->first +":0.1";
	for (int i = 0; i < npop; i++)	toreturn+= "):0.1";


	toreturn = toreturn.substr(0, toreturn.size()-4);
	toreturn+= ";";
	return toreturn;
}


vector<string> CountData::list_pops(){
	vector<string> toreturn;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++) {
		if ( params->restrict_pop == true ){
			if (it->second < params->pops2use ) toreturn.push_back( it->first);
		}
		else toreturn.push_back( it->first);
	}
	return toreturn;
}


void CountData::read_counts(string infile){
	//cout << "here?\n"; cout.flush();
    allele_counts.clear();
    pop2id.clear();
    id2pop.clear();
    npop = 0;
    nsnp = 0;
    rss.clear();
    chr.clear();
    pos.clear();
    a1.clear();
    a2.clear();
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
   // cout << "here2\n"; cout.flush();
    while (ss>> buf){
    	line.push_back(buf);
     }

    /*
     * make map from header, number populations according to order
     */
    int start = 0;
    //cout << "here4\n"; cout.flush();
    for(int i = start; i < line.size(); i++) {
    	//cout << line.size() << "\n"; cout.flush();
    	pop2id.insert(make_pair(line[i], i-start));
    	id2pop.insert(make_pair(i-start, line[i]));
    	npop ++;
    }
   // cout << "here3\n"; cout.flush();
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
            for ( int i = start; i < line.size(); i++){
            	//cout <<  line[i] << "\n";
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

pair< vector<string>, vector<double> > CountData::get_freqs(int i){
	pair<vector<string>, vector<double> > toreturn;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		toreturn.first.push_back(it->first);
		int j = it->second;
		double f = gsl_matrix_get(alfreqs, i, j);
		toreturn.second.push_back(f);
	}
	return toreturn;


}


pair< vector<string>, vector<double> > CountData::get_centered_freqs(int i){
	pair<vector<string>, vector<double> > toreturn;
	double mean = 0;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		int j = it->second;
		double f = gsl_matrix_get(alfreqs, i, j);
		if (!isnan(f)) mean+= f;
		//toreturn.second.push_back(f);
	}
	mean = mean/ (double) npop;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		toreturn.first.push_back(it->first);
		int j = it->second;
		double f = gsl_matrix_get(alfreqs, i, j) - mean;
		toreturn.second.push_back(f);
	}
	return toreturn;
}

void CountData::set_alfreqs(){
	mean_ninds.clear();
	mean_hzy.clear();
	id2nsnp.clear();

	for (int i = 0; i < npop; i++){
		mean_ninds.insert(make_pair(i, 0.0));
		mean_hzy.insert(make_pair(i, 0.0));
		id2nsnp.insert(make_pair(i, 0));
	}
	for (int i = 0; i < nsnp; i++){
		for (int j = 0; j < npop; j++){
			int c1 = allele_counts[i][j].first;
			int c2 = allele_counts[i][j].second;
			double f = (double) c1 / ( (double) c1 + (double) c2 );
			if ( c1+c2 < 1){
				cerr << "Warning: no counts at SNP "<< i << " population "<< j <<"\n";
				gsl_matrix_set(alfreqs, i, j, f);
				continue;
			}
			gsl_matrix_set(alfreqs, i, j, f);
			mean_ninds[j] += ((double) c1+ (double) c2)/2.0;
			double tmp2 = (double) c2 / ((double) c1+ (double) c2 - 1.0);
			double tmphzy = 2* f * tmp2;
			if (c1+c2 < 2){
				tmphzy = 2*f*(1-f);
			}
			//if (id2pop[j] == "San"){
			//	cout << i << " "<< c1 << " "<< c2 << " "<< f << " "<< tmp2 << " "<< tmphzy << " "<< mean_hzy[j]<< "\n";
		//	}
			mean_hzy[j] += tmphzy; //2*f*(1-f);
			id2nsnp[j]++;
		}
	}
	for (int i = 0; i < npop; i++){
		mean_ninds[i] = mean_ninds[i]/ id2nsnp[i];
		mean_hzy[i] = mean_hzy[i]/ id2nsnp[i];
		//cout << id2pop[i] << " "<< mean_hzy[i] << "\n";
	}
}


void CountData::scale_alfreqs(){
	for (int i = 0; i < nsnp; i++){
		double total = 0;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			double scaled;
			if (params->alfreq_scaling == 1) scaled = asin(sqrt(f));
			else scaled = f;
			total = total+scaled;
			gsl_matrix_set(alfreqs, i, j, scaled);
		}

		double m = total/ (double) npop;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			if (params->micro) gsl_matrix_set(alfreqs, i, j, f-m);
			else{
				double f = gsl_matrix_get(alfreqs, i, j);
				if (params->alfreq_scaling == 3) {
					gsl_matrix_set(alfreqs, i, j, (f-m)/sqrt(m *(1-m)) );
					if (m < 1e-8) gsl_matrix_set(alfreqs, i, j, 0);
				}
				else if (params->alfreq_scaling == 4) gsl_matrix_set(alfreqs, i, j, f);
				else gsl_matrix_set(alfreqs, i, j, f-m);
			}
		}
	}
}

void CountData::print_alfreqs(string outfile){
	ogzstream out(outfile.c_str());
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++)	out << it->first << " ";
	out << "\n";
	for (int i = 0; i < nsnp ; i++){
		for(map<string, int>::iterator it2 = pop2id.begin(); it2 != pop2id.end(); it2++)	 out << " "<< gsl_matrix_get(alfreqs, i, it2->second);
		out << "\n";
	}
}


string CountData::get_pop_in_index(int index){
	string toreturn;
	map<string , int>::iterator it = pop2id.begin();
	while (it != pop2id.end()){
		if (it->second == index) return it->first;
		it++;
	}
	if (it == pop2id.end()) {
		cerr << "Trying to get index "<< index << " in CountData, none found\n";
		exit(1);
	}
	return toreturn;
}

