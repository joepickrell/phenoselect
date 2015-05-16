/*
 * PhyloPop_params.h
 *
 *  Created on: Jul 1, 2011
 *      Author: pickrell
 */

#ifndef PHYLOPOP_PARAMS_H_
#define PHYLOPOP_PARAMS_H_

#include "Settings.hpp"

class PhyloPop_params{
public:
	PhyloPop_params();
	vector<string> pops;
	int nburn;
	int sfreq;

};

#endif /* PHYLOPOP_PARAMS_H_ */
