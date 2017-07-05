#include "SimpleShells/include/SimpleShellsControlArchitecture.h"

#include <algorithm>
#include <set>
#include <math.h>       /* log */



/**
 * If defined, the agent will log its position to stdout.
 * NOTE: DO NOT USE FOR REGULAR RUNS!! Creates vast amounts of output!!
 */
//#define TRACK_AGENTS

bool genomeSortByFitness(const Genome & genomeOne, const Genome & genomeTwo) {
	return genomeOne.fitness > genomeTwo.fitness;
}

/// Fitness should be positive
Genome SimpleShellsControlArchitecture::selectWeighted(std::vector<Genome> & genomes, bool onFitness) {
	std::vector<Genome>::const_iterator a, b, c;
	// std::cout  << "Iteration: " << gWorld->getIterations() << "; Size: " << genomes.size();

	if (genomes.size() < 2)
	{
		// std::cout << std::endl;
		return *(genomes.begin());
	}

	// else

	// NOTE : expect sorted genomes on fitnesses.
	a = genomes.begin();
	b = genomes.begin() + genomes.size() / 2;
	c = genomes.begin() + genomes.size() - 1;

	// if >= 2 then do tournament, else do rank-based roulette
	if (_tournamentSize >= 2) { // true in default

		std::set<size_t> tournament;
		size_t tournamentSize( std::min(_tournamentSize, static_cast<size_t>(genomes.size())) );
		size_t selected(genomes.size());

		// randomly select tournamentSize many genomes from pool with replacement
		for (unsigned i = 0; i < tournamentSize; ++i) {
			size_t candidate(static_cast<size_t>( Rand::randint(0, genomes.size()) ));
			if (tournament.find(candidate) != tournament.end()) // try again if duplicate candidate drawn
				--i;
			else {
				tournament.insert(candidate);
				selected = std::min(candidate, selected);
			}
		}

		// std::cout  << "; tournament: "  << selected <<  std::endl;

		// Note: genomes expected to be sorted!
		// returns the winner = the most front randomly selected genome (selected is index offset)
		return *(genomes.begin() + selected);
	}

	double n(genomes.size() - 1);
	double pos(0.0);

	// Note: 1 <= _selectionPressure <= 2

	double fitnessTotal = .0;
	double rankedTotal = .0;

	if(onFitness) {
		// Loop reversed and replace fitness w. rank-based variant
		for (std::vector<Genome>::reverse_iterator genome = genomes.rbegin(); genome != genomes.rend(); genome++) {
			fitnessTotal += genome->fitness;
			genome->fitness = 2.0 - _selectionPressure + 2.0 * (_selectionPressure - 1.0) * pos/n;
			rankedTotal += genome->fitness;
			pos++;
		}

		double random = Rand::randouble() * rankedTotal;
		// std::cout  << "; Fitness: " << a->fitness << "; " << b->fitness << "; " << c->fitness << "; Fitness avg: " << fitnessTotal / genomes.size() << "; roulette wheel: "  << random <<  std::endl;
		for (std::vector<Genome>::const_iterator it = genomes.begin(); it < genomes.end(); it++) {
			random -= it->fitness;
			if (random <= 0)
				return *it;
		}
	} else { // on novelty
		// Loop reversed and replace fitness w. rank-based variant
		for (std::vector<Genome>::reverse_iterator genome = genomes.rbegin(); genome != genomes.rend(); genome++) {
			fitnessTotal += genome->novelty;
			genome->novelty = 2.0 - _selectionPressure + 2.0 * (_selectionPressure - 1.0) * pos/n;
			rankedTotal += genome->novelty;
			pos++;
		}

		double random = Rand::randouble() * rankedTotal;
		// std::cout  << "; Fitness: " << a->novelty << "; " << b->novelty << "; " << c->novelty << "; Fitness avg: " << fitnessTotal / genomes.size() << "; roulette wheel: "  << random <<  std::endl;
		for (std::vector<Genome>::const_iterator it = genomes.begin(); it < genomes.end(); it++) {
			random -= it->novelty;
			if (random <= 0)
				return *it;
		}
	}

	// Should never get here
	return *(genomes.begin());
}


void mutate(std::vector<double> & parameters, double delta) {

	for (std::vector<double>::iterator it = parameters.begin(); it < parameters.end(); it++)
	{
		//*it += (ranf() - 0.5) * delta;
		*it += Rand::gaussianf(0.0, delta);
	}
}


SimpleShellsControlArchitecture::SimpleShellsControlArchitecture( RobotAgentWorldModel *__wm ) : BehaviorControlArchitecture ( __wm ) {
	_wm = (SimpleShellsAgentWorldModel*)__wm;

	int tmpInt = 0;
	gProperties.checkAndGetPropertyValue("gMaxGenePool", &tmpInt, true);
	gMaxGenePool = tmpInt;
	gProperties.checkAndGetPropertyValue("gMaxLifetimeGathering", &tmpInt, true);
	_maxLifetime[PHASE_GATHERING] = tmpInt;
	gProperties.checkAndGetPropertyValue("gMaxLifetimeMating", &tmpInt, true);
	_maxLifetime[PHASE_MATING] = tmpInt;
	gProperties.checkAndGetPropertyValue("gHiddenNeuronCount", &_hiddenNeuronCount, true);

	gProperties.checkAndGetPropertyValue("gTournamentSize", &tmpInt, 0);
	_tournamentSize = static_cast<unsigned>(tmpInt);
	// std::cout << "Tournament size: " << _tournamentSize << std::endl;

	_randomSelection = false;
	gProperties.checkAndGetPropertyValue("gRandomSelection", &_randomSelection, false);
	_useMarket = true;
	gProperties.checkAndGetPropertyValue("gUseMarket", &_useMarket, true);
	_useSpecBonus = false;
	gProperties.checkAndGetPropertyValue("gUseSpecBonus", &_useSpecBonus, false);
	_task1Premium = 1.0;
	gProperties.checkAndGetPropertyValue("gTask1Premium", &_task1Premium, 1.0);
	_selectionPressure = 1.5;
	gProperties.checkAndGetPropertyValue("gSelectionPressure", &_selectionPressure, true);

	//manuel: here the number of parameters is calculated
	if (_hiddenNeuronCount > 0) {
		_parameterCount = ((_wm->_sensorCount * (gPuckColors + 1) + 1 + 2) * _hiddenNeuronCount + (_hiddenNeuronCount + 1) * 2 ); //manuel: 1 more for direction change
		_response.assign(_hiddenNeuronCount, .0);
	} else {
		_parameterCount = 2 + (_wm->_sensorCount * (gPuckColors + 1) + 1 + 2) * 2;
		// 2 + (8 * (2 + 1) + 1 + 2) * 2 = 
		// 2 + 27 * 2 = 56
	}

	_wm->_genePool.reserve(gMaxGenePool);
	_nearbyGenomes.reserve(gNbOfAgents); // Agent counter is unpredictable at this time

	//manuel: i guess here all parameters are set to 0
	_activeGenome.parameters.assign(_parameterCount, .0);
	//manuel: genomes are initalized to some random between 0 and 1
	mutate(_activeGenome.parameters, 1.0);
	
	_wm->_puckCounters = &(_activeGenome.pucks);
    _age = 0;

    // NOV
    // _activeGenome.archive.reserve(gMaxArchiveSize);
    _fitnessThreshold = 1.5;
	gProperties.checkAndGetPropertyValue("gfitnessThreshold", &_fitnessThreshold, true);
 	_archiveThreshold = 3.0;
	gProperties.checkAndGetPropertyValue("gArchiveThresholdStart", &_archiveThreshold, true);
	_lastAddedToArchive = 0;
	_activeGenome.novelty = 0.0;
    calcBehaviourAndId();
}


SimpleShellsControlArchitecture::~SimpleShellsControlArchitecture() {
	// nothing to do.
}


void SimpleShellsControlArchitecture::reset() {

	logStats();

	_wm->reset(_maxLifetime); // Genepool is cleaned here
    _age = 0;
	// Set all puck counters to zero
	_activeGenome.pucks.assign(gPuckColors, 0);
	_activeGenome._distance = 0.0;
	_activeGenome._lastXPos = _wm->getXReal();
	_activeGenome._lastYPos = _wm->getYReal();
	// Finally, update the active genome.
	_activeGenome.id = gNextGenomeId++;

	// NOV
	_activeGenome.age = 0;
    _activeGenome.novelty = 0.0;

	// std::cout << "[descendant] " <<  _wm->_world->getIterations() << ' ' << _wm->_winnerId << ' ' << _activeGenome.id << '\n';
}

void SimpleShellsControlArchitecture::logStats()
{
	//
        // Dump lifetime stats: time, ID, amount of pucks gathered, distance covered, lifetime
        //
        // std::cout << "[gathered] " << _wm->_world->getIterations() << ' ' << _activeGenome.id;
        // for (int i = 0; i < _activeGenome.pucks.size(); ++i) {
        //         if (i == _wm->_energyPuckId && _wm->_excludeEnergyPucks)
        //                 // std::cout << ' ' << _wm->_nrBoosts;
        //         else
        //                 // std::cout << ' ' << _activeGenome.pucks[i];
        // }

        // std::cout << ' ' << _activeGenome._distance;
        // std::cout << ' ' << _age;
        // std::cout << '\n';
}


void SimpleShellsControlArchitecture::prepareShutdown() {
	logStats();
}
    

/* obsolete */
double getSpecDeg(std::vector<int> & counts) {
	// Sum counts
	double totalCount = gPuckColors;
	for (uint i = 0; i < counts.size(); i++) {
		totalCount += counts[i];
	}
	// Translate counts to fractions, treating count of zero as a very small fraction
	std::vector<double> fractions;
	fractions.resize(counts.size());
	for (uint i = 0; i < counts.size(); i++) {
		fractions[i] = (counts[i] + 1) / totalCount;
	}
	// Calculate entropy of the given fractions
	double entropy = .0;
	for (uint i = 0; i < counts.size(); i++) {
		entropy -= fractions[i] * std::log(fractions[i]);
	}
	// Calculate base entropy (assuming all counts are same)
	double baseFraction = 1.0 / gPuckColors;
	double baseEntropy = -gPuckColors * baseFraction * std::log(baseFraction);
	// Calculate how much our case is different than this base = degree of specialization.
	// Will be around 1.0 for perfect specialist, and around .0 for perfect generalist.
	// Given that puck counts are skewed, all will be seen as specialists of small degree
	return 1.0 - entropy / baseEntropy;
}

void SimpleShellsControlArchitecture::assignFitness(std::vector<Genome> &genomes) {
	// In the default properties file I'm using
	// _excludeEnergyPucks = false, so all pucks are counted!
	// Fitness = 0.005 + # all pucks collected by genome

	// Calculate how many pucks of each color and how many puck in total were gathered.
	std::vector<int> puckTotals(gPuckColors, 0);
	int cumTotal = 0;

	// Count different colors over received genomes
	for (std::vector<Genome>::iterator itGenome = genomes.begin(); itGenome < genomes.end(); itGenome++) {
		for (int i = 0; i < gPuckColors; i++) {

			if (i != _wm->_energyPuckId  || !_wm->_excludeEnergyPucks) { // true in default
				puckTotals[i] += itGenome->pucks[i];
			}
		}
	}

	// Calc cumulative total
	for (int i = 0; i < gPuckColors; i++)
	{
		if (i != _wm->_energyPuckId  || !_wm->_excludeEnergyPucks) { // true in default
			cumTotal += puckTotals[i];
			//std::cout << puckTotals[i] << ' ';
		}
	}
	//std::cout << "(tot: " << cumTotal << ") ";

	// Calculate the puck prices (depending on puck fraction)
	std::vector<double> puckPrices;
	puckPrices.resize(gPuckColors);

	for (int i = 0; i < gPuckColors; i++) {
		if (i != _wm->_energyPuckId || !_wm->_excludeEnergyPucks) {

			if (_useMarket) { // false in default
				puckPrices[i] = (puckTotals[i] == 0) ? 0.0 : double(cumTotal) / double(puckTotals[i]);
			} else {
				puckPrices[i] = 1.0;
			}
			//std::cout << puckPrices[i] << ' ';
		}
	}
	//std::cout << std::endl;

	// TODO: better premium calculation for pucks
	puckPrices[0] *= _task1Premium; // = 1 in default

	// Finally, use the prices to calculate fitness of a genome (and store it within)
	for (std::vector<Genome>::iterator itGenome = genomes.begin(); itGenome < genomes.end(); itGenome++) {
		// We are using roulette selection, so genomes which gathered no pucks will still receive a small fitness.
		double rating = .005;

		for (int i = 0; i < gPuckColors; i++) {
			if (i != _wm->_energyPuckId || !_wm->_excludeEnergyPucks) { // true in default
				rating += double(itGenome->pucks[i]) * puckPrices[i];
				// std::cout << itGenome->pucks[i] << ' ';
			}
		}
		
		if (_useSpecBonus) { // false by default
			std::cout << "here" << std::endl;
			// Calculate specialization degree (value in [0, 1] with 1 corresponding to ideal specialist)
			double degree = getSpecDeg(itGenome->pucks);
			// Use inverse logistic function to bring the specialization degree value from typical 0.3 .. 0.7 region into 0 .. 30 region, in such a way, that both specialists and generalists ends will be a bit amplified.
			double value = std::log(1.01/(degree+0.01) - 0.99) * 10.0 + 15.0;
			if (value < 0) value = 0; else if (value > 30) value = 30;
			itGenome->fitness = rating * (1.0 + value);
		} else {
			itGenome->fitness = rating;
		}
		// std::cout << itGenome->fitness << std::endl;
	}

	// Assign fitness to _activeGenome for MCNS comparison
	double rating = .005;

	for (int i = 0; i < gPuckColors; i++) {
		if (i != _wm->_energyPuckId || !_wm->_excludeEnergyPucks) {
			rating += double(_activeGenome.pucks[i]) * puckPrices[i];
		}
	}
	_activeGenome.fitness = rating;
}

// Just gather all nearby genomes. Duplicates will be handled elsewhere.
void SimpleShellsControlArchitecture::gatherGenomes(std::vector<Genome*> & genePool, int commDistSquared) {
	genePool.clear();
	std::vector<RobotAgentPtr>* agents = gWorld->listAgents();
	std::vector<RobotAgentPtr>::iterator itAgent;
	if (commDistSquared < 0) {
		for (itAgent = agents->begin(); itAgent < agents->end(); itAgent++) {
			RobotAgentPtr& agent = *itAgent;
			SimpleShellsControlArchitecture* controller = static_cast<SimpleShellsControlArchitecture*>(agent->getBehavior());
			genePool.push_back(controller->getGenome()); // only a ref is returned here, so no copy!
		}
	} else {
		for (itAgent = agents->begin(); itAgent < agents->end(); itAgent++) {

			RobotAgentPtr& agent = *itAgent;
			if (isRadioConnection(agent->getWorldModel()->_agentId,_wm->_agentId)) {
				SimpleShellsControlArchitecture* controller = static_cast<SimpleShellsControlArchitecture*>(agent->getBehavior());
				genePool.push_back(controller->getGenome());
			}
		}
	}
}


// check if newly gathered genomes are the same as previously gathered and update their puckcount if so
void SimpleShellsControlArchitecture::updateGenomes() {
//std::cout << "Gathering Genomes...";
	gatherGenomes(_nearbyGenomes, gMaxCommDistSquared);

//std::cout << _nearbyGenomes.size() << " found\n";
	std::vector<Genome*>::iterator itNewGenomes;
	std::vector<Genome>::iterator itGatheredGenomes;
	for (itNewGenomes = _nearbyGenomes.begin(); itNewGenomes < _nearbyGenomes.end(); itNewGenomes++) {
		bool sameDNA = false;
		// Check if new genome has actually same vector as one of the gathered.
		for (itGatheredGenomes = _wm->_genePool.begin(); itGatheredGenomes < _wm->_genePool.end(); itGatheredGenomes++) {
			sameDNA = (itGatheredGenomes->id == (*itNewGenomes)->id);
			if (sameDNA) break;
		}
		// If parameter vector is the same, update puck counts within the genome (which are increasing monotonously).
		if (sameDNA) {
			for (int i = 0; i < gPuckColors; i++) itGatheredGenomes->pucks[i] = (*itNewGenomes)->pucks[i];
		} else {
			Genome newGene = **itNewGenomes;
			if (!gWithSync) {
				newGene.archive.clear();
			}
			_wm->_genePool.push_back(newGene); // This makes a full deep copy 

		if (_wm->_genePool.size() == gMaxGenePool) {
				// stop gathering genomes once all individuals are in _genePool.
				// FIXME: ugly hack depends on interpretation by worldmodel. Calling nextAction will skip selection phase. Update state machine implementation to be more robust.
				_wm->_lifetime[_wm->_phase] = 0;
				break;
			}
		}
	}
}

// manuel: making the move: desired values calcucated with use of nn weights
void SimpleShellsControlArchitecture::updateActuators() {
	std::vector<double> & parameters = _activeGenome.parameters;
	int geneToUse = 0;

	double trans = .0, rotor = .0;

	if (_hiddenNeuronCount > 0) {
		// Calculating response of the hidden layer.
		// Initiate with biases.
		for (int n = 0; n < _hiddenNeuronCount; n++) {
			_response[n] = 1.0 * parameters[geneToUse++];
			// manuel: making them random may not be useful to me
			// Using those in hopes to initiate feedback loop, which will make trajectories more random, and hence more likely to encounter pucks
			_response[n] += parameters[geneToUse++] * _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
			_response[n] += parameters[geneToUse++] * _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
		}

		// Distance to obstacle along each of the sensor rays.
		// Distances to pucks of various colors along each of the sensor rays.
		for (int s = 0; s < _wm->_sensorCount; s++) {
			double normObstacleRange = _wm->_rangeSensors[s]->getNormObstacleRange();
			for (int n = 0; n < _hiddenNeuronCount; n++) {
				_response[n] += parameters[geneToUse++] * normObstacleRange;
			}
			for (int c = 0; c < gPuckColors; c++) {
				double normPuckRange = _wm->_rangeSensors[s]->getNormPuckRange(c);
				for (int n = 0; n < _hiddenNeuronCount; n++) {
					_response[n] += parameters[geneToUse++] * normPuckRange;
				}
			}
		}

		// Activation
		for (int n = 0; n < _hiddenNeuronCount; n++) {
			_response[n] = tanh(_response[n]);
		}

		for (int n = 0; n < _hiddenNeuronCount; n++) {
			trans += _response[n] * parameters[geneToUse++];
			rotor += _response[n] * parameters[geneToUse++];
		}
	} else if (false) {
		// for collecting sample sensor inputs
		trans += parameters[geneToUse++] * _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
		trans += parameters[geneToUse++] * _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
		rotor += parameters[geneToUse++] * _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
		rotor += parameters[geneToUse++] * _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
		for (int s = 0; s < _wm->_sensorCount; s++) {
			double normObstacleRange = _wm->_rangeSensors[s]->getNormObstacleRange();
			trans += normObstacleRange * parameters[geneToUse++];
			rotor += normObstacleRange * parameters[geneToUse++];
			for (int c = 0; c < gPuckColors; c++) {
				double normPuckRange = _wm->_rangeSensors[s]->getNormPuckRange(c);
				trans += normPuckRange * parameters[geneToUse++];
				rotor += normPuckRange * parameters[geneToUse++];
			}
		}
	} else {
		// we don't use hidden, so the calculuations are done here
		trans += parameters[geneToUse++] * _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
		trans += parameters[geneToUse++] * _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
		rotor += parameters[geneToUse++] * _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
		rotor += parameters[geneToUse++] * _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
		for (int s = 0; s < _wm->_sensorCount; s++) {
			double normObstacleRange = _wm->_rangeSensors[s]->getNormObstacleRange();
			trans += normObstacleRange * parameters[geneToUse++];
			rotor += normObstacleRange * parameters[geneToUse++];
			for (int c = 0; c < gPuckColors; c++) {
				double normPuckRange = _wm->_rangeSensors[s]->getNormPuckRange(c);
				trans += normPuckRange * parameters[geneToUse++];
				rotor += normPuckRange * parameters[geneToUse++];
			}
		}
	}

	// Output layer
	// Initiate with biases.
	trans += 1.0 * parameters[geneToUse++];
	rotor += 1.0 * parameters[geneToUse++];	

	trans = tanh(trans);
	rotor = tanh(rotor);

	_wm->_desiredTranslationalValue = trans * gMaxTranslationalSpeed * _wm->_speedPenalty ;
	_wm->_desiredRotationalVelocity = rotor * gMaxRotationalSpeed * _wm->_speedPenalty;
}

// *********************************** NOV ***********************************

void SimpleShellsControlArchitecture::calcBehaviourAndId() {
	// Run a genome on sample sensor inputs, to simulate and record behavior for novelty comparision
	// Nr sensors = 8
	// Nr pucktyes = 2
	// So in loop: 1 obstacle-range, 2 puck-range, 1 ob, 2 puck ... = 8 x 3 = 24 inputs from sensors

	int sensorToUse, geneToUse;
	double trans, rotor;
	std::vector<double> &parameters = _activeGenome.parameters;
	std::hash<std::string> hasher;
	std::string behaviour2str = "";

	_activeGenome.simSensorResults.clear();

	for(unsigned int i = 0; i < simSensorInput.size(); i++) {
		std::vector<int8_t> result;
		std::vector<double> sensorRanges = simSensorInput[i];
		sensorToUse = 0;
		geneToUse = 0;
		trans = 0.0;
		rotor = 0.0;

		trans += parameters[geneToUse++] * _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
		trans += parameters[geneToUse++] * _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
		rotor += parameters[geneToUse++] * _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
		rotor += parameters[geneToUse++] * _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
		
		for (int s = 0; s < _wm->_sensorCount; s++) {
			double normObstacleRange = sensorRanges[sensorToUse++];

			trans += normObstacleRange * parameters[geneToUse++];
			rotor += normObstacleRange * parameters[geneToUse++];
			for (int c = 0; c < gPuckColors; c++) {
				double normPuckRange = sensorRanges[sensorToUse++];

				trans += normPuckRange * parameters[geneToUse++];
				rotor += normPuckRange * parameters[geneToUse++];
			}
		}	
		// Output layer
		// Initiate with biases.
		trans += 1.0 * parameters[geneToUse++];
		rotor += 1.0 * parameters[geneToUse++];	

		// Binarize resulting behaviour
		int trans2 = binarize(trans);
		int rotor2 = binarize(rotor);
		// ... and save back to genome
		result.push_back(trans2);
		result.push_back(rotor2);
		_activeGenome.simSensorResults.push_back(result);
		// Add for Id generation
		behaviour2str += std::to_string(trans2);
		behaviour2str += std::to_string(rotor2);

	}	
   	unsigned long hashed = hasher(behaviour2str);
   	_activeGenome.behaviourId = hashed;
}

int8_t SimpleShellsControlArchitecture::binarize(double d) {
	// Binarize for simpler comparision
	return (d > 0) ? 1 : 0;
}

bool genomeSortByNovelty(const Genome &g1, const Genome &g2) {
	return g1.novelty > g2.novelty;
}

bool genomeSortByDistance(const double dist1, const double dist2) {
	return dist1 < dist2;
}

bool genomeSortById(const Genome &g1, const Genome &g2) {
	return g1.behaviourId > g2.behaviourId;
}

int SimpleShellsControlArchitecture::calcPucksCollected(const Genome &g) {
		int result = 0;
		for (int i = 0; i < gPuckColors; i++) {
			if (i != _wm->_energyPuckId || !_wm->_excludeEnergyPucks) {
				result += double(g.pucks[i]);
			}
		}
		return result;
	}

void SimpleShellsControlArchitecture::select() {
	// If the egg has gathered enough genomes or it's time is up, reactivate it with one of the gathered genomes.

	if(gFitnessMode == 1) {
		selectClassic();
	} else if (gFitnessMode == 2) {
		selectNovelty();
	} else if (gFitnessMode == 3) {
		//nothing for now
	} else if (gFitnessMode == 4) {
		selectRandom();
	}
}

void SimpleShellsControlArchitecture::selectClassic() {
	// Assign fitness based on puckcount
	assignFitness(_wm->_genePool);
	// Do roulette selection (need to sort genomes vector first, because we are calculating medians for logging)
	std::sort(_wm->_genePool.begin(), _wm->_genePool.end(), genomeSortByFitness);
	Genome winner = chooseParent(true);
	// std::vector<Genome>::const_iterator winner = selectWeighted(_wm->_genePool);
	_wm->_winnerId = winner.id;
	// Reactivate with the winning genome
	for (int i = 0; i < _parameterCount; i++) {
		_activeGenome.parameters[i] = winner.parameters[i];
	}	
	// Mutate a bit
	mutate(_activeGenome.parameters, 0.1);	
	// Do some logging
	calcBehaviourAndId();

	_wm->dumpGenePoolStats();
}

void SimpleShellsControlArchitecture::selectNovelty() {
	// Assign fitness based on novelty

	//log poolsize
	poolSizes.push_back(_wm->_genePool.size());

	// Calc differences between bahaviors
	assignNovelty(_wm->_genePool, false, false);
	// Find winner
	std::sort(_wm->_genePool.begin(), _wm->_genePool.end(), genomeSortByNovelty);
	Genome winner = chooseParent(false);
	_wm->_winnerId = winner.id;
	gWinners.push_back(winner.id);
	// Reactivate with the winning genome
	for (int i = 0; i < _parameterCount; i++) {
		_activeGenome.parameters[i] = winner.parameters[i];
	}
	// Mutate a bit
	mutate(_activeGenome.parameters, 0.1);
	// Sync with winner archive
	if (gWithSync) {
		syncArchives(winner);	
		winner.archive.clear();
	}
	// .. and add to archive?
	evaluateIndividual(_wm->_genePool, winner, true, false);
	// Calc behaviour and Id for new genome
	calcBehaviourAndId();
	// Do some logging
	logPool(_wm->_genePool, winner);
	_wm->dumpGenePoolStats();
}

void SimpleShellsControlArchitecture::selectRandom() {
	// Select a random genome (for baseline testing)

	Genome winner = _wm->_genePool.at(Rand::randint(0, _wm->_genePool.size()));
	_wm->_winnerId = winner.id;
	// Reactivate with the winning genome
	for (int i = 0; i < _parameterCount; i++) {
		_activeGenome.parameters[i] = winner.parameters[i];
	}	
	// Mutate a bit
	mutate(_activeGenome.parameters, 0.1);
	// Do some logging
	_wm->dumpGenePoolStats();
}

Genome SimpleShellsControlArchitecture::chooseParent(bool onFitness) {
	return selectWeighted(_wm->_genePool, onFitness);
}

bool SimpleShellsControlArchitecture::aboveFitnessThreshold(Genome &winner) {
	bool bigger = winner.fitness >= (_activeGenome.fitness * _fitnessThreshold);
	
	// TODO dynamic threshold
	// if (bigger) {
	// 	_fitnessThreshold -= 0.1; // take all the fitness you can get
	// } else {
	// 	_fitnessThreshold += 0.3; // try explore more
	// }
	return bigger;
}

void SimpleShellsControlArchitecture::syncArchives(Genome &winner) {
	int i, j;
	std::vector< std::pair<Genome, int> > newArchive;
	// Backwards loop both archives
	i = _activeGenome.archive.size() - 1;
	j = winner.archive.size() - 1;
	while ((int)newArchive.size() < gMaxArchiveSize && j >= 0 && i >= 0) {
		if (_activeGenome.archive[i].second > winner.archive[j].second) {
			newArchive.push_back(_activeGenome.archive[i]);
			i--;
		} else if (_activeGenome.archive[i].second == winner.archive[j].second) {
			if (_activeGenome.archive[i].first.behaviourId == winner.archive[j].first.behaviourId) {
				newArchive.push_back(_activeGenome.archive[i]);
			} else {
				newArchive.push_back(_activeGenome.archive[i]);
				newArchive.push_back(winner.archive[j]);
			}
				i--;
				j--;
		} else if (_activeGenome.archive[i].second < winner.archive[j].second) {
			newArchive.push_back(winner.archive[j]);
			j--;
		}
	}
	
	std::reverse(newArchive.begin(), newArchive.end());
	_activeGenome.archive = newArchive;
}

void SimpleShellsControlArchitecture::copyGenome(Genome &g) {
	g.id = _activeGenome.id;
	g.parameters = _activeGenome.parameters;
	g.pucks = _activeGenome.pucks;
	g.fitness = _activeGenome.fitness;
	g._distance = _activeGenome._distance;
	g._lastXPos = _activeGenome._lastXPos;
	g._lastYPos = _activeGenome._lastYPos;
	g.simSensorResults = _activeGenome.simSensorResults;
}

bool SimpleShellsControlArchitecture::aboveArchiveThreshold(double result) {
	// Adjust threshold if necessary 
	int diff;

	diff = gWorld->getIterations() - _lastAddedToArchive;
	// numbers taken from novelty paper
	if(diff > 2500) { 
		_archiveThreshold *= 0.95;
	} else if (diff < 500) {
		_archiveThreshold *= 1.2;
	}
	return result >= _archiveThreshold;
}

bool SimpleShellsControlArchitecture::genomeInArchive(const Genome &ind) {
	for (auto &g : _activeGenome.archive) {
		if (ind.behaviourId == g.first.behaviourId) {
			return true;
		}
	}
	return false;
}

void SimpleShellsControlArchitecture::addGenomeToArchive(Genome ind) {
	if (!genomeInArchive(ind)) {
		addedToArchCount++;
		if((int)_activeGenome.archive.size() >= gMaxArchiveSize) {
			_activeGenome.archive.erase(_activeGenome.archive.begin()); // delete the oldest
		}
		// Add genome and timestamp
		_activeGenome.archive.push_back(std::make_pair(ind, gWorld->getIterations()));
		_lastAddedToArchive = gWorld->getIterations();
	}
}

void SimpleShellsControlArchitecture::assignNovelty(std::vector<Genome> &genomes, bool addToArchive, bool logging) {
	// Evaluate the novelty of the gathered genepool
	// From https://github.com/jal278/novelty-search-cpp

	for(auto &ind : genomes) {
		evaluateIndividual(genomes, ind, addToArchive, logging);
	}
}

void SimpleShellsControlArchitecture::evaluateIndividual(std::vector<Genome> &genomes, Genome &ind, bool addToArchive, bool logging) {
	// Evaluate the novelty of a single individual
	// From https://github.com/jal278/novelty-search-cpp

	double result;
	// Assign fitness according to average novelty
	if(addToArchive) { // Consider adding a genome to the archive
		int r = rand() % 1000; // Add any Genome with 0.1% chance
		int chance = rand() % 1000; // Add any Genome with 1% chance

		result = noveltyAvgNN(genomes, ind, addToArchive, logging);
		if(aboveArchiveThreshold(result) || _activeGenome.archive.size() < 1 || r == chance) {
			addGenomeToArchive(ind);
		} 
	} else { 
		result = noveltyAvgNN(genomes, ind, addToArchive, logging); // add genepool in kNN function
		ind.novelty = result;
	}
}

double SimpleShellsControlArchitecture::noveltyAvgNN(std::vector<Genome> &genomes, Genome &ind, bool addToArchive, bool logging) {
		// From https://github.com/jal278/novelty-search-cpp
		// 2 types of calls: with pop or without
		// called on 1 individual A to assess its novelty
		
		// 1. Make tuple of all in pop and archive: for each the novelty in relatio to A (calculated) and A itself
		std::vector<double> novelties;

		if(logging) {
			novelties = mapNoveltyPop(genomes, ind);
		} else {
			if (addToArchive) {
				// Map the novelty metric across the archive
				novelties = mapNoveltyArchive(ind);
			} else {
				// Map the novelty metric across the archive + current population
				novelties = mapNovelty(genomes, ind);
			}	
		}
		std::sort(novelties.begin(), novelties.end(), genomeSortByDistance);
		
		// 2. Init parameters
		int neigh = (!addToArchive) ? gKNNeighbours : 1; // Set k for kNN
		double density = 0.0;
		double sum = 0.0;
		double weight = 0.0;
		int i = 0;
		
		while (weight < neigh && i < (int)novelties.size()) {
			double term = novelties[i];
			double w = 1.0;

			// Sum of all distances to A
			sum += term * w;
			weight += w;		
			i++;				
		}
		if (weight != 0) {
			// Sum of distances / # neighbors (max = k = 15)
			density = sum/weight; 
		}
		return density;
	}

	double SimpleShellsControlArchitecture::calcBehaviorDistance(Genome &g1, Genome &g2) {
		// Calculate distance between 2 gnome's behaviour

		double diff = 0.0;

		for(unsigned int i = 0; i < g1.simSensorResults.size(); i++) {
			diff += hammingDistance(g1.simSensorResults[i], g2.simSensorResults[i]);
		}
		return diff;
	}

	double SimpleShellsControlArchitecture::hammingDistance(std::vector<int8_t> &v1, std::vector<int8_t> &v2) {
		// Calc eudlidean distance btw 2 vectors
		// dist((x, y), (a, b)) = √(x - a)² + (y - b)²

		double d1, d2;

		d1 = v1[0] - v2[0];
		d2 = v1[1] - v2[1];

		return std::abs(d1 + d2);
	}

	std::vector<double> SimpleShellsControlArchitecture::mapNovelty(std::vector<Genome> &genomes, Genome &ind) {
	// Map the novelty metric across the archive + current population
		std::vector<double> novelties;
		for(auto &g : _activeGenome.archive) {
			novelties.push_back(calcBehaviorDistance(g.first, ind));
		}
		for (auto &g : genomes) {
			novelties.push_back(calcBehaviorDistance(g, ind));
		}
		return novelties;
	}

	std::vector<double> SimpleShellsControlArchitecture::mapNoveltyArchive(Genome &ind) {
	// Map the novelty metric across the archive 
		std::vector<double> novelties;
		for(auto &g : _activeGenome.archive) {
			novelties.push_back(calcBehaviorDistance(g.first, ind));
		}
		return novelties;
	}

	std::vector<double> SimpleShellsControlArchitecture::mapNoveltyPop(std::vector<Genome> &genomes, Genome &ind) {
	// Map the novelty metric across some population of genomes, for diversity calc
		std::vector<double> novelties;
		for(auto &g : genomes) {
			novelties.push_back(calcBehaviorDistance(g, ind));
		}
		return novelties;
	}

	void SimpleShellsControlArchitecture::printAvgFitness() {
		if (_wm->_agentId == 5 && gWorld->getIterations()%500 == 0) {
			std::vector<Genome> genomes;
			double avg = 0;

			for (int i = 0; i < gWorld->getNbOfAgent(); i++) {
				SimpleShellsControlArchitecture *agent = (SimpleShellsControlArchitecture *)((gWorld->getAgent(i))->getBehavior());
				genomes.push_back(agent->_activeGenome);
			}
			assignFitness(genomes);
			for (int i = 0; i < (int)genomes.size(); i++) {
				avg += genomes[i].fitness;
			}
			std::cout << std::fixed << std::setprecision(3) <<  (avg / genomes.size()) << std::endl; 
		}	
	}

	void SimpleShellsControlArchitecture::printPucksPer1000() {
		if (_wm->_agentId == 5 && gWorld->getIterations()%1000 == 0) {
			std::cout << gPucksCollected << std::endl;
			std::cout << gPucksCollected1 + gPucksCollected2 << std::endl;
	
			
		}
	}

	void SimpleShellsControlArchitecture::collectAnalysis() {
		if (_wm->_agentId == 5 && gWorld->getIterations()%1000 == 0) {
			std::vector<Genome> genomes = logDiversity();
			logAvgsAndPucks(genomes);
			logWinners();	
		}
	}

	void SimpleShellsControlArchitecture::logWinners() {
		std::ofstream outfile;
		// std::string filepath = "nov-experiments/samplesize/winners_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		std::string filepath = "nov-experiments/archivesize/winners_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		// std::string filepath = "nov-experiments/fit/winners_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		outfile.open(filepath, std::ios_base::app);

		for (int a : gWinners) {
			outfile << std::to_string(a) + ","; 
		}
		outfile << "\n";

		gWinners.clear();

	}

	void SimpleShellsControlArchitecture::logAvgsAndPucks(std::vector<Genome> &genomes) {
		std::ofstream outfile;
		std::ofstream outfileIndCounters;
		// std::string filepath = "nov-experiments/samplesize/genstats_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		// std::string filepathIndCounters = "nov-experiments/samplesize/indstats_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		// std::string filepath = "nov-experiments/fit/genstats_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		// std::string filepathIndCounters = "nov-experiments/fit/indstats_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		std::string filepath = "nov-experiments/archivesize/genstats_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		std::string filepathIndCounters = "nov-experiments/archivesize/indstats_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		outfile.open(filepath, std::ios_base::app);
		outfileIndCounters.open(filepathIndCounters, std::ios_base::app);

		double avgDist = 0.0;
		double avgLife = 0.0;
		int avgPool = 0;

		//calc avg pool
    	avgPool =  std::accumulate(poolSizes.begin(), poolSizes.end(), 0) / poolSizes.size();

		for (int i = 0; i < gWorld->getNbOfAgent(); i++) {
			avgLife += ((SimpleShellsControlArchitecture *)((gWorld->getAgent(i))->getBehavior()))->_wm->_lifetime[PHASE_GATHERING];
		}
		avgLife /= gWorld->getNbOfAgent();

		for (auto &g : genomes) {
			avgDist += g._distance;
			outfileIndCounters << std::to_string(g.id) + "," + std::to_string(g._distance) + "," + std::to_string(g.pucks.at(0) + g.pucks.at(1)) + "," + std::to_string(g.novelty) + "," + std::to_string(g.archive.size()) + "\n";
		}
		avgDist /= genomes.size();
		outfile << std::to_string(avgDist) + "," + std::to_string(avgLife) + "," + std::to_string(gPucksCollected1)  + "," + std::to_string(gPucksCollected2) + "," + std::to_string(avgPool) + "," + std::to_string(addedToArchCount) + "\n";
		
		gPucksCollected = 0;
		gPucksCollected1 = 0;
		gPucksCollected2 = 0;
		addedToArchCount = 0;
		poolSizes.clear();

	}

	void SimpleShellsControlArchitecture::logPool(std::vector<Genome> pool, Genome winner) {
		std::ofstream outfile;
		std::string filepath = "nov-experiments/archivesize/poolstats_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		outfile.open(filepath, std::ios_base::app);

		int isWinner = 0;
		for (auto &g : pool) {
			isWinner = (g.id == winner.id) ? 1:0;
			outfile << std::to_string(gWorld->getIterations()) + "," + std::to_string(g.id) + "," + std::to_string(g._distance) + "," + std::to_string(g.pucks.at(0) + g.pucks.at(1)) + "," + std::to_string(g.novelty) + "," + std::to_string(isWinner) + ";";
			isWinner = 0;
		}
		outfile << "\n";
	}

	std::vector<Genome> SimpleShellsControlArchitecture::logDiversity() {
		std::ofstream outfile;
		// std::string filepath = "nov-experiments/samplesize/div_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		// std::string filepath = "nov-experiments/fit/div_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		std::string filepath = "nov-experiments/archivesize/div_" + std::to_string(gMaxCommDistSquared) + "_" + std::to_string(gMaxArchiveSize) + ".csv";
		// std::string filepath = "nov-experiments/samplesize/div_" + std::to_string(gSampleSize) + ".csv";
		outfile.open(filepath, std::ios_base::app);
		std::vector<Genome> genomes;

		for (int i = 0; i < gWorld->getNbOfAgent(); i++) {
			genomes.push_back(((SimpleShellsControlArchitecture *)((gWorld->getAgent(i))->getBehavior()))->_activeGenome);
		}
		assignNovelty(genomes, false, true);
		for (auto &g : genomes) {
			outfile << std::to_string(g.novelty) + ","; 
		}
		outfile << "\n";
		return genomes; 
	}

// *********************************** NOV ***********************************

// manuel: making the move: at some case...
void SimpleShellsControlArchitecture::step() {

	// Print some results
	if (_wm->_agentId == 5 && gWorld->getIterations()%10000 == 0) {
		std::cout << "size " <<_activeGenome.archive.size() << std::endl;			
	}	
	// collectAnalysis();

	// Random robot death
	double randomDeathProb(1.5e-4);
	if (_wm->_action == ACTION_GATHER && Rand::randouble() < randomDeathProb) {
		_wm->_lifetime[PHASE_GATHERING] = 1;
	}

	bool done = false;
	while (!done) {

		// In this order
		// ACTION_ACTIVATE, ACTION_GATHER, ACTION_BRAKE, ACTION_MATE, ACTION_SELECT, ACTION_SIZE
		
		switch (_wm->_action) {
		case ACTION_SELECT: // 5
			select();
			_wm->advance();
			break;
		case ACTION_ACTIVATE : // 1 // this is the intial action of a robot when sim starts
			reset();
			_wm->_phase = PHASE_GATHERING;
			_wm->setRobotLED_colorValues(34, 139, 34);
			_wm->advance();
			break;
		case ACTION_GATHER : // 2
			_age++;
			_activeGenome.age = _age;
			updateActuators();
			_wm->advance();
			// Last turn of gathering can't be combined with the braking turn.
			done = true;
			break;
		case ACTION_BRAKE : // 3
			_wm->_desiredTranslationalValue = 0;
			_wm->_desiredRotationalVelocity = 0;
			_wm->_phase = PHASE_MATING;
			_wm->setRobotLED_colorValues(255, 69, 0);
			_wm->advance();
			break;
		case ACTION_MATE : // 4
			updateGenomes();
			// Last turn of mating can be combined with selection, activation and first turn of gathering
			if (_wm->advance() > 0) done = true;
			break;
		}
	}
	if (gWorld->getIterations() % 10 == 0)
	{
		double dx = _activeGenome._lastXPos - _wm->getXReal();
		double dy = _activeGenome._lastYPos - _wm->getYReal();
		_activeGenome._distance += sqrt(dx*dx + dy*dy);
		_activeGenome._lastXPos = _wm->getXReal();
		_activeGenome._lastYPos = _wm->getYReal();
	}

	#ifdef TRACK_AGENTS
	if (gWorld->getIterations() % 50 == 0)
		// std::cout << "[position] " << gWorld->getIterations() << ' ' << _activeGenome.id << ' ' <<_wm->getXReal() << ' ' << _wm->getYReal() << '\n';
	#endif
}
