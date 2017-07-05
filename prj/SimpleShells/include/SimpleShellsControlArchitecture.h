#ifndef SIMPLESHELLSCONTROLARCHITECTURE_H
#define SIMPLESHELLSCONTROLARCHITECTURE_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"

#include "BehaviorControlArchitectures/BehaviorControlArchitecture.h"
#include "SimpleShells/include/SimpleShellsAgentWorldModel.h"

#include "World/World.h"
#include "Agents/RobotAgent.h"

#include "Utilities/Misc.h"
#include "Utils/Rand.h"

#include <cmath>
#include <functional> // std::hash
#include <cstdio>
#include <numeric>




class RangeSensor;

class SimpleShellsControlArchitecture : public BehaviorControlArchitecture {
        // manuel: here is where the each robots saves all its stats
        // use getBehavior to get this, then access as usual

        SimpleShellsAgentWorldModel* _wm;
        int _maxLifetime[2];
	int _age;
        int _hiddenNeuronCount;
        int _parameterCount;
        size_t _tournamentSize; /// << tournamentSize smaller than 2 implies rank-based roulettewheel selection
        
        bool _useMarket;
        bool _useSpecBonus;
        bool _randomSelection;
	double _task1Premium;
        double _selectionPressure;
        double _archiveThreshold;
        double _fitnessThreshold;
        int _lastAddedToArchive;

        std::vector<double> _response;

        void assignFitness(std::vector<Genome> & genomes);
        void gatherGenomes(std::vector<Genome*> & genePool, int commDistSquared);

        std::vector<Genome*> _nearbyGenomes;

        Genome _activeGenome;

        Genome selectWeighted(std::vector<Genome> & genomes, bool onFitness);

	void logStats();

	public:
		SimpleShellsControlArchitecture(RobotAgentWorldModel *__wm);
		~SimpleShellsControlArchitecture();
		void step();
        	void reset();
                void prepareShutdown();
		void updateActuators();
		void updateGenomes();
                void select();
                void selectClassic();
                // void selectFRankendN();
                // void selectNRankendF();
                void selectFitness();
                void selectNovelty();
        	void selectRandom();
                void selectProChoice();
                void selectMCNS();
                void selectMCNS2();
		Genome* getGenome() {return &_activeGenome;}
                void copyGenome(Genome &g);

        private:
                // NOV
                void calcBehaviourAndId();
                int8_t binarize(double d);
                double calcBehaviorDistance(Genome &g1, Genome &g2);
                double hammingDistance(std::vector<int8_t> &v1, std::vector<int8_t> &v2);
                void changeEnv();
                void changeEnv2();
                void assignNovelty(std::vector<Genome> &genomes, bool addToArchive, bool logging);
                void evaluateIndividual(std::vector<Genome> &genomes, Genome &ind, bool addToArchive, bool logging);
                double noveltyAvgNN(std::vector<Genome> &genomes, Genome &ind, bool addToArchive, bool logging);
                std::vector<double> mapNovelty(std::vector<Genome> &genomes, Genome &ind);
                std::vector<double> mapNoveltyArchive(Genome &ind);
                std::vector<double> mapNoveltyPop(std::vector<Genome> &genomes, Genome &ind);
                int calcPucksCollected(const Genome &g);
                bool aboveArchiveThreshold(double result);
                void addGenomeToArchive(Genome ind);
                bool genomeInArchive(const Genome &ind);
                void syncArchives(Genome &winner);
                Genome chooseParent(bool onFitness);
                Genome chooseParentRoulette(bool onFitness);
                Genome chooseParentTournament(bool onFitness);
                void printAvgFitness();
                void printAvgAge();
                void printPucksPer500();
                void printPucksPer1000();
                void printChoicesPer500();
                bool aboveFitnessThreshold(Genome &winner);
                void collectAnalysis();
                std::vector<Genome> logDiversity();
                void logAvgsAndPucks(std::vector<Genome> &genomes);
                void logWinners();    
                void logPool(std::vector<Genome> pool, Genome winner);            

};


#endif

