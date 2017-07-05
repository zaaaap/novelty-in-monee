/*
 *  roborobo.h
 *  roborobo
 *
 *  Created by Nicolas on 16/12/08.
 *
 */

#ifndef ROBOROBO_H
#define ROBOROBO_H

#include "RoboroboMain/common.h"
#include "Utilities/ExtendedProperties.h"

#include <string.h>
#include <vector>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <iostream>
#include <string>
#include <array>
#include <set>
#include <queue>






class RobotAgent;
typedef boost::shared_ptr<RobotAgent> RobotAgentPtr;

class Puck;
class PuckGen;

// Pre-processor options


//#include "World/World.h"

// project-wide Global data

// *********************************** NOV ***********************************


struct Genome {
    int id;
    std::vector<double> parameters;
    std::vector<int> pucks;
    double fitness;
    double _distance, _lastXPos, _lastYPos;
    // NOV added
    int age;
    unsigned long behaviourId;
    double novelty;
    std::vector< std::vector<int8_t> > simSensorResults;
    std::vector< std::pair<Genome, int> > archive;
};

extern std::vector< std::vector<double> > simSensorInput; // collected samples of sensor input for calc genome behaviour
extern int gFitnessMode; // What fitness function to use: 1 = classic, 2 = novelty, 3 = random
extern int gChangeEnvOn; // In % of iterations. For now only the energyPuck is swapped
extern int gKNNeighbours; // How many neighbours to consider for kNN in novelty search
extern int gPucksCollected;
extern int gPucksCollected1;
extern int gPucksCollected2;
extern int gMaxArchiveSize;
extern bool gWithSync;
extern int gFitChoice;
extern int gNovChoice;
extern int gSelectionFunction;
extern int gMaxRepeatSelection;
extern int gSampleSize;
extern std::vector<int> gWinners;
extern std::vector<int> poolSizes;
extern int addedToArchCount;

// *********************************** NOV ***********************************

#define PAUSE_COMMAND 200 //400 // Pause between each user keyboard command (GUI mode only)

extern ExtendedProperties gProperties;

extern std::string gPropertiesFilename; // default is default.properties
extern std::string gCurrentBuildInfo;  // display through "-v[ersion]" command line option -- check/set value in config.h
extern std::string gStartTime; // a string containing a readable form of the time when the current program was started (useful for naming log files).

extern std::string gCompileTime;
extern std::string gCompileDate;

extern bool gVideoRecording;

// files

extern std::string gLogCommentText; // user comment that will be written in the log file (e.g. description of experimental setup)
extern std::string gLogFilename;
extern std::string gCollisionLogFilename;
extern std::ofstream gLogFile;
extern std::ofstream gCollisionLogFile;

extern std::string gAgentMaskImageFilename;
extern std::string gAgentSpecsImageFilename;
extern std::string gForegroundImageFilename;
extern std::string gEnvironmentImageFilename;
extern std::string gBackgroundImageFilename;
extern std::string gZoneImageFilename;
extern std::string gZoneCaptionPrefixFilename;
extern std::string gConfigFilename;

// general purpose

extern int gRandomSeed;				// random seed. Default value (="-1") means time based.

extern bool gVerbose;				// trigger console verbose on/off
extern bool gBatchMode;				// trigger no display / no verbose mode.
extern int  gDisplayMode;			// [0]: 60fps standard display speed, [1]: ignore some refresh (see gFastDisplayModeSpeed value), [2]: fastest, but no display)
									// much faster without display, in mode [0] (default); display speed == world refresh speed
extern int  gFastDisplayModeSpeed;  // #images ignored in fast display mode (the higher, the faster)
extern bool gInspectAgent;			// display currently selected robot input/output value on console
extern long int gVersion;

// global data

class World;
extern World *gWorld;				// pointer to the World

extern int	gNbOfAgents;			// no. of agents in envt.
extern int	gAgentCounter;			// current number of existing agents
extern int	gAgentIndexFocus;		// focused agent id.


extern bool gRadioNetwork;			// update radio network. note that complexity is : N(N-1)/2 (with N, the nb of agents). *DO* impact speed, use with care.
extern int gMaxRadioDistanceToSquare;		// maximum distance between agents below which there is a network connection. STORED TO THE SQUARE (ie. d^2) while properties file must contain the "real" value (ie. nb of pixels of radius)

extern bool gDisplayZoneCaption;
extern bool gDisplaySensors; // display sensor rays on screen
extern bool gPauseMode;				// pause on/off (default is on)
extern bool gInspectorMode;			// defines if manual control is robot-centered or "god"-centered (ie. envt viewpoint)
extern bool gInspectorAgent;		// defines if there *actually* is an inspector agent (note: if using batch mode for, e.g. evolution, no inspector agent is needed)

extern int gMaxIt;

/*
extern bool gEnergyPoints_alwaysRender; // if active, energy points are re-display at each iterations - this does not impact accuracy, but guarantee display is correct at more or less important computational cost. re-display is only performed if simulation is at normal speed.
extern bool gEnergyMode; // define if agents will use multiple energy points. If it's the case, other following parameters have to be defined, and energy should be managed in a robot observer class
extern int gMaxEnergyPoints; // define the number of energy point present in the map
extern float gEnergyPointRadius; // define the "size" of energy points, that is to say, the distance from which a robot can take the energy from an energy point
extern float gEnergyPointValue; // define the value of an erngy point, i.e. how much energy the robot will earn by taking an energy point
extern int gEnergyPointRespawnLagMaxValue; // used to know how long will it take to an energy point to reapear
extern std::vector<EnergyPoint> gEnergyPoints; //list of the energy points currently active
*/
extern bool gEnergyMode;
extern int gPuckCount;
extern int gPuckColors;
extern std::vector<Puck> gPucks;
extern std::vector<PuckGen> gPuckGens;
extern std::vector< std::vector< int > > gPuckMap;

extern std::vector<Genome> gGenomes;

extern bool gUseDitchSensors;



//Screen attributes
extern int gScreenWidth;
extern int gScreenHeight;
extern int gScreenBPP; // bits per pixel color (32 bits)
extern int gFramesPerSecond; // frame rate

extern bool gAgentRegistration;

//The agent dimensions
extern int gAgentWidth;
extern int gAgentHeight;

extern int gMoveStepWidth; // inspector agent X speed
extern int gMoveStepHeight;// inspector agent Y speed

extern int gMaxSpeedOnXaxis;
extern int gMaxSpeedOnYaxis;

extern int gInspectorAgentXStart;
extern int gInspectorAgentYStart;

// Rapport entre la vitesse de scrolling du parallax et du premier plan
extern int gParallaxFactor;

extern bool gKeyUp, gKeyDown, gKeyRight, gKeyLeft;
extern bool gFullScreenMode;
extern bool gNiceRendering;

extern bool gUserCommandMode; // agent tele-operation (true) or auto-control (false) ?

//The dimensions of the level

extern int gAreaWidth;
extern int gAreaHeight;

// SDL Objects

//camera
extern SDL_Rect gCameraParallax;
extern SDL_Rect gCamera;


//image surfaces
extern SDL_Surface *gScreen;			// surface for rendering
extern SDL_Surface *gForegroundImage;	// "foreground" image, mostly decorative -- intended to be a "nice" representation of egEnvironmentImage
extern SDL_Surface *gEnvironmentImage;	// the real environment: all obstacle assume color #000001, should be close to gForegroundImage
										// note: for the sake of simplicity, one may load the same image into both gForegroundImage and gEnvironmentImage
extern SDL_Surface *gBackgroundImage;	// decorative background image -- see gParallaxFactor for setting parallax (imply background image should be larger than envt and foregd images)
extern SDL_Surface *gAgentMaskImage;	// mask used for collision
extern SDL_Surface *gAgentSpecsImage;	// image contains specs wrt. proximity sensor (location and orientation)


//stories
// note: key stories (first story stories of each decimal) must be visited prior to story in the decimal to be displayed (eg. story11 is displayed only if story10 has been displayed -- exception if story10 does not exist)
// lexicon: key stories (e.g. 0, 10, 20, 100) and sub-stories (e.g. 11, 144, 1, 255)
extern SDL_Surface *gZoneImage;  // gZoneImage , gZoneCaptions, gZoneStatus
extern SDL_Surface *gZoneCaptionImage[256]; // all zone captions. total depends on max color component value (ie. defined in 0...255)
extern bool gZoneStatus[256]; // told stories

// general purpose
extern SDL_Event gEvent; //The event structure

// agent dynamics and structure
extern float gMaxTranslationalSpeed; // pixels/sec.
extern float gMaxTranslationalDeltaValue; // default: = gMaxTranslationSpeed (however, may favor bang-bang)
extern float gMaxRotationalSpeed; // degrees/sec.
extern int gSensorRange;

extern int gLocomotionMode; // locomotion mode dynamics
							// =0: agent crash but can still rotate - ie. rotate while crash (no obstacle friction)
							// =1: agent crash and can rotate only if not crashing again - ie. rotate if no translation twd obstacle (full obstacle friction)

// Organism related globals
extern bool gUseOrganisms; // agents can aggregate into organisms, default := false
                                              // if set to true agents that are part of an organism are moved at the same time
extern bool gRiverbed; // use the implementation for a stream in organisms
extern int gStreamDirection; // direction of the stream in degrees
extern double gStreamStrength; // strength of the stream
extern int gRiverbedZone; // id of riverbed zone
extern std::string gAgentPositiveMaskImageFilename;
extern std::string gAgentNegativeMaskImageFilename;
extern std::string gOrganismLogFilename;
extern std::string gOrganismSizesLogFilename;
extern SDL_Surface *gAgentPositiveMaskImage;	// mask used for showing that an agent wants to connect
extern SDL_Surface *gAgentNegativeMaskImage;	// mask used for showing that an agent wants to connect
extern int gConnectionGap; // the gap allowed between agents in pixels
extern int gOrganismSampleFrequency; // the frequency of logging samples in iterations. Higher means less frequent

extern int gSpawnWidth; // Maximum x coordinate where a robot can be placed initially
extern int gSpawnHeight; // Maximum y coordinate where a robot can be placed initially

extern bool gRenderRobotId;
extern std::string gFontName;



extern Uint32 gNextGenomeId;


extern Uint16 gMaxGenePool;

extern int gApproximateDiameterSquared;
extern std::vector<int> gDistSquared;
extern int gMaxCommDistSquared;


extern std::vector<bool> gRadioNetworkArray; // a 2D adjacent matrix stating network connections and signal strength. By convention, 0 means no signal.



extern Uint32 G_COLOR_WHITE;




inline bool isRadioConnection(int & idOne, int & idTwo) {
    int index = (idOne < idTwo) ? idOne * gAgentCounter + idTwo : idTwo * gAgentCounter + idOne;
    return gRadioNetworkArray[index];
}

inline bool isCollisionPotential(int & idOne, int & idTwo) {
	int index = (idOne < idTwo) ? idOne * gAgentCounter + idTwo : idTwo * gAgentCounter + idOne;
    return gDistSquared[index] < gApproximateDiameterSquared * 9;
}



// function hearders

void clean_up();
void displayHelp();
bool checkQuitEvent();
bool handleKeyEvent(Uint8 *keyboardStates);
void updateDisplay();
void initLogging();
void stopLogging();
void initMonitor( bool __isInspectorAgent = false );
void updateMonitor(Uint8* __keyboardStates);
bool loadProperties( std::string __propertiesFilename );

void initRoborobo();
bool runRoborobo( int __maxIt = -1);
void closeRoborobo();
void resetRoboroboWorld();

int launchRoborobo(); // the main entry point


#endif // ROBOROBO_H

