#include "ofxGALib.h"

#include <ga/GASimpleGA.h>
#include <ga/GABin2DecGenome.h>
#include <ga/std_stream.h>

float Objective(GAGenome& g);

void ofxGALib::setup(const vector<RangeInfo>& ranges, int repeat, int popsize, int ngen, float pmut, float pcross)
{
    mRepeat = repeat;

    // Create a phenotype then fill it with the phenotypes we will need to map to
    // the values we read from the file.  The arguments to the add() method of a
    // Bin2Dec phenotype are (1) number of bits, (2) min value, and (3) max value.
    // The phenotype maps a floating-point number onto the number of bits that
    // you designate.  Here we just make everything use 8 bits and use the max and
    // min that were used to generate the target values.  You can experiment with
    // the number of bits and max/min values in order to make the GA work better
    // or worse.

    GABin2DecPhenotype map;

    for (int i = 0; i < ranges.size() * mRepeat; ++i) {
        const RangeInfo& ri = ranges[i % ranges.size()];
        map.add(ri.mSize, ri.mMin, ri.mMax);
    }

    mInter.resize(ranges.size() * mRepeat);
    mOut.resize(ranges.size() * mRepeat);

    // Create the template genome using the phenotype map we just made.  The
    // GA will use this genome to clone the population that it uses to do the
    // evolution.  We pass the objective function to create the genome.  We 
    // also use the user data function in the genome to keep track of our
    // target values.

    if (mGenome) delete mGenome;

    mGenome = new GABin2DecGenome(map, Objective, (void *)this);

    // Generate a sequence of random numbers using the values in the min and max
    // arrays.  We also set one of them to integer value to show how you can get
    // explicit integer representations by choosing your number of bits
    // appropriately.

	unsigned int seed = ofGetUnixTime();

	GARandomSeed(seed);

	// Now create the GA using the genome, set the parameters, and run it.
	if (ga) delete ga;
	ga = new GASimpleGA(*mGenome);
	ga->populationSize(popsize);
	ga->nGenerations(ngen);
	ga->pMutation(pmut);
	ga->pCrossover(pcross);
	ga->scoreFilename("bog.dat");
	ga->flushFrequency(0);	// dump scores to disk every 50th generation
	ga->initialize(seed); 

	started = true;
}

float ofxGALib::evaluate( const vector<float>& values )
{
    cout << "No" << endl;
    return 0;   
}

ofxGALib::ofxGALib() : mFunc(0), mGenome(0), ga(0)
{
    setFitness(this, &ofxGALib::evaluate);
	started = false;
}

ofxGALib::~ofxGALib()
{
	if (mGenome) delete mGenome;
	if (ga) delete ga;
}

float ofxGALib::run(int times)
{
	for (int i = 0; i < times && !ga->done(); ++i) {
		ga->step();
	} 

    *mGenome = ga->statistics().bestIndividual();
    for(int i = 0; i < mGenome->nPhenotypes(); i++){
        mOut[i] = mGenome->phenotype(i);
    }
    return mGenome->score();
}

bool ofxGALib::done()
{
	return ga->done();
}

float Objective(GAGenome& g)
{
    GABin2DecGenome& genome = (GABin2DecGenome &)g;
    ofxGALib* gaLib = (ofxGALib*)g.userData();

    assert(genome.nPhenotypes() == gaLib->mInter.size());

    for(int i = 0; i < genome.nPhenotypes(); i++)
        gaLib->mInter[i] = genome.phenotype(i);

    return gaLib->mFunc->call(gaLib->mInter);
}
