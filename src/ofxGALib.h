#pragma once
#include "ofMain.h"
#include <functional>

class GABin2DecGenome;
class GASimpleGA;

struct RangeInfo 
{
	RangeInfo() {};
    RangeInfo(int size, float mn, float mx) : 
        mSize(size), mMin(mn), mMax(mx)
    {};
    int mSize;
    float mMin;
    float mMax;
};

template <typename R, typename P > 
class GenericFunctor
{
public:
    virtual R call(P param) = 0;        
    virtual void* getObjPtr() = 0;
};

template <typename R, typename T, typename P> 
class TGenericFunctor : public GenericFunctor<R, P>
{
private:
    T* m_pObj;
    R (T::*m_pFunc)(P param);

public:
    TGenericFunctor(T* pObj, R(T::*pFunc)(P)) { m_pObj = pObj;  m_pFunc = pFunc; };
    virtual R call(P param) { return (*m_pObj.*m_pFunc)(param);};             
    virtual void* getObjPtr() { return (void*)m_pObj; };
    bool operator==(const TGenericFunctor& other) const { return m_pObj == other.m_pObj && m_pFunc == other.m_pFunc; };
};

class ofxGALib
{
public:
    ofxGALib();
	virtual ~ofxGALib();
	void setup(const vector<RangeInfo>& ranges, int repeat, int popsize = 25, int ngen = 50, float pmut = 0.05, float pcross = 0.4);
    float run(int times);
	bool done();

    vector<float> mInter;
    vector<float> mOut;

	void init();

	GenericFunctor<float, const vector<float>&>* mFunc;

    template <typename T>
    void setFitness(T* pObj, float(T::*pFunc)(const vector<float>&)) { if (mFunc) delete mFunc; mFunc = new TGenericFunctor<float, T, const vector<float>&>(pObj, pFunc); }

    float evaluate(const vector<float>& values);

    GABin2DecGenome* mGenome;
	GASimpleGA* ga;
    int mRepeat;
	bool started;
};
