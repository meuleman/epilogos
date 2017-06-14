#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <utility> // for pair()
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <string>

using namespace std;

enum measurementType {KL, KLs, KLss}; // corresponding to measurements using KL, KL*, or KL**

bool PairFloatInt_LT(const pair<float,int>& a, const pair<float,int>& b);
bool PairFloatInt_LT(const pair<float,int>& a, const pair<float,int>& b)
{
  return a.first < b.first;
}

class Model {
public:
  virtual bool init(const char *pObsFname, const char *pQcatFname, const char *pNullFname,
		    const string& chrom, const int& firstBegPos, const int& regionWidth) = 0;
  virtual int  size(void) const = 0;
  virtual bool getQcontrib(ifstream& infile, const char *pFilename, const int& Nsites) = 0;
  virtual bool processInputValue(const int& val) = 0;
  virtual void computeAndWriteMetric(void) = 0;
};

class KLModel : public Model {
public:
  KLModel() {};
  bool init(const char *pObsFname, const char *pQcatFname, const char *pNullFname,
	    const string& chrom, const int& firstBegPos, const int& regionWidth);
  int  size(void) const { return m_size; }
  bool getQcontrib(ifstream& infile, const char *pFilename, const int& Nsites);
  bool processInputValue(const int& val);
  void computeAndWriteMetric(void);
protected:
  int m_numStates;
  int m_size; // number of values required on each line of input
  int m_group1size, m_group2size;
  int m_numValsProcessedForGroup1, m_numValsProcessedForGroup2;
  bool m_writeNullMetric;
  ofstream m_ofsObs, m_ofsNullValues, m_ofsQcat;
  string m_chrom;
  int m_curBegPos, m_curEndPos, m_regionWidth;
private:
  KLModel(const KLModel&); // we have no need for a copy constructor, so disable it
  vector<int> m_P1numerators, m_P2numerators; 
  vector<float> m_Q1contrib, m_Q2contrib;
  vector<float> m_logsOfObservationTallies;
};

class KLsModel : public KLModel {
public:
  KLsModel() {};
  bool getQcontrib(ifstream& infile, const char *pFilename, const int& Nsites);
  bool processInputValue(const int& val);
  void computeAndWriteMetric(void);
private:
  KLsModel(const KLsModel&); // we have no need for a copy constructor, so disable it
  vector<int> m_Ps1numerators, m_Ps2numerators;
  vector<float> m_Qs1contrib, m_Qs2contrib;
  vector<float> m_logsOfObservationTallies;
  vector<pair<int, int> > m_unorderedStatePairDecompositions;
};

class KLssModel : public KLModel {
public:
  KLssModel() {};
  bool getQcontrib(ifstream& infile, const char *pFilename, const int& Nsites);
  bool processInputValue(const int& val);
  void computeAndWriteMetric(void);
private:
  KLssModel(const KLssModel&); // we have no need for a copy constructor, so disable it
  map<int, map<int,float> > m_Qss1contrib, m_Qss2contrib;
  map<int, pair<map<int, set<int> >, map<int, set<int> > > > m_statePairGroupObservationsAtThisSite;
};


bool KLModel::init(const char *pObsFname, const char *pQcatFname, const char *pNullFname,
		   const string& chrom, const int& firstBegPos, const int& regionWidth)
{
  if (pObsFname != NULL)
    {
      m_ofsObs.open(pObsFname);
      if (!m_ofsObs)
	{
	  cerr << "Error:  Unable to open file \"" << pObsFname << "\" for writing." << endl << endl;
	  return false;
	}
    }
  if (pQcatFname != NULL)
    {
      m_ofsQcat.open(pQcatFname);
      if (!m_ofsQcat)
	{
	  cerr << "Error:  Unable to open file \"" << pQcatFname << "\" for writing." << endl << endl;
	  return false;
	}
    }
    if (pNullFname != NULL)
    {
      m_ofsNullValues.open(pNullFname);
      if (!m_ofsNullValues)
	{
	  cerr << "Error:  Unable to open file \"" << pNullFname << "\" for writing." << endl << endl;
	  return false;
	}
    }
  m_chrom = chrom;
  m_curBegPos = firstBegPos;
  m_regionWidth = regionWidth;
  m_curEndPos = m_curBegPos + m_regionWidth;
  m_numValsProcessedForGroup1 = m_numValsProcessedForGroup2 = 0;
  m_group1size = m_group2size = 0;
  m_writeNullMetric = m_ofsNullValues.is_open();
  m_numStates = m_size = 0; // these will be set by getQcontrib()
  return true;
}

bool KLModel::getQcontrib(ifstream& infile, const char *pFilename, const int& Nsites)
{
  const int BUFSIZE(100000);
  char buf[BUFSIZE], *p;
  const float LOG_Nsites(log(static_cast<float>(Nsites)));
  unsigned long thisNumTallies, totalNumTallies(0);

  if (!infile.getline(buf, BUFSIZE))
    {
      cerr << "Error:  File " << pFilename << " is empty." << endl << endl;
      return false;
    }
  p = strtok(buf, "\t");
  thisNumTallies = atol(p);
  totalNumTallies += thisNumTallies;
  if (0 == m_group1size)
    m_Q1contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
  else
    m_Q2contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
  while (p = strtok(NULL, "\t"))
    {
      thisNumTallies = atol(p);
      totalNumTallies += thisNumTallies;
      if (0 == m_group1size)
	m_Q1contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
      else
	m_Q2contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
    }
  // infile should only contain 1 line of data
  if (infile.getline(buf, BUFSIZE))
    {
      cerr << "Error:  File " << pFilename << " contains multiple lines of data; "
	   << "it should contain a single line of tab-delimited state-pair tallies."
	   << endl << endl;
      return false;
    }

  // The number of unique (unordered) state pairs is numStates*(numStates + 1)/2
  // (yes, +1, not -1).  Thus if x is the number of elements in the Q* vector,
  // numStates satisfies numStates^2 + numStates - 2*x == 0.
  if (0 == m_group1size)
    m_numStates = m_Q1contrib.size();
  else
    {
      // Perform a sanity check for safety's sake.
      int thisNumStates = m_Q2contrib.size();
      if (thisNumStates != m_numStates)
	{
	  cerr << "Error:  The file containing tallies for Q for group 1 implies there are "
	       << m_numStates << " possible states,\n"
	       << "but file " << pFilename << " (containing Q for group 2) implies there are "
	       << thisNumStates << " possible states." << endl << endl;
	  return false;
	}
    }

  // The sum of the tallies in this Q* equals the total number of sites
  // times the number of epigenomes
  if (0 == m_group1size)
    {
      m_group1size = static_cast<int>(floor(static_cast<float>(totalNumTallies)/static_cast<float>(Nsites) + 0.01));
      m_P1numerators.assign(m_Q1contrib.size(), 0);
    }
  else
    {
      m_group2size = static_cast<int>(floor(static_cast<float>(totalNumTallies)/static_cast<float>(Nsites) + 0.01));
      m_P2numerators.assign(m_Q2contrib.size(), 0);
    }
  // Note:  The factor of m_group1size in each element of Q
  // is not stored in m_Qcontrib, because each term of Q enters into the metric
  // only via the ratio P/Q, and this factor cancels out within this ratio.

  // At any site, the number of times any state is observed
  // is a number between 0 and numEpigenomes, inclusive.
  // At each site, we'll need the log2 of one or more of these tallies.
  // We compute them here to avoid needlessly computing logs of the same numbers thousands or millions of times.
  if (0 == m_group2size)
    {
      m_logsOfObservationTallies.push_back(0); // unused
      for (int i = 1; i <= m_group1size; i++)
	m_logsOfObservationTallies.push_back(log(static_cast<float>(i)));
    }
  else
    {
      if (m_group2size > m_logsOfObservationTallies.size() - 1)
	for (int i = m_logsOfObservationTallies.size(); i <= m_group2size; i++)
	  m_logsOfObservationTallies.push_back(log(static_cast<float>(i)));
    }

  m_size = m_P1numerators.size() + m_P2numerators.size();
  return true;
}

bool KLModel::processInputValue(const int& thisTally)
{
  bool processingGroup1(true);

  // Check these variables, in case excess columns appear in this line of input.
  // This is also how we determine, in the case of two groups of epigenomes being compared,
  // whether the input value is for group 1 or group 2.
  if (m_numValsProcessedForGroup1 == m_numStates)
    {
      if (0 == m_group2size)
	{
	  cerr << "Error:  Found excess columns in a line of input; expected "
	       << m_numStates << "." << endl;
	  return false;
	}
      else
	{
	  if (m_numValsProcessedForGroup2 == m_numStates)
	    {
	      cerr << "Error:  Found excess columns in a line of input; expected "
		   << m_numStates * 2 << "." << endl;
	      return false;
	    }
	  else
	    processingGroup1 = false;
	}
    }

  if (processingGroup1)
    m_P1numerators[m_numValsProcessedForGroup1++] = thisTally;
  else
    m_P2numerators[m_numValsProcessedForGroup2++] = thisTally;

  return true;
}

void KLModel::computeAndWriteMetric(void)
{
  static const float LOG2(0.6931471806);
  static const float denom1 = LOG2 * static_cast<float>(m_group1size);
  static const float denom2 = LOG2 * static_cast<float>(m_group2size);
  static vector<pair<float, int> > contribOfEachState(m_numStates, make_pair(0,0));
  float retVal(0);
  
  if (!m_writeNullMetric)
    {
      // Initialize the vector of per-state contributions. This is needed for the "qcat" file.
      for (int i = 0; i < m_numStates; i++)
	{
	  contribOfEachState[i].first = 0;
	  contribOfEachState[i].second = i+1;
	}
    }

  for (int i = 0; i < m_P1numerators.size(); i++)
    {
      float term(0);
      if (m_P1numerators[i] != 0)
	{
	  if (m_Q1contrib[i] < -999.)
	    term = m_Q1contrib[i];
	  else
	    term += (static_cast<float>(m_P1numerators[i])/denom1) *
	      (m_logsOfObservationTallies[m_P1numerators[i]] + m_Q1contrib[i]);
	}
      if (m_group2size != 0)
	{
	  if (m_P2numerators[i] != 0)
	    {
	      if (m_Q2contrib[i] < -999.)
		term = -m_Q2contrib[i];
	      else
		term -= (static_cast<float>(m_P2numerators[i])/denom2) *
		  (m_logsOfObservationTallies[m_P2numerators[i]] + m_Q2contrib[i]);
	    }
	}
      if (!m_writeNullMetric) // no need to break down the metric by state if we're solely tasked with writing null metric values
	contribOfEachState[i].first = term;

      retVal += (0 == m_group2size ? term : fabs(term));
    }

  if (!m_writeNullMetric)
    {
      float contribOfMaxTerm;
      int stateWithMaxTerm_1based;
      bool signOfMaxTerm;
      char formattedQcatFloat[10];
      sort(contribOfEachState.begin(), contribOfEachState.end(), PairFloatInt_LT);
      if (fabs(contribOfEachState.front().first) > contribOfEachState.back().first)
	{
	  contribOfMaxTerm = contribOfEachState.front().first;
	  stateWithMaxTerm_1based = contribOfEachState.front().second;
	}
      else
	{
	  contribOfMaxTerm = contribOfEachState.back().first;
	  stateWithMaxTerm_1based = contribOfEachState.back().second;
	}

      signOfMaxTerm = contribOfMaxTerm > 0 ? true : false;
      contribOfMaxTerm = fabs(contribOfMaxTerm);
      m_ofsObs << m_chrom << '\t' << m_curBegPos << '\t' << m_curEndPos << '\t'
	       << stateWithMaxTerm_1based << '\t'
	       << contribOfMaxTerm << '\t'
	       << (signOfMaxTerm ? "1" : "-1") << '\t'
	       << retVal << endl;
      sprintf(formattedQcatFloat, "%.4g", contribOfEachState[0].first);
      m_ofsQcat << m_chrom << '\t' << m_curBegPos << '\t' << m_curEndPos
		<< "\t[ [" << formattedQcatFloat << ',' << contribOfEachState[0].second << ']';
      for (int i = 1; i < contribOfEachState.size(); i++)
	{
	  sprintf(formattedQcatFloat, "%.4g", contribOfEachState[i].first);	  
	  m_ofsQcat << ", [" << formattedQcatFloat << ',' << contribOfEachState[i].second << ']';
	}
      m_ofsQcat << " ]" << endl;
      m_curBegPos = m_curEndPos;
      m_curEndPos += m_regionWidth;
    }
  else
    m_ofsNullValues << retVal << endl;
  
  // reset the counting variables and the "P* numerator" (m_Ps1numerators, m_Ps2numerators) tallies
  m_numValsProcessedForGroup1 = m_numValsProcessedForGroup2 = 0;
  m_P1numerators.assign(m_P1numerators.size(), 0);
  if (!m_P2numerators.empty())
    m_P2numerators.assign(m_P2numerators.size(), 0);
}


bool KLsModel::getQcontrib(ifstream& infile, const char *pFilename, const int& Nsites)
{
  const int BUFSIZE(100000);
  char buf[BUFSIZE], *p;
  const float LOG_Nsites(log(static_cast<float>(Nsites)));
  unsigned long thisNumTallies, totalNumTallies(0);

  if (!infile.getline(buf, BUFSIZE))
    {
      cerr << "Error:  File " << pFilename << " is empty." << endl << endl;
      return false;
    }
  p = strtok(buf, "\t");
  thisNumTallies = atol(p);
  totalNumTallies += thisNumTallies;
  if (0 == m_group1size)
    m_Qs1contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
  else
    m_Qs2contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
  while (p = strtok(NULL, "\t"))
    {
      thisNumTallies = atol(p);
      totalNumTallies += thisNumTallies;
      if (0 == m_group1size)
	m_Qs1contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
      else
	m_Qs2contrib.push_back(0 == thisNumTallies ? -999999. : LOG_Nsites - log(static_cast<float>(thisNumTallies)));
    }
  // infile should only contain 1 line of data
  if (infile.getline(buf, BUFSIZE))
    {
      cerr << "Error:  File " << pFilename << " contains multiple lines of data; "
	   << "it should contain a single line of tab-delimited state-pair tallies."
	   << endl << endl;
      return false;
    }

  // The number of unique (unordered) state pairs is numStates*(numStates + 1)/2
  // (yes, +1, not -1).  Thus if x is the number of elements in the Q* vector,
  // numStates satisfies numStates^2 + numStates - 2*x == 0.
  if (0 == m_group1size)
    m_numStates = static_cast<int>(floor((sqrt(1. + 8.*static_cast<float>(m_Qs1contrib.size())) - 1.)/2. + 0.01));
  else
    {
      // Perform a sanity check for safety's sake.
      int thisNumStates = static_cast<int>(floor((sqrt(1. + 8.*static_cast<float>(m_Qs2contrib.size())) - 1.)/2. + 0.01));
      if (thisNumStates != m_numStates)
	{
	  cerr << "Error:  The file containing tallies for Q* for group 1 implies there are "
	       << m_numStates << " possible states,\n"
	       << "but file " << pFilename << " (containing Q* for group 2) implies there are "
	       << thisNumStates << " possible states." << endl << endl;
	  return false;
	}
    }

  // The sum of the tallies in this Q* equals the total number of sites
  // times the number of unique epigenome pairs (the latter of which equals numEpigenomes*(numEpigenomes-1)/2).
  if (0 == m_group1size)
    {
      m_group1size = static_cast<int>(floor((sqrt(1. + 8.*static_cast<float>(totalNumTallies)/static_cast<float>(Nsites)) + 1.)/2. + 0.01));
      m_Ps1numerators.assign(m_Qs1contrib.size(), 0);
    }
  else
    {
      m_group2size = static_cast<int>(floor((sqrt(1. + 8.*static_cast<float>(totalNumTallies)/static_cast<float>(Nsites)) + 1.)/2. + 0.01));
      m_Ps2numerators.assign(m_Qs2contrib.size(), 0);
    }
  // Note:  The factor of numEpiPairs*(numEpiPairs-1)/2 in each nonzero element of Q*
  // is not stored in m_QsContrib, because each term of Q* enters into the metric
  // only via the ratio P*/Q*, and this factor cancels out within this ratio.

  // At any site, the number of times any unique state pair is observed
  // is a number between 0 and numEpigenomes*(numEpigenomes-1)/2, inclusive.
  // At each site, we'll need the log2 of one or more of these tallies.
  // We compute them here to avoid needlessly computing logs of the same numbers thousands or millions of times.
  if (0 == m_group2size)
    {
      m_logsOfObservationTallies.push_back(0); // unused
      for (int i = 1; i <= m_group1size*(m_group1size-1)/2; i++)
	m_logsOfObservationTallies.push_back(log(static_cast<float>(i)));
    }
  else
    {
      if (m_group2size*(m_group2size-1)/2 > m_logsOfObservationTallies.size() - 1)
	for (int i = m_logsOfObservationTallies.size(); i <= m_group2size*(m_group2size-1)/2; i++)
	  m_logsOfObservationTallies.push_back(log(static_cast<float>(i)));
    }

  // Lastly, it's helpful to have a lookup table to map from unique unordered state pairs 0, 1, 2, ...
  // to the states that form those state pairs (i.e., 0 --> (1,1), 1 --> (1,2), ...,
  // and 119 --> (15,15) if m_numStates = 15).
  // Imagine assigning the unique unordered state pairs to the upper triangle (including the diagonal)
  // of a matrix of m_numStates*m_numStates values, and imagine starting from the lower right corner element
  // (maxUniqueStatePairID, which maps to (m_numStates,m_numStates)) and subtracting increments from it
  // (delta = 0, 1, 2, ...) and deriving from those deltas the amount that needs to be subtracted
  // from each row index of the matrix (delta_row) and each column index of the matrix (delta_column)
  // to transform the state pair (m_numStates,m_numStates) into
  // (row index, column index) = (state1, state2).
  int maxUniqueStatePairID(m_numStates*(m_numStates+1)/2 - 1);
  m_unorderedStatePairDecompositions.assign(maxUniqueStatePairID+1, make_pair(0,0));
  for (int delta = 0; delta <= maxUniqueStatePairID; delta++)
    {
      int delta_row, delta_column;
      delta_row = static_cast<int>(floor((-1. + sqrt(1. + 8.*static_cast<float>(delta)))/2. + 0.01));
      if (delta <= 2)
	{
	  // Edge cases for delta = 0, 1, 2
	  if (2 == delta)
	    delta_column = 1;
	  else
	    delta_column = 0;
	}
      else
	delta_column = delta - delta_row*(delta_row + 1)/2;
      m_unorderedStatePairDecompositions[maxUniqueStatePairID - delta].first = m_numStates - delta_row;
      m_unorderedStatePairDecompositions[maxUniqueStatePairID - delta].second = m_numStates - delta_column;
    }
  
  m_size = m_Ps1numerators.size() + m_Ps2numerators.size();
  return true;
}

bool KLsModel::processInputValue(const int& thisTally)
{
  bool processingGroup1(true);

  // Check these variables, in case excess columns appear in this line of input.
  // This is also how we determine, in the case of two groups of epigenomes being compared,
  // whether the input value is for group 1 or group 2.
  if (m_numValsProcessedForGroup1 == m_Ps1numerators.size())
    {
      if (m_numValsProcessedForGroup2 == m_Ps2numerators.size())
	{
	  cerr << "Error:  Found excess columns in a line of input; expected "
	       << m_group1size + m_group2size << "." << endl;
	  return false;
	}
      else
	processingGroup1 = false;
    }

  if (processingGroup1)
    m_Ps1numerators[m_numValsProcessedForGroup1++] = thisTally;
  else
    m_Ps2numerators[m_numValsProcessedForGroup2++] = thisTally;

  return true;
}

void KLsModel::computeAndWriteMetric(void)
{
  static const float LOG2(0.6931471806);
  static const float denom1 = LOG2 * static_cast<float>(m_group1size)*static_cast<float>(m_group1size - 1)/2.;
  static const float denom2 = LOG2 * static_cast<float>(m_group2size)*static_cast<float>(m_group2size - 1)/2.;
  static vector<pair<float, int> > contribOfEachState(m_numStates, make_pair(0,0));
  float retVal(0), contribOfMaxStatePairTerm(0), contribOfMaxStateTerm(0);
  int statePairWithMaxTerm_1based, stateWithMaxTerm_1based;
  bool signOfMaxStatePairTerm, signOfMaxStateTerm;

  if (!m_writeNullMetric)
    {
      // Initialize the vector of per-state contributions. This is needed for the "qcat" file.
      for (int i = 0; i < m_numStates; i++)
	{
	  contribOfEachState[i].first = 0;
	  contribOfEachState[i].second = i+1;
	}
    }

  for (int uniqueStatePairID = 0; uniqueStatePairID < m_Ps1numerators.size(); uniqueStatePairID++)
    {
      float term = 0;
      float absTerm; // |term|
      int stateOfEpi1_1based, stateOfEpi2_1based;

      if (!m_writeNullMetric) // no need to break down the metric by state if we're solely tasked with writing null metric values
	{
	  stateOfEpi1_1based = m_unorderedStatePairDecompositions[uniqueStatePairID].first;
	  stateOfEpi2_1based = m_unorderedStatePairDecompositions[uniqueStatePairID].second;
	}
      if (m_Ps1numerators[uniqueStatePairID] != 0)
	{
	  if (m_Qs1contrib[uniqueStatePairID] < -999.)
	    term = m_Qs1contrib[uniqueStatePairID];
	  else
	    term += (static_cast<float>(m_Ps1numerators[uniqueStatePairID])/denom1) *
	      (m_logsOfObservationTallies[m_Ps1numerators[uniqueStatePairID]] + m_Qs1contrib[uniqueStatePairID]);
	}
      if (m_group2size != 0)
	{
	  if (m_Ps2numerators[uniqueStatePairID] != 0)
	    {
	      if (m_Qs2contrib[uniqueStatePairID] < -999.)
		term = -m_Qs2contrib[uniqueStatePairID];
	      else
		term -= (static_cast<float>(m_Ps2numerators[uniqueStatePairID])/denom2) *
		  (m_logsOfObservationTallies[m_Ps2numerators[uniqueStatePairID]] + m_Qs2contrib[uniqueStatePairID]);
	    }
	}
      absTerm = fabs(term);

      if (!m_writeNullMetric) // no need to break down the metric by state if we're solely tasked with writing null metric values
	{
	  if (absTerm > fabs(contribOfMaxStatePairTerm))
	    {
	      contribOfMaxStatePairTerm = term;
	      statePairWithMaxTerm_1based = uniqueStatePairID + 1;
	    }
	  contribOfEachState[stateOfEpi1_1based - 1].first += 0.5*term;
	  contribOfEachState[stateOfEpi2_1based - 1].first += 0.5*term;
	}

      retVal += (0 == m_group2size ? term : absTerm);
    }

  if (!m_writeNullMetric)
    {
      // extract the states (s1,s2) from the unique unordered state pair
      // that contributed the most to the metric
      int s1 = m_unorderedStatePairDecompositions[statePairWithMaxTerm_1based - 1].first,
	s2 = m_unorderedStatePairDecompositions[statePairWithMaxTerm_1based - 1].second;
      char formattedQcatFloat[10];
      sort(contribOfEachState.begin(), contribOfEachState.end(), PairFloatInt_LT);
      if (fabs(contribOfEachState.front().first) > contribOfEachState.back().first)
	{
	  contribOfMaxStateTerm = -1*contribOfEachState.front().first;
	  stateWithMaxTerm_1based = contribOfEachState.front().second;
	  signOfMaxStateTerm = false;
	}
      else
	{
	  contribOfMaxStateTerm = contribOfEachState.back().first;
	  stateWithMaxTerm_1based = contribOfEachState.back().second;
	  signOfMaxStateTerm = true;
	}
      signOfMaxStatePairTerm = contribOfMaxStatePairTerm > 0 ? true : false;      
      contribOfMaxStatePairTerm = fabs(contribOfMaxStatePairTerm);
      m_ofsObs << m_chrom << '\t' << m_curBegPos << '\t' << m_curEndPos << '\t'
	       << stateWithMaxTerm_1based << '\t' << contribOfMaxStateTerm << '\t'
	       << (signOfMaxStateTerm ? "1" : "-1") << '\t'
	       << '(' << s1 << ',' << s2 << ')' << '\t' // << statePairWithMaxTerm_1based << '\t'
	       << contribOfMaxStatePairTerm << '\t'
	       << (signOfMaxStatePairTerm ? "1" : "-1") << '\t'
	       << retVal << endl;
      sprintf(formattedQcatFloat, "%.4g", contribOfEachState[0].first);
      m_ofsQcat << m_chrom << '\t' << m_curBegPos << '\t' << m_curEndPos
		<< "\t[ [" << formattedQcatFloat << ',' << contribOfEachState[0].second << ']';
      for (int i = 1; i < contribOfEachState.size(); i++)
	{
	  sprintf(formattedQcatFloat, "%.4g", contribOfEachState[i].first);
	  m_ofsQcat << ", [" << contribOfEachState[i].first << ',' << contribOfEachState[i].second << ']';
	}
      m_ofsQcat << " ]" << endl;
      m_curBegPos = m_curEndPos;
      m_curEndPos += m_regionWidth;
    }
  else
    m_ofsNullValues << retVal << endl;
  
  // reset the counting variables and the "P* numerator" (m_Ps1numerators, m_Ps2numerators) tallies
  m_numValsProcessedForGroup1 = m_numValsProcessedForGroup2 = 0;
  m_Ps1numerators.assign(m_Ps1numerators.size(), 0);
  if (!m_Ps2numerators.empty())
    m_Ps2numerators.assign(m_Ps2numerators.size(), 0);
}


bool KLssModel::getQcontrib(ifstream& infile, const char *pFilename, const int& Nsites)
{
  const int BUFSIZE(100000);
  char buf[BUFSIZE], *p;
  const float LOG2(0.6931471806), LOG_Nsites(log(static_cast<float>(Nsites)));
  float denom;
  vector<int> row, transposeRow;
  vector<vector<int> > tallyMatrix, transposedTallyMatrix;
  map<int, map<int,float> > tempMap;
  int linenum(0), fieldnum, numCols(0);

  while (infile.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 0;
      if (!(p = strtok(buf, "\t")))
	{
	  cerr << "Error:  Failed to parse line " << linenum << " of file " << pFilename << '.'
	       << endl << endl;
	  return false;
	}
      if (1 == linenum)
	{
	  row.push_back(atoi(p));
	  numCols++;
	  while (p = strtok(NULL, "\t"))
	    {
	      row.push_back(atoi(p));
	      numCols++;
	    }
	  m_numStates = static_cast<int>(floor(sqrt(static_cast<float>(numCols)) + 0.01));
	}
      else
	{
	  row[fieldnum++] = atoi(p);
	  while ((p = strtok(NULL, "\t")) && fieldnum < row.size())
	    row[fieldnum++] = atoi(p);
	  if (fieldnum != numCols)
	    {
	      cerr << "Error:  Found " << numCols << " columns on line 1 of " << pFilename
		   << " but only " << fieldnum << " columns on line " << linenum
		   << ".\nEach row must have the same number of columns; the # of columns "
		   << "must equal the square of the number of possible states\n"
		   << "(i.e., it must equal the number of possible state pairs)." << endl << endl;
	      return false;
	    }
	  else
	    {
	      if (p = strtok(NULL, "\t"))
		{
		  cerr << "Error:  Found " << numCols << " columns on line 1 of " << pFilename
		       << " but at least " << numCols+1 << " columns on line " << linenum
		       << ".\nEach row must have the same number of columns; the # of columns "
		       << "must equal the square of the number of possible states\n"
		       << "(i.e., it must equal the number of possible state pairs)." << endl << endl;
		  return false;
		}
	    }
	}
      tallyMatrix.push_back(row);
    }

  // The number of rows (linenum) equals numEpigenomes*(numEpigenomes - 1)/2,
  // so the number of epigenomes satisfies numEpigenomes^2 - numEpigenomes - 2*linenum == 0.
  if (0 == m_group1size)
    m_group1size = static_cast<int>(floor(1. + (sqrt(1. + 8.*static_cast<float>(linenum)))/2. + 0.001));
  else
    m_group2size = static_cast<int>(floor(1. + (sqrt(1. + 8.*static_cast<float>(linenum)))/2. + 0.001));
  
  transposeRow.assign(tallyMatrix.size(), 0);
  for (int c = 0; c < numCols; c++)
    {
      for (int r = 0; r < tallyMatrix.size(); r++)
	transposeRow[r] = tallyMatrix[r][c];
      transposedTallyMatrix.push_back(transposeRow);
    }
  
  denom = LOG2 * static_cast<float>(tallyMatrix.size());

  for (int statePairID = 0; statePairID < transposedTallyMatrix.size(); statePairID++)
    {
      map<int,float> epigenomePairID_to_QssContrib;
      for (int epigenomePairID = 0; epigenomePairID < transposedTallyMatrix[0].size(); epigenomePairID++)
	{
	  if (transposedTallyMatrix[statePairID][epigenomePairID] != 0)
	    epigenomePairID_to_QssContrib[epigenomePairID] =
	      (LOG_Nsites - log(static_cast<float>(transposedTallyMatrix[statePairID][epigenomePairID]))) / denom;
	  else
	    epigenomePairID_to_QssContrib[epigenomePairID] = 999999.;	    
	}
      tempMap[statePairID + 1] = epigenomePairID_to_QssContrib;
    }

  if (m_Qss1contrib.empty())
    m_Qss1contrib = tempMap;
  else
    m_Qss2contrib = tempMap;

  m_size = m_group1size*(m_group1size - 1)/2 + m_group2size*(m_group2size - 1)/2;
  return true;
}

bool KLssModel::processInputValue(const int& statePairID)
{
  static map<int, set<int> > emptyMap;
  static pair<map<int, set<int> >, map<int, set<int> > > pairOfEmptyMaps(emptyMap, emptyMap);
  bool processingGroup1(true);
  int remainder = statePairID % m_numStates, statePairGroupID;

  // Check these variables, in case excess columns appear in this line of input.
  // This is also how we determine, in the case of two groups of epigenomes being compared,
  // whether the input value is for group 1 or group 2.
  if (m_numValsProcessedForGroup1 == m_group1size*(m_group1size-1)/2)
    {
      if (m_numValsProcessedForGroup2 == m_group2size*(m_group2size-1)/2)
	{
	  cerr << "Error:  Found excess columns in a line of input; expected "
	       << m_group1size + m_group2size << "." << endl;
	  return false;
	}
      else
	processingGroup1 = false;
    }

  if (remainder != 0)
    {
      int quotient = statePairID / m_numStates;
      if (quotient + 1 > remainder) // reflect statePairID across the matrix diagonal, from the lower triangular matrix to the upper one
	statePairGroupID = m_numStates*(remainder - 1) + (quotient + 1);
      else
	statePairGroupID = statePairID;
    }
  else
    statePairGroupID = statePairID;
  map<int, pair<map<int, set<int> >, map<int, set<int> > > >::iterator it =
    m_statePairGroupObservationsAtThisSite.find(statePairGroupID);
  if (m_statePairGroupObservationsAtThisSite.end() == it)
    {
      set<int> thisEpigenomePairID;      
      m_statePairGroupObservationsAtThisSite[statePairGroupID] = pairOfEmptyMaps;
      if (processingGroup1)
	{
	  thisEpigenomePairID.insert(m_numValsProcessedForGroup1++);
	  m_statePairGroupObservationsAtThisSite[statePairGroupID].first[statePairID] = thisEpigenomePairID;
	}
      else
	{
	  thisEpigenomePairID.insert(m_numValsProcessedForGroup2++);
	  m_statePairGroupObservationsAtThisSite[statePairGroupID].second[statePairID] = thisEpigenomePairID;
	}
    }
  else
    {
      map<int, set<int> >::iterator it2;
      if (processingGroup1)
	{
	  it2 = it->second.first.find(statePairID);
	  if (it->second.first.end() == it2)
	    {
	      set<int> thisEpigenomePairID;
	      thisEpigenomePairID.insert(m_numValsProcessedForGroup1++);
	      it->second.first[statePairID] = thisEpigenomePairID;
	    }
	  else
	    it2->second.insert(m_numValsProcessedForGroup1++);
	}
      else
	{
	  it2 = it->second.second.find(statePairID);
	  if (it->second.second.end() == it2)
	    {
	      set<int> thisEpigenomePairID;
	      thisEpigenomePairID.insert(m_numValsProcessedForGroup2++);
	      it->second.second[statePairID] = thisEpigenomePairID;
	    }
	  else
	    it2->second.insert(m_numValsProcessedForGroup2++);
	}
    }

  return true;
}

void KLssModel::computeAndWriteMetric(void)
{
  float retVal(0);
  static vector<pair<float, int> > contribOfEachState(m_numStates, make_pair(0,0));
  float contribOfMaxStatePairGroupTerm(0), contribOfMaxStateTerm(0);
  int statePairGroupWithMaxTerm_1based, stateWithMaxTerm_1based;
  bool signOfMaxStatePairGroupTerm, signOfMaxStateTerm;

  if (!m_writeNullMetric)
    {
      // Initialize the vector of per-state contributions. This is needed for the "qcat" file.
      for (int i = 0; i < m_numStates; i++)
	{
	  contribOfEachState[i].first = 0;
	  contribOfEachState[i].second = i+1;
	}
    }
  
  for (map<int, pair<map<int, set<int> >, map<int, set<int> > > >::const_iterator statePairGroup_it =
	 m_statePairGroupObservationsAtThisSite.begin();
       statePairGroup_it != m_statePairGroupObservationsAtThisSite.end(); statePairGroup_it++)
    {
      const int &statePairGroupID = statePairGroup_it->first;
      float term = 0; // The contribution to D_KL from each state pair group.  Each encompasses state pairs (a,b) and (b,a), or (a,a) alone.
      float absTerm; // |term|
      int row, column; // 1-based row and column numbers of the upper triangular matrix of state pair groups;
                       // these are the states of the two epigenomes within an epigenome pair
      if (!m_writeNullMetric) // no need to break down the metric by state if we're solely tasked with writing null metric values
	{
	  column = statePairGroupID % m_numStates;
	  row = statePairGroupID / m_numStates + 1;
	  if (0 == column)
	    {
	      column = m_numStates;
	      row -= 1;
	    }
	}

      for (map<int, set<int> >::const_iterator group1_it = statePairGroup_it->second.first.begin();
	   group1_it != statePairGroup_it->second.first.end(); group1_it++)
	{
	  const int &statePairID = group1_it->first;
	  for (set<int>::const_iterator epigenomePair_it = group1_it->second.begin();
	       epigenomePair_it != group1_it->second.end(); epigenomePair_it++)
	    {
	      const int &epigenomePairID = *epigenomePair_it;
	      term += (m_Qss1contrib.find(statePairID)->second).find(epigenomePairID)->second;
	    }
	}
      // The following loop will only be executed if a comparison between groups is being made.
      for (map<int, set<int> >::const_iterator group2_it = statePairGroup_it->second.second.begin();
	   group2_it != statePairGroup_it->second.second.end(); group2_it++)
	{
	  const int &statePairID = group2_it->first;
	  for (set<int>::const_iterator epigenomePair_it = group2_it->second.begin();
	       epigenomePair_it != group2_it->second.end(); epigenomePair_it++)
	    {
	      const int &epigenomePairID = *epigenomePair_it;
	      term -= (m_Qss2contrib.find(statePairID)->second).find(epigenomePairID)->second;
	    }
	}

      absTerm = fabs(term);
      if (!m_writeNullMetric) // no need to break down the metric by state if we're solely tasked with writing null metric values
	{
	  if (absTerm > fabs(contribOfMaxStatePairGroupTerm))
	    {
	      contribOfMaxStatePairGroupTerm = term;
	      statePairGroupWithMaxTerm_1based = statePairGroupID;
	    }
	  contribOfEachState[row - 1].first += 0.5*term;    // row = state of epigenome 1 (1-based)
	  contribOfEachState[column - 1].first += 0.5*term; // column = state of epigenome 2 (1-based)
	}
      retVal += (0 == m_group2size ? term : absTerm);
    }

  if (!m_writeNullMetric)
    {
      int s1 = statePairGroupWithMaxTerm_1based / m_numStates + 1,
	s2 = statePairGroupWithMaxTerm_1based % m_numStates; // statePairGroupID = (s1,s2)
      char formattedQcatFloat[10];
      if (0 == s2)
	{
	  s2 = m_numStates;
	  s1 -= 1;
	}
      sort(contribOfEachState.begin(), contribOfEachState.end(), PairFloatInt_LT);
      if (fabs(contribOfEachState.front().first) > contribOfEachState.back().first)
	{
	  contribOfMaxStateTerm = -1*contribOfEachState.front().first;
	  stateWithMaxTerm_1based = contribOfEachState.front().second;
	  signOfMaxStateTerm = false;
	}
      else
	{
	  contribOfMaxStateTerm = contribOfEachState.back().first;
	  stateWithMaxTerm_1based = contribOfEachState.back().second;
	  signOfMaxStateTerm = true;
	}
      signOfMaxStatePairGroupTerm = contribOfMaxStatePairGroupTerm > 0 ? true : false;
      contribOfMaxStatePairGroupTerm = fabs(contribOfMaxStatePairGroupTerm);
      m_ofsObs << m_chrom << '\t' << m_curBegPos << '\t' << m_curEndPos << '\t'
	       << stateWithMaxTerm_1based << '\t' << contribOfMaxStateTerm << '\t'
	       << (signOfMaxStateTerm ? "1" : "-1") << '\t'
	       << '(' << s1 << ',' << s2 << ')' << '\t' // << statePairWithMaxTerm_1based << '\t'
	       << contribOfMaxStatePairGroupTerm << '\t'
	       << (signOfMaxStatePairGroupTerm ? "1" : "-1") << '\t'
	       << retVal << endl;
      sprintf(formattedQcatFloat, "%.4g", contribOfEachState[0].first);
      m_ofsQcat << m_chrom << '\t' << m_curBegPos << '\t' << m_curEndPos
		<< "\t[ [" << formattedQcatFloat << ',' << contribOfEachState[0].second << ']';
      for (int i = 1; i < contribOfEachState.size(); i++)
	{
	  sprintf(formattedQcatFloat, "%.4g", contribOfEachState[i].first);
	  m_ofsQcat << ", [" << formattedQcatFloat << ',' << contribOfEachState[i].second << ']';
	}
      m_ofsQcat << " ]" << endl;
      m_curBegPos = m_curEndPos;
      m_curEndPos += m_regionWidth;
    }
  else
    m_ofsNullValues << retVal << endl;
  
  // reset counting variables and the map
  m_numValsProcessedForGroup1 = m_numValsProcessedForGroup2 = 0;
  m_statePairGroupObservationsAtThisSite.clear();
}


bool parseInputWriteOutput(ifstream& ifs, const char *pFilename, Model* pModel);
bool parseInputWriteOutput(ifstream& ifs, const char *pFilename, Model* pModel)
{
  const int BUFSIZE(100000);
  char buf[BUFSIZE], *p;
  int linenum(0), numColsProcessed;

  while (ifs.getline(buf,BUFSIZE))
    {
      linenum++;
      numColsProcessed = 0;

      p = strtok(buf, "\t");
      do {
	int inputValue(atoi(p));
	numColsProcessed++;
	if (!pModel->processInputValue(inputValue))
	  {
	    cerr << "The error was detected in column " << numColsProcessed
		 << " of line " << linenum << " of file " << pFilename
		 << "." << endl << endl;
	    return false;
	  }
      } while ((p = strtok(NULL, "\t")));

      if (numColsProcessed != pModel->size())
	{
	  cerr << "Error:  Expected to find " << pModel->size() << " columns of integers "
	       << "on line " << linenum << " of " << pFilename
	       << ", but only found " << numColsProcessed << "." << endl << endl;
	  return false;
	}

      pModel->computeAndWriteMetric();
    }
  return true;
}


int main(int argc, const char* argv[])
{
  if (10 != argc && 11 != argc && 7 != argc)
    {
    Usage:
      cerr << "Usage type #1:  " << argv[0] << " infile KLtype NsitesGenomewide infileQ outfileObs outfileQcat chr firstBegPos regionWidth [infileQ2]\n"
	   << "where\n"
	   << "* infile holds the tab-delimited state or state pair IDs observed in the epigenomes or pairs of epigenomes, one line per genomic segment\n"
	   << "* KLtype is either 0 (for KL), 1 (for KL*), or 2 (for KL**)\n"
	   << "  KL compares states, KL* compares tallies of state pairs, and KL** compares state pairs of individual epigenome pairs\n"
	   << "* NsitesGenomewide is the total number of sites observed genome-wide\n"
	   << "* infileQ contains the Q, Q*, or Q** tally matrix (also see below)\n"
	   << "* outfileObs will receive genomic coordinates (regions on chromosome \"chr\" of width regionWidth, starting at firstBegPos),\n"
	   << "  the state (or state pair) making the largest contribution to the metric,\n"
	   << "  the magnitude of that contribution, and the total value of the metric.\n"
	   << "  If two groups are specified (see below), it will also include a column containing +/-1,\n"
	   << "  specifying whether the first group (+1) or the second (-1) contributes more to the overall metric.\n"
	   << "* outfileQcat will be in \"qcat\" format, uncompressed\n"
	   << "* Optional additional argument infileQ2 can be used to specify Q, Q*, or Q** for a 2nd group of epigenomes,\n"
	   << "  in which case the metric quantifies the difference (distance) between them.\n"
	   << "\n"
	   << "Usage type #2:  " << argv[0] << " infile KLtype NsitesGenomewide infileQ1 infileQ2 outfileNulls\n"
	   << "where\n"
	   << "* infile contains random permutations of states (or state pairs) observed in the initial input data\n"
	   << "* outfileNulls will receive the total difference metric for each line of permuted states (or state pairs)\n"
	   << "* the remaining arguments are the same as described above\n"
	   << "This second \"usage type\" is used to generate a null distribution, for estimating significance\n"
	   << "of the metric values calculated via \"usage type 1.\""
	   << endl << endl;
      return -1;
    }

  const char *pOutfileObsFilename(NULL), *pOutfileQcatFilename(NULL), *pOutfileNullValsFilename(NULL),
    *pInfilename(argv[1]), *pQ1filename(argv[4]), *pQ2filename(NULL);
  ifstream infile(pInfilename), infileQ1(pQ1filename), infileQ2;
  ofstream outfileObs, outfileQcat, outfileNulls;
  const int measurementTypeInt(atoi(argv[2])), Nsites(atoi(argv[3]));
  int firstBegPos, regionWidth;
  string chrom;
  KLModel mKL;
  KLsModel mKLs;
  KLssModel mKLss;
  Model *pM;
  
  if (KL != measurementTypeInt && KLs != measurementTypeInt && KLss != measurementTypeInt)
    {
      cerr << "Error:  Invalid \"measurementType\" received (2nd argument, \"" << argv[2] << "\").\n"
	   << "The valid options are 0 (to use KL), 1 (to use KL*), and 2 (to use KL**)." << endl << endl;
      goto Usage;
    }
  if (!infile)
    {
      cerr << "Error:  Unable to open file \"" << pInfilename << "\" for reading." << endl << endl;
      goto Usage;
    }
  if (!infileQ1)
    {
      cerr << "Error:  Unable to open file \"" << pQ1filename << "\" for reading." << endl << endl;
      goto Usage;
    }  
  if (7 == argc)
    {
      pQ2filename = argv[5];
      pOutfileNullValsFilename = argv[6];
    }
  else
    {
      pOutfileObsFilename = argv[5];
      pOutfileQcatFilename = argv[6];
      chrom = string(argv[7]);
      firstBegPos = atoi(argv[8]);
      regionWidth = atoi(argv[9]);
      if (11 == argc)
	pQ2filename = argv[10];
    }
  if (pQ2filename != NULL)
    {
      infileQ2.open(pQ2filename);
      if (!infileQ2)
	{
	  cerr << "Error:  Unable to open file \"" << pQ2filename << "\" for reading." << endl << endl;
	  goto Usage;
	}
    }

  if (KL == measurementTypeInt)
    pM = &mKL;
  else
    {
      if (KLs == measurementTypeInt)
	pM = &mKLs;
      else
	pM = &mKLss;
    }

  if (!pM->init(pOutfileObsFilename, pOutfileQcatFilename, pOutfileNullValsFilename,
		chrom, firstBegPos, regionWidth))
    return -1;
  if (!pM->getQcontrib(infileQ1, pQ1filename, Nsites))
    return -1;
  if (infileQ2.is_open() && !pM->getQcontrib(infileQ2, pQ2filename, Nsites))
    return -1;

  if (!parseInputWriteOutput(infile, pInfilename, pM))
    return -1;

  return 0;
}
