#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <utility> // for pair
#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

enum measurementType {KL, KLs, KLss}; // corresponding to measurements using KL, KL*, or KL**

bool parseOneSetOfColumnSpecs(char* pString, set<int>& colsOut);
bool parseOneSetOfColumnSpecs(char* pString, set<int>& colsOut)
{
  string origString(pString);
  char *p1 = NULL, *p2 = pString;
  int first(-1), last(-1); // first and last column numbers in a range, e.g. 1 and 10 in the range 1-10
  while (true)
    {
      if (',' == *p2 || '\0' == *p2)
	{
	  if (p1 != NULL)
	    {
	      if (',' == *p2)
		*p2++ = '\0';
	      last = atoi(p1);
	      if (first < 0)
		first = last;
	      for (int i = first; i <= last; i++)
		{
		  set<int>::iterator it = colsOut.find(i);
		  if (it != colsOut.end())
		    {
		      cerr << "Warning:  Value " << i << " specified multiple times in \""
			   << origString << "\"" << endl;
		      continue;
		    }
		  colsOut.insert(i);
		}
	      if ('\0' == *p2)
		break;
	    }
	  else
	    cerr << "Warning:  Range string \"" << origString << "\" begins with \',\'; anything intended to precede it "
		 << "will necessarily be excluded." << endl << endl;
	  p1 = p2;
	  first = last = -1;
	  continue;
	}
      else
	{
	  if ('-' == *p2)
	    {
	      if (p1 != NULL)
		{
		  *p2++ = '\0';
		  first = atoi(p1);
		  p1 = p2;
		  if ('\0' == *p1)
		    {
		      cerr << "Error:  Range \"" << origString << "\" ends with \'-\'. The final value in the range must be explicitly specified."
			   << endl << endl;
		      return false;
		    }
		  continue;
		}
	      else
		{
		  cerr << "Error:  Range \"" << origString << "\" begins with \'-\'. Must begin with a valid column number, 1 or greater."
		       << endl << endl;
		  return false;
		}
	    }
	  else
	    {
	      if (!isdigit(*p2))
		{
		  cerr << "Error:  Invalid character \'" << *p2 << "\' in range \"" << origString << "\"."
		       << endl << endl;
		  return false;
		}
	      if (NULL == p1)
		p1 = p2;
	      p2++;
	      continue;
	    }
	}
    }
  return true;
}

// The format of the input file is chromosome, beg position, end position,
// state of epigenome 1 at that site/region, state of epigenome 2 there, ....
// group1 and group2 define the columns of input data that should be assigned to groups 1 and 2
// (input column 4 = epigenome 1, column 5 = epigenome 2, etc., since columns 1-3 = chr:beg-end).
// group1 and group2 are guaranteed to have no entries in common;
// we only need to ensure there's actually a column for each entry.
// group2 may be empty, in which case we're measuring the properties of a single group of epigenomes.
// If group2 is not empty, then we're comparing the properties of two groups of epigenomes,
// and at each site, we shuffle the states observed in the union of epigenomes from the two groups
// and write additional results to ofsRand, to be used later to estimate P-values for the observations.
// If group2 is not empty, tallies contributing to Q1, Q1*, or Q1** (measurement types KL, KLs, KLss
// respectively) are written to output file ofsQ and tallies contributing to Q2, Q2*, or Q2**
// (tallies for group 2) are written to output file ofsQ2.
// If group2 is empty, then tallies contributing to Q, Q*, or Q** are written to ofsQ
// and nothing is written to ofsQ2 (which is not an open ofstream in this case).
// The total number of sites (i.e., the number of lines in input file ifs) is written to ofsNsites.

bool onePassThroughData(ifstream& ifs, const measurementType& KLtype, const set<int>& group1, const set<int>& group2,
			const int& numStates, ofstream& ofsP, ofstream& ofsQ,
			ofstream& ofsQ2, ofstream& ofsRand, ofstream& ofsNsites);
bool onePassThroughData(ifstream& ifs, const measurementType& KLtype, const set<int>& group1, const set<int>& group2,
			const int& numStates, ofstream& ofsP, ofstream& ofsQ,
			ofstream& ofsQ2, ofstream& ofsRand, ofstream& ofsNsites)
{
  const int BUFSIZE(10000);
  char buf[BUFSIZE], *p;
  const bool comparisonOfGroups(group2.empty() ? false : true);
  const int numStatePairs(numStates * numStates), numUniqueStatePairs(numStates*(numStates+1)/2);
  int numEpiPairs1, numEpiPairs2;
  vector<int> allStatesAtThisSite, zeroes, shuffledStatesInThe2groupsAtThisSite;
  vector<int> P1, P2, Ps1, Ps2, Q1, Q2, Qs1, Qs2, randP1, randP2, randPs1, randPs2;
  vector<vector<int> > Qss1, Qss2;
  int linenum(0), fieldnum, numFieldsOnLineOne, k;

  // Set up the variables which will hold the parsed state data.
  
  if (KLtype != KL)
    {
      if (KLss == KLtype)
	{
	  // We track specific state pairs in specific epigenome pairs.
	  numEpiPairs1 = (group1.size() * (group1.size()-1)) / 2;
	  zeroes.assign(numStatePairs, 0);
	  for (int i = 0; i < numEpiPairs1; i++)
	    Qss1.push_back(zeroes);
	}
      else
	{
	  // We tally observances of state pairs across epigenome pairs.
	  // For numStates states, there are numStates*numStates ordered state pairs,
	  // and (numStates*numStates - numStates)/2 + numStates = numStates*(numStates+1)/2
	  // unique (unordered) state pairs.
	  Ps1.assign(numUniqueStatePairs, 0);
	  Qs1.assign(numUniqueStatePairs, 0);
	}
    }
  else
    {
      // We tally observances of states.
      P1.assign(numStates, 0);
      Q1.assign(numStates, 0);
    }
  
  if (comparisonOfGroups)
    {
      shuffledStatesInThe2groupsAtThisSite.assign(group1.size() + group2.size(), 0);
      if (KLtype != KL)
	{
	  if (KLss == KLtype)
	    {
	      numEpiPairs2 = (group2.size() * (group2.size()-1)) / 2;
	      for (int i = 0; i < numEpiPairs2; i++)
		Qss2.push_back(zeroes);
	    }
	  else
	    {
	      Ps2.assign(numUniqueStatePairs, 0);
	      Qs2.assign(numUniqueStatePairs, 0);
	      randPs1.assign(numUniqueStatePairs, 0);
	      randPs2.assign(numUniqueStatePairs, 0);
	    }
	}
      else
	{
	  P2.assign(numStates, 0);
	  Q2.assign(numStates, 0);
	  randP1.assign(numStates, 0);
	  randP2.assign(numStates, 0);
	}
    }

  // One line at a time, read in the states observed in all epignomes,
  // including all epigenomes _not_ being analyzed.
  // Then select the states observed in the epigenomes of interest,
  // possibly in two groups of epigenomes that will be compared later.
  // If there are two groups, shuffle the observations between them
  // and write those random observations to ofsRand.
  
  while (ifs.getline(buf,BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      p = strtok(buf, "\t");
      // field 1:  chromosome, ignore
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
	{
	MissingField:
	  cerr << "Error:  Failed to find field " << fieldnum
	       << " on line " << linenum << " of the input file."
	       << endl << endl;
	  return false;
	}
      // field 2:  begin site, ignore
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
	goto MissingField;
      // field 3:  end site, ignore

      if (1 == linenum)
	{
	  while (p = strtok(NULL, "\t"))
	    {
	      const int thisState(atoi(p));
	      if (thisState > numStates || thisState < 1)
		{
		  cerr << "Error:  Illegal state (" << thisState
		       << ") detected in field " << ++fieldnum << " on line "
		       << linenum << ".  Re-specify the correct number of possible states,\n"
		       << "and/or ensure each state label is a positive (nonzero) integer."
		       << endl << endl;
		  return false;
		}
	      allStatesAtThisSite.push_back(thisState);
	      fieldnum++;
	    }
	  numFieldsOnLineOne = fieldnum;
	}
      else
	{
	  while (p = strtok(NULL, "\t"))
	    {
	      const int thisState(atoi(p));
	      if (thisState > numStates || thisState < 1)
		{
		  cerr << "Error:  Illegal state (" << thisState
		       << ") detected in field " << ++fieldnum << " on line "
		       << linenum << ".  Re-specify the correct number of possible states,\n"
		       << "and/or ensure each state label is a positive (nonzero) integer."
		       << endl << endl;
		  return false;
		}
	      if (fieldnum - 3 < allStatesAtThisSite.size())
		allStatesAtThisSite[fieldnum - 3] = thisState;
	      else
		{
		  cerr << "Error:  Expected to find " << numFieldsOnLineOne
		       << " fields of data on line " << linenum
		       << ", to match the # found on line 1; found at least "
		       << ++fieldnum << " instead." << endl << endl;
		  return false;
		}
	      fieldnum++;
	    }
	  if (fieldnum != numFieldsOnLineOne)
	    {
	      cerr << "Error:  Expected to find " << numFieldsOnLineOne
		   << " fields of data on line " << linenum
		   << ", to match the # found on line 1; but only found "
		   << fieldnum << "." << endl << endl;
	      return false;
	    }
	}

      if (comparisonOfGroups)
	{
	  // Shuffle the subset of states observed in the two groups at this site.
	  k = 0;
	  for (set<int>::const_iterator it = group1.begin(); it != group1.end(); it++)
	    shuffledStatesInThe2groupsAtThisSite[k++] = allStatesAtThisSite[*it - 1];
	  for (set<int>::const_iterator it = group2.begin(); it != group2.end(); it++)
	    shuffledStatesInThe2groupsAtThisSite[k++] = allStatesAtThisSite[*it - 1];
	  random_shuffle(shuffledStatesInThe2groupsAtThisSite.begin(), shuffledStatesInThe2groupsAtThisSite.end());
	}

      if (KL == KLtype)
	{
	  P1.assign(P1.size(), 0);
	  if (comparisonOfGroups)
	    {
	      randP1.assign(P1.size(), 0);
	      k = 0;
	    }
	  for (set<int>::const_iterator it = group1.begin(); it != group1.end(); it++)
	    {
	      P1[allStatesAtThisSite[*it - 1] - 1]++;
	      Q1[allStatesAtThisSite[*it - 1] - 1]++;
	      if (comparisonOfGroups)
		randP1[shuffledStatesInThe2groupsAtThisSite[k++] - 1]++;
	    }
	  ofsP << P1[0];
	  if (comparisonOfGroups)
	    ofsRand << randP1[0];
	  for (int i = 1; i < P1.size(); i++)
	    {
	      ofsP << '\t' << P1[i];
	      if (comparisonOfGroups)
		ofsRand << '\t' << randP1[i];
	    }
	  if (comparisonOfGroups)
	    {
	      P2.assign(P2.size(), 0);
	      randP2.assign(P1.size(), 0);
	      for (set<int>::const_iterator it = group2.begin(); it != group2.end(); it++)
		{
		  P2[allStatesAtThisSite[*it - 1] - 1]++;
		  Q2[allStatesAtThisSite[*it - 1] - 1]++;
		  randP2[shuffledStatesInThe2groupsAtThisSite[k++] - 1]++;
		}
	      for (int i = 0; i < P2.size(); i++)
		{
		  ofsP << '\t' << P2[i];
		  ofsRand << '\t' << randP2[i];
		}
	    }
	}
      else
	{
	  int k2, j(0);
	  if (KLs == KLtype)
	    {
	      Ps1.assign(numUniqueStatePairs, 0);
	      if (comparisonOfGroups)
		{
		  Ps2.assign(numUniqueStatePairs, 0);
		  randPs1.assign(numUniqueStatePairs, 0);
		  randPs2.assign(numUniqueStatePairs, 0);
		}
	    }
	  // Loop over all ordered state pairs in all ordered epigenome IDs,
	  // recording each such observation in P* or P**,
	  // and recording the tally of such observations in Q* or Q**.
	  // (Technically, these aren't P*, P**, Q*, or Q**,
	  // but rather the main components of them; they'll be "completed"
	  // during subsequent processing by another program.)
	  k = 0;
	  k2 = k + 1;
	  for (set<int>::const_iterator it1 = group1.begin(); it1 != group1.end(); it1++)
	    {
	      set<int>::const_iterator it2 = it1;
	      it2++;
	      for (; it2 != group1.end(); it2++)
		{
		  // Example with 15 states:  Ordered state pairs (1,1), (1,2), ..., (1,15) map to 1, 2, ..., 15;
		  // ordered state pairs (2,1), (2,2), ..., (2,15) map to 16, 17, ..., 30;
		  // state pair (15,15) maps to 225.
		  int stateOfEpi1 = allStatesAtThisSite[*it1 - 1], stateOfEpi2 = allStatesAtThisSite[*it2 - 1];
		  if (KLss == KLtype)
		    {
		      int thisStatePairID = (stateOfEpi1 - 1)*numStates + stateOfEpi2;
		      if (0 == j)
			ofsP << thisStatePairID;
		      else
			ofsP << '\t' << thisStatePairID;
		      Qss1[j++][thisStatePairID - 1]++;
		    }
		  else // KLs == KLtype
		    {
		      int thisUniqueStatePairID;
		      if (stateOfEpi1 > stateOfEpi2)
			{
			  // swap these, to make the calculation simpler
			  int temp = stateOfEpi1;
			  stateOfEpi1 = stateOfEpi2;
			  stateOfEpi2 = temp;
			}
		      // Upper triangular matrix, entries = 0, 1, ..., numStates-1 in row 1; numStates, ..., 2*numStates - 2 in row 2;
		      // entry for the bottom right corner is numStates*(numStates+1)/2 - 1.
		      thisUniqueStatePairID = (stateOfEpi1 - 1)*numStates - (stateOfEpi1-2)*(stateOfEpi1-1)/2 + (stateOfEpi2-stateOfEpi1);
		      Ps1[thisUniqueStatePairID]++;
		      Qs1[thisUniqueStatePairID]++;
		    }
		  if (comparisonOfGroups)
		    {
		      stateOfEpi1 = shuffledStatesInThe2groupsAtThisSite[k];
		      stateOfEpi2 = shuffledStatesInThe2groupsAtThisSite[k2];
		      if (KLss == KLtype)
			{
			  int thisStatePairID = (stateOfEpi1 - 1)*numStates + stateOfEpi2;
			  if (0 == k && 1 == k2)
			    ofsRand << thisStatePairID;
			  else
			    ofsRand << '\t' << thisStatePairID;
			}
		      else // KLs == KLtype
			{
			  int thisUniqueStatePairID;
			  if (stateOfEpi1 > stateOfEpi2)
			    {
			      // swap these, to make the calculation simpler
			      int temp = stateOfEpi1;
			      stateOfEpi1 = stateOfEpi2;
			      stateOfEpi2 = temp;
			    }
			  thisUniqueStatePairID = (stateOfEpi1 - 1)*numStates - (stateOfEpi1-2)*(stateOfEpi1-1)/2 + (stateOfEpi2-stateOfEpi1);
			  randPs1[thisUniqueStatePairID]++;
			}
		      if (group1.size() == ++k2)
			{
			  k++;
			  k2 = k + 1;
			}
		    }
		}
	    }
	  if (comparisonOfGroups)
	    {
	      k = j = 0;
	      k2 = k + 1;
	      for (set<int>::const_iterator it1 = group2.begin(); it1 != group2.end(); it1++)
		{
		  set<int>::const_iterator it2 = it1;
		  it2++;
		  for (; it2 != group2.end(); it2++)
		    {
		      // See the comments a few lines above for an explanation of the mapping
		      // from state pair to a unique ID (integer) representing the ordered (or unordered, unique) state pair.
		      int stateOfEpi1 = allStatesAtThisSite[*it1 - 1], stateOfEpi2 = allStatesAtThisSite[*it2 - 1];
		      if (KLss == KLtype)
			{
			  int thisStatePairID = (stateOfEpi1 - 1)*numStates + stateOfEpi2;
			  ofsP << '\t' << thisStatePairID;
			  Qss2[j++][thisStatePairID - 1]++;
			}
		      else // KLs == KLtype
			{
			  int thisUniqueStatePairID;
			  if (stateOfEpi1 > stateOfEpi2)
			    {
			      int temp = stateOfEpi1;
			      stateOfEpi1 = stateOfEpi2;
			      stateOfEpi2 = temp;
			    }
			  thisUniqueStatePairID = (stateOfEpi1 - 1)*numStates - (stateOfEpi1-2)*(stateOfEpi1-1)/2 + (stateOfEpi2-stateOfEpi1);
			  Ps2[thisUniqueStatePairID]++;
			  Qs2[thisUniqueStatePairID]++;
			}
		      stateOfEpi1 = shuffledStatesInThe2groupsAtThisSite[k + group1.size()];
		      stateOfEpi2 = shuffledStatesInThe2groupsAtThisSite[k2 + group1.size()];
		      if (KLss == KLtype)
			{
			  int thisStatePairID = (stateOfEpi1 - 1)*numStates + stateOfEpi2;			  
			  ofsRand << '\t' << thisStatePairID;
			}
		      else
			{
			  int thisUniqueStatePairID;
			  if (stateOfEpi1 > stateOfEpi2)
			    {
			      int temp = stateOfEpi1;
			      stateOfEpi1 = stateOfEpi2;
			      stateOfEpi2 = temp;
			    }
			  thisUniqueStatePairID = (stateOfEpi1 - 1)*numStates - (stateOfEpi1-2)*(stateOfEpi1-1)/2 + (stateOfEpi2-stateOfEpi1);
			  randPs2[thisUniqueStatePairID]++;
			}
		      if (++k2 + group1.size() == shuffledStatesInThe2groupsAtThisSite.size())
			{
			  k++;
			  k2 = k + 1;
			}
		    }
		}
	    }
	  if (KLs == KLtype) // then we haven't yet written the P* contributions to disk
	    {
	      ofsP << Ps1[0];
	      if (comparisonOfGroups)
		ofsRand << randPs1[0];
	      for (int i = 1; i < Ps1.size(); i++)
		{
		  ofsP << '\t' << Ps1[i];
		  if (comparisonOfGroups)
		    ofsRand << '\t' << randPs1[i];
		}
	      if (comparisonOfGroups)
		{
		  for (int i = 0; i < Ps2.size(); i++)
		    {
		      ofsP << '\t' << Ps2[i];
		      ofsRand << '\t' << randPs2[i];
		    }
		}
	    }
	} // end of if/else to handle choice of KL, KL*, or KL**

      ofsP << endl;
      ofsRand << endl;
    } // end of loop for reading and processing all input data

  // Write out the tallies over sites, for eventual use in Q, Q*, or Q**.

  switch (KLtype) {

  case KL:
    ofsQ << Q1[0];
    for (int i = 1; i < Q1.size(); i++)
      ofsQ << '\t' << Q1[i];
    if (comparisonOfGroups)
      {
	ofsQ2 << Q2[0];
	for (int i = 1; i < Q2.size(); i++)
	  ofsQ2 << '\t' << Q2[i];
	ofsQ2 << endl;
      }
    ofsQ << endl;
    break;

  case KLs:
    ofsQ << Qs1[0];
    for (int i = 1; i < Qs1.size(); i++)
      ofsQ << '\t' << Qs1[i];
    if (comparisonOfGroups)
      {
	ofsQ2 << Qs2[0];
	for (int i = 1; i < Qs2.size(); i++)
	  ofsQ2 << '\t' << Qs2[i];
	ofsQ2 << endl;
      }
    ofsQ << endl;
    break;
    
  case KLss:  
    for (int i = 0; i < Qss1.size(); i++)
      {
	ofsQ << Qss1[i][0];
	for (int j = 1; j < Qss1[0].size(); j++)
	  ofsQ << '\t' << Qss1[i][j];
	ofsQ << endl;
      }
    if (comparisonOfGroups)
      {
	for (int i = 0; i < Qss2.size(); i++)
	  {
	    ofsQ2 << Qss2[i][0];
	    for (int j = 1; j < Qss2[0].size(); j++)
	      ofsQ2 << '\t' << Qss2[i][j];
	    ofsQ2 << endl;
	  }
      }
    break;
  }

  ofsNsites << linenum << endl;
  
  return true;
}

int main(int argc, char* argv[])
{
  if (8 != argc && 11 != argc && 2 != argc && 3 != argc)
    {
    Usage:
      cerr << "Usage flavor 1:  " << argv[0] << " infile measurementType numStates outfileP outfileQ outfileNsites groupSpec [group2spec outfileRandP outfileQ2]\n"
	   << "where\n"
	   << "* infile is tab-delimited: chrom, start, stop, state of epigenome1, state of epigenome2, ...\n"
	   << "* measurementType is either 0 (to use KL), 1 (KL*), or 2 (KL**)\n"
	   << "* numStates is the number of possible states (e.g. 15) that can be observed in \"infile\"\n"
	   << "* outfileP, in 1-to-1 correspondence with \"infile,\" will contain state tallies or state-pair tallies\n"
	   << "* outfileQ will contain chromosome-wide state tallies or state-pair tallies\n"
	   << "* outfileNsites will contain the total number of lines in \"infile\"\n"
	   << "* groupSpec (group definition) is a comma- and/or dash-delimited range of integers enclosed in double-quotes, such as \"1,4-6,8,9\"\n"
	   << "  (1 = epigenome in column 4, 2 = epigenome in column 5, ...)\n"
	   << "Calling " << argv[0] << " with the above 7 arguments specifies measurement of a single group of epigenomes.\n"
	   << "Specification of optional additional arguments group2spec, outfileRandP, and outfileQ2  will instigate a comparison between specified groups 1 and 2;\n"
	   << "outfileP will contain tallies for both groups, outfileQ will contain tallies for group1, outfileQ2 will contain tallies for group2,\n"
	   << "and outfileRandP will contain tallies obtained after randomly permuting the states observed in group 1 and group 2 among the union of all epigenomes.\n"
	   << "\n"
	   << "Usage flavor 2:  " << argv[0] << " groupSpec [group2spec]\n"
	   << "where groupSpec (and optional group2spec) are defined as above.\n"
	   << "In this case, the group definition is parsed, and its size (number of epigenomes) is written to standard output;\n"
	   << "if two group specifications are provided, the numbers of epigenomes in the two groups are written to standard output, separated by a tab character.\n"
	   << "The program does nothing further in this scenario."
	   << endl << endl;
      return -1;
    }

  if (2 == argc || 3 == argc)
    {
      set<int> grp1;
      if (!parseOneSetOfColumnSpecs(argv[1], grp1))
	return -1;
      if (3 == argc)
	{
	  set<int> grp2;
	  if (!parseOneSetOfColumnSpecs(argv[2], grp2))
	    return -1;
	  for (set<int>::const_iterator it2 = grp2.begin(); it2 != grp2.end(); it2++)
	    {
	      set<int>::const_iterator it1 = grp1.find(*it2);
	      if (it1 != grp1.end())
		{
		  cerr << "Error:  Value " << *it1 << " found in both group specifications." << endl << endl;
		  return -1;
		}
	    }
	  cout << grp1.size() << '\t' << grp2.size() << endl;
	}
      else
	cout << grp1.size() << endl;
      return 0;
    }
  
  ifstream infile(argv[1]);
  const int measurementTypeInt(atoi(argv[2])), numStates(atoi(argv[3]));
  ofstream outfileP(argv[4]), outfileQ(argv[5]), outfileNsites(argv[6]), outfileRand, outfileQ2;

  set<int> group1, group2;

  if (KL != measurementTypeInt && KLs != measurementTypeInt && KLss != measurementTypeInt)
    {
      cerr << "Error:  Invalid \"measurementType\" received (2nd argument, \"" << argv[2] << "\").\n"
	   << "The valid options are 0 (to use KL), 1 (to use KL*), and 2 (to use KL**)." << endl << endl;
      goto Usage;
    }
  if (!infile)
    {
      cerr << "Error:  Unable to open input file \"" << argv[1] << "\" for read." << endl << endl;
      return -1;
    }
  if (!outfileP)
    {
      cerr << "Error:  Unable to open output file \"" << argv[4] << "\" for write." << endl << endl;
      return -1;
    }
  if (!outfileQ)
    {
      cerr << "Error:  Unable to open output file \"" << argv[5] << "\" for write." << endl << endl;
      return -1;
    }
  if (!outfileNsites)
    {
      cerr << "Error:  Unable to open output file \"" << argv[6] << "\" for write." << endl << endl;
      return -1;
    }
  if (!parseOneSetOfColumnSpecs(argv[7], group1))
    return -1;

  if (11 == argc)
    {
      if (!parseOneSetOfColumnSpecs(argv[8], group2))
	return -1;
      for (set<int>::const_iterator it2 = group2.begin(); it2 != group2.end(); it2++)
	{
	  set<int>::const_iterator it1 = group1.find(*it2);
	  if (it1 != group1.end())
	    {
	      cerr << "Error:  Value " << *it1 << " was found in both group specifications." << endl << endl;
	      return -1;
	    }
	}
      outfileRand.open(argv[9]);
      if (!outfileRand)
	{
	  cerr << "Error:  Unable to open output file \"" << argv[9] << "\" for write." << endl << endl;
	  return -1;
	}
      outfileQ2.open(argv[10]);
      if (!outfileQ2)
	{
	  cerr << "Error:  Unable to open output file \"" << argv[10] << "\" for write." << endl << endl;
	  return -1;
	}
    }

  if (!onePassThroughData(infile, static_cast<measurementType>(measurementTypeInt), group1, group2, numStates,
			  outfileP, outfileQ, outfileQ2, outfileRand, outfileNsites))
    return -1;

  return 0;
}
