#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <string>
#include <utility> // for pair()

using namespace std;

const float g_changeOfScale(1.0e+7); // for converting metrics to integers

struct NullData {
  long metricAsInt;
  int numOccs;
  float pvalue;
};

struct OutputValue {
  string inputEntry;
  long metricAsInt;
  long outputOrder;
  float pvalue;
};

bool Metric_GT(const OutputValue& a, const OutputValue& b);
bool Metric_GT(const OutputValue& a, const OutputValue& b)
{
  return a.metricAsInt > b.metricAsInt;
}

bool OutputOrder_LT(const OutputValue& a, const OutputValue& b);
bool OutputOrder_LT(const OutputValue& a, const OutputValue& b)
{
  return a.outputOrder < b.outputOrder;
}

// The input file is assumed to contain an arbitrary number of columns of data,
// the last of which is the metric for which a p-value will be estimated.
// Columns 1-3 are assumed to be genomic coordinates (chr, begin, end).
// The output file will be same as the input file, but with p-value estimates appended.
// FDR estimates will need to be made for the p-values by another program/procedure.

bool loadDataAndReport(ifstream& ifs, ofstream& ofs, const vector<NullData>& nullDistn);
bool loadDataAndReport(ifstream& ifs, ofstream& ofs, const vector<NullData>& nullDistn)
{
  const int BUFSIZE(10000);
  char buf[BUFSIZE], *p;
  string prevColumnStr;
  OutputValue ov;
  vector<OutputValue> outvec;
  long linenum(0);
  int fieldnum, expectedFinalFieldNum(-1);
  
  while (ifs.getline(buf,BUFSIZE))
    {
      linenum++;
      fieldnum = 0;
      ov.inputEntry = string(buf);
      ov.outputOrder = linenum;
      prevColumnStr.clear();
      p = strtok(buf, "\t");
      while (p != NULL)
	{
	  fieldnum++;
	  prevColumnStr = string(p);
	  p = strtok(NULL, "\t");
	}
      if (1 == linenum)
	expectedFinalFieldNum = fieldnum;
      else
	{
	  if (fieldnum != expectedFinalFieldNum)
	    {
	      cerr << "Error:  Detected " << expectedFinalFieldNum << " column(s) of data on line 1 of the input file,\n"
		   << "but detected " << fieldnum << " column(s) of data on line " << linenum << '.' << endl << endl;
	      return false;
	    }
	}
      ov.metricAsInt = static_cast<long>(floor(atof(prevColumnStr.c_str())*g_changeOfScale + 0.5));
      outvec.push_back(ov);
    }

  sort(outvec.begin(), outvec.end(), Metric_GT);
  float pval(0);
  unsigned int idxNull(0);
  for (unsigned int idxObs = 0; idxObs < outvec.size(); idxObs++)
    {
      while (idxNull < nullDistn.size() && nullDistn[idxNull].metricAsInt > outvec[idxObs].metricAsInt)
	idxNull++;
      if (nullDistn.size() == idxNull || (idxNull != 0 && nullDistn[idxNull].metricAsInt < outvec[idxObs].metricAsInt))
	idxNull--;
      if (!(0 == idxNull && nullDistn[idxNull].metricAsInt < outvec[idxObs].metricAsInt))
	pval = nullDistn[idxNull].pvalue;
      outvec[idxObs].pvalue = pval;
    }
  sort(outvec.begin(), outvec.end(), OutputOrder_LT);
  for (unsigned int idxObs = 0; idxObs < outvec.size(); idxObs++)
    ofs << outvec[idxObs].inputEntry << '\t' << outvec[idxObs].pvalue << endl;

  return true;
}

void loadNullDistn(ifstream& ifs, vector<NullData>& ndistn);
void loadNullDistn(ifstream& ifs, vector<NullData>& ndistn)
{
  const int BUFSIZE(100);
  char buf[BUFSIZE];
  int numNullValues(0);
  map<int,int> tempMap; // This map provides an efficient means for tallying occurrences of values.
  map<int,int>::iterator it;
  pair<int,int> mapElementToInsert;
  mapElementToInsert.second = 1;

  while (ifs.getline(buf,BUFSIZE))
    {
      int value = static_cast<int>(floor(atof(buf)*g_changeOfScale + 0.5));
      // We use lower_bound() rather than find() because for values that haven't yet been inserted
      // into the map, the former will give us a "hint" of where the new value should be inserted.
      it = tempMap.lower_bound(value); // >=
      if (tempMap.end() == it)
	{
	  // no lower bound found in the map; every map element is smaller, or the map is empty
	  mapElementToInsert.first = value;
	  if (tempMap.begin() != it)
	    tempMap.insert(--it, mapElementToInsert);
	  else // map is empty
	    tempMap.insert(mapElementToInsert);
	}
      else
	{
	  if (it->first == value)
	    it->second++;
	  else
	    {
	      if (it != tempMap.begin())
		it--;
	      mapElementToInsert.first = value;
	      tempMap.insert(it, mapElementToInsert);
	    }
	}
      numNullValues++;
    }

  if (0 == numNullValues)
    {
      cerr << "Error:  Received an empty file of null values." << endl << endl;
      exit(2);
    }

  float N(static_cast<float>(numNullValues));
  int runningTallyOfOccurrences(0);
  NullData ndata;
  it = tempMap.end();
  it--;
  while (it != tempMap.begin())
    {
      ndata.metricAsInt = it->first;
      ndata.numOccs = it->second;
      runningTallyOfOccurrences += ndata.numOccs;
      ndata.pvalue = static_cast<float>(runningTallyOfOccurrences) / N;
      ndistn.push_back(ndata);
      it--;
    }
  ndata.metricAsInt = it->first;
  ndata.numOccs = it->second;
  ndata.pvalue = 1.;
  ndistn.push_back(ndata);

}

int main(int argc, const char* argv[])
{
  if (4 != argc)
    {
      cerr << "Usage:  " << argv[0] << " infile nullDistnFile outfile\n"
	   << "where \"nullDistnFile\" contains random values that constitute a null distribution,\n"
	   << "and the values in the final column of \"infile\" are to be compared with the null values\n"
	   << "to obtain p-value estimates.\n"
	   << "The contents of \"infile,\", with p-values appended, are written to \"outfile.\"\n"
	   << "FDR estimates will need to be made for the p-values by another program/procedure."
	   << endl << endl;
      return -1;
    }

  ifstream infile(argv[1]);
  if (!infile)
    {
      cerr << "Error:  Failed to open file \"" << argv[1] << "\" for read." << endl << endl;
      return -1;
    }
  ifstream nullDistnFile(argv[2]);
  if (!nullDistnFile)
    {
      cerr << "Error:  Failed to open file \"" << argv[2] << "\" for read." << endl << endl;
      return -1;
    }
  ofstream outfile(argv[3]);
  if (!outfile)
    {
      cerr << "Error:  Failed to open file \"" << argv[3] << "\" for write." << endl << endl;
      return -1;
    }
  vector<NullData> nullDistn;

  loadNullDistn(nullDistnFile, nullDistn);
  if (!loadDataAndReport(infile, outfile, nullDistn))
    return -1;

  return 0;
}

