#ifndef DAIRYMODEL_HPP
#define DAIRYMODEL_HPP

#include <vector>
#include <deque>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

//#include "R.h"

#include "bodyweight.hh"
#include "basicdt.hh"


//-----------------------------------------------------------------------------

/** Class for creating the binary files of the HMDP with daily yield
measurements used in the JDS paper.

The model is generated in 2 steps:

1)
*/
class DairyModelDaily {
 public:
public:
    /**
     * \param binName Vector containing the names of the binary files.
     * \param vS1 Vector of state numbers at level 1 (vS1[n] = number of states
     *        at stage n (int)).
     * \param vS2 Vector of state numbers at level 2 (vS2[n] = number of states
     *        at stage n (int)).
     * \param matLab1 2-dim vector of labels at level 1 (matLab1[n][s] = label
     *        of state s at stage n (string)).
     * \param matLab2 2-dim vector of labels at level 2 (matLab2[n][s] = label
     *        of state s at stage n (string)).
     * \param sDry Vector of possible dry weeks (int).
     * \param prIC 3-dim vector of IC probabilities (prIC[lac][dfc][preg] = pr
     *        at lactation lac and day dfc given preganancy preg (flt)).
     * \param prICDry Vector of IC pr in dry period (prICDry[lac] = pr in
     *        lactation lac (flt)).
     * \param 3-dim vector of vectors (prM prM[lac][dfc][iM] = vector containing
     *        pairs (iM,pr) denoting the pr to iM at next stage (vector<flt>)).
     * \param prM2A 2-dim vector of vectors (prM2A[lac][iM] = vector containing
     *        pairs (iA,pr) (vector<flt>)).
     * \param prPregT 2-dim vector of pr of positive pregnancy test (prPreg[lac][dfc]).
     * \param yield yield[lac][dfc][iM]
     */
    DairyModelDaily(vector<string> binName,
        int dryPeriodLth,
        int gestLth,
        flt probCalf,
        flt probBullCalf,
        flt priceCalf,
        flt priceBullCalf,
        flt priceHeifer,
        flt priceSFU,
        flt priceECM,
        flt priceCarcassKg,
        int insemStart, int insemFinish,
        int pregTestLth,
        int iAZero);

    /** Set vector from array. */
    void SetVS1(intPtr p, int size);

    /** Set vector from array. */
    void SetVS2(intPtr p, int size);

    /** Set vector from array. */
    void SetVA2M(intPtr p, int size);

    /** Set vector from array. */
    /*void SetSDry(intPtr p, int size) {
        sDry.assign(p, p+size);
    }*/

    /** Set vector from array. */
    void SetPrICDry(fltPtr p, int size) {
        prICDry.assign(p, p+size);
    }

    /** Get the pointer. */
    vector< vector<flt> > * GetPrM2APtr() {
        //Rprintf("tjk: %d %d\n",&prM2A,prM2A.size());
        return &prM2A;
    }

    /** Get the pointer. */
    vector< vector< vector<flt> > > * GetPrICPtr() {
        return &prIC;
    }

    /** Get the pointer. */
    vector< vector< vector< vector<flt> > > > * GetPrMPtr() {
        return &prM;
    }

    /** Get the pointer. */
    vector< vector< vector<flt> > > * GetYieldPtr() {
        return &yield;
    }

    /** Get the pointer. */
    vector< vector< flt> > * GetPrPregTPtr() {
        return &prPregT;
    }

    /** Set the labels of level 1 stored in matLab1.
     * \pre Assume that vS1 and have been set.
     */
    string LabelS1(int d1, int s1);

    string LabelS2(int d2,int s2,bool last);

    /** Get the pointer to the string matrix. */
    /*vector< vector<string> > * GetMatLab1Ptr() {
        vector< vector<string> > * p = &matLab1;
        return p;
    }*/

    /** Get the pointer to the string matrix. */
    /*vector< vector<string> > * GetMatLab2Ptr() {
        vector< vector<string> > * p = &matLab2;
        return p;
    }*/

    ~DairyModelDaily() {}

    /** Generate all binary state files for the model */
    void GenerateStatesBinary();

    /** Generate all binary action files for the model */
    void GenerateActionsBinary();


 private:

    /** Create action weight labels */
    void AddActionWeightLabels(const string & filename);

    /** Add founder action.
     \pre The state/action variables must have been set.
     */
    void AddFounderAction(const int & sId, int & aId, string & label,
        vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
        ofstream & fA, ofstream & fALbl, ofstream & fTransP,
        ofstream & fAWeight);


    /** Add level 1 action.
     \pre The state/action variables must have been set.
     */
    void AddLev1Action(const int & sId, const vector<int> & sIdx, int & aId,
        string & label,
        vector<int> & aScpIdx,
        vector<flt> & pr,
        vector<flt>& weights,
        ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight);


    /** Create actions for a specific state which are added to the binary files.
     \pre The state/action variables must have been set.
     */
    void AddLev2Actions(const int & sId, const vector<int> & sIdx, int & aId,
        string & label, vector<int> & aScpIdx, vector<flt> & pr,
        vector<flt>& weights,
        ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight);


    /** Replace action at level 2.
     */
    void ReplaceAction(int & d1, int & d2, int & iDry,
        const int & sId, int & aId, string & label,
        vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
        ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight);


    /** Dry action at level 2.
     */
    void DryAction(int & d1, int & d2, int & iM, int & iDry,
        const int & sId, int & aId, string & label,
        vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
        ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight);


    /** Keep action at level 2.
     */
    void KeepAction(int & d1, int & d2, int & s2, int & iM, int & iDry,
        const int & sId, int & aId, string & label,
        vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
        ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight);


    /** Reward of the carcass. */
    flt RewardCarcass(int lac, int dfc, int con) {
        //Rprintf("CarcassW:%f price:%f\n",oBW.CarcassWeight(lac, dfc, con),priceCarcassKg);
        return oBW.CarcassWeight(lac, dfc, con)*priceCarcassKg;
    }

    /** Add the action to the binary files.
     \pre Assume that aScpIdx, weights, pr and label have been set.
     */
    void AddActionBinary(const int & sId, int & aId, const string & label,
        const vector<int> & aScpIdx,
        const vector<flt> & pr,
        const vector<flt> & weights,
        ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight);

    /** Add the state to the binary files.
     \pre Assume that sIdx and label have been set.
     */
    void AddStateBinary(int & sId, const string & label,
        const vector<int> & sIdx, ofstream & fS, ofstream & fSLbl);

    /** Write to binary file a vector of integers. */
    void WriteBinary(ofstream & out, vector<int> content) {
        out.write(reinterpret_cast<const char*>(&*content.begin()), content.size()*sizeof(int));
    }

    /** Write to binary file a vector of floats. */
    void WriteBinary(ofstream & out, vector<flt> content) {
        out.write(reinterpret_cast<const char*>(&*content.begin()), content.size()*sizeof(flt));
    }

    /** Write to binary file a string. */
    void WriteBinary(ofstream & out, string label) {
        //Rprintf("write:%s ",label.c_str());
        out.write(label.c_str(), label.size()+1);   // add the null terminator also (\0)
    }

    /** Write to binary file an integer. */
    void WriteBinary(ofstream & out, int id) {
        out.write(reinterpret_cast<char*>(&id), sizeof(id));
    }

    /** Write to binary file a float. */
    void WriteBinary(ofstream & out, flt id) {
        out.write(reinterpret_cast<char*>(&id), sizeof(id));
    }


private:


	/**
	 * Returns the index of the shortest possible week where dry at this stage.
	 * Do not consider the unknown state.
	 * @param n The stage number (dfc).
	 * @return The index of the shortest possible week where dry.
	 */
    int GetMinDryWeekIdx(int n) {
		if (n < insemStart+pregTestLth) // status still unknown. Return last idx
			return sizeDry-1;
		if (n > insemStart + gestLth - dryPeriodLth)
			return Dfc2Week(n)-Dfc2Week(insemStart+gestLth-dryPeriodLth);	// idx of current week
		return 0;
	}

	/**
	 * Returns the shortest possible week where dry at this stage.
	 * Do not consider the unknown state.
	 * @param n The stage number (dfc).
	 * @return The shortest possible week where dry.
	 */
    int GetMinDryWeek(int n) {
		if (n < insemStart+pregTestLth) // status still unknown.
			return -1;
		if (n <= insemStart + gestLth - dryPeriodLth)
			return Dfc2Week(insemStart+gestLth-dryPeriodLth);
		if (n > insemStart + gestLth - dryPeriodLth)
			return Dfc2Week(n);
		return -200;
	}

	/**
	 * Returns the longest possible calving interval to be defined
	 * at this stage. Assume that the preg test is taken prior such that the
	 * observation at day insemStart+pregTestLth could be that the cow i pregnant!
	 * @param n The stage number.
	 * @return The longest possible calving interval.
	 */
    int GetMaxDryWeek(int n) {
    	if (n < insemStart+pregTestLth) // status still unknown.
			return -1;
		if (n < insemFinish+pregTestLth)
			return Dfc2DryWeek(n);
    	return Dfc2Week(insemFinish+gestLth-dryPeriodLth);	// last possible week
	}

	/**
	 * Returns the index of the longest possible calving interval to be defined
	 * at this stage. Assume that the preg test is taken prior such that the
	 * observation at day insemStart+pregTestLth could be that the cow i pregnant!
	 * @param n The stage number.
	 * @return The index of the longest possible calving interval.
	 */
    int GetMaxDryWeekIdx(int n) {
    	if (n < insemStart+pregTestLth) // status still unknown. Return last idx
			return sizeDry-1;
		if (n < insemFinish+pregTestLth)
			return Dfc2DryWeek(n)-Dfc2DryWeek(insemStart+pregTestLth);
    	return Dfc2Week(insemFinish+gestLth-dryPeriodLth)-Dfc2Week(insemStart+gestLth-dryPeriodLth);	// last possible week idx
	}

	/**
	 * Returns the index of the calving interval corresponding to a specified state at a
	 * given stage.
	 * @param n The stage number.
	 * @param i The state number.
	 * @return The corresponding index of the dry week (numbered from zero).
	 */
    int GetDryWeekIdx(int n, int i) {
    	if (i >= GetLev2States(n)-1-sizeM)	// preg unknown state
    		return sizeDry-1;
		int minI = GetMinDryWeekIdx(n);
		return (i / sizeM) + minI;
	}

    /**
     * Returns the day of conception.
     * @param idxDry
     * @return Day of conception.
     */
    int GetDOC(int idxDry) {
    	if (idxDry==sizeDry-1) return maxDfc + 1000;
    	return Week2Dfc(GetDryWeek(idxDry))+dryPeriodLth-gestLth;
    }

	/**
	 * Returns the dry week corresponding to idxDry. Note idxDry must not correspond to the unknown state!
	 * @param idxDry
	 * @return The dry week.
	 */
    int GetDryWeek(int idxDry) {
		return firstDryWeek+idxDry;
	}

	/**
	 * Returns the index of the latent milk yield corresponding to a specified state at a
	 * given stage.
	 * @param i The state number.
	 * @return The corresponding index of the latent milk yield (numbered from zero).
	 */
    int GetMIdx(int i) {
		return (i % sizeM);
	}

    /**
     * Return the index of the state at stage n.
     * @param n Stage number (dfc).
     * @param idxM Index of latent mean.
     * @param idxDry Index of dry week.
     * @return The index of the state at stage n.
     */
    int GetStateIdx(int n, int idxM, int idxDry) {
    	if (idxDry==sizeDry-1) // preg unknown
    		return (GetLev2States(n)-1-sizeM) + idxM;
    	int minI = GetMinDryWeekIdx(n);
		return sizeM * (idxDry-minI) + idxM;
	}

	/**
	 * Returns the number of states at a specified stage of a process at level 2.
	 * @param n The stage number.
	 * @return The number of states.
	 */
    int GetLev2States(int n) {
    	if (n>maxDfc) return 0;
    	if (this->GetMaxDryWeekIdx(n)==sizeDry-1) {	// pregnancy unknown
    		return sizeM + 1;
    	}
		return (this->GetMaxDryWeekIdx(n) - this->GetMinDryWeekIdx(n) + 2)*sizeM + 1;
	}

	/** Return the probability of M (prM[lac][dfc][iM] is a vector with (idx, pr, idx, pr, ...))*/
	flt GetPrM(int & lac, int & dfc, int & iM, int & i) {
        if (lac>=(int)prM.size()) return prM[(int)prM.size()][dfc][iM][i];    // assume that prM equal for lac's higher
        return prM[lac][dfc][iM][i];
    }

    /**
     * Return week number given that we are at day dfc.
     * @param dfc Days from calving.
     * @return Week number.
     */
    int Dfc2Week(int dfc) {
		return (dfc-1)/7+1;
	}



    /**
     * Return week number were dry given a positive preg test at day dfc.
     * @param dfc Days from calving.
     * @return Dry week number.
     */
    int Dfc2DryWeek(int dfc) {
		return Dfc2Week(dfc-pregTestLth+gestLth-dryPeriodLth);
	}

    /**
     * Return idxDry were dry given a positive preg test at day dfc.
     * @param dfc Days from calving.
     * @return Index of Dry week number.
     */
    int Dfc2DryWeekIdx(int dfc) {
		return Dfc2DryWeek(dfc)-firstDryWeek;
	}

    /**
     * Return last day of week.
     * @param w Week in lactation.
     * @return Last day of week.
     */
    int Week2Dfc(int w) {
		return w*7;
	}

	/*
	 * Read the prM csv file and store it in prM.
	 * @param file File name of csv.
	 */
	/*void ReadPrM(string file) {
        ifstream csvFile(file);
        string line;
        //int linenum = 0;
        while (getline (csvFile, line))
        {
            //linenum++;
            //cout << "\nLine #" << linenum << ":" << endl;
            istringstream linestream(line);

            int lac,dfc,iM;
            //int itemnum = 0;
            getline (linestream, lac, ',')
            getline (linestream, dfc, ',')
            getline (linestream, iM, ',')
            flt item;
            while (getline (linestream, item, ','))
            {
                //itemnum++;
                //cout << "Item #" << itemnum << ": " << item << endl;
                prM[lac][dfc][iM].push_back(item)
            }
        }

        return 0;
	}*/



    vector< vector<int> > sIdS2;    ///< Matrix storing the sId's used for shared child process linking (stage 1 in the child process). sIdS2[d1][s2] = sId of state s2 (corresponding to state s2=iM) at level 2 and stage 1 given level one at stage d1. The dim of the matrix will be (maxLac+1,sizeM).

    vector<int> vS1;    ///< Vector containing the number of level 1 states (size=maxLac+2). vS1[i] = number of states at stage i.
	vector<int> vS2;    ///< Vector containing the number of level 2 states (size=maxDfc+1). vS2[i] = number of states at stage i. Note vS2[1] = sizeM.
	vector<int> vA2M;   ///< Vector containing the mapping of iA to iM (vA2M[iA] = iM).
	int iAZero;         ///< The index iA to the interval containing A=0.
	//vector< vector<string> > matLab1;    ///< Matrix containing the labels of level 1.  matLab1[i][j] = label of state j at stage i.
	//vector< vector<string> > matLab2;    ///< Matrix containing the labels of level 2. matLab2[i][j] = label of state j at stage i.
    //vector<int> sDry;   ///< Possible dry weeks sDry[i] = dry week for iDry = i.
    vector< vector< vector<flt> > > prIC;       ///< Array of pr for IC (prIC[lac][dfc][preg=T/F]).
    vector<flt> prICDry;    ///< Array of pr for IC in dry period (prICDry[lac]).
    vector< vector< vector< vector<flt> > > > prM; ///< Array of pr M (prM[lac][dfc][iM] is a vector with (idx, pr, idx, pr, ...))
    vector< vector<flt> > prM2A; ///< Array of pr M to A (prM2A[iM] is a vector with (idx, pr, idx, pr, ...))
    vector< vector< vector<flt> > > yield; ///< Array of daily yield in ECM kg (yield[lac][dfc][iM]).
    vector< vector<flt> > prPregT; ///< Array of pr of positive pregnancy test (prPregT[lac][dfc]).


    // model variables
    uSInt maxLac;             ///< Max number of lactations.
    uSInt maxDfc;             ///< Max dfc.
    flt priceHeifer;
    uSInt sizeDry;
    uSInt firstDryWeek;       ///< First week where possible to dry.
    uSInt sizeM;
    uSInt sizeA;
    uSInt dryPeriodLth;
    uSInt gestLth;
    flt probCalf;
    flt priceCalf;
    flt probBullCalf;
    flt priceBullCalf;
    flt priceSFU;
    flt priceECM;
	flt priceCarcassKg;
    uSInt insemStart, insemFinish;
    uSInt pregTestLth;

    // binary files
    string fSFile;        ///< Binary file for state idx.
    string fSLblFile;     ///< Binary file for state labels.
    string fAFile;        ///< Binary file for action idx.
    string fALblFile;         ///< Binary file for action labels.
    string fAWeightFile;      ///< Binary file for action costs/weights.
    string fTransPFile;       ///< Binary file for action transition probabilities.

    BodyWeight oBW;     ///< Body weight object.
};

//-----------------------------------------------------------------------------

#endif
