#include "dairymodel.hh"

// ----------------------------------------------------------------------------

void DairyModelDaily::GenerateStatesBinary() {
    ofstream fS,fSLbl;
    fS.open(fSFile.c_str(), ios::out | ios::binary);
    fSLbl.open(fSLblFile.c_str(), ios::out | ios::binary);
    // state variables
    int sId;            ///< Id of the state we consider.
    vector<int> sIdx;   ///< Vector of index (d0,s0,a0,d1,s1,a1,d2,s2) (or less in at level 0 or 1).
    string label;       ///< Label of the action or state.

    sIdx.clear();
    sId = 0;
    label = "Dummy";
    sIdx.assign(2,0);   // set to dummy state at level 0
    AddStateBinary(sId, label, sIdx, fS, fSLbl);   // Add founder state
    sIdx.push_back(0);  // a0 action idx, i.e. sIdx = (0,0,0)
    for (int d1=0;d1<=maxLac+1;++d1) {
        sIdx.push_back(d1);
        for (int s1=0; s1<vS1[d1]; ++s1) {
            sIdx.push_back(s1);
            label = LabelS1(d1, s1);
            //Rprintf("sIdx: %s - label: %s\n",VectorStr<int>(sIdx).c_str(),label.c_str());
            AddStateBinary(sId, label, sIdx, fS, fSLbl);
            if (s1==0 & d1>0 & d1<maxLac+1) { // if subprocess
                for (int a1=0; a1<1; ++a1) {
                    sIdx.push_back(a1);
                    for (int d2=0; d2<=maxDfc; ++d2) {
                        sIdx.push_back(d2);
                        for (int s2=0; s2<vS2[d2]; ++s2) {
                            sIdx.push_back(s2);
                            label = LabelS2(d2,s2,s2==vS2[d2]-1);
                            //Rprintf("  sIdx: %s - label: %s\n",VectorStr<int>(sIdx).c_str(),label.c_str());
                            if (d2==1) {
                                sIdS2[d1][s2] = sId; // store id for shared child linking.
                                //Rprintf("sIdS2,sId:%d\n",sId);
                            }
                            AddStateBinary(sId, label, sIdx, fS, fSLbl);
                            sIdx.pop_back();
                        }
                        sIdx.pop_back();
                    }
                    sIdx.pop_back();
                }
            }
            sIdx.pop_back();
        }
        sIdx.pop_back();
    }
    fS.close();
    fSLbl.close();
    //Rprintf("sIdS2:%s\n",VectorStr<int>(sIdS2[1]).c_str());
}

// ----------------------------------------------------------------------------

DairyModelDaily::DairyModelDaily(vector<string> binName,
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
    int iAZero)
{
    fSFile = binName[0];
    fSLblFile = binName[1];
    fAFile = binName[2];
    fALblFile = binName[3];
    fAWeightFile = binName[4];
    AddActionWeightLabels(binName[5]);
    fTransPFile = binName[6];

    this->dryPeriodLth = dryPeriodLth;
    this->gestLth = gestLth;
    this->probCalf = probCalf;
    this->probBullCalf = probBullCalf;
    this->priceCalf = priceCalf;
    this->priceBullCalf = priceBullCalf;
    this->priceHeifer = priceHeifer;
    this->priceSFU = priceSFU;
    this->priceECM = priceECM;
    this->priceCarcassKg = priceCarcassKg;
    this->insemStart = insemStart;
    this->insemFinish = insemFinish;
    this->pregTestLth = pregTestLth;
    this->iAZero = iAZero;

    firstDryWeek = Dfc2Week(insemStart+gestLth-dryPeriodLth);
    sizeDry = Dfc2Week(insemFinish+gestLth-dryPeriodLth)-Dfc2Week(insemStart+gestLth-dryPeriodLth)+2;	// size of sDry
}

// ---------------------------------------------------------------------------

void DairyModelDaily::SetVS1(intPtr p, int size) {
    vS1.assign(p, p+size);
    maxLac = vS1.size()-2;
    sizeA = vS1[1]-1;
    /*Rprintf("vS1: ");
    for (int i=0;i<vS1.size();++i) Rprintf("i,n = %d,%d,  ",i,vS1[i]);
    Rprintf("\n");*/
}

// ---------------------------------------------------------------------------

void DairyModelDaily::SetVS2(intPtr p, int size) {
    vS2.assign(p, p+size);
    maxDfc = vS2.size()-1;
    sizeM = vS2[1]-1;
    //Rprintf("sizeM:%d\n",sizeM);
    // set size of sIdS2
    sIdS2.resize(maxLac+1);
    for (int d1=0; d1<maxLac+1; ++d1) sIdS2[d1].resize(sizeM);
    /*Rprintf("vS2: ");
    for (int i=0;i<vS2.size();++i) Rprintf("i,n = %d,%d,  ",i,vS2[i]);
    Rprintf("\n");*/
}

// ---------------------------------------------------------------------------

void DairyModelDaily::SetVA2M(intPtr p, int size) {
    vA2M.assign(p, p+size);
}

// ---------------------------------------------------------------------------

void DairyModelDaily::GenerateActionsBinary() {
    ofstream fA,fALbl,fAWeight,fTransP;
    fA.open(fAFile.c_str(), ios::out | ios::binary);
    fALbl.open(fALblFile.c_str(), ios::out | ios::binary);
    fAWeight.open(fAWeightFile.c_str(), ios::out | ios::binary);
    fTransP.open(fTransPFile.c_str(), ios::out | ios::binary);

    // state variables
    int sId;            ///< Id of the state we consider.
    vector<int> sIdx;   ///< Vector of index (d0,s0,a0,d1,s1,a1,d2,s2) (or less in at level 0 or 1).
    string label;       ///< Label of the action or state.
    // action variables (label above is also an action variable)
    int aId;            ///< Id of the action we consider.
    vector<int> aScpIdx;    ///< Used to store (scp, idx, scp, idx, ...).
    vector<flt> weights;    ///< Used to store action weights.
    vector<flt> pr;         ///< Used to store probabilities.

    sIdx.clear();
    sId = aId = 0;
    sIdx.assign(2,0);   // set to dummy state at level 0
    AddFounderAction(sId, aId, label, aScpIdx, pr, weights, fA, fALbl, fTransP, fAWeight);
    ++sId;
    sIdx.push_back(0);  // a0 action idx
    for (int d1=0;d1<=maxLac+1;++d1) {
        sIdx.push_back(d1);
        for (int s1=0; s1<vS1[d1]; ++s1) {
            sIdx.push_back(s1);
            AddLev1Action(sId, sIdx, aId, label, aScpIdx, pr, weights, fA, fALbl, fTransP, fAWeight);
            ++sId;
            if (s1==0 & d1>0 & d1<maxLac+1) { // if subprocess
                for (int a1=0; a1<1; ++a1) {
                    sIdx.push_back(a1);
                    for (int d2=0; d2<=maxDfc; ++d2) {
                        sIdx.push_back(d2);
                        //Rprintf("vS2[%d]:%d\n",d2,vS2[d2]);
                        for (int s2=0; s2<vS2[d2]; ++s2) {
                            sIdx.push_back(s2);
                            AddLev2Actions(sId, sIdx, aId, label, aScpIdx, pr, weights, fA, fALbl, fTransP, fAWeight);
                            ++sId;
                            sIdx.pop_back();
                        }
                        sIdx.pop_back();
                    }
                    sIdx.pop_back();
                }
            }
            sIdx.pop_back();
        }
        sIdx.pop_back();
    }
    fA.close();
    fALbl.close();
    fAWeight.close();
    fTransP.close();
}

// ---------------------------------------------------------------------------

void DairyModelDaily::AddActionWeightLabels(const string & filename) {
    ofstream fAWLbl;    ///< Binary file for action weight labels.
    fAWLbl.open(filename.c_str(), ios::out | ios::binary);
    WriteBinary(fAWLbl,string("Time"));
    WriteBinary(fAWLbl,string("Reward"));
    WriteBinary(fAWLbl,string("Yield"));
    fAWLbl.close();
}

// ---------------------------------------------------------------------------

void DairyModelDaily::AddActionBinary(const int & sId, int & aId,
    const string & label,
    const vector<int> & aScpIdx,
    const vector<flt> & pr,
    const vector<flt> & weights,
    ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight)
{
    int sepInt = -1;
    flt sepFlt = -1;
    string aIdStr = ToString<int>(aId);
    //Rprintf("aId:%d sId:%d aScpIdx:%s ",aId,sId,VectorStr<int>(aScpIdx).c_str());

    WriteBinary(fA,sId);
    WriteBinary(fA,aScpIdx);
    WriteBinary(fA,sepInt);

    //Rprintf("label:%s ",label.c_str());
    WriteBinary(fALbl,aIdStr);
    WriteBinary(fALbl,label);

    //Rprintf("pr:%s ",VectorStr<flt>(pr).c_str());
    WriteBinary(fTransP,pr);
    WriteBinary(fTransP,sepFlt);

    //Rprintf("w:%s \n",VectorStr<flt>(weights).c_str());
    WriteBinary(fAWeight,weights);
    ++aId;
}

// ---------------------------------------------------------------------------

void DairyModelDaily::AddStateBinary(int & sId, const string & label,
    const vector<int> & sIdx, ofstream & fS, ofstream & fSLbl)
{
    int sepInt = -1;
    string sIdStr = ToString<int>(sId);

    WriteBinary(fS,sIdx);
    WriteBinary(fS,sepInt);

    WriteBinary(fSLbl,sIdStr);
    WriteBinary(fSLbl,label);
    ++sId;
}

// ---------------------------------------------------------------------------

void DairyModelDaily::AddLev1Action(const int & sId, const vector<int> & sIdx, int & aId,
    string & label,
    vector<int> & aScpIdx,
    vector<flt> & pr,
    vector<flt>& weights,
    ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight)
{
    int d1 = sIdx[3];
    int s1 = sIdx[4];
    pr.clear(); aScpIdx.clear();
    label = "Dummy";
    weights.assign(3,0);
    pr.push_back(1);
    if (d1==0) {
        label = "Buy cow";
        weights[0]=0; weights[1] = -priceHeifer; weights[2]=0;
        aScpIdx.push_back(1);
        aScpIdx.push_back(iAZero);
        AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
        return;
    }
    if (d1==maxLac+1) {     // TODO (LRE#1#): Note have switched the two states so replaced last like in other stages!!
        if (s1==0) {    // Replace it state
            label = "Replace";
            weights[0]=0;
            weights[1] = RewardCarcass(maxLac+1, 1, maxDfc);
            weights[2]=0;
            aScpIdx.assign(2,0);

        }
        else {  // Replaced state
            aScpIdx.assign(2,0);
        }
        AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
        return;
    }
    // otherwise lactation stage
    if (s1==vS1[d1]-1) {   // replaced state
        aScpIdx.push_back(0);
        aScpIdx.push_back(0);
        AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
        return;
    }
    if (s1==0) {    // state where define child process
        aScpIdx.push_back(2);
        aScpIdx.push_back(0);
    }
    else {    // just make a hyperarc to stage 1 in the child process defined above (shared child process)
        aScpIdx.push_back(3);
        int iM=vA2M[s1];
        aScpIdx.push_back(sIdS2[d1][iM]);
    }
    AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
}

// ---------------------------------------------------------------------------

void DairyModelDaily::AddLev2Actions(const int & sId, const vector<int> & sIdx,
    int & aId, string & label, vector<int> & aScpIdx, vector<flt> & pr,
    vector<flt>& weights,
    ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight)
{
    int d1 = sIdx[3];
    int d2 = sIdx[6];
    int s2 = sIdx[7];
    //Rprintf("lac:%d dfc:%d s2:%d ",d1,d2,s2);
    pr.clear(); aScpIdx.clear();
    label = "Dummy";
    weights.assign(3,0);
    pr.push_back(1);

    if (d2==0) {   // action to the iM corresponding to iA=0
        aScpIdx.push_back(1);
        aScpIdx.push_back(vA2M[0]);
        AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
        return;
    }
    if (s2==vS2[d2]-1) {    // IC state
        aScpIdx.push_back(0);
        aScpIdx.push_back(vS1[d1+1]-1);
        AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
        return;
    }
    //label = LabelS2(d2,s2,false);
    int iM = GetMIdx(s2);
    int iDry = GetDryWeekIdx(d2, s2);
    //Rprintf("label:%s iM:%d iDry:%d ",label.c_str(),iM,iDry);
    if (d2==maxDfc) { // last stage
        if (iDry==sizeDry-1) {  // if not pregnant
            ReplaceAction(d1,d2,iDry,sId,aId,label,aScpIdx,pr,weights, fA, fALbl,
                fTransP, fAWeight);
        }
        else {  // dry or replace
            DryAction(d1,d2,iM,iDry,sId,aId,label,aScpIdx,pr,weights, fA, fALbl,
                fTransP, fAWeight);
            ReplaceAction(d1,d2,iDry,sId,aId,label,aScpIdx,pr,weights, fA, fALbl,
                fTransP, fAWeight);
        }
        return;
    }
    if (iDry==sizeDry-1) {    // if cow not pregnant
        KeepAction(d1,d2,s2,iM,iDry,sId,aId,label,aScpIdx,pr,weights,fA,fALbl,
            fTransP,fAWeight);
        ReplaceAction(d1,d2,iDry,sId,aId,label,aScpIdx,pr,weights,fA,fALbl,
            fTransP,fAWeight);
        return;
    }
    if (d2==Week2Dfc(GetDryWeek(iDry))) { // dry the cow
        DryAction(d1,d2,iM,iDry,sId,aId,label,aScpIdx,pr,weights, fA, fALbl,
            fTransP, fAWeight);
        ReplaceAction(d1,d2,iDry,sId,aId,label,aScpIdx,pr,weights, fA, fALbl,
            fTransP, fAWeight);
        return;
    }
    // else
    KeepAction(d1,d2,s2,iM,iDry,sId,aId,label,aScpIdx,pr,weights,fA,fALbl,
        fTransP,fAWeight);
    ReplaceAction(d1,d2,iDry,sId,aId,label,aScpIdx,pr,weights,fA,fALbl,
        fTransP,fAWeight);

}

// ---------------------------------------------------------------------------

void DairyModelDaily::ReplaceAction(int & d1, int & d2, int & iDry,
    const int & sId, int & aId, string & label,
    vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
    ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight)
{
    label = "Replace";
    weights[0] = 0;
    int con = GetDOC(iDry);
    weights[1] = RewardCarcass(d1, d2, con);
    weights[2] = 0;

    aScpIdx.clear();
    aScpIdx.push_back(0);
    aScpIdx.push_back(vS1[d1+1]-1);

    pr.clear();
    pr.push_back(1);

    AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
}

// ---------------------------------------------------------------------------

void DairyModelDaily::DryAction(int & d1, int & d2, int & iM, int & iDry,
    const int & sId, int & aId, string & label,
    vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
    ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight)
{
    label = "Dry";
    weights[0] = dryPeriodLth;
    int con = GetDOC(iDry);
    flt rewRep = prICDry[d1] * RewardCarcass(d1, (int)(d2+floor(dryPeriodLth/2)), con);
    flt rewKeep = (1-prICDry[d1]) * (probCalf*priceCalf + probBullCalf*priceBullCalf - priceSFU*oBW.sfuDry*dryPeriodLth);
    weights[1] =rewRep+rewKeep;
    weights[2] = 0;

    aScpIdx.clear(); pr.clear();
    if (d1==maxLac) {   // Replace it state
        aScpIdx.push_back(0);
        aScpIdx.push_back(0);
        pr.push_back(1);
        AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
        return;
    }

    for (idx i=0;i<prM2A[iM].size();++i) {  // add transpr from m to A
        if (i % 2 == 0) {  // even
            aScpIdx.push_back(0);
            aScpIdx.push_back((int)prM2A[iM][i]);
            continue;
        }
        pr.push_back(prM2A[iM][i]*(1-prICDry[d1]));
    }
    aScpIdx.push_back(0);   // add trans pr for IC
    //if (d2!=maxDfc)
    //if (d2==483) Rprintf("iDry=%d, iM=%d, s=%d\n ",iDry,iM,vS1[d1+1]-1);
    aScpIdx.push_back(vS1[d1+1]-1);
    //else aScpIdx.push_back(1);
    pr.push_back(prICDry[d1]);

    AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
        fTransP, fAWeight);
}

// ---------------------------------------------------------------------------

void DairyModelDaily::KeepAction(int & d1, int & d2, int & s2, int & iM, int & iDry,
    const int & sId, int & aId, string & label,
    vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
    ofstream & fA, ofstream & fALbl, ofstream & fTransP, ofstream & fAWeight)
{
    //Rprintf("\nKeep Action\n");
    label = "Keep";
    flt icPr;
    if (iDry==sizeDry-1) icPr = prIC[d1][d2][0];
    else icPr = prIC[d1][d2][1];    // pregnant
    int dfc = d2;
    int lac = d1;

    weights[0] = 1;
    int con = GetDOC(iDry);
    flt rewKeep = (1-icPr) * (yield[d1][d2][iM]*priceECM - priceSFU*oBW.feedingSFU(d1, d2, con, yield[d1][d2][iM], gestLth));
    flt rewRep = icPr * RewardCarcass(d1, d2, con);
    weights[1] =rewRep+rewKeep;
    weights[2] = yield[d1][d2][iM];
    //Rprintf("icPr:%f yield:%f sfu:%f\n",icPr,yield[d1][d2][iM],oBW.feedingSFU(d1, d2, con, yield[d1][d2][iM], gestLth));

    aScpIdx.clear(); pr.clear();
    if (iDry==sizeDry-1 & dfc<insemFinish+pregTestLth
        & dfc>=insemStart+pregTestLth-1)    // i.e. pregnancy unknown and possible to run preg test next day
    {
        flt prPT = prPregT[lac][dfc+1];
        // if continue open
        for (idx i=0;i<prM[d1][d2][iM].size();++i) {
            if (i % 2 == 0) {  // even
                aScpIdx.push_back(1);
                aScpIdx.push_back(GetStateIdx(d2+1, (int)prM[d1][d2][iM][i], iDry));
                continue;
            }
            pr.push_back(prM[d1][d2][iM][i]*(1-prPT)*(1-icPr));
        }
        // if become pregnant
        int iDryW = Dfc2DryWeekIdx(dfc+1);
        //if (dfc==74) Rprintf("iDryW=%d, iM=%d, prPT=%f\n",iDryW,iM,prPT);
        for (idx i=0;i<prM[d1][d2][iM].size();++i) {
            if (i % 2 == 0) {  // even
                aScpIdx.push_back(1);
                aScpIdx.push_back(GetStateIdx(d2+1, (int)prM[d1][d2][iM][i], iDryW));
                //if (dfc==273) Rprintf("sIdx=%d\n",GetStateIdx(d2+1, (int)prM[d1][d2][iM][i], iDryW));
                continue;
            }
            pr.push_back(prM[d1][d2][iM][i]*prPT*(1-icPr));
        }
    }
    else {  // already pregnant
        for (idx i=0;i<prM[d1][d2][iM].size();++i) {
            if (i % 2 == 0) {  // even
                aScpIdx.push_back(1);
                aScpIdx.push_back(GetStateIdx(d2+1, (int)prM[d1][d2][iM][i], iDry));
                continue;
            }
            pr.push_back(prM[d1][d2][iM][i]*(1-icPr));
        }
    }
    aScpIdx.push_back(1);
    aScpIdx.push_back(vS2[d2+1]-1);
    pr.push_back(icPr);
    AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
            fTransP, fAWeight);
}

// ---------------------------------------------------------------------------

string DairyModelDaily::LabelS1(int d1, int s1) {
    if (d1==0) return "Dummy";
    if (d1==maxLac+1 & s1==1) return "Replaced";
    if (d1==maxLac+1 & s1==0) return "Replace it";
    if (s1==sizeA) return "Replaced";
    string str = "idxA=";
    str += ToString<int>(s1);
    return str;
}

// ---------------------------------------------------------------------------

string DairyModelDaily::LabelS2(int d2,int s2,bool last) {
    if (d2==0) return "Dummy";
    if (last) return "Replaced due to IC";
    int iM = GetMIdx(s2);
    int iDry = GetDryWeekIdx(d2, s2);
    string str = ToString<int>(iM);
    str += ",";
    str += ToString<int>(iDry);
    return str;
}

// ---------------------------------------------------------------------------

void DairyModelDaily::AddFounderAction(const int & sId, int & aId, string & label,
    vector<int> & aScpIdx, vector<flt> & pr, vector<flt>& weights,
    ofstream & fA, ofstream & fALbl, ofstream & fTransP,
    ofstream & fAWeight)
{
    aScpIdx.clear();
    aScpIdx.push_back(2);
    aScpIdx.push_back(0);
    weights.assign(3,0);
    pr.assign(1,1);
    label = "Dummy";
    AddActionBinary(sId, aId, label, aScpIdx, pr, weights, fA, fALbl,
        fTransP, fAWeight);
}

// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
