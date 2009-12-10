#include <vector>
#include <string>
#include <fstream>
#include "basicdt.hh"
#include "dairymodel.hh"
using namespace std;

#include "Rdefines.h"


#ifdef __cplusplus
extern "C" {
#endif

static SEXP type_tag;

typedef class DairyModelDaily* DairyModelDailyPtr;

/* macro to check if ptr valid */
#define CHECK_PTR(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || \
        R_ExternalPtrTag(s) != type_tag) \
        error("External pointer not valid!"); \
} while (0)

/* Install the type tag */
SEXP DAIRY_init(void)
{
    type_tag = install("TYPE_TAG");
    return R_NilValue;
}

/* Create an DairyModelDaily object. */
SEXP DAIRY_newDairyModelDaily(SEXP rBinNames,
    SEXP rDryPeriodLth,
    SEXP rGestLth,
    SEXP rProbCalf,
    SEXP rProbBullCalf,
    SEXP rPriceCalf,
    SEXP rPriceBullCalf,
    SEXP rPriceHeifer,
    SEXP rPriceSFU,
    SEXP rPriceECM,
    SEXP rPriceCarcassKg,
    SEXP rInsemStart,
    SEXP rInsemFinish,
    SEXP rPregTestLth,
    SEXP rVS1,
    SEXP rVS2,
    //SEXP rSDry,
    SEXP rPrICDry,
    SEXP rPrIC,
    SEXP rPrM,
    SEXP rPrM2A,
    SEXP rPrPregT,
    SEXP rYield,
    SEXP rVA2M,
    SEXP rIAZero,
    SEXP rFun)
{
    // set binary names
    vector<string> binName;
    for(int i=0; i<7; ++i) binName.push_back((string)CHAR(STRING_ELT(rBinNames, i)));

    DairyModelDaily * pMod = new DairyModelDaily(
        binName,
        INTEGER_POINTER(rDryPeriodLth)[0],
        INTEGER_POINTER(rGestLth)[0],
        NUMERIC_POINTER(rProbCalf)[0],
        NUMERIC_POINTER(rProbBullCalf)[0],
        NUMERIC_POINTER(rPriceCalf)[0],
        NUMERIC_POINTER(rPriceBullCalf)[0],
        NUMERIC_POINTER(rPriceHeifer)[0],
        NUMERIC_POINTER(rPriceSFU)[0],
        NUMERIC_POINTER(rPriceECM)[0],
        NUMERIC_POINTER(rPriceCarcassKg)[0],
        INTEGER_POINTER(rInsemStart)[0],
        INTEGER_POINTER(rInsemFinish)[0],
        INTEGER_POINTER(rPregTestLth)[0],
        INTEGER_POINTER(rIAZero)[0]
    );

    // assign vectors ------------------
    pMod->SetVS1(INTEGER_POINTER(rVS1),GET_LENGTH(rVS1));
    pMod->SetVS2(INTEGER_POINTER(rVS2),GET_LENGTH(rVS2));
    pMod->SetVA2M(INTEGER_POINTER(rVA2M),GET_LENGTH(rVA2M));
    //pMod.SetSDry(INTEGER_POINTER(rSDry),GET_LENGTH(rSDry));
    pMod->SetPrICDry(NUMERIC_POINTER(rPrICDry),GET_LENGTH(rPrICDry));

    // assign prIC (enclose to be sure)
    {SEXP dim;
    int q,l,i,j,k;
    flt *p;
    PROTECT(dim = GET_DIM(rPrIC));
    int * pDim = INTEGER_POINTER(AS_INTEGER(dim));
    vector< vector< vector< flt > > > * pPrIC = pMod->GetPrICPtr();
    pPrIC->resize(pDim[0]);     // alloc mem for vector
    for (i=0;i<pDim[0];i++){
        (*pPrIC)[i].resize(pDim[1]);
        for (j=0;j<pDim[1];j++) (*pPrIC)[i][j].resize(pDim[2]);
    }
    //Rprintf("prIC : Dim %d %d %d\n",pDim[0],pDim[1],pDim[2]);
    for(k=0;k<pDim[2];++k) {
        for(j=0;j<pDim[1];++j) {
            for(i=0;i<pDim[0];++i) {
                q = i + j*pDim[0] + k*(pDim[0]*pDim[1]);
                l = GET_LENGTH(VECTOR_ELT(rPrIC, q));
                p = NUMERIC_POINTER(VECTOR_ELT(rPrIC, q));
                //Rprintf("q: %d ",q);Rprintf("lth: %d ",l);Rprintf("v = %f\n",p[0]);
                (*pPrIC)[i][j][k] = p[0];
            }
        }
    }
    UNPROTECT(1);}

    // assign prM (enclose to be sure)
    {SEXP dim;
    int q,l,i,j,k,r;
    flt *p;
    PROTECT(dim = GET_DIM(rPrM));
    int * pDim = INTEGER_POINTER(AS_INTEGER(dim));
    vector< vector< vector< vector<flt> > > > * pPrM = pMod->GetPrMPtr();
    pPrM->resize(pDim[0]);
    for (i=0;i<pDim[0];i++){
        (*pPrM)[i].resize(pDim[1]);
        for (j=0;j<pDim[1];j++) (*pPrM)[i][j].resize(pDim[2]);
    }
    //Rprintf("prM : Dim %d %d %d\n",pDim[0],pDim[1],pDim[2]);
    for(k=0;k<pDim[2];++k) {
        for(j=0;j<pDim[1];++j) {
            for(i=0;i<pDim[0];++i) {
                q = i + j*pDim[0] + k*(pDim[0]*pDim[1]);
                l = GET_LENGTH(VECTOR_ELT(rPrM, q));
                p = NUMERIC_POINTER(VECTOR_ELT(rPrM, q));
                //Rprintf("q: %d ",q);Rprintf("lth: %d v = ",l);
                (*pPrM)[i][j][k].resize(l);
                for (r=0;r<l;++r) {
                    (*pPrM)[i][j][k][r] = p[r];
                    //Rprintf("%f ",p[r]);
                }
                //Rprintf("\n");
            }
        }
    }
    UNPROTECT(1);}

    // assign prM2A (enclose to be sure)
    {
    int i,l,r,lth;
    flt *p;
    lth = GET_LENGTH(rPrM2A);
    vector< vector<flt> >  * pPrM2A = pMod->GetPrM2APtr();
    //Rprintf("prM2A : Dim %d\n",pPrM2A->size());
    pPrM2A->resize(lth);
    //Rprintf("prM2A : Dim %d\n",lth);
    for(i=0;i<lth;++i) {
        l = GET_LENGTH(VECTOR_ELT(rPrM2A, i));
        p = NUMERIC_POINTER(VECTOR_ELT(rPrM2A, i));
        //Rprintf("i: %d ",i);Rprintf("lth: %d v = ",l);
        (*pPrM2A)[i].resize(l);
        for (r=0;r<l;++r) {
            (*pPrM2A)[i][r] = p[r];
            //Rprintf("%f ",p[r]);
        }
        //Rprintf("\n");
    }
    }

    // assign prPregT (enclose to be sure)
    {SEXP dim;
    int q,l,i,j;
    flt *p;
    PROTECT(dim = GET_DIM(rPrPregT));
    int * pDim = INTEGER_POINTER(AS_INTEGER(dim));
    vector< vector<flt> > * pPrPregT = pMod->GetPrPregTPtr();
    pPrPregT->resize(pDim[0]);     // alloc mem for vector
    for (i=0;i<pDim[0];i++){
        (*pPrPregT)[i].resize(pDim[1]);
    }
    //Rprintf("prPregT : Dim %d %d\n",pDim[0],pDim[1]);
    for(j=0;j<pDim[1];++j) {
        for(i=0;i<pDim[0];++i) {
            q = i + j*pDim[0];
            l = GET_LENGTH(VECTOR_ELT(rPrPregT, q));
            p = NUMERIC_POINTER(VECTOR_ELT(rPrPregT, q));
            //Rprintf("q: %d ",q);Rprintf("lth: %d ",l);Rprintf("v = %f\n",p[0]);
            (*pPrPregT)[i][j] = p[0];
        }
    }
    UNPROTECT(1);
    //Rprintf("prPreg[1,75]=%f\n",(*pPrPregT)[1][75]);
    }

    // assign yield (enclose to be sure)
    {SEXP dim;
    int q,l,i,j,k;
    flt *p;
    PROTECT(dim = GET_DIM(rYield));
    int * pDim = INTEGER_POINTER(AS_INTEGER(dim));
    vector< vector< vector< flt > > > * pYield = pMod->GetYieldPtr();
    pYield->resize(pDim[0]);     // alloc mem for vector
    for (i=0;i<pDim[0];i++){
        (*pYield)[i].resize(pDim[1]);
        for (j=0;j<pDim[1];j++) (*pYield)[i][j].resize(pDim[2]);
    }
    //Rprintf("yield : Dim %d %d %d\n",pDim[0],pDim[1],pDim[2]);
    for(k=0;k<pDim[2];++k) {
        for(j=0;j<pDim[1];++j) {
            for(i=0;i<pDim[0];++i) {
                q = i + j*pDim[0] + k*(pDim[0]*pDim[1]);
                l = GET_LENGTH(VECTOR_ELT(rYield, q));
                p = NUMERIC_POINTER(VECTOR_ELT(rYield, q));
                //Rprintf("q: %d ",q);Rprintf("lth: %d ",l);Rprintf("v = %f\n",p[0]);
                (*pYield)[i][j][k] = p[0];
            }
        }
    }
    UNPROTECT(1);}



	/*vector< vector<string> > * pMatLab1 = pMod->GetMatLab1Ptr();    ///< Matrix containing the labels of level 1.  matLab1[i][j] = label of state j at stage i.
	int * matLab1Dims = INTEGER_POINTER(AS_INTEGER(GET_DIM(rMatLab1)));
	pMatLab1->assign(matLab1Dims[0], vector<string>() );  // allocate mem for d1
	for (int d1=0; d1<matLab1Dims[0]; ++d1) {
	    (*pMatLab1)[d1].resize(matLab1Dims[1]);
        for (int s1=0; s1<matLab1Dims[1]; ++s1) (*pMatLab1)[d1][s1] = (string)CHAR(STRING_ELT(rMatLab1, d1+matLab1Dims[0]*s1));
    }

	vector< vector<string> > * pMatLab2 = pMod->GetMatLab2Ptr();    ///< Matrix containing the labels of level 2. matLab2[i][j] = label of state j at stage i.
	int * matLab2Dims = INTEGER_POINTER(AS_INTEGER(GET_DIM(rMatLab2)));
	pMatLab2->assign(matLab2Dims[0], vector<string>() );
	for (int d2=0; d2<matLab2Dims[0]; ++d2) {
	    (*pMatLab2)[d2].resize(matLab2Dims[1]);
        for (int s2=0; s2<matLab2Dims[1]; ++s2) (*pMatLab2)[d2][s2] = (string)CHAR(STRING_ELT(rMatLab2, d2+matLab2Dims[0]*s2));
    }*/
    // End assign vectors ------------------

    /*

    pMod.SetPrIC((string)CHAR(STRING_ELT(rPrICFile, 0)));
    pMod.SetPrM((string)CHAR(STRING_ELT(rPrMFile, 0)));
    pMod.SetYield((string)CHAR(STRING_ELT(rYieldFile, 0)));
    pMod.SetPrPregT((string)CHAR(STRING_ELT(rPrPregTFile, 0)));

*/
    if (pMod == NULL)
        return R_NilValue;
    else {
        SEXP val = R_MakeExternalPtr(pMod, type_tag, R_NilValue);
        R_RegisterFinalizer(val, rFun);
        return val;
    }
}

/* Remove the DairyModelDaily object */
SEXP DAIRY_deleteDairyModelDaily(SEXP s)
{
    CHECK_PTR(s);
    cout << "test -1\n" << flush;
    DairyModelDailyPtr p = (DairyModelDailyPtr)R_ExternalPtrAddr(s);
    if (p != NULL) {
        cout << "test 0\n" << flush;
        delete p;
        cout << "test 1\n" << flush;
        R_ClearExternalPtr(s);
        cout << "test 2\n" << flush;
    }
    cout << "test 3\n" << flush;
    return R_NilValue;
}

/* Create binary files for the states. */
SEXP DAIRY_genStates(SEXP ptr)
{
	CHECK_PTR(ptr);
	DairyModelDailyPtr p = (DairyModelDailyPtr)R_ExternalPtrAddr(ptr);
	if (p == NULL) error("pointer is NULL");
	p->GenerateStatesBinary();
	return R_NilValue;
}

/* Create binary files for the actions. */
SEXP DAIRY_genActions(SEXP ptr)
{
	CHECK_PTR(ptr);
	DairyModelDailyPtr p = (DairyModelDailyPtr)R_ExternalPtrAddr(ptr);
	if (p == NULL) error("pointer is NULL");
	p->GenerateActionsBinary();
	return R_NilValue;
}






























/** Write state idx to csv. */
/*void state(ofstream & outf, const vector<int> & idx, const string & label, int & sId) {
	outf << sId << ",";
	for (int i=0;i<idx.size();++i) {
		if (idx[i]>=0) outf << idx[i];
		outf << ",";
	}
	if (label.compare("") != 0) outf << '"' << label << '"';
	outf << endl;
	sId++;
}*/

/** Create all states of the HMDP.
	\param vS1 Vector containing the number of level 1 states (size=maxLac+2).
	\param vS2 Vector containing the number of level 2 states (size=maxDfc+1).
	\param matL1 Matrix containing the labels of level 1.
	\param matL2 Matrix containing the labels of level 2.
 */
/*SEXP DAIRY_GenStates2Csv(SEXP vS1, SEXP vS2, SEXP matL1, SEXP matL2,
	SEXP fileName)
{
	int maxLac = GET_LENGTH(vS1)-2;
	int maxDfc = GET_LENGTH(vS2)-1;
	string label;
	int lastS2;
	int sId=0;
	ofstream outf(CHAR(STRING_ELT(fileName, 0)));
	vector<int> idx(5,-1); // containing the stage, state or action idx's (-1 represent NA)
	int * pVS1=INTEGER_POINTER(vS1);
	int * pVS2=INTEGER_POINTER(vS2);
	int * matL1Dims = INTEGER_POINTER(AS_INTEGER(GET_DIM(matL1)));
	int * matL2Dims = INTEGER_POINTER(AS_INTEGER(GET_DIM(matL2)));
	int nRMatL1,nRMatL2;
	nRMatL1 = matL1Dims[0];
	nRMatL2 = matL2Dims[0];
	PROTECT(matL1 = AS_CHARACTER(matL1));
	PROTECT(matL2 = AS_CHARACTER(matL2));
	outf << "sId,d1,s1,a1,d2,s2,label\n";
	state(outf, idx,"Dummy",sId);  // founder state
	for (int d1=0;d1<=maxLac+1;++d1) {
		idx[0]=d1;
		for (int s1=0; s1<pVS1[d1]; ++s1) {
			idx[1]=s1;
			label = (string)CHAR(STRING_ELT(matL1, d1+nRMatL1*s1));
			state(outf, idx, label,sId);
			//Rprintf("d1 s1 states %d %d %d\n",d1,s1,pVS1[d1]);
			if (s1==0 & d1>0 & d1<maxLac+1) { // if subprocess
				for (int a1=0; a1<1; ++a1) {
					idx[2]=a1;
					for (int d2=0; d2<=maxDfc; ++d2) {
						idx[3]=d2;
						lastS2=pVS2[d2]-1;
						for (int s2=0; s2<=lastS2; ++s2) {
							//Rprintf("d2 s2 states %d %d %d\n",d2,s2,pVS2[d2]);
							idx[4]=s2;
							label = (string)CHAR(STRING_ELT(matL2, d2+nRMatL2*s2));
							state(outf, idx, label,sId);
						}
						idx[4]=-1;
					}
					idx[3]=-1;
				}
			}
			idx[2]=-1;
		}
		idx[1]=-1;
	}
	UNPROTECT(2);
	outf.close();
	return R_NilValue;
}*/



// --------------------------
// Test functions

/** Load a 3-dim array 'a' containing a vector for each entry and store it in a
stl 4-dim vector.

Let d = dim(a) and consider entry a[[i,j,k]] which corresponds to entry q in the
list, i.e. a[[q]]. Then we have that:

q = i + j*d[1] + k*(d[1]*d[2]) - (d[1] + d[1]*d[2]) (index from 1 (R))

q = i + j*d[1] + k*(d[1]*d[2]) (index from 0 (C++))
*/
SEXP DAIRY_PrintValue(SEXP ary) {
    SEXP dim;
    int q,l,i,j,k,r;
    flt *p;
    PROTECT(dim = GET_DIM(ary));
    int * pDim = INTEGER_POINTER(AS_INTEGER(dim));
    vector< vector< vector< vector<flt> > > > v;
    v.resize(pDim[0]);
    for (i=0;i<pDim[0];i++){
        v[i].resize(pDim[1]);
        for (j=0;j<pDim[1];j++) v[i][j].resize(pDim[2]);
    }
    Rprintf("Array dim: %d %d %d\n",pDim[0],pDim[1],pDim[2]);
    for(k=0;k<pDim[2];++k) {
        for(j=0;j<pDim[1];++j) {
            for(i=0;i<pDim[0];++i) {
                q = i + j*pDim[0] + k*(pDim[0]*pDim[1]);
                l = GET_LENGTH(VECTOR_ELT(ary, q));
                p = NUMERIC_POINTER(VECTOR_ELT(ary, q));
                Rprintf("q: %d ",q);Rprintf("lth: %d v = ",l);
                v[i][j][k].resize(l);
                for (r=0;r<l;++r) {
                    v[i][j][k][r] = p[r];
                    Rprintf("%f ",p[r]);
                }
                Rprintf("\n");
            }
        }
    }
    UNPROTECT(1);
    return R_NilValue;
}

#ifdef __cplusplus
}
#endif


