#include "bodyweight.hh"

// ----------------------------------------------------------------------------

BodyWeight::BodyWeight() {
    udvegt=680;
    gompM=2.5483;
    gompN=0.00314;
    foster=0.02;
    fostervaegt=0.1133;
    L0=3.0;  // start bcs
    T=70;
    LT=2.6;  // minimum bcs
    Lnext=3.0; // end bcs
    gest=282;
    maxLLoss=0.0324;
    kgBCSUnitPrStdBW=0.090;
    cowAgeA=493.1604;
    cowAgeB=379.3949;
    sfuFetus=0.03647;
    sfuECM=0.4;
    feStdBWG=4.0;
    feKgBcsGainA=0.48570;
    feKgBcsGainB=1.38570;
    feKgBcsLossA=0.47140;
    feKgBcsLossB=1.08570;
    sfuDry=7;
}

// ---------------------------------------------------------------------------

BodyWeight::BodyWeight(flt udvegt, flt gompM, flt gompN,
        flt foster, flt fostervaegt,
        flt L0, flt T, flt LT, flt Lnext,
        flt maxLLoss, flt kgBCSUnitPrStdBW, flt gest,
        flt cowAgeA, flt cowAgeB,
        flt sfuFetus,
        flt sfuECM,
        flt feStdBWG,
        flt feKgBcsGainA,
        flt feKgBcsGainB,
        flt feKgBcsLossA,
        flt feKgBcsLossB
        ) {
    this->udvegt=udvegt;
    this->gompM=gompM;
    this->gompN=gompN;
    this->foster=foster;
    this->fostervaegt=fostervaegt;
    this->L0=L0;  // start bcs
    this->T=T;
    this->LT=LT;  // minimum bcs
    this->Lnext=Lnext; // end bcs
    this->maxLLoss=maxLLoss;
    this->kgBCSUnitPrStdBW=kgBCSUnitPrStdBW;
    this->gest=gest;
    this->cowAgeA=cowAgeA;
    this->cowAgeB=cowAgeB;
    this->sfuFetus=sfuFetus;
    this->sfuECM=sfuECM;
    this->feStdBWG=feStdBWG;
    this->feKgBcsGainA=feKgBcsGainA;
    this->feKgBcsGainB=feKgBcsGainB;
    this->feKgBcsLossA=feKgBcsLossA;
    this->feKgBcsLossB=feKgBcsLossB;
}

// ---------------------------------------------------------------------------

flt BodyWeight::GetBCS(int dfc, int con) {
    Rprintf("BCS: dfc:%d, con:%d\n",dfc,con);
    flt TT = T;
    flt a,b,c;
    if (con>=TT) {
        b=2*(LT-L0)/TT;
        if (b< -maxLLoss) {
            b= -maxLLoss;
            TT=2*(LT-L0)/-maxLLoss;
        }
        if (dfc<TT) return (L0-LT)/(TT*TT)*(dfc*dfc)+2*(LT-L0)/TT*dfc+L0;
        if (dfc<=con) return(LT);
        return (Lnext-LT)/(gest*gest)*((dfc*dfc)-2*con*dfc+(con*con))+LT;
    } else {
        flt dLT=2*(TT-con)*(Lnext-LT)/(gest*gest);
        flt LTarGet=LT+(dLT*(TT-con))/2;
        flt dL0=2*(LTarGet-L0)/TT-dLT;
        if (dL0< -maxLLoss) {
            dL0= -maxLLoss;
            TT=2*(LTarGet-L0)/(dL0+dLT);
        }
        if (dfc<TT) {
            a=(dLT-dL0)/TT;
            b=dL0;
            c=L0;
        } else {
            flt dLnext=2*(Lnext-LTarGet)/gest;
            a=(dLnext-dLT)/(con+gest-TT);
            b=dLT-a*TT;
            c= LTarGet - (a/2*(TT*TT)+b*TT);
        }
        return 0.5*a*(dfc*dfc)+b*dfc+c;
    }
}

// ---------------------------------------------------------------------------

flt BodyWeight::feedingSFU(int lac, int dfc, int con, flt milk, int gest) {
    flt age = CowAge(lac);
    flt eFetusW = 0;
    if (dfc>con) eFetusW = this->energyFetus(GetFetusW(con+gest-dfc));
    Rprintf("main:%f bcs:%f bwg:%f fetus:%f milk:%f\n",
        energyMain(GetBWStd(age+dfc)),
        energyBCS(DeltaKgBCSAge((int)age, dfc, con), GetBCS(dfc,con)),
        energyBWG(GetBWStd(age+dfc+1)-GetBWStd(age+dfc)),
        eFetusW,
        energyMilk(milk));
    return
        energyMain(GetBWStd(age+dfc)) +
        energyBCS(DeltaKgBCSAge((int)age, dfc, con), GetBCS(dfc,con)) +
        energyBWG(GetBWStd(age+dfc+1)-GetBWStd(age+dfc)) +
        eFetusW +
        energyMilk(milk);
}

// ---------------------------------------------------------------------------

flt BodyWeight::energyBCS(flt dKgBCS, flt bcs) {
    if (dKgBCS<0) return dKgBCS * (feKgBcsLossA * bcs + feKgBcsLossB);
    return dKgBCS * (feKgBcsGainA * bcs + feKgBcsGainB);
}

// ---------------------------------------------------------------------------

flt BodyWeight::feedingSFU(flt bwStd, flt dKgBCS, flt bcs, flt dKgBWG,
        flt fetusW, flt milkECM) {
    return (this->energyMain(bwStd) + this->energyBCS(dKgBCS, bcs) +
    this->energyBWG(dKgBWG) + this->energyFetus(fetusW) +
    this->energyMilk(milkECM));
}

// ---------------------------------------------------------------------------

flt BodyWeight::energyBWG(flt dKgBWG) {
    return feStdBWG*dKgBWG;
}

// ---------------------------------------------------------------------------

flt BodyWeight::energyMilk(flt milkECM) {
    return milkECM*sfuECM;
}

// ---------------------------------------------------------------------------

flt BodyWeight::energyFetus(flt fetusW) {
    return fetusW*sfuFetus;
}

// ---------------------------------------------------------------------------

flt BodyWeight::energyMain(flt bwStd) {
    return 1.1*(bwStd/200 + 1.5);
    // TODO Fandt en fejl har brugt 1.1*bwStd/200 + 1.5 i artiklen!!
}

// ---------------------------------------------------------------------------

flt BodyWeight::CarcassWeight(int lac, int dfc, int con) {
    return GetBW(lac, dfc, con)*0.5;
}

// ---------------------------------------------------------------------------

flt BodyWeight::GetBWStd(int lac, int dfc) {
    return this->GetBWStd(CowAge(lac)+dfc);
}

// ---------------------------------------------------------------------------

flt BodyWeight::GetBW(int lac, int dfc, int con) {
    flt bw = GetBWStd(CowAge(lac)+dfc);
    flt bcs = GetBCS(dfc, con);
    //Rprintf("CowAge:%f bWStd:%f BCS:%f\n",CowAge(lac),bw,bcs);
    return bw + Bcs2Kg(bw, 1)*(bcs-3);
}

// ---------------------------------------------------------------------------

flt BodyWeight::DeltaKgBCSAge(int age, int dfc, int con) {
    flt dBCS = GetBCS(dfc+1, con)-GetBCS(dfc, con);
    return Bcs2Kg(GetBWStd(age+dfc), dBCS);
}

// ---------------------------------------------------------------------------

flt BodyWeight::DeltaKgBCS(int lac, int dfc, int con) {
    flt dBCS = GetBCS(dfc+1, con)-GetBCS(dfc, con);
    return Bcs2Kg(GetBWStd(lac, dfc), dBCS);
}

// ---------------------------------------------------------------------------

flt BodyWeight::CowAge(int lac) {
    return floor(cowAgeA + cowAgeB*(flt)lac);
}

// ---------------------------------------------------------------------------

flt BodyWeight::GetBWStd(flt age) {
    //Rprintf("Age (bwStd):%f ",age);
    return udvegt*exp(-gompM*exp(-gompN*age));
}

// ---------------------------------------------------------------------------

flt BodyWeight::GetFetusW(int dtc) {
    return exp(-foster*dtc)*fostervaegt*udvegt;
}

// ---------------------------------------------------------------------------

flt BodyWeight::Bcs2Kg(flt bwStd,flt bcs) {
    return (bwStd*kgBCSUnitPrStdBW)*bcs;
}

// ---------------------------------------------------------------------------




