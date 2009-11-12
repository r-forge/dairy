#ifndef BODYWEIGHT_HPP
#define BODYWEIGHT_HPP

#include <cmath>
using namespace std;
#include "Rdefines.h"
//#include "basicdt.hh"
typedef double flt;                 ///< A floating number datatype.
//-----------------------------------------------------------------------------

/**
 * Methods related to feeding and the body weight of the cow.
 * @author Lars Relund
 */
class BodyWeight
{
public:
    friend class DairyModelDaily;     // so can access private members

	BodyWeight();

	BodyWeight(flt udvegt, flt gompM, flt gompN,
        flt foster, flt fostervaegt,
        flt L0, flt T, flt LT, flt Lnext, flt maxLLoss, flt kgBCSUnitPrStdBW,
        flt gest,
        flt cowAgeA, flt cowAgeB,
        flt sfuFetus, flt sfuECM,
        flt feStdBWG, flt feKgBcsGainA, flt feKgBcsGainB,
        flt feKgBcsLossA, flt feKgBcsLossB);

	/**
	 *
	 * @param dtc Days to calving
	 * @return Fetus weight.
	 */
	flt GetFetusW(int dtc);


	/**
	 * Gompertz curve (body weight at BCS=3, fetus not included)
	 * @param age Age of the cow in days
	 * @return
	 */
	flt GetBWStd(flt age);

	/**
	 * BCS (1-5 scale).
	 * @param dfc
	 * @param con Day of conception.
	 * @return
	 */
	flt GetBCS(int dfc, int con);

	flt CowAge(int lac);

	/**
	 * Daily change in BCS converted to kg.
	 * @param lac
	 * @param dfc
	 * @param con Day of conception
	 * @return
	 */
	flt DeltaKgBCS(int lac, int dfc, int con);

	/**
	 * Daily change in BCS converted to kg.
	 * @param age Age at start of parity
	 * @param dfc
	 * @param con Day of conception
	 * @return
	 */
	flt DeltaKgBCSAge(int age, int dfc, int con);

	/**
	 * Return the BW at day dfc assuming conception at day con. Note that con = 0 indicate
	 * that the cow still not pregnant.
	 * @param dfc Days from calving
	 * @param con Day of conception
	 * @return The BW.
	 */
	flt GetBW(int lac, int dfc, int con);

	/**
	 * Return the std BW at day dfc.
	 * @param lac Lactation number.
	 * @param dfc Days from calving
	 * @return The standard body weight.
	 */
	flt GetBWStd(int lac, int dfc);

	/**
	 * Return the price of the cow when culled at day dfc
	 * @param dfc Days from calving
	 * @param con Day of conception
	 * @return The price.
	 */
	flt CarcassWeight(int lac, int dfc, int con);

	/**
	 *
	 * @param bwStd Standard BW in kg.
	 * @return
	 */
	flt energyMain(flt bwStd);

	/**
	 *
	 * @param fetusW Weight of fetus.
	 * @return
	 */
	flt energyFetus(flt fetusW);
	/**
	 *
	 * @param milkECM ECM milk in Kg.
	 * @return
	 */
	flt energyMilk(flt milkECM);

	/**
	 * Energy needs for body weight gain.
	 * @param dKgBWG Daily change in BWG in Kg.
	 * @return
	 */
	flt energyBWG(flt dKgBWG);

	/**
	 *
	 * @param dKgBCS Daily change in BCS in Kg.
	 * @param bcs BCS at start of day.
	 * @return
	 */
	flt energyBCS(flt dKgBCS, flt bcs);

	/**
	 * Total daily energy needs in SFU
	 * @param bwStd Standard BW in kg.
	 * @param dKgBCS Daily change in BCS in Kg.
	 * @param bcs BCS at start of day.
	 * @param dKgBWG Daily change in BWG in Kg.
	 * @param fetusW fetusW Weight of fetus.
	 * @param milkECM ECM milk in Kg.
	 * @return
	 */
	flt feedingSFU(flt bwStd, flt dKgBCS, flt bcs, flt dKgBWG,
			flt fetusW, flt milkECM);

	/**
	 * Total daily energy needs in SFU
	 * @param lac Lactation.
	 * @param dfc Days from calving.
	 * @param con Day of conception (if greater than dfc assume not pregnant).
	 * @param milk Yild in kg ECM.
	 * @param gest Length of gestration.
	 * @return
	 */
	flt feedingSFU(int lac, int dfc, int con, flt milk, int gest);

private:

	/**
	 * Covert BCS to Kg.
	 * @param bwStd Std. body weight.
	 * @param bcs BCS of the cow.
	 * @return
	 */
	flt Bcs2Kg(flt bwStd,flt bcs);

private:

	flt udvegt;
	flt gompM;
	flt gompN;
	flt foster;
	flt fostervaegt;

	flt L0;  // start bcs
	flt T;
	flt LT;  // minimum bcs
	flt Lnext; // end bcs
	flt gest;
	flt maxLLoss;
	flt kgBCSUnitPrStdBW;

	flt cowAgeA, cowAgeB;	// parameters in cowAge method

	flt sfuFetus;
	flt sfuECM;
	flt feStdBWG;
	flt feKgBcsGainA;
	flt feKgBcsGainB;
	flt feKgBcsLossA;
	flt feKgBcsLossB;
	flt sfuDry;
};

//-----------------------------------------------------------------------------

#endif
