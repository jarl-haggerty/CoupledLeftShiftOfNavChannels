#include "HodgkinHuxley.hpp"
#include <math.h>
#include <algorithm>
#include <sstream>

using namespace std;

namespace Jarl {
  const Scalar potassiumDissociationConstant = 3.5;
  const Scalar sodiumDissociationConstant = 10;
  const Scalar gasConstant = 8.314472;
  const Scalar faradayConstant = 96485.3399;

  Scalar alphaN(const Scalar potential) {
    if(potential == -55) {
      return .1;
    } else {
      return .01 * (potential + 55) / (1 - exp(-(potential + 55) / 10));
    }
  }

  Scalar betaN(const Scalar potential) {
    return .125 * exp(-(potential + 65) / 80);
  }

  Scalar infinityN(const Scalar potential) {
    return alphaN(potential) / (alphaN(potential) + betaN(potential));
  }

  Scalar derivativeN(const Scalar potential, const Scalar n) {
    return alphaN(potential) * (1 - n) - betaN(potential) * n;
  }

  Scalar alphaM(const Scalar potential) {
    if(potential == -40) {
      return 1;
    } else {
      return .1 * (potential + 40) / (1 - exp(-(potential + 40) / 10));
    }
  }

  Scalar betaM(const Scalar potential) {
    return 4 * exp(-(potential + 65) / 18);
  }

  Scalar infinityM(const Scalar potential) {
    return alphaM(potential) / (alphaM(potential) + betaM(potential));
  }

  Scalar derivativeM(const Scalar potential, const Scalar m) {
    return alphaM(potential) * (1 - m) - betaM(potential) * m;
  }

  Scalar alphaH(const Scalar potential) {
    return .07 * exp(-(potential + 65) / 20);
  }

  Scalar betaH(const Scalar potential) {
    return 1 / (1 + exp(-(potential + 35) / 10));
  }

  Scalar infinityH(const Scalar potential) {
    return alphaH(potential) / (alphaH(potential) + betaH(potential));
  }

  Scalar derivativeH(const Scalar potential, const Scalar h) {
    return alphaH(potential) * (1 - h) - betaH(potential) * h;
  }

  HodgkinHuxley::Settings::Settings() {
    maxPumpCurrent = 0;
    potassiumLeakConductance = 0;
    sodiumLeakConductance = 0;
    blebbing = 0;
    leftShift = 0;
    stimulation = 0;
  }

  class HodgkinHuxley::Private {
  public:
    Scalar potential;
    Scalar capacitance;
    Scalar leakConductance;
    Scalar leakReversalPotential;
    Scalar potassiumConductance;
    Scalar potassiumReversalPotential;
    Scalar sodiumConductance;
    Scalar sodiumReversalPotential;
    Scalar maxPumpCurrent;
    Scalar innerPotassiumConcentration;
    Scalar outerPotassiumConcentration;
    Scalar innerSodiumConcentration;
    Scalar outerSodiumConcentration;
    Scalar surfaceArea;
    Scalar innerVolume;
    Scalar outerVolume;
    Scalar temperature;
    Scalar threshold;
    Scalar potassiumLeakConductance;
    Scalar sodiumLeakConductance;
    Scalar blebbing;
    Scalar leftShift;
    Scalar n;
    Scalar m;
    Scalar h;
    Scalar blebbedN;
    Scalar blebbedM;
    Scalar blebbedH;
    Scalar lastPotential;
    Scalar stimulation;

    Private(const Settings& settings) {
      potential = settings.potential;
      capacitance = settings.capacitance;
      leakConductance = settings.leakConductance;
      leakReversalPotential = settings.leakReversalPotential;
      potassiumConductance = settings.potassiumConductance;
      potassiumReversalPotential = settings.potassiumReversalPotential;
      sodiumConductance = settings.sodiumConductance;
      sodiumReversalPotential = settings.sodiumReversalPotential;
      maxPumpCurrent = settings.maxPumpCurrent;
      innerPotassiumConcentration = settings.innerPotassiumConcentration;
      outerPotassiumConcentration = settings.outerPotassiumConcentration;
      innerSodiumConcentration = settings.innerSodiumConcentration;
      outerSodiumConcentration = settings.outerSodiumConcentration;
      surfaceArea = settings.surfaceArea;
      innerVolume = settings.innerVolume;
      outerVolume = settings.outerVolume;
      temperature = settings.temperature;
      threshold = settings.threshold;
      potassiumLeakConductance = settings.potassiumLeakConductance;
      sodiumLeakConductance = settings.sodiumLeakConductance;
      blebbing = settings.blebbing;
      leftShift = settings.leftShift;
      stimulation = settings.stimulation;

      n = infinityN(potential);
      m = infinityM(potential);
      h = infinityH(potential);
      blebbedN = infinityN(potential + leftShift);
      blebbedM = infinityM(potential + leftShift);
      blebbedH = infinityH(potential + leftShift);
    }

    Scalar getLeakCurrent(const Scalar potential) {
      return leakConductance * (potential - leakReversalPotential);
    }


    Scalar getPotassiumCurrent(const Scalar potential, const Scalar n,
			       Scalar reversalPotential) {
      return potassiumConductance * n * n * n * n
	* (potential - reversalPotential);
    }

    Scalar getPotassiumPumpLeakCurrent(const Scalar potential,
				       Scalar reversalPotential) {
      return potassiumLeakConductance * (potential - reversalPotential);
    }

    Scalar getTotalPotassiumCurrent(const Scalar potential, const Scalar n,
				    Scalar reversalPotential, const Scalar pumpBaseCurrent) {
      return getPotassiumCurrent(potential, n, reversalPotential) - 2
	* pumpBaseCurrent
	+ getPotassiumPumpLeakCurrent(potential, reversalPotential);
    }

    Scalar getSodiumCurrent(const Scalar potential, const Scalar m, const Scalar h, const Scalar blebbedM, const Scalar blebbedH,
			    Scalar reversalPotential) {
      return sodiumConductance * (m * m * m * h * (1 - blebbing) + blebbedM * blebbedM * blebbedM * blebbedH * blebbing)
	* (potential - reversalPotential);
    }

    Scalar getSodiumPumpLeakCurrent(const Scalar potential,
				    Scalar reversalPotential) {
      return sodiumLeakConductance * (potential - reversalPotential);
    }

    Scalar getTotalSodiumCurrent(const Scalar potential, const Scalar m, const Scalar h, const Scalar blebbedM, const Scalar blebbedH,
				 Scalar reversalPotential, const Scalar pumpBaseCurrent) {
      return getSodiumCurrent(potential, m, h, blebbedM, blebbedH, reversalPotential) + 3
	* pumpBaseCurrent
	+ getSodiumPumpLeakCurrent(potential, reversalPotential);
    }

    Scalar getPumpBaseCurrent(const Scalar outerPotassiumConcentration,
			      Scalar innerSodiumConcentration) {
      Scalar potassiumTemp = 1 + potassiumDissociationConstant/outerPotassiumConcentration;
      Scalar sodiumTemp = 1 + sodiumDissociationConstant/innerSodiumConcentration;
      return maxPumpCurrent/(potassiumTemp*potassiumTemp*sodiumTemp*sodiumTemp*sodiumTemp);
    }

    Scalar derivativeInnerConcentration(const Scalar current) {
      return -1e-6 * current * surfaceArea / faradayConstant / innerVolume;
    }
	
    Scalar derivativeOuterConcentration(const Scalar current) {
      return 1e-6 * current * surfaceArea / faradayConstant / outerVolume;
    }

    Scalar calculateReversalPotential(const Scalar innerConcentration, const Scalar outerConcentration) {
      return -gasConstant * temperature / faradayConstant * 1000 * log(innerConcentration / outerConcentration);
    }
  };

  HodgkinHuxley::HodgkinHuxley(const Settings& settings) : priv(new Private(settings)) {}

  HodgkinHuxley::~HodgkinHuxley() {
    delete priv;
  }

  Scalar HodgkinHuxley::simulate(const Scalar time, const Scalar limit, const bool force) {
    Scalar pumpBaseCurrent;
    Scalar potassiumReversalPotential, sodiumReversalPotential;
    Scalar potassiumCurrent, sodiumCurrent, leakCurrent, totalCurrent;
    Scalar nK1, nK2, nK3, nK4;
    Scalar mK1, mK2, mK3, mK4;
    Scalar hK1, hK2, hK3, hK4;
    Scalar blebbedNK1, blebbedNK2, blebbedNK3, blebbedNK4;
    Scalar blebbedMK1, blebbedMK2, blebbedMK3, blebbedMK4;
    Scalar blebbedHK1, blebbedHK2, blebbedHK3, blebbedHK4;
    Scalar innerPotassiumConcentrationK1, innerPotassiumConcentrationK2, innerPotassiumConcentrationK3, innerPotassiumConcentrationK4;
    Scalar outerPotassiumConcentrationK1, outerPotassiumConcentrationK2, outerPotassiumConcentrationK3, outerPotassiumConcentrationK4;
    Scalar innerSodiumConcentrationK1, innerSodiumConcentrationK2, innerSodiumConcentrationK3, innerSodiumConcentrationK4;
    Scalar outerSodiumConcentrationK1, outerSodiumConcentrationK2, outerSodiumConcentrationK3, outerSodiumConcentrationK4;
    Scalar potentialK1, potentialK2, potentialK3, potentialK4;

    potassiumReversalPotential = priv->calculateReversalPotential(priv->innerPotassiumConcentration,
								  priv->outerPotassiumConcentration);
    sodiumReversalPotential = priv->calculateReversalPotential(priv->innerSodiumConcentration, 
							       priv->outerSodiumConcentration);
    pumpBaseCurrent = priv->getPumpBaseCurrent(priv->outerPotassiumConcentration,
					       priv->innerSodiumConcentration);
    potassiumCurrent = priv->getTotalPotassiumCurrent(priv->potential, priv->n,
						      potassiumReversalPotential, pumpBaseCurrent);
    sodiumCurrent = priv->getTotalSodiumCurrent(priv->potential, priv->m, priv->h, priv->blebbedM, priv->blebbedH,
						sodiumReversalPotential, pumpBaseCurrent);
    leakCurrent = priv->getLeakCurrent(priv->potential);
    totalCurrent = potassiumCurrent + sodiumCurrent + leakCurrent + priv->stimulation;
    nK1 = time * derivativeN(priv->potential, priv->n);
    mK1 = time * derivativeM(priv->potential, priv->m);
    hK1 = time * derivativeH(priv->potential, priv->h);
    //blebbedNK1 = time * derivativeN(priv->potential + priv->leftShift, priv->blebbedN);
    blebbedMK1 = time * derivativeM(priv->potential + priv->leftShift, priv->blebbedM);
    blebbedHK1 = time * derivativeH(priv->potential + priv->leftShift, priv->blebbedH);
    innerPotassiumConcentrationK1 = time
      * priv->derivativeInnerConcentration(potassiumCurrent);
    outerPotassiumConcentrationK1 = time
      * priv->derivativeOuterConcentration(potassiumCurrent);
    innerSodiumConcentrationK1 = time
      * priv->derivativeInnerConcentration(sodiumCurrent);
    outerSodiumConcentrationK1 = time
      * priv->derivativeOuterConcentration(sodiumCurrent);
    potentialK1 = time
      * (-totalCurrent / priv->capacitance);
		
    potassiumReversalPotential = priv->calculateReversalPotential(priv->innerPotassiumConcentration + innerPotassiumConcentrationK1/2,
								  priv->outerPotassiumConcentration + outerPotassiumConcentrationK1/2);
    sodiumReversalPotential = priv->calculateReversalPotential(priv->innerSodiumConcentration + innerSodiumConcentrationK1/2, 
							       priv->outerSodiumConcentration + outerSodiumConcentrationK1/2);
    pumpBaseCurrent = priv->getPumpBaseCurrent(priv->outerPotassiumConcentration + outerPotassiumConcentrationK1/2,
					       priv->innerSodiumConcentration + innerSodiumConcentrationK1/2);
    potassiumCurrent = priv->getTotalPotassiumCurrent(priv->potential + potentialK1/2, priv->n + nK1/2,
						      potassiumReversalPotential, pumpBaseCurrent);
    sodiumCurrent = priv->getTotalSodiumCurrent(priv->potential + potentialK1/2, 
						priv->m + mK1/2, priv->h + hK1/2, 
						priv->blebbedM + blebbedMK1/2, priv->blebbedH + blebbedHK1/2,
						sodiumReversalPotential, pumpBaseCurrent);
    leakCurrent = priv->getLeakCurrent(priv->potential + potentialK1/2);
    totalCurrent = potassiumCurrent + sodiumCurrent + leakCurrent + priv->stimulation;
    nK2 = time * derivativeN(priv->potential + potentialK1/2, priv->n + nK1/2);
    mK2 = time * derivativeM(priv->potential + potentialK1/2, priv->m + mK1/2);
    hK2 = time * derivativeH(priv->potential + potentialK1/2, priv->h + hK1/2);
    //blebbedNK2 = time * derivativeN(priv->potential + priv->leftShift + potentialK1/2, priv->blebbedN + blebbedNK1/2);
    blebbedMK2 = time * derivativeM(priv->potential + priv->leftShift + potentialK1/2, priv->blebbedM + blebbedMK1/2);
    blebbedHK2 = time * derivativeH(priv->potential + priv->leftShift + potentialK1/2, priv->blebbedH + blebbedHK1/2);
    innerPotassiumConcentrationK2 = time
      * priv->derivativeInnerConcentration(potassiumCurrent);
    outerPotassiumConcentrationK2 = time
      * priv->derivativeOuterConcentration(potassiumCurrent);
    innerSodiumConcentrationK2 = time
      * priv->derivativeInnerConcentration(sodiumCurrent);
    outerSodiumConcentrationK2 = time
      * priv->derivativeOuterConcentration(sodiumCurrent);
    potentialK2 = time
      * (-totalCurrent / priv->capacitance);

    potassiumReversalPotential = priv->calculateReversalPotential(priv->innerPotassiumConcentration + innerPotassiumConcentrationK2/2,
								  priv->outerPotassiumConcentration + outerPotassiumConcentrationK2/2);
    sodiumReversalPotential = priv->calculateReversalPotential(priv->innerSodiumConcentration + innerSodiumConcentrationK2/2, 
							       priv->outerSodiumConcentration + outerSodiumConcentrationK2/2);
    pumpBaseCurrent = priv->getPumpBaseCurrent(priv->outerPotassiumConcentration + outerPotassiumConcentrationK2/2,
					       priv->innerSodiumConcentration + innerSodiumConcentrationK2/2);
    potassiumCurrent = priv->getTotalPotassiumCurrent(priv->potential + potentialK2/2, priv->n + nK2/2,
						      potassiumReversalPotential, pumpBaseCurrent);
    sodiumCurrent = priv->getTotalSodiumCurrent(priv->potential + potentialK2/2, 
						priv->m + mK2/2, priv->h + hK2/2, 
						priv->blebbedM + blebbedMK2/2, priv->blebbedH + blebbedHK2/2,
						sodiumReversalPotential, pumpBaseCurrent);
    leakCurrent = priv->getLeakCurrent(priv->potential + potentialK2/2);
    totalCurrent = potassiumCurrent + sodiumCurrent + leakCurrent + priv->stimulation;
    nK3 = time * derivativeN(priv->potential + potentialK2/2, priv->n + nK2/2);
    mK3 = time * derivativeM(priv->potential + potentialK2/2, priv->m + mK2/2);
    hK3 = time * derivativeH(priv->potential + potentialK2/2, priv->h + hK2/2);
    //blebbedNK3 = time * derivativeN(priv->potential + priv->leftShift + potentialK2/2, priv->blebbedN + blebbedNK2/2);
    blebbedMK3 = time * derivativeM(priv->potential + priv->leftShift + potentialK2/2, priv->blebbedM + blebbedMK2/2);
    blebbedHK3 = time * derivativeH(priv->potential + priv->leftShift + potentialK2/2, priv->blebbedH + blebbedHK2/2);
    innerPotassiumConcentrationK3 = time
      * priv->derivativeInnerConcentration(potassiumCurrent);
    outerPotassiumConcentrationK3 = time
      * priv->derivativeOuterConcentration(potassiumCurrent);
    innerSodiumConcentrationK3 = time
      * priv->derivativeInnerConcentration(sodiumCurrent);
    outerSodiumConcentrationK3 = time
      * priv->derivativeOuterConcentration(sodiumCurrent);
    potentialK3 = time
      * (-totalCurrent / priv->capacitance);
		
    potassiumReversalPotential = priv->calculateReversalPotential(priv->innerPotassiumConcentration + innerPotassiumConcentrationK3,
								  priv->outerPotassiumConcentration + outerPotassiumConcentrationK3);
    sodiumReversalPotential = priv->calculateReversalPotential(priv->innerSodiumConcentration + innerSodiumConcentrationK3, 
							       priv->outerSodiumConcentration + outerSodiumConcentrationK3);
    pumpBaseCurrent = priv->getPumpBaseCurrent(priv->outerPotassiumConcentration + outerPotassiumConcentrationK3,
					       priv->innerSodiumConcentration + innerSodiumConcentrationK3);
    potassiumCurrent = priv->getTotalPotassiumCurrent(priv->potential + potentialK3, priv->n + nK3,
						      potassiumReversalPotential, pumpBaseCurrent);
    sodiumCurrent = priv->getTotalSodiumCurrent(priv->potential + potentialK3, 
						priv->m + mK3, priv->h + hK3, 
						priv->blebbedM + blebbedMK3, priv->blebbedH + blebbedHK3,
						sodiumReversalPotential, pumpBaseCurrent);
    leakCurrent = priv->getLeakCurrent(priv->potential + potentialK3);
    totalCurrent = potassiumCurrent + sodiumCurrent + leakCurrent + priv->stimulation;
    nK4 = time * derivativeN(priv->potential + potentialK3, priv->n + nK3);
    mK4 = time * derivativeM(priv->potential + potentialK3, priv->m + mK3);
    hK4 = time * derivativeH(priv->potential + potentialK3, priv->h + hK3);
    //blebbedNK4 = time * derivativeN(priv->potential + priv->leftShift + potentialK3, priv->blebbedN + blebbedNK3);
    blebbedMK4 = time * derivativeM(priv->potential + priv->leftShift + potentialK3, priv->blebbedM + blebbedMK3);
    blebbedHK4 = time * derivativeH(priv->potential + priv->leftShift + potentialK3, priv->blebbedH + blebbedHK3);
    innerPotassiumConcentrationK4 = time
      * priv->derivativeInnerConcentration(potassiumCurrent);
    outerPotassiumConcentrationK4 = time
      * priv->derivativeOuterConcentration(potassiumCurrent);
    innerSodiumConcentrationK4 = time
      * priv->derivativeInnerConcentration(sodiumCurrent);
    outerSodiumConcentrationK4 = time
      * priv->derivativeOuterConcentration(sodiumCurrent);
    potentialK4 = time
      * (-totalCurrent / priv->capacitance);

    Scalar change = 0, nTemp, mTemp, hTemp, blebbedMTemp, blebbedHTemp, innerPotassiumConcentrationTemp, outerPotassiumConcentrationTemp,
      innerSodiumConcentrationTemp, outerSodiumConcentrationTemp, potentialTemp;

    nTemp = nK1 + nK2 + nK3 + nK4;
    mTemp = mK1 + mK2 + mK3 + mK4;
    hTemp = hK1 + hK2 + hK3 + hK4;
    //Temp = blebbedNK1 + blebbedNK2 + blebbedNK3 + blebbedNK4blebbedN;
    blebbedMTemp = blebbedMK1 + blebbedMK2 + blebbedMK3 + blebbedMK4;
    blebbedHTemp = blebbedHK1 + blebbedHK2 + blebbedHK3 + blebbedHK4;
    potentialTemp = potentialK1 + potentialK2 + potentialK3 + potentialK4;

    change = max<Scalar>(change, fabs(nTemp/priv->n));
    change = max<Scalar>(change, fabs(mTemp/priv->m));
    change = max<Scalar>(change, fabs(hTemp/priv->h));
    change = max<Scalar>(change, fabs(blebbedMTemp/priv->blebbedM));
    change = max<Scalar>(change, fabs(blebbedHTemp/priv->blebbedH));
    change = max<Scalar>(change, fabs(potentialTemp/priv->potential));
    change /= 6;

    if(change < limit || force) {
      innerPotassiumConcentrationTemp = innerPotassiumConcentrationK1 + innerPotassiumConcentrationK2 + innerPotassiumConcentrationK3 + innerPotassiumConcentrationK4;
      outerPotassiumConcentrationTemp = outerPotassiumConcentrationK1 + outerPotassiumConcentrationK2 + outerPotassiumConcentrationK3 + outerPotassiumConcentrationK4;
      innerSodiumConcentrationTemp = innerSodiumConcentrationK1 + innerSodiumConcentrationK2 + innerSodiumConcentrationK3 + innerSodiumConcentrationK4;
      outerSodiumConcentrationTemp = outerSodiumConcentrationK1 + outerSodiumConcentrationK2 + outerSodiumConcentrationK3 + outerSodiumConcentrationK4;

      priv->n += nTemp/6;
      priv->m += mTemp/6;
      priv->h += hTemp/6;
      //priv->blebbedN += (blebbedNK1 + blebbedNK2 + blebbedNK3 + blebbedNK4)/6;
      priv->blebbedM += blebbedMTemp/6;
      priv->blebbedH += blebbedHTemp/6;
      priv->innerPotassiumConcentration += innerPotassiumConcentrationTemp/6;
      priv->outerPotassiumConcentration += outerPotassiumConcentrationTemp/6;
      priv->innerSodiumConcentration += innerSodiumConcentrationTemp/6;
      priv->outerSodiumConcentration += outerSodiumConcentrationTemp/6;
      priv->lastPotential = priv->potential;
      priv->potential += potentialTemp/6;
      priv->potassiumReversalPotential = priv->calculateReversalPotential(priv->innerPotassiumConcentration, priv->outerPotassiumConcentration);
      priv->sodiumReversalPotential = priv->calculateReversalPotential(priv->innerSodiumConcentration, priv->outerSodiumConcentration);
    }
    return change;
  }

  void HodgkinHuxley::setStimulation(const Scalar stimulation) {
    priv->stimulation = stimulation;
  }

  void HodgkinHuxley::setBlebbing(const Scalar blebbing) {
    priv->blebbing = blebbing;
  }

  void HodgkinHuxley::setLeftShift(const Scalar leftShift) {
    priv->leftShift = leftShift;
  }

  bool HodgkinHuxley::isSpiked() const {
    return priv->potential > priv->threshold && priv->lastPotential <= priv->threshold;
  }

  Scalar HodgkinHuxley::getPotential() const {
    return priv->potential;
  }

  Scalar HodgkinHuxley::getPotassiumReversalPotential() const {
    return priv->potassiumReversalPotential;
  }

  Scalar HodgkinHuxley::getSodiumReversalPotential() const {
    return priv->sodiumReversalPotential;
  }

  string HodgkinHuxley::toString() const {
    stringstream stream;
    stream << priv->potential << "\t" << 
      priv->n << "\t" << priv->m << "\t" << priv->h << "\t" << priv->stimulation << "\t" << priv->blebbing << "\t" << priv->leftShift << endl;
    return stream.str();
  }

  ostream &operator<<(ostream &stream, HodgkinHuxley neuron) {
    stream << neuron.priv->potential << "\t" << neuron.priv->n << "\t" << neuron.priv->m << "\t" << neuron.priv->h << "\t" << endl;
  }
}


