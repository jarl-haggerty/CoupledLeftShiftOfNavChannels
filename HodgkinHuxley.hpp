#ifndef HODGKIN_HUXLEY_HPP
#define HODGKIN_HUXLEY_HPP

#include <iostream>
#include <string>

namespace Jarl {
  typedef double Scalar;
  
  extern const Scalar potassiumDissociationConstant;
  extern const Scalar sodiumDissociationConstant;
  extern const Scalar gasConstant;
  extern const Scalar faradayConstant;

  class HodgkinHuxley {
  public:
    class Settings {
    public:
      Settings();
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
      Scalar stimulation;
    };
    HodgkinHuxley(const Settings& settings);
    ~HodgkinHuxley();
    Scalar simulate(const Scalar time, const Scalar limit, const bool force);
    void setStimulation(const Scalar stimulation);
    void setBlebbing(const Scalar blebbing);
    void setLeftShift(const Scalar leftShift);
    bool isSpiked() const;
    Scalar getPotential() const;
    Scalar getPotassiumReversalPotential() const;
    Scalar getSodiumReversalPotential() const;
    std::string toString() const;
    friend std::ostream &operator<<(std::ostream &stream, HodgkinHuxley neuron);
  private:
    class Private;
    Private* const priv;
  };
}

#endif
