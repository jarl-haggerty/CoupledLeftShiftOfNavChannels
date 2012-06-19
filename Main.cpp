#include "HodgkinHuxley.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <time.h>

using namespace std;
using namespace Jarl;

HodgkinHuxley::Settings settings;
Scalar resolution;
int writeResolution, what;
/*
void stimulationControlExperiment() {
  HodgkinHuxley neuron(settings);

  ofstream output("stimulation_control_experiment.tsv");
  output << "time (ms)\tpotential (mV)\tpotassium reversal potential (mv)\tsodium reversal potential (mv)" << endl;

  Scalar stimulationStart = 250, stimulationDuration = 333, stimulation = -12, currentTime;
  int limit = 5000/resolution;
  time_t start = time(NULL);
  for (int i = 0; i < limit; i++) {
    currentTime = i*resolution;
    if (i % writeResolution == 0) {
      output << currentTime << "\t" << neuron.getPotential() << "\t" << neuron.getPotassiumReversalPotential() << "\t" << neuron.getSodiumReversalPotential() << endl;
    }
    if (currentTime > stimulationStart
	&& currentTime < stimulationStart + stimulationDuration) {
      neuron.setStimulation(stimulation);
    } else {
      neuron.setStimulation(0);
    }
    neuron.simulate(resolution);
  }
  cout << difftime(time(NULL), start) << endl;

  output.close();
}

void spontaneousBlebbingExperiment() {
  HodgkinHuxley neuron(settings);

  ofstream output("spontaneous_blebbing_experiment.tsv");
  output << "time (ms)\tpotential (mV)\tpotassium reversal potential (mv)\tsodium reversal potential (mv)" << endl;

  Scalar blebTime = 250, blebbing = 1, leftShift = 19, stimulationDuration = 333, limit = 20000/resolution, currentTime, stimulation = -12;
  time_t start = time(NULL);
  for (int i = 0; i < limit; i++) {
    currentTime = i * resolution;
    if (i % writeResolution == 0) {
      output << currentTime << "\t" << neuron.getPotential() << "\t" << neuron.getPotassiumReversalPotential() << "\t" << neuron.getSodiumReversalPotential() << endl;
    }
    if (currentTime > blebTime) {
      neuron.setBlebbing(blebbing);
      neuron.setLeftShift(leftShift);
    } else {
      neuron.setBlebbing(0);
      neuron.setLeftShift(0);
    }
    neuron.simulate(resolution);
  }
  cout << difftime(time(NULL), start) << endl;
  output.close();
}
*/
void rateExperiment(Scalar blebbing, Scalar leftShift, Scalar stimulation) {
  HodgkinHuxley neuron(settings);
  //cout << "    " << neuron.getPotential() << endl;

  stringstream nameStream;
  nameStream << "experiments/blebbing_" << blebbing << "_left_shift_" << leftShift << "_stimulation_" << stimulation << ".tsv";
  ofstream output(nameStream.str().c_str());

  output << "time (ms)\tpotential (mV)\tpotassium reversal potential (mv)\tsodium reversal potential (mv)" << endl;
  Scalar blebbingStart = 100, stimulationStart = 500, duration = 200500, currentTime, change;

  time_t start = time(NULL);

  int limit = duration/resolution, steps = 1, writes = 0, second = 0;
  what = 0;
  for (int i = 0; i < limit;) {
    currentTime = i * resolution;
    if (i > writeResolution*writes) {
      writes++;
      //cout << steps << endl;
      output << currentTime << "\t" << neuron.getPotential() << "\t" << neuron.getPotassiumReversalPotential() << "\t" << neuron.getSodiumReversalPotential() << endl;
    }
    //cout << currentTime << "\t" << steps << "\t" << change << endl;
    if(difftime(time(NULL), start) >= second) {
      second++;
      cout << currentTime << "\t" << steps << "\t" << change << endl;
    }
    //if(neuron.isSpiked() && currentTime < 300) {
    //  cout << currentTime << endl;
    //}
    //if(what == 0) {
    //  cout << neuron.toString() << endl;
    //  what = 1;
    //}
    if(currentTime > blebbingStart) {
      neuron.setBlebbing(blebbing);
      neuron.setLeftShift(leftShift);
      if (currentTime > stimulationStart) {
	neuron.setStimulation(stimulation);
      } else {
	neuron.setStimulation(0);
      }
    } else {
      neuron.setBlebbing(0);
      neuron.setLeftShift(0);
    }
    change = neuron.simulate(steps*resolution, 1e-3, steps == 1);
    while(true) {
      if(change > 1e-3) {
	if(steps > 1) {
	  steps /= 10;
	  change = neuron.simulate(steps*resolution, 1e-3, steps == 1);
	} else {
	  i += steps;
	  break;
	}
      } else if(change < 1e-4) {
	i += steps;
	if(steps < 100) {
	  steps *= 10;
	}
	break;
      } else {
	i += steps;
	break;
      }
    }
  }
  what  = 1;

  cout << blebbing << " " << leftShift << " " << stimulation << " " << difftime(time(NULL), start) << endl;
  output.close();
}

int main(int argc, char* argv[]) {
  resolution = 1e-3;
  writeResolution = (int)(.1/resolution);
  what = 0;

  settings = HodgkinHuxley::Settings();
  settings.potential = -59.9;
  settings.threshold = -15;
  settings.capacitance = 1;
  settings.leakConductance = .5;
  settings.leakReversalPotential = -59.9;
  settings.potassiumConductance = 36;
  settings.potassiumReversalPotential = -81.3;
  settings.sodiumConductance = 120;
  settings.sodiumReversalPotential = 51.5;
  settings.maxPumpCurrent = 90.9;
  settings.potassiumLeakConductance = .1;
  settings.sodiumLeakConductance = .25;
  settings.innerPotassiumConcentration = 150;
  settings.outerPotassiumConcentration = 6;
  settings.innerSodiumConcentration = 20;
  settings.outerSodiumConcentration = 154;
  settings.temperature = 293.15;
  settings.innerVolume = 3e-15;
  settings.outerVolume = 3e-15;
  settings.surfaceArea = 6e-8;

  /*
  for(Scalar blebbing = 0;blebbing < 1;blebbing += .01) {
    for(Scalar leftShift = 0;leftShift < 40;leftShift += .4) {
      rateExperiment(blebbing, leftShift, 0);
    }
  }*/

  rateExperiment(1, 2, 0);
  
  /*
  for(Scalar blebbing = 0;blebbing < 1;blebbing += .01) {
    for(Scalar leftShift = 0;leftShift < 40;leftShift += .4) {
      rateExperiment(blebbing, leftShift, -12);
    }
  }
  */
}
