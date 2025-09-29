// PYTHIA script to produce a number of plots, including:
// - dN/dy
// - dpT/dy
// - pdfs of delta y, for regular vs. joining step hadrons
// - hadron species ratios, for regular vs. joining step hadrons

// Author: Jade Abidi
// Created: 25/09/2025

#include "Pythia8/Pythia.h"
#include "Pythia8/StringFragmentation.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>

using namespace std;
using namespace Pythia8;

namespace Pythia8 {
  extern map<int, int> getFlavCountsReg();
  extern map<int, int> getFlavCountsJoin();
}

int main() {
  // Specify string CM energy (GeV).
  double cme = 500;

  // Specify id of quark for q-qbar hadronisation.
  int qId = 1;

  // Set up generator.
  Pythia pythia;
  Event& event = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Read in settings from file.
  pythia.readFile("plots.cmnd");

  // Retrieve number of events and subruns.
  int nSubruns = pythia.mode("Main:numberOfSubruns");
  int nEvents = pythia.mode("Main:numberOfEvents");

  // Iterate over subruns.
  for (int iRun = 1; iRun <= nSubruns; ++iRun) {
    // Initialise.
    pythia.readFile("plots.cmnd", iRun);
    string runName = pythia.word("Main:spareWord1");
    cout << "Initialising PYTHIA for q-qbar hadronisation, string mass = "
	 << cme << endl;
    if (!pythia.init()) return 1;

    // Book histograms.
    Hist dNdy("dN/dy distribution of all hadrons", 100, -10., 10., false, true);
    Hist deltayReg("delta y pdf for regular hadrons", 100, -5., 5., false,
		   true);
    Hist deltayJoinNeg("delta y pdf for left joining neighbour",
			     100, -5., 5., false, true);
    Hist deltayJoinPos("delta y pdf for right joining neighbour", 100, -5., 5.,
		       false, true);
    Hist deltayJoinBetween("delta y pdf between joining hadron", 100, -5., 5.,
			   false, true);

    // Initialise maps that will store counts of each primary hadron ID, for
    // joining step and regular hadrons.
    map<int, int> hadronCountsReg;
    map<int, int> hadronCountsJoin;
    int totalReg = 0;
    int totalJoin = 0;

    // Event loop.
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
      // Reset event record, add q-qbar pair.
      event.reset();
      double mm = pdt.m0(qId);
      double ee = cme / 2;
      double pp = sqrtpos(pow2(ee) - pow2(mm));
      event.append(qId, 23, 101, 0, 0., 0., pp, ee, mm);
      event.append(-qId, 23, 0, 101, 0., 0., -pp, ee, mm);

      // Generate event.
      bool eventSuccess = pythia.next();
      if (!eventSuccess) {
	cout << "Error: Event generation failed." << endl;
	break;
      }
       
      // Get indices of primary hadrons in event record.
      // Also add all rapidities to dN/dy histograms.
      vector<int> primary;
      for (int i = 0; i < event.size(); ++i) {
	int status = event[i].status();
	if (status == 1216 || (status > 80 && status < 90)) {
	  primary.push_back(i);
	  dNdy.fill(event[i].y());
	  if (hadronCountsReg.find(event[i].id()) == hadronCountsReg.end())
            hadronCountsReg[event[i].id()] = 0;
	  if (hadronCountsJoin.find(event[i].id()) == hadronCountsJoin.end())
            hadronCountsJoin[event[i].id()] = 0;
	}
      }

      // Step through string, adding rapidity spacings to appropriate histograms.
      // Also populate hadron counts.
      for (size_t i = 0; i < primary.size() - 1; ++i) {
	int status = event[primary[i]].status();
	double deltay = event[primary[i]].y() - event[primary[i + 1]].y();
	if (status == 1216 && event[primary[i + 1]].status() == 1216) {
	  // Between joining hadrons.
	  deltayJoinBetween.fill(deltay);
	  hadronCountsJoin[event[primary[i]].id()] =
	    hadronCountsJoin[event[primary[i]].id()] + 1;
	  totalJoin += 1;
	} else if (status == 1216) {
	  // Neighbour to join, and itself a joining hadron.
	  deltayJoinPos.fill(deltay);
	  hadronCountsJoin[event[primary[i]].id()] =
	    hadronCountsJoin[event[primary[i]].id()] + 1;
	  totalJoin += 1;
	} else if (event[primary[i + 1]].status() == 1216) {
	  // Neighbour to join.
          deltayJoinNeg.fill(deltay);
	} else if (i == 0) {
	  // First hadron.
	} else if (i == primary.size() - 2) {
	  // Don't count rapidity spacing to last hadron, but count as regular.
	  hadronCountsReg[event[primary[i]].id()] =
	    hadronCountsReg[event[primary[i]].id()] + 1;
	  totalReg += 1;
	} else {
	  // Regular hadron.
	  deltayReg.fill(deltay);
	  hadronCountsReg[event[primary[i]].id()] =
	    hadronCountsReg[event[primary[i]].id()] + 1;
	  totalReg += 1;
	}
      }
    }

    // Calculate ratios of hadron counts between regular and joining step.
    // Print each out.
    for (const auto& cpair : hadronCountsJoin) {
      int hadId = cpair.first;
      double regProb = (double)hadronCountsReg[hadId] / (double)totalReg;
      double joinProb = (double)cpair.second / (double)totalJoin;
      double ratio = joinProb / regProb;
      string hadName = pdt.name(hadId);
      double hadMass = pdt.m0(hadId);
      cout << "Joining / regular ratio for " << hadName << " (" << hadId
	   << ", m=" << hadMass << "):    " << ratio << endl;
    }

    // Calculate ratios of flavour counts between regular and joining step.
    // Print each out.
    map<int, int> flavCountsJoin = getFlavCountsJoin();
    map<int, int> flavCountsReg = getFlavCountsReg();
    int totalFlavJoin = 0;
    int totalFlavReg = 0;
    for (const auto& cpair : flavCountsJoin)
      totalFlavJoin += cpair.second;
    for (const auto& cpair : flavCountsReg)
      totalFlavReg += cpair.second;
    
    for (const auto& cpair : flavCountsJoin) {
      int flavCombId = cpair.first;
      double regProb = (double)flavCountsReg[flavCombId] / (double)totalFlavReg;
      double joinProb = (double)cpair.second / (double)totalFlavJoin;
      double ratio = joinProb / regProb;
      cout << "Joining / regular ratio for " << flavCombId
	   << ":    " << ratio << endl;
    }

    // Normalise histograms.
    dNdy.normalizeSpectrum(nEvents);
    deltayReg.normalizeIntegral();
    deltayJoinPos.normalizeIntegral();
    deltayJoinNeg.normalizeIntegral();
    deltayJoinBetween.normalizeIntegral();

    // Print histograms.
    pythia.stat();
    cout << dNdy << deltayReg << deltayJoinPos << deltayJoinNeg << deltayJoinBetween;
  }

  // Finalise.
  return 0;
}
