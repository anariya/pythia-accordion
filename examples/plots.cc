// PYTHIA script to produce a number of plots, including:
// - dN/dy
// - dpT/dy
// - pdfs of delta y, for regular vs. joining step hadrons
// - hadron species ratios, for regular vs. joining step hadrons

// Author: Jade Abidi
// Created: 25/09/2025

#include "Pythia8/Pythia.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>

using namespace std;
using namespace Pythia8;

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

    Hist deltayReg("delta y pdf for regular hadrons", 100, -5., 5., false, true);
    Hist deltayJoin("delta y pdf for joining hadrons", 100, -5., 5., false, true);
    Hist massReg("mass of regular hadrons", 500, 0., 2., false, true);
    Hist massJoin("mass of joining hadrons", 500, 0., 2., false, true);

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
      // Also add all rapidities to dN/dy histogram.
      vector<int> primary;
      for (int i = 0; i < event.size(); ++i) {
	int status = event[i].status();
	if (status == 1216 || (status > 80 && status < 90)) {
	  primary.push_back(i);
	  dNdy.fill(event[i].y());
	  if (status == 1216)
	    massJoin.fill(event[i].m());
	  else
	    massReg.fill(event[i].m());
	}
      }

      // Step through string, adding rapidity spacings to appropriate histograms.
      for (size_t i = 0; i < primary.size() - 1; ++i) {
	int status = event[primary[i]].status();
	double deltay = event[primary[i]].y() - event[primary[i + 1]].y();
	if (status == 1216 || event[primary[i + 1]].status() == 1216) {
	  // Joining step.
	  deltayJoin.fill(deltay);
	} else if (i == 0) {
	  // First hadron.
	} else if (i == primary.size() - 2) {
	  // Last hadron.
	} else {
	  deltayReg.fill(deltay);
	}
      }
    }

    // Normalise histograms.
    dNdy.normalizeSpectrum(nEvents);
    deltayReg.normalizeIntegral();
    deltayJoin.normalizeIntegral();
    massReg.normalizeIntegral();
    massJoin.normalizeIntegral();

    // Print histograms.
    pythia.stat();
    cout << dNdy << deltayReg << deltayJoin << massReg << massJoin;
  }

  // Finalise.
  return 0;
}
