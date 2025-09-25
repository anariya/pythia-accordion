// PYTHIA script to generate histograms of rapidities of primary hadrons
// produced by fragmentation. Simulation is done for a single q-qbar string
// and only the hadronisation process is considered, with parton shower and
// other effects disabled. The invariant mass of the string can be varied.
// Other histograms such as dE/dy and dpT/dy are also generated.

// Author: Jade Abidi
// Created: 30/07/2025

#include "Pythia8/Pythia.h"
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

using namespace Pythia8;
using namespace std;

int main() {
  // Specify invariant string CM energy (GeV).
  double stringMass = 500;

  // Specify number of events to simulate.
  int nEvent = 1000000;

  // Specify id of quark.
  // 1 - down. 2 - up. 3 - strange. 4 - charm. 5 - bottom. 6 - tnop.
  int qid = 1;

  // Option for massless quarks.
  bool masslessQuarks = false;

  // Set up generator.
  Pythia pythia;
  Event& event = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Disable parton shower and hard process since q-qbar will be manually
  // input.
  pythia.readString("ProcessLevel:all = off");

  // Optional: set tune.
  pythia.readString("Tune:ee = 1");

  // Use the bug-fixed version of the aExtraDiquark parameter. 
  pythia.readString("StringZ:useOldAExtra = off");
    
  // Disable hadron decay.
  pythia.readString("HadronLevel:Decay = off");

  // Optional: Disable transverse momentum (enforce 1+1 dimensions)
  //pythia.readString("StringPT:sigma = 0");

  // Customise output to be more readable and less cluttered.
  pythia.readString("Next:numberCount = 100000");

  // Customise jet joining top mass
  pythia.readString("StringFragmentation:stopMass = 0.8");
  // pythia.readString("StringFragmentation:stopNewFlav = 1.5");

  // Initialise.
  cout << "Initialising PYTHIA for q-qbar hadronisation, string mass = " \
       << stringMass << endl;
  if (!pythia.init()) return 1;
    
  // Set up histograms.
  Hist dndy("Rapidity distribution dn/dy of primary hadrons", 100,
		   -10., 10., false, true);
  Hist dpTdy("Distribution of total transverse momentum over rapidity dpT/dy",
	     100, -10., 10., false, true);
  Hist histZ("z+ distribution of primary hadrons", 100, 0., 1., false, true);
  Hist histZ1("z+ distribution of 1st-rank primary hadron", 100, 0., 1., false, true);
  Hist histZ2("z+ distribution of 2nd-rank primary hadron", 100, 0., 1., false, true);
  Hist histZ3("z+ distribution of 3rd-rank primary hadron", 100, 0., 1., false, true);
  Hist histZ4("z+ distribution of 4th-rank primary hadron", 100, 0., 1., false, true);
  Hist histZ5("z+ distribution of 5th-rank primary hadron", 100, 0., 1., false, true);
  Hist histZ6("z+ distribution of 6th-rank primary hadron", 100, 0., 1., false, true);
  Hist histZLast("z+ distribution of last-rank primary hadron", 100, 0., 1., false, true);
  Hist histZMid("z+ distribution of mid-rank primary hadrons", 100, 0., 1., false, true);

  // Event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    // Reset event record and add q-qbar pair.
    event.reset();
    double mm = masslessQuarks ? 0 : pdt.m0(qid);
    double ee = stringMass / 2;
    double pp = sqrtpos(ee*ee - mm*mm);
    event.append(qid, 23, 101, 0, 0., 0., pp, ee, mm);
    event.append(-qid, 23, 0, 101, 0., 0., -pp, ee, mm);

    // Generate event.
    if (!pythia.next()) {
      cout << "Error: Event generation failed." << endl;
      break;
    }

    // Loop over particles.
    double pPosTot = stringMass;
    int iRank = 0;
    for (int i = 0; i < event.size(); ++i) {
      // Add primary hadron rapidities to histograms.
      int status = event[i].statusAbs();
      double y = event[i].y();
      if (status > 80 && status < 90) {
        if ( event[i].status() == 83 ) {
          iRank += 1;
          double pPos = event[i].e() + event[i].pz();
          double zPos = pPos / pPosTot; 
          histZ.fill(zPos);
          if ( event[i+1].statusAbs() == 1216 ) histZLast.fill(zPos);         
          else {
            if ( iRank != 1 ) histZMid.fill(zPos);
            if ( iRank == 1 ) histZ1.fill( zPos ); 
            else if ( iRank == 2 ) histZ2.fill( zPos ); 
            else if ( iRank == 3 ) histZ3.fill( zPos ); 
            else if ( iRank == 4 ) histZ4.fill( zPos ); 
            else if ( iRank == 5 ) histZ5.fill( zPos ); 
            else if ( iRank == 6 ) histZ6.fill( zPos ); 
          }
          // Update total remaining pPlus. 
          pPosTot -= pPos; 
        }
        // Fill dn/dy histogram. 
	dndy.fill(y);
	// Fill dpT/dy histogram.
	dpTdy.fill(y, event[i].pT());
      }
    }
  }

  // Rescale histograms to show dn/dy, dpT/dy.
  dndy.normalizeSpectrum(nEvent);
  dpTdy.normalizeSpectrum(nEvent);
  histZ.normalizeSpectrum(nEvent);
  histZ1.normalizeSpectrum(nEvent);
  histZ2.normalizeSpectrum(nEvent);
  histZ3.normalizeSpectrum(nEvent);
  histZ4.normalizeSpectrum(nEvent);
  histZ5.normalizeSpectrum(nEvent);
  histZ6.normalizeSpectrum(nEvent);
  histZLast.normalizeSpectrum(nEvent);
  histZMid.normalizeSpectrum(nEvent);
    
  // Print statistics and histograms.
  pythia.stat();
  cout << dndy << dpTdy;
  cout << histZ << histZ1 << histZ2 << histZ3 << histZ4 << histZ5 << histZ6
       << histZLast << histZMid;

  // Matplotlib output.
  HistPlot hpl("rapidityplot");
  hpl.frame("dndy_latest", "dn/dy for quark-antiquark pair at 500 GeV",
	    "y", "dn/dy");
  hpl.add(dndy, "-");
  hpl.plot();
  hpl.frame("dpTdy_latest", "dpT/dy for quark-antiquark pair at 500 GeV",
	    "y", "dpT/dy");
  hpl.add(dpTdy, "-");
  hpl.plot();

  // Output histograms to pyplot table.
  dndy.pyplotTable("dndy_latest.csv", false);
  dpTdy.pyplotTable("dpTdy_latest.csv", false);
  
  return 0;
}
