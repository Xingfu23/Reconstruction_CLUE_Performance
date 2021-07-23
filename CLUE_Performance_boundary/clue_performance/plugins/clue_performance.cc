// -*- C++ -*-
//
// Package:    CLUE_Performance_boundary/clue_performance
// Class:      clue_performance
//
/**\class clue_performance clue_performance.cc CLUE_Performance_boundary/clue_performance/plugins/clue_performance.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Xing-Fu Su
//         Created:  Fri, 23 Jul 2021 03:17:24 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"

#include <TTree.h>
//
// class declaration
//

const int kmax = 20000;
using reco::TrackCollection;
using reco::PFCandidateCollection;
using reco::PhotonCollection;

class clue_performance : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit clue_performance(const edm::ParameterSet&);
  ~clue_performance();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  //-----------physics variables----------------------
  int evt_no = 0;
  hgcal::RecHitTools tool;
  TTree* tree = new TTree("ntuple", "ntuple");
  int layercluster_number;
  int thickness_type[kmax];
  int layercluster_layer[kmax];
  double layercluster_energy[kmax];
  double layercluster_phi[kmax];
  double layercluster_eta[kmax];
  double layercluster_z[kmax];
  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  
  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<reco::CaloClusterCollection> layerclustersToken_;
  edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster>> multisToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfcandsticlToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfcandsToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersToken_;
  // ----------member data ---------------------------

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
clue_performance::clue_performance(const edm::ParameterSet& iConfig)
    //: tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {
    : layerclustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("layerclusters"))),
      photonsToken_(consumes<reco::PhotonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("photons"))),
      multisToken_(consumes<std::vector<reco::HGCalMultiCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("multis"))),
      pfcandsticlToken_(consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfcandsticl"))),
      pfcandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfcands"))),
      trackstersToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksters")))
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
          setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

clue_performance::~clue_performance() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void clue_performance::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  evt_no += 1;
  std::cout << "Event Number: " << evt_no << std::endl;
  Handle<reco::CaloClusterCollection> layerclusters;
  iEvent.getByToken(layerclustersToken_, layerclusters);

  layercluster_number = 0;
  for (reco::CaloClusterCollection::const_iterator cluster_idx = layerclusters->begin(); cluster_idx != layerclusters->end(); cluster_idx++) {
    if (TMath::Abs(cluster_idx->z())>367.699 && tool.getLayer(cluster_idx->seed())<28){
      layercluster_layer[layercluster_number] = tool.getLayer(cluster_idx->seed())+28; 
    }
    else{
      layercluster_layer[layercluster_number] = tool.getLayer(cluster_idx->seed());
    }
    thickness_type[layercluster_number] = tool.getSiThickIndex(cluster_idx->seed());
    layercluster_energy[layercluster_number] = cluster_idx->energy();
    layercluster_phi[layercluster_number] = cluster_idx->phi();
    layercluster_eta[layercluster_number] = cluster_idx->eta();
    layercluster_z[layercluster_number] = cluster_idx->z();
    layercluster_number += 1;
  }

  tree->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void clue_performance::beginJob() {
  // please remove this method if not needed

  tree->Branch("layercluster_number", &layercluster_number, "layercluster_number/I");
  tree->Branch("layercluster_layer", &layercluster_layer, "layercluster_layer[layercluster_number]/I");
  tree->Branch("thickness_type",  &thickness_type, "thickness_type[layercluster_number]/I");
  tree->Branch("layercluster_energy", &layercluster_energy, "layercluster_energy[layercluster_number]/D");
  tree->Branch("layercluster_z", &layercluster_z, "layercluster_z[layercluster_number]/D");
  tree->Branch("layercluster_phi", &layercluster_phi, "layercluster_phi[layercluster_number]/D");
  tree->Branch("layercluster_eta", &layercluster_eta, "layercluster_eta[layercluster_number]/D");
}

// ------------ method called once each job just after ending the event loop  ------------
void clue_performance::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void clue_performance::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(clue_performance);
