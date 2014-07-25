// -*- C++ -*-
//
// Package:    Analysis/MiniAODjetTester
// Class:      MiniAODjetTester
// 
/**\class MiniAODjetTester MiniAODjetTester.cc Analysis/MiniAODjetTester/plugins/MiniAODjetTester.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Dolen
//         Created:  Mon, 07 Jul 2014 15:53:35 GMT
//
//


// system include files
#include <memory>

// core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormats
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h" 
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"


#include "DataFormats/Math/interface/LorentzVector.h"

// utilities
#include "DataFormats/Math/interface/deltaR.h"

// TFile
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// root
#include "TH1.h"

// fastjet
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
// #include "Nsubjettiness.hh"
// #include "Njettiness.hh"

//
// class declaration
//

class MiniAODjetTester : public edm::EDAnalyzer {
   public:
      explicit MiniAODjetTester(const edm::ParameterSet&);
      ~MiniAODjetTester();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

      // ----------root ---------------------------

      TH1F * LeadingJetMassHasTop          ;
      TH1F * LeadingJetMassNoTop           ;
};


//
// constructors and destructor
//
MiniAODjetTester::MiniAODjetTester(const edm::ParameterSet& iConfig):
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGen"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedGen")))
{
    edm::Service<TFileService> fs;
    LeadingJetMassHasTop        = fs->make<TH1F>("LeadingJetMassHasTop",            "LeadingJetMassHasTop",         100,  0,  400 );
    LeadingJetMassNoTop        = fs->make<TH1F>("LeadingJetMassNoTop",            "LeadingJetMassNoTop",         100,  0,  400 );

}


MiniAODjetTester::~MiniAODjetTester()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAODjetTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  using namespace fastjet;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  edm::Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(fatjetToken_, fatjets);


  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);

  double counttop = 0;
  for(size_t i=0; i<pruned->size();i++){
    int  id = (*pruned)[i].pdgId();
    int  status = (*pruned)[i].status();
    if (status<30 && status>=20 && fabs(id)==6) counttop++;

  }
  if (counttop<1){
    // cout<<"NO TOPS!!!!!"<<endl;
    // cout<<"NO TOPS!!!!!"<<endl;
    // cout<<"NO TOPS!!!!!"<<endl;
    // cout<<"NO TOPS!!!!!"<<endl;
    // cout<<"NO TOPS!!!!!"<<endl;
    for(size_t i=0; i<pruned->size();i++){
      // int  id = (*pruned)[i].pdgId();
      // int  status = (*pruned)[i].status();
      // if (status<30 && status>=20) cout<<"  id "<<id<<" status "<<status<<endl;

  //     if(abs((*pruned)[i].pdgId()) > 500 && abs((*pruned)[i].pdgId()) <600){
  //       const Candidate * bMeson = &(*pruned)[i];
  //       std::cout << "PdgID: " << bMeson->pdgId() << " pt " << bMeson->pt() << " eta: " << bMeson->eta() << " phi: " << bMeson->phi() << std::endl;
  //       std::cout << "  found daugthers: " << std::endl;
  //       for(size_t j=0; j<packed->size();j++){
  // //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
  //         const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
  //         if(motherInPrunedCollection != nullptr && isAncestor( bMeson , motherInPrunedCollection)){
  //                 std::cout << "     PdgID: " << (*packed)[j].pdgId() << " pt " << (*packed)[j].pt() << " eta: " << (*packed)[j].eta() << " phi: " << (*packed)[j].phi() << std::endl;
  //         }
  //       }
  //     }
    }
  }
  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);




  //AK4 jet loop
  bool verboseAK4=true;

  for (const pat::Jet &j :  *jets) {
    if (j.pt() < 40 || fabs(j.eta()) > 2.4) continue;

      if (verboseAK4){
        cout<<"AK4 jet - pt "<<j.pt()<<" mass "<<j.mass()<<"  nconst "<<j.numberOfDaughters()<<j.numberOfDaughters()<<" area "<<j.jetArea()
        <<" pileupID "<<j.userFloat("pileupJetId:fullDiscriminant")
        <<" vtxMass "<<j.userFloat("vtxMass")
        <<" vtxNtracks "<<j.userFloat("vtxNtracks")
        <<" vtx3DVal "<<j.userFloat("vtx3DVal")
        <<" vtx3DSig "<<j.userFloat("vtx3DSig")
        <<endl;
      }
  }

  //AK8 jet loop
  
  bool verboseAK8=true;

  for (const pat::Jet &j :  *fatjets) {
    if (j.pt() < 100 || fabs(j.eta()) > 2.4) continue;
     
    std::vector<fastjet::PseudoJet> jetparticles;

    for (unsigned int id = 0, nd = j.numberOfDaughters(); id < nd; ++id) {
        const pat::PackedCandidate &dau = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(id));  
        jetparticles.push_back( fastjet::PseudoJet( dau.px(), dau.py(), dau.pz(),dau.energy() ));
    }
    fastjet::PseudoJet combJet = fastjet::join(jetparticles);

    if (counttop>0) LeadingJetMassHasTop        ->Fill( j.userFloat("ak8PFJetsCHSPrunedLinks"));
    if (counttop<1) LeadingJetMassNoTop         ->Fill( j.userFloat("ak8PFJetsCHSPrunedLinks"));


    if (verboseAK8){ 
      cout<<"FAT jet - pt "<< j.pt()<<"  mass "<<j.mass()<<"  nconst "<<j.numberOfDaughters()<<" area "<<j.jetArea()
      // <<" pileupID "<<j.userFloat("pileupJetId:fullDiscriminant")
      // <<" vtxMass "<<j.userFloat("vtxMass")
      // <<" vtxNtracks "<<j.userFloat("vtxNtracks")
      // <<" vtx3DVal "<<j.userFloat("vtx3DVal")
      // <<" vtx3DSig "<<j.userFloat("vtx3DSig")
      <<" pruned mass "<<j.userFloat("ak8PFJetsCHSPrunedLinks")
      <<" trimmed mass "<<j.userFloat("ak8PFJetsCHSTrimmedLinks")
      <<" filtered mass "<<j.userFloat("ak8PFJetsCHSFilteredLinks")
      <<" cms mass "<<j.userFloat("cmsTopTagPFJetsCHSLinksAK8")
      <<" join mass "<<combJet.m()
      <<" uncorrected mass "<<j.correctedP4(0).mass()
      <<endl;
    }
  }
  
  // PF candidate loop
  std::vector<fastjet::PseudoJet> allparticles;
  std::vector<fastjet::PseudoJet> mychsparticles;
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate &pf = (*pfs)[i];
    allparticles.push_back( fastjet::PseudoJet( pf.px(), pf.py(), pf.pz(),pf.energy() ));
    if (pf.charge()==0 || pf.fromPV() >= 1) mychsparticles.push_back( fastjet::PseudoJet( pf.px(), pf.py(), pf.pz(), pf.energy() ));
  }

  // Cluster jets by hand
  JetDefinition jet_def_ak8 (antikt_algorithm, 0.8);
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  
  ClusterSequenceArea cs_pf_ak8      (allparticles    , jet_def_ak8  , area_def);
  ClusterSequenceArea cs_chs_ak8     (mychsparticles  , jet_def_ak8  , area_def);

  Selector selectorquick = SelectorNHardest(3); 
  vector<PseudoJet> pfJets_ak8       = selectorquick(sorted_by_pt(cs_pf_ak8    .inclusive_jets()));
  vector<PseudoJet> chsJets_ak8      = selectorquick(sorted_by_pt(cs_chs_ak8   .inclusive_jets()));

  bool verboseByHand = false;

  if (verboseByHand)
  {
    for (unsigned int i=0; i<pfJets_ak8.size(); i++ )
    {
      if (pfJets_ak8[i].pt() < 100 || fabs(pfJets_ak8[i].eta()) > 2.4) continue;
      cout<<"Clustered by hand pf jet pt - "<<pfJets_ak8[i].pt()<<" mass "<<pfJets_ak8[i].m()<<endl;
    }
    for (unsigned int i=0; i<chsJets_ak8.size(); i++ )
    {
      if (chsJets_ak8[i].pt() < 100 || fabs(chsJets_ak8[i].eta()) > 2.4) continue;
      cout<<"Clustered by hand chs jet pt - "<<chsJets_ak8[i].pt()<<" mass "<<chsJets_ak8[i].m()<<endl;
    }
  }


//++recoPFJets "ak15PFJets" "" "Demo" (productId = 5:600)
  edm::Handle<reco::PFJetCollection> ak15recojets;
  iEvent.getByLabel("ak15PFJets", ak15recojets);
  for ( size_t iJet = 0; iJet < ak15recojets->size(); ++iJet ) 
  {
    reco::PFJetRef jet(ak15recojets, iJet);

    // if ( !(jet->pt() > 30) ) continue;
        cout<<" ak15 pt "<<jet->pt()<<" mass "<<jet->mass() <<" area "<< jet->jetArea() <<" nConst "<< jet->nConstituents ()  <<endl;

  }
  // Collections not in miniaod
    // patJetsAK15PF
  // edm::Handle<reco::PFJetCollection> AK15JETS;
  // iEvent.getByLabel("patJetsAK15PF", AK15JETS);
 


  edm::Handle<std::vector<pat::Jet> > AK15JETS;
    iEvent.getByLabel( "patJetsAK15PF", AK15JETS );

    int jet_number=0;
  for ( std::vector<pat::Jet>::const_iterator jetBegin = AK15JETS->begin(), jetEnd = AK15JETS->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
    {

      double pt             = ijet->pt();
      //double phi            = ijet->phi();
      //double eta            = ijet->eta();
      //double rapidity       = ijet->rapidity();
      double mass           = ijet->mass();
      int    nconst       = ijet->numberOfDaughters();
      double csvbdisc        = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
      double    area       = ijet->jetArea();

      double ntracks = ijet->associatedTracks().size();

      reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
      //double flavour = ijet->partonFlavour();

      //double costheta = 10;
     
      //std::abs(jet->partonFlavour())

     

      // bool jet_btagged_csvm = csvmva > 0.679 ;
      // int btag =0;
      // if (jet_btagged_csvm) btag=1;
     
      cout<<"UnFilt " <<jet_number<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" csvbdisc "<<csvbdisc<< " area "<<area<<" JEC L1 "<<ijet->jecFactor("L1FastJet")<<" L2 "<< ijet->jecFactor("L2Relative")<<" L3 "<< ijet->jecFactor("L3Absolute") << " current "<<ijet->currentJECLevel()<<" uncor pt "<<uncorrJet.pt()<<endl;
      jet_number++;

    }

    edm::Handle<std::vector<pat::Jet> > AK15JETSFILTERED;
    iEvent.getByLabel( "patJetsAK15PFfilteredPacked", AK15JETSFILTERED );

    edm::Handle<std::vector<pat::Jet> > ak15filtsubjets;
    iEvent.getByLabel( "patJetsAK15PFfilteredSubjets", ak15filtsubjets );


    int jet_number2=0;
  for ( std::vector<pat::Jet>::const_iterator jetBegin = AK15JETSFILTERED->begin(), jetEnd = AK15JETSFILTERED->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
    {
      double pt             = ijet->pt();
      //double phi            = ijet->phi();
      //double eta            = ijet->eta();
      //double rapidity       = ijet->rapidity();
      double mass           = ijet->mass();
      int    nsubjets       = ijet->numberOfDaughters();
      double csvbdisc        = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
      double    area       = ijet->jetArea();


      double ntracks = ijet->associatedTracks().size();
      //double flavour = ijet->partonFlavour();

      //double costheta = 10;
     
      //std::abs(jet->partonFlavour())

     

      // bool jet_btagged_csvm = csvmva > 0.679 ;
      // int btag =0;
      // if (jet_btagged_csvm) btag=1;
     


       cout<<"Filt "<<jet_number2<<"  pt "<<pt<<" mass "<<mass<<" nsubjets "<<nsubjets<<" ntracks "<<ntracks<<" csvbdisc "<<csvbdisc<< " area "<<area<<endl;

      // //std::vector<edm::Ptr<reco::PFCandidate> > pfCands = ijet->getPFConstituents();
      // std::vector< edm::Ptr<reco::Candidate> > nextSubjets;
      
      for (int i = 0; i < nsubjets; i++ ) {

        reco::Candidate const * subjet =  ijet->daughter(i);
      //   // cout<<"subjet " <<subjet->mass()<<endl;
      //   // reco::Jet const * recosubjet = dynamic_cast<reco::Jet const *>(subjet);  
      //   // cout<<"recosubjet " <<recosubjet->mass()<<endl;
      //   // pat::Jet const * patsubjet = dynamic_cast<pat::Jet const *>(recosubjet);
      //   // cout<<"patsubjet " <<patsubjet->mass()<<endl;

        pat::Jet const * patsubjet = dynamic_cast<pat::Jet const *>(subjet);

      //   //const pat::Jet* this_subjet = dynamic_cast<const pat::Jet*>((ijet->daughter(i)));
        
           cout<<"  subjet pt "<<patsubjet->pt()<<"  mass "<<patsubjet->mass()<<"  bdisc "<<patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags") <<" area "<< patsubjet->jetArea() <<endl;
        }

      jet_number2++;

}

}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAODjetTester::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAODjetTester::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAODjetTester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODjetTester);
