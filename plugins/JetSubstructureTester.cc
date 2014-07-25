// -*- C++ -*-
//
// Package:    Analysis/JetSubstructureTester
// Class:      JetSubstructureTester
// 
/**\class JetSubstructureTester JetSubstructureTester.cc Analysis/JetSubstructureTester/plugins/JetSubstructureTester.cc

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

// CMS Top Tagger
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"
#include "RecoJets/JetAlgorithms/interface/CATopJetHelper.h"

// utilities
#include "DataFormats/Math/interface/deltaR.h"

// TFile
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// root
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

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

class JetSubstructureTester : public edm::EDAnalyzer {
   public:
      explicit JetSubstructureTester(const edm::ParameterSet&);
      ~JetSubstructureTester();

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

      TH1F * CA8PF_PT           ; 
      TH1F * CA8PF_PHI          ; 
      TH1F * CA8PF_ETA          ; 
      TH1F * CA8PF_RAP          ; 
      TH1F * CA8PF_MASS         ; 
      TH1F * CA8PF_NCONST       ; 
      TH1F * CA8PF_BDISC        ; 
      TH1F * CA8PF_AREA         ; 
      TH1F * CA8PF_NTRACKS      ; 
      TH1F * CA8PF_FLAVOUR      ; 
      TH1F * CA8PF_MASS_UNCORR  ; 
      TH1F * CA8PF_CH_MULT      ; 
      TH1F * CA8PF_NE_MULT      ; 
      TH1F * CA8PF_CHEF         ; 
      TH1F * CA8PF_NHEF         ; 
      TH1F * CA8PF_CEEF         ; 
      TH1F * CA8PF_NEEF         ; 
      TH1F * CA8PF_CMEF         ; 

      TH1F * CA8CHS_PT           ; 
      TH1F * CA8CHS_PHI          ; 
      TH1F * CA8CHS_ETA          ; 
      TH1F * CA8CHS_RAP          ; 
      TH1F * CA8CHS_MASS         ; 
      TH1F * CA8CHS_NCONST       ; 
      TH1F * CA8CHS_BDISC        ; 
      TH1F * CA8CHS_AREA         ; 
      TH1F * CA8CHS_NTRACKS      ; 
      TH1F * CA8CHS_FLAVOUR      ; 
      TH1F * CA8CHS_MASS_UNCORR  ; 
      TH1F * CA8CHS_CH_MULT      ; 
      TH1F * CA8CHS_NE_MULT      ; 
      TH1F * CA8CHS_CHEF         ; 
      TH1F * CA8CHS_NHEF         ; 
      TH1F * CA8CHS_CEEF         ; 
      TH1F * CA8CHS_NEEF         ; 
      TH1F * CA8CHS_CMEF         ;

      TH1F * AK15_PT                      ;
      TH1F * AK15_PHI                     ;
      TH1F * AK15_ETA                     ;
      TH1F * AK15_RAP                     ;
      TH1F * AK15_MASS                    ;
      TH1F * AK15_NCONST                  ;
      TH1F * AK15_BDISC                   ;
      TH1F * AK15_AREA                    ;
      TH1F * AK15_NTRACKS                 ;
      TH1F * AK15_FLAVOUR                 ;
      TH1F * AK15_MASS_UNCORR             ;

      TH1F * AK15FILT_PT                    ; 
      TH1F * AK15FILT_PHI                   ; 
      TH1F * AK15FILT_ETA                   ; 
      TH1F * AK15FILT_RAP                   ; 
      TH1F * AK15FILT_MASS                  ; 
      TH1F * AK15FILT_NCONST                ; 
      TH1F * AK15FILT_BDISC                 ; 
      TH1F * AK15FILT_AREA                  ; 
      TH1F * AK15FILT_NTRACKS               ; 
      TH1F * AK15FILT_FLAVOUR               ; 
      TH1F * AK15FILT_MASS_UNCORR           ; 
      TH1F * AK15FILT_MAXSUBJETBDISC        ; 
      TH1F * AK15FILT_MAXSUBJETBDISCFLAVOUR ; 

      TH1F * CA8PRUNE_PT                    ; 
      TH1F * CA8PRUNE_PHI                   ; 
      TH1F * CA8PRUNE_ETA                   ; 
      TH1F * CA8PRUNE_RAP                   ; 
      TH1F * CA8PRUNE_MASS                  ; 
      TH1F * CA8PRUNE_NCONST                ; 
      TH1F * CA8PRUNE_BDISC                 ; 
      TH1F * CA8PRUNE_AREA                  ; 
      TH1F * CA8PRUNE_NTRACKS               ; 
      TH1F * CA8PRUNE_FLAVOUR               ; 
      TH1F * CA8PRUNE_MASS_UNCORR           ; 
      TH1F * CA8PRUNE_MAXSUBJETBDISC        ; 
      TH1F * CA8PRUNE_MAXSUBJETBDISCFLAVOUR ; 

      TH1F * CA8PRUNE_CH_MULT    ; 
      TH1F * CA8PRUNE_NE_MULT    ; 
      TH1F * CA8PRUNE_CHEF       ; 
      TH1F * CA8PRUNE_NHEF       ; 
      TH1F * CA8PRUNE_CEEF       ; 
      TH1F * CA8PRUNE_NEEF       ; 
      TH1F * CA8PRUNE_CMEF       ; 

      TH1F * CMS_PT                    ; 
      TH1F * CMS_PHI                   ; 
      TH1F * CMS_ETA                   ; 
      TH1F * CMS_RAP                   ; 
      TH1F * CMS_MASS                  ; 
      TH1F * CMS_NCONST                ; 
      TH1F * CMS_BDISC                 ; 
      TH1F * CMS_AREA                  ; 
      TH1F * CMS_NTRACKS               ; 
      TH1F * CMS_FLAVOUR               ; 
      TH1F * CMS_MASS_UNCORR           ; 
      TH1F * CMS_MAXSUBJETBDISC        ; 
      TH1F * CMS_MAXSUBJETBDISCFLAVOUR ; 
      TH1F * CMS_GROOMMASS             ;
      TH1F * CMS_MINMASS               ;

      TH1F * CMSFJ_PT                    ; 
      TH1F * CMSFJ_PHI                   ; 
      TH1F * CMSFJ_ETA                   ; 
      TH1F * CMSFJ_RAP                   ; 
      TH1F * CMSFJ_MASS                  ; 
      TH1F * CMSFJ_NCONST                ; 
      TH1F * CMSFJ_BDISC                 ; 
      TH1F * CMSFJ_AREA                  ; 
      TH1F * CMSFJ_NTRACKS               ; 
      TH1F * CMSFJ_FLAVOUR               ; 
      TH1F * CMSFJ_MASS_UNCORR           ; 
      TH1F * CMSFJ_MAXSUBJETBDISC        ; 
      TH1F * CMSFJ_MAXSUBJETBDISCFLAVOUR ; 
      TH1F * CMSFJ_GROOMMASS             ;
      TH1F * CMSFJ_MINMASS               ;


      TH1F * HEP_PT                       ; 
      TH1F * HEP_PHI                      ; 
      TH1F * HEP_ETA                      ; 
      TH1F * HEP_RAP                      ; 
      TH1F * HEP_MASS                     ; 
      TH1F * HEP_NCONST                   ; 
      TH1F * HEP_BDISC                    ; 
      TH1F * HEP_AREA                     ; 
      TH1F * HEP_NTRACKS                  ; 
      TH1F * HEP_FLAVOUR                  ; 
      TH1F * HEP_MASS_UNCORR              ; 
      TH1F * HEP_MAXSUBJETBDISC           ; 
      TH1F * HEP_MAXSUBJETBDISCFLAVOUR    ; 
      TH1F * HEP_M123                     ;
      TH1F * HEP_M12                      ;
      TH1F * HEP_M13                      ;
      TH1F * HEP_M23                      ;
      TH1F * HEP_M23M123                  ;
      TH1F * HEP_ATANM13M12               ;
      TH2F * HEP_ATANM13M12_M23M123       ;


      TTree *AK15FiltJetTree;
      Float_t Mass;
      Float_t Pt;
      Float_t Eta;
      Float_t Rapidity;
      Float_t Phi;
      Float_t Bdisc;
      Float_t SubjetBdisc;


};


//
// constructors and destructor
//
JetSubstructureTester::JetSubstructureTester(const edm::ParameterSet& iConfig):
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGen"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedGen")))
{
  edm::Service<TFileService> fs;

  CA8PF_PT         = fs->make<TH1F>("CA8PF_PT",         "",100,0,4000   );
  CA8PF_PHI        = fs->make<TH1F>("CA8PF_PHI",        "",100,-3.2,3.2 );
  CA8PF_ETA        = fs->make<TH1F>("CA8PF_ETA",        "",100,-4,4     );
  CA8PF_RAP        = fs->make<TH1F>("CA8PF_RAP",        "",100,-4,4     );
  CA8PF_MASS       = fs->make<TH1F>("CA8PF_MASS",       "",100,0,500   );
  CA8PF_NCONST     = fs->make<TH1F>("CA8PF_NCONST",     "",1000,0,1000  );
  CA8PF_BDISC      = fs->make<TH1F>("CA8PF_BDISC",      "",100,0,1      );
  CA8PF_AREA       = fs->make<TH1F>("CA8PF_AREA",       "",100,0,10     );
  CA8PF_NTRACKS    = fs->make<TH1F>("CA8PF_NTRACKS",    "",200,0,200    );
  CA8PF_FLAVOUR    = fs->make<TH1F>("CA8PF_FLAVOUR",    "",24,0,24      );
  CA8PF_MASS_UNCORR= fs->make<TH1F>("CA8PF_MASS_UNCORR","",100,0,500   );
  CA8PF_CH_MULT    = fs->make<TH1F>("CA8PF_CH_MULT",    "",500,0,500   );    
  CA8PF_NE_MULT    = fs->make<TH1F>("CA8PF_NE_MULT",    "",500,0,500   );    
  CA8PF_CHEF       = fs->make<TH1F>("CA8PF_CHEF"   ,    "",100,0,1     );    
  CA8PF_NHEF       = fs->make<TH1F>("CA8PF_NHEF"   ,    "",100,0,1     );    
  CA8PF_CEEF       = fs->make<TH1F>("CA8PF_CEEF"   ,    "",100,0,1     );    
  CA8PF_NEEF       = fs->make<TH1F>("CA8PF_NEEF"   ,    "",100,0,1     );    
  CA8PF_CMEF       = fs->make<TH1F>("CA8PF_CMEF"   ,    "",100,0,1     );    

  CA8CHS_PT         = fs->make<TH1F>("CA8CHS_PT",         "",100,0,4000   );
  CA8CHS_PHI        = fs->make<TH1F>("CA8CHS_PHI",        "",100,-3.2,3.2 );
  CA8CHS_ETA        = fs->make<TH1F>("CA8CHS_ETA",        "",100,-4,4     );
  CA8CHS_RAP        = fs->make<TH1F>("CA8CHS_RAP",        "",100,-4,4     );
  CA8CHS_MASS       = fs->make<TH1F>("CA8CHS_MASS",       "",100,0,500   );
  CA8CHS_NCONST     = fs->make<TH1F>("CA8CHS_NCONST",     "",1000,0,1000  );
  CA8CHS_BDISC      = fs->make<TH1F>("CA8CHS_BDISC",      "",100,0,1      );
  CA8CHS_AREA       = fs->make<TH1F>("CA8CHS_AREA",       "",100,0,10     );
  CA8CHS_NTRACKS    = fs->make<TH1F>("CA8CHS_NTRACKS",    "",200,0,200    );
  CA8CHS_FLAVOUR    = fs->make<TH1F>("CA8CHS_FLAVOUR",    "",24,0,24      );
  CA8CHS_MASS_UNCORR= fs->make<TH1F>("CA8CHS_MASS_UNCORR","",100,0,500   );
  CA8CHS_CH_MULT    = fs->make<TH1F>("CA8CHS_CH_MULT",    "",500,0,500   );    
  CA8CHS_NE_MULT    = fs->make<TH1F>("CA8CHS_NE_MULT",    "",500,0,500   );    
  CA8CHS_CHEF       = fs->make<TH1F>("CA8CHS_CHEF"   ,    "",100,0,1     );    
  CA8CHS_NHEF       = fs->make<TH1F>("CA8CHS_NHEF"   ,    "",100,0,1     );    
  CA8CHS_CEEF       = fs->make<TH1F>("CA8CHS_CEEF"   ,    "",100,0,1     );    
  CA8CHS_NEEF       = fs->make<TH1F>("CA8CHS_NEEF"   ,    "",100,0,1     );    
  CA8CHS_CMEF       = fs->make<TH1F>("CA8CHS_CMEF"   ,    "",100,0,1     );    

  AK15_PT         = fs->make<TH1F>("AK15_PT",         "",100,0,4000   );
  AK15_PHI        = fs->make<TH1F>("AK15_PHI",        "",100,-3.2,3.2 );
  AK15_ETA        = fs->make<TH1F>("AK15_ETA",        "",100,-4,4     );
  AK15_RAP        = fs->make<TH1F>("AK15_RAP",        "",100,-4,4     );
  AK15_MASS       = fs->make<TH1F>("AK15_MASS",       "",100,0,1000   );
  AK15_NCONST     = fs->make<TH1F>("AK15_NCONST",     "",1000,0,1000  );
  AK15_BDISC      = fs->make<TH1F>("AK15_BDISC",      "",100,0,1      );
  AK15_AREA       = fs->make<TH1F>("AK15_AREA",       "",100,0,10     );
  AK15_NTRACKS    = fs->make<TH1F>("AK15_NTRACKS",    "",200,0,200    );
  AK15_FLAVOUR    = fs->make<TH1F>("AK15_FLAVOUR",    "",24,0,24      );
  AK15_MASS_UNCORR= fs->make<TH1F>("AK15_MASS_UNCORR","",100,0,1000   );

  AK15FILT_PT                   = fs->make<TH1F>("AK15FILT_PT",                   "",100,0,4000);
  AK15FILT_PHI                  = fs->make<TH1F>("AK15FILT_PHI",                  "",100,-3.2,3.2  );
  AK15FILT_ETA                  = fs->make<TH1F>("AK15FILT_ETA",                  "",100,-4,4  );
  AK15FILT_RAP                  = fs->make<TH1F>("AK15FILT_RAP",                  "",100,-4,4  );
  AK15FILT_MASS                 = fs->make<TH1F>("AK15FILT_MASS",                 "",100,0,400 );
  AK15FILT_NCONST               = fs->make<TH1F>("AK15FILT_NCONST",               "",10,0,10 );
  AK15FILT_BDISC                = fs->make<TH1F>("AK15FILT_BDISC",                "",100,0,1   );
  AK15FILT_AREA                 = fs->make<TH1F>("AK15FILT_AREA",                 "",100,0,10  );
  AK15FILT_NTRACKS              = fs->make<TH1F>("AK15FILT_NTRACKS",              "",200,0,200 );
  AK15FILT_FLAVOUR              = fs->make<TH1F>("AK15FILT_FLAVOUR",              "",24,0,24     );
  AK15FILT_MASS_UNCORR          = fs->make<TH1F>("AK15FILT_MASS_UNCORR",          "",100,0,400 );
  AK15FILT_MAXSUBJETBDISC       = fs->make<TH1F>("AK15FILT_MAXSUBJETBDISC",       "",100,0,1 );
  AK15FILT_MAXSUBJETBDISCFLAVOUR= fs->make<TH1F>("AK15FILT_MAXSUBJETBDISCFLAVOUR","",24,0,24 );

  CA8PRUNE_PT                     = fs->make<TH1F>("CA8PRUNE_PT",                   "",100,0,4000   );
  CA8PRUNE_PHI                    = fs->make<TH1F>("CA8PRUNE_PHI",                  "",100,-3.2,3.2 );
  CA8PRUNE_ETA                    = fs->make<TH1F>("CA8PRUNE_ETA",                  "",100,-4,4     );
  CA8PRUNE_RAP                    = fs->make<TH1F>("CA8PRUNE_RAP",                  "",100,-4,4     );
  CA8PRUNE_MASS                   = fs->make<TH1F>("CA8PRUNE_MASS",                 "",100,0,500   );
  CA8PRUNE_NCONST                 = fs->make<TH1F>("CA8PRUNE_NCONST",               "",1000,0,1000  );
  CA8PRUNE_BDISC                  = fs->make<TH1F>("CA8PRUNE_BDISC",                "",100,0,1      );
  CA8PRUNE_AREA                   = fs->make<TH1F>("CA8PRUNE_AREA",                 "",100,0,10     );
  CA8PRUNE_NTRACKS                = fs->make<TH1F>("CA8PRUNE_NTRACKS",              "",200,0,200    );
  CA8PRUNE_FLAVOUR                = fs->make<TH1F>("CA8PRUNE_FLAVOUR",              "",24,0,24      );
  CA8PRUNE_MASS_UNCORR            = fs->make<TH1F>("CA8PRUNE_MASS_UNCORR",          "",100,0,500   );
  CA8PRUNE_MAXSUBJETBDISC         = fs->make<TH1F>("CA8PRUNE_MAXSUBJETBDISC",       "",100,0,1      );
  CA8PRUNE_MAXSUBJETBDISCFLAVOUR  = fs->make<TH1F>("CA8PRUNE_MAXSUBJETBDISCFLAVOUR","",24,0,24 );
  
  CA8PRUNE_CH_MULT  = fs->make<TH1F>("CA8PRUNE_CH_MULT",            "",500,0,500   );    
  CA8PRUNE_NE_MULT  = fs->make<TH1F>("CA8PRUNE_NE_MULT",            "",500,0,500   );    
  CA8PRUNE_CHEF     = fs->make<TH1F>("CA8PRUNE_CHEF"   ,            "",100,0,1     );    
  CA8PRUNE_NHEF     = fs->make<TH1F>("CA8PRUNE_NHEF"   ,            "",100,0,1     );    
  CA8PRUNE_CEEF     = fs->make<TH1F>("CA8PRUNE_CEEF"   ,            "",100,0,1     );    
  CA8PRUNE_NEEF     = fs->make<TH1F>("CA8PRUNE_NEEF"   ,            "",100,0,1     );    
  CA8PRUNE_CMEF     = fs->make<TH1F>("CA8PRUNE_CMEF"   ,            "",100,0,1     );    

  CMS_PT                     = fs->make<TH1F>("CMS_PT",                   "",100,0,4000   );
  CMS_PHI                    = fs->make<TH1F>("CMS_PHI",                  "",100,-3.2,3.2 );
  CMS_ETA                    = fs->make<TH1F>("CMS_ETA",                  "",100,-4,4     );
  CMS_RAP                    = fs->make<TH1F>("CMS_RAP",                  "",100,-4,4     );
  CMS_MASS                   = fs->make<TH1F>("CMS_MASS",                 "",100,0,500   );
  CMS_NCONST                 = fs->make<TH1F>("CMS_NCONST",               "",1000,0,1000  );
  CMS_BDISC                  = fs->make<TH1F>("CMS_BDISC",                "",100,0,1      );
  CMS_AREA                   = fs->make<TH1F>("CMS_AREA",                 "",100,0,10     );
  CMS_NTRACKS                = fs->make<TH1F>("CMS_NTRACKS",              "",200,0,200    );
  CMS_FLAVOUR                = fs->make<TH1F>("CMS_FLAVOUR",              "",24,0,24      );
  CMS_MASS_UNCORR            = fs->make<TH1F>("CMS_MASS_UNCORR",          "",100,0,500   );
  CMS_MAXSUBJETBDISC         = fs->make<TH1F>("CMS_MAXSUBJETBDISC",       "",100,0,1      );
  CMS_MAXSUBJETBDISCFLAVOUR  = fs->make<TH1F>("CMS_MAXSUBJETBDISCFLAVOUR","",24,0,24 );
  CMS_GROOMMASS              = fs->make<TH1F>("CMS_GROOMMASS",            "",100,0,400    );
  CMS_MINMASS                = fs->make<TH1F>("CMS_MINMASS",              "",100,0,200    );
        
  CMSFJ_PT                     = fs->make<TH1F>("CMSFJ_PT",                   "",100,0,4000   );
  CMSFJ_PHI                    = fs->make<TH1F>("CMSFJ_PHI",                  "",100,-3.2,3.2 );
  CMSFJ_ETA                    = fs->make<TH1F>("CMSFJ_ETA",                  "",100,-4,4     );
  CMSFJ_RAP                    = fs->make<TH1F>("CMSFJ_RAP",                  "",100,-4,4     );
  CMSFJ_MASS                   = fs->make<TH1F>("CMSFJ_MASS",                 "",100,0,500   );
  CMSFJ_NCONST                 = fs->make<TH1F>("CMSFJ_NCONST",               "",1000,0,1000  );
  CMSFJ_BDISC                  = fs->make<TH1F>("CMSFJ_BDISC",                "",100,0,1      );
  CMSFJ_AREA                   = fs->make<TH1F>("CMSFJ_AREA",                 "",100,0,10     );
  CMSFJ_NTRACKS                = fs->make<TH1F>("CMSFJ_NTRACKS",              "",200,0,200    );
  CMSFJ_FLAVOUR                = fs->make<TH1F>("CMSFJ_FLAVOUR",              "",24,0,24      );
  CMSFJ_MASS_UNCORR            = fs->make<TH1F>("CMSFJ_MASS_UNCORR",          "",100,0,500   );
  CMSFJ_MAXSUBJETBDISC         = fs->make<TH1F>("CMSFJ_MAXSUBJETBDISC",       "",100,0,1      );
  CMSFJ_MAXSUBJETBDISCFLAVOUR  = fs->make<TH1F>("CMSFJ_MAXSUBJETBDISCFLAVOUR","",24,0,24 );
  CMSFJ_GROOMMASS              = fs->make<TH1F>("CMSFJ_GROOMMASS",            "",100,0,400    );
  CMSFJ_MINMASS                = fs->make<TH1F>("CMSFJ_MINMASS",              "",100,0,200    );
  
  HEP_PT                      = fs->make<TH1F>("HEP_PT",                   "",100,0,4000   );                              
  HEP_PHI                     = fs->make<TH1F>("HEP_PHI",                  "",100,-3.2,3.2 );                              
  HEP_ETA                     = fs->make<TH1F>("HEP_ETA",                  "",100,-4,4     );                              
  HEP_RAP                     = fs->make<TH1F>("HEP_RAP",                  "",100,-4,4     );                              
  HEP_MASS                    = fs->make<TH1F>("HEP_MASS",                 "",100,0,500   );                              
  HEP_NCONST                  = fs->make<TH1F>("HEP_NCONST",               "",1000,0,1000  );                              
  HEP_BDISC                   = fs->make<TH1F>("HEP_BDISC",                "",100,0,1      );                              
  HEP_AREA                    = fs->make<TH1F>("HEP_AREA",                 "",100,0,10     );                                
  HEP_NTRACKS                 = fs->make<TH1F>("HEP_NTRACKS",              "",200,0,200    );                                
  HEP_FLAVOUR                 = fs->make<TH1F>("HEP_FLAVOUR",              "",24,0,24      );                                
  HEP_MASS_UNCORR             = fs->make<TH1F>("HEP_MASS_UNCORR",          "",100,0,500   );                                
  HEP_MAXSUBJETBDISC          = fs->make<TH1F>("HEP_MAXSUBJETBDISC",       "",100,0,1      );                                
  HEP_MAXSUBJETBDISCFLAVOUR   = fs->make<TH1F>("HEP_MAXSUBJETBDISCFLAVOUR","",100,0,1      );                                
  HEP_M123                    = fs->make<TH1F>("HEP_M123",                 "",100,0,400);                                
  HEP_M12                     = fs->make<TH1F>("HEP_M12",                  "",100,0,200);                                
  HEP_M13                     = fs->make<TH1F>("HEP_M13",                  "",100,0,200);                                
  HEP_M23                     = fs->make<TH1F>("HEP_M23",                  "",100,0,200);                                
  HEP_M23M123                 = fs->make<TH1F>("HEP_M23M123",              "",100,0,2);                                
  HEP_ATANM13M12              = fs->make<TH1F>("HEP_ATANM13M12",           "",100,0,2);                                
  HEP_ATANM13M12_M23M123      = fs->make<TH2F>("HEP_ATANM13M12_M23M123",   "",100,0,2,100,0,2);

                 
                /////////////////////////////////////////////
  //TTree
  /////////////////////////////////////////////
  TFileDirectory subDir3 = fs->mkdir("Trees");
  AK15FiltJetTree = new TTree("AK15FiltJetTree","Tree for saving info of events with two tags");
  AK15FiltJetTree->Branch("Mass",             & Mass,           "Mass/F");
  AK15FiltJetTree->Branch("Pt",               & Pt,             "Pt/F");
  AK15FiltJetTree->Branch("Eta",              & Eta,            "Eta/F");
  AK15FiltJetTree->Branch("Rapidity",         & Rapidity,       "Rapidity/F");
  AK15FiltJetTree->Branch("Phi",              & Phi,            "Phi/F");
  AK15FiltJetTree->Branch("Bdisc",            & Bdisc,          "Bdisc/F");
  AK15FiltJetTree->Branch("SubjetBdisc",      & SubjetBdisc,    "SubjetBdisc/F");



                 
     

}


JetSubstructureTester::~JetSubstructureTester()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
JetSubstructureTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  using namespace fastjet;

  bool verbose = false;

  int run = iEvent.id().run();
  int event = iEvent.id().event();
  int lumi = iEvent.id().luminosityBlock();
  
  if (verbose)cout<<"\n\nAnalyze event "<<event<<" run "<<run <<" lumi "<<lumi<<endl;
  

  // CA8PF ungroomed
  edm::Handle<std::vector<pat::Jet> > CA8PF;
  iEvent.getByLabel( "patJetsCA8PF", CA8PF );
  int count_CA8PF = 0;
  for ( std::vector<pat::Jet>::const_iterator jetBegin = CA8PF->begin(), jetEnd = CA8PF->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
  {
    if (count_CA8PF > 2) break;
    double pt            = ijet->pt();
    double phi           = ijet->phi();
    double eta           = ijet->eta();
    double rapidity      = ijet->rapidity();
    double mass          = ijet->mass();
    int    nconst        = ijet->numberOfDaughters();
    double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
    double area          = ijet->jetArea();
    double ntracks       = ijet->associatedTracks().size();
    double flavour       = abs(ijet->partonFlavour());

    reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
    double massuncorr = uncorrJet.mass();
    
    count_CA8PF++;

    if (pt<400) continue;
    CA8PF_PT          ->Fill( pt);              
    CA8PF_PHI         ->Fill( phi);              
    CA8PF_ETA         ->Fill( eta);              
    CA8PF_RAP         ->Fill( rapidity);              
    CA8PF_MASS        ->Fill( mass);              
    CA8PF_NCONST      ->Fill( nconst);              
    CA8PF_BDISC       ->Fill( bdisc);              
    CA8PF_AREA        ->Fill( area);              
    CA8PF_NTRACKS     ->Fill( ntracks);              
    CA8PF_FLAVOUR     ->Fill( flavour);              
    CA8PF_MASS_UNCORR ->Fill( massuncorr);              
    
    CA8PF_CH_MULT     ->Fill( ijet->chargedMultiplicity() );              
    CA8PF_NE_MULT     ->Fill( ijet->neutralMultiplicity() );              
    CA8PF_CHEF        ->Fill( ijet->chargedHadronEnergyFraction() );              
    CA8PF_NHEF        ->Fill( ijet->neutralHadronEnergyFraction() );              
    CA8PF_CEEF        ->Fill( ijet->chargedEmEnergyFraction() );              
    CA8PF_NEEF        ->Fill( ijet->neutralEmEnergyFraction() );              
    CA8PF_CMEF        ->Fill( ijet->chargedMuEnergyFraction() );              
  
    if (verbose) cout<<"CA8PF " <<count_CA8PF<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<" JEC L1 "<<ijet->jecFactor("L1FastJet")<<" L2 "<< ijet->jecFactor("L2Relative")<<" L3 "<< ijet->jecFactor("L3Absolute") << " current "<<ijet->currentJECLevel()<<" uncor pt "<<uncorrJet.pt()<<endl;
  }


  // CA8CHS ungroomed
  edm::Handle<std::vector<pat::Jet> > CA8CHS;
  iEvent.getByLabel( "patJetsCA8CHS", CA8CHS );
  int count_CA8CHS = 0;
  for ( std::vector<pat::Jet>::const_iterator jetBegin = CA8CHS->begin(), jetEnd = CA8CHS->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
  {
    if (count_CA8CHS > 2) break;
    double pt            = ijet->pt();
    double phi           = ijet->phi();
    double eta           = ijet->eta();
    double rapidity      = ijet->rapidity();
    double mass          = ijet->mass();
    int    nconst        = ijet->numberOfDaughters();
    double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
    double area          = ijet->jetArea();
    double ntracks       = ijet->associatedTracks().size();
    double flavour       = abs(ijet->partonFlavour());

    reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
    double massuncorr = uncorrJet.mass();
    
    count_CA8CHS++;

    if (pt<400) continue;
    CA8CHS_PT          ->Fill( pt);              
    CA8CHS_PHI         ->Fill( phi);              
    CA8CHS_ETA         ->Fill( eta);              
    CA8CHS_RAP         ->Fill( rapidity);              
    CA8CHS_MASS        ->Fill( mass);              
    CA8CHS_NCONST      ->Fill( nconst);              
    CA8CHS_BDISC       ->Fill( bdisc);              
    CA8CHS_AREA        ->Fill( area);              
    CA8CHS_NTRACKS     ->Fill( ntracks);              
    CA8CHS_FLAVOUR     ->Fill( flavour);              
    CA8CHS_MASS_UNCORR ->Fill( massuncorr);              
    
    CA8CHS_CH_MULT     ->Fill( ijet->chargedMultiplicity() );              
    CA8CHS_NE_MULT     ->Fill( ijet->neutralMultiplicity() );              
    CA8CHS_CHEF        ->Fill( ijet->chargedHadronEnergyFraction() );              
    CA8CHS_NHEF        ->Fill( ijet->neutralHadronEnergyFraction() );              
    CA8CHS_CEEF        ->Fill( ijet->chargedEmEnergyFraction() );              
    CA8CHS_NEEF        ->Fill( ijet->neutralEmEnergyFraction() );              
    CA8CHS_CMEF        ->Fill( ijet->chargedMuEnergyFraction() );              
  
    if (verbose) cout<<"CA8CHS " <<count_CA8CHS<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<" JEC L1 "<<ijet->jecFactor("L1FastJet")<<" L2 "<< ijet->jecFactor("L2Relative")<<" L3 "<< ijet->jecFactor("L3Absolute") << " current "<<ijet->currentJECLevel()<<" uncor pt "<<uncorrJet.pt()<<endl;
  }


  // AK15 ungroomed
  edm::Handle<std::vector<pat::Jet> > AK15JETS;
  iEvent.getByLabel( "patJetsAK15PF", AK15JETS );
  int count_AK15JETS = 0;
  for ( std::vector<pat::Jet>::const_iterator jetBegin = AK15JETS->begin(), jetEnd = AK15JETS->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
  {
    if (count_AK15JETS > 2) break;
    double pt            = ijet->pt();
    double phi           = ijet->phi();
    double eta           = ijet->eta();
    double rapidity      = ijet->rapidity();
    double mass          = ijet->mass();
    int    nconst        = ijet->numberOfDaughters();
    double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
    double area          = ijet->jetArea();
    double ntracks       = ijet->associatedTracks().size();
    double flavour       = abs(ijet->partonFlavour());

    reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
    double massuncorr = uncorrJet.mass();
    
    count_AK15JETS++;

    if (pt<200) continue;
    AK15_PT          ->Fill( pt);              
    AK15_PHI         ->Fill( phi);              
    AK15_ETA         ->Fill( eta);              
    AK15_RAP         ->Fill( rapidity);              
    AK15_MASS        ->Fill( mass);              
    AK15_NCONST      ->Fill( nconst);              
    AK15_BDISC       ->Fill( bdisc);              
    AK15_AREA        ->Fill( area);              
    AK15_NTRACKS     ->Fill( ntracks);              
    AK15_FLAVOUR     ->Fill( flavour);              
    AK15_MASS_UNCORR ->Fill( massuncorr);              
     
    if (verbose) cout<<"UnFilt " <<count_AK15JETS<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<" JEC L1 "<<ijet->jecFactor("L1FastJet")<<" L2 "<< ijet->jecFactor("L2Relative")<<" L3 "<< ijet->jecFactor("L3Absolute") << " current "<<ijet->currentJECLevel()<<" uncor pt "<<uncorrJet.pt()<<endl;
  }

  // AK15 filtered
  edm::Handle<std::vector<pat::Jet> > AK15JETSFILTERED;
  iEvent.getByLabel( "patJetsAK15PFfilteredPacked", AK15JETSFILTERED );
  edm::Handle<std::vector<pat::Jet> > ak15filtsubjets;
  iEvent.getByLabel( "patJetsAK15PFfilteredSubjets", ak15filtsubjets );

  int count_AK15JETSFILTERED =0;
  for ( std::vector<pat::Jet>::const_iterator jetBegin = AK15JETSFILTERED->begin(), jetEnd = AK15JETSFILTERED->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
  {
    if (count_AK15JETSFILTERED >2 ) break;
    double pt            = ijet->pt();
    double phi           = ijet->phi();
    double eta           = ijet->eta();
    double rapidity      = ijet->rapidity();
    double mass          = ijet->mass();
    int    nconst      = ijet->numberOfDaughters();
    double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
    double area          = ijet->jetArea();
    double ntracks       = ijet->associatedTracks().size();
    double flavour       = abs(ijet->partonFlavour());

    reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
    double massuncorr = uncorrJet.mass();
     
    if (verbose) cout<<"Filt "<<count_AK15JETSFILTERED<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<endl;

    double maxbdisc = -20;
    double maxbdiscpartonflavour = -999;
    for (int i = 0; i < nconst; i++ ) {
      reco::Candidate const * subjet =  ijet->daughter(i);
      pat::Jet const * patsubjet = dynamic_cast<pat::Jet const *>(subjet);
      double sbdisc = patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags");
      double spartonflavour = abs(patsubjet->partonFlavour());
      if (sbdisc  > maxbdisc ) {maxbdisc = sbdisc; maxbdiscpartonflavour =spartonflavour; }
      if (verbose) cout<<"  subjet pt "<<patsubjet->pt()<<"  mass "<<patsubjet->mass()<<"  bdisc "<<patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags") <<" area "<< patsubjet->jetArea() <<endl;
    }
    count_AK15JETSFILTERED++;

    if (pt<200) continue;

    AK15FILT_PT                    ->Fill( pt);              
    AK15FILT_PHI                   ->Fill( phi);              
    AK15FILT_ETA                   ->Fill( eta);              
    AK15FILT_RAP                   ->Fill( rapidity);              
    AK15FILT_MASS                  ->Fill( mass);              
    AK15FILT_NCONST                ->Fill( nconst);              
    AK15FILT_BDISC                 ->Fill( bdisc);              
    AK15FILT_AREA                  ->Fill( area);              
    AK15FILT_NTRACKS               ->Fill( ntracks);              
    AK15FILT_FLAVOUR               ->Fill( flavour);              
    AK15FILT_MASS_UNCORR           ->Fill( massuncorr);
    AK15FILT_MAXSUBJETBDISC        ->Fill( maxbdisc);
    AK15FILT_MAXSUBJETBDISCFLAVOUR ->Fill( maxbdiscpartonflavour);

    Mass = mass;
    Pt = pt;
    Eta = eta;
    Rapidity = rapidity;
    Phi = phi;
    Bdisc = bdisc;
    SubjetBdisc = maxbdisc;
    AK15FiltJetTree->Fill();

  }


    edm::Handle<std::vector<pat::Jet> > CA8JETSPRUNED;
    iEvent.getByLabel( "patJetsCA8CHSprunedPacked", CA8JETSPRUNED );
    int count_CA8JETSPRUNED=0;
    for ( std::vector<pat::Jet>::const_iterator jetBegin = CA8JETSPRUNED->begin(), jetEnd = CA8JETSPRUNED->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
    {
      if (count_CA8JETSPRUNED >2 ) break;
      double pt            = ijet->pt();
      double phi           = ijet->phi();
      double eta           = ijet->eta();
      double rapidity      = ijet->rapidity();
      double mass          = ijet->mass();
      int    nconst        = ijet->numberOfDaughters();
      double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
      double area          = ijet->jetArea();
      double ntracks       = ijet->associatedTracks().size();
      double flavour       = abs(ijet->partonFlavour());

      reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
      double massuncorr = uncorrJet.mass();

      if (verbose) cout<<"Prune "<<count_CA8JETSPRUNED<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<endl;

      double maxbdisc = -20;  
      double maxbdiscpartonflavour = -999;
   
      for (int i = 0; i < nconst; i++ ) {
        reco::Candidate const * subjet =  ijet->daughter(i);
        pat::Jet const * patsubjet = dynamic_cast<pat::Jet const *>(subjet);        
                   
        double sbdisc = patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags");
        double spartonflavour = abs(patsubjet->partonFlavour());
        if (sbdisc  > maxbdisc ) {maxbdisc = sbdisc; maxbdiscpartonflavour =spartonflavour; }
        if (verbose) cout<<"  subjet pt "<<patsubjet->pt()<<"  mass "<<patsubjet->mass()<<"  bdisc "<<patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags") <<" area "<< patsubjet->jetArea() <<endl;
      }
      count_CA8JETSPRUNED++;
      cout<<"beforept"<<endl;
      if (pt<400) continue;

cout<<" CA8PRUNE_PT                    "<<endl;   CA8PRUNE_PT                    ->Fill( pt);              
cout<<" CA8PRUNE_PHI                   "<<endl;   CA8PRUNE_PHI                   ->Fill( phi);              
cout<<" CA8PRUNE_ETA                   "<<endl;   CA8PRUNE_ETA                   ->Fill( eta);              
cout<<" CA8PRUNE_RAP                   "<<endl;   CA8PRUNE_RAP                   ->Fill( rapidity);              
cout<<" CA8PRUNE_MASS                  "<<endl;   CA8PRUNE_MASS                  ->Fill( mass);              
cout<<" CA8PRUNE_NCONST                "<<endl;   CA8PRUNE_NCONST                ->Fill( nconst);              
cout<<" CA8PRUNE_BDISC                 "<<endl;   CA8PRUNE_BDISC                 ->Fill( bdisc);              
cout<<" CA8PRUNE_AREA                  "<<endl;   CA8PRUNE_AREA                  ->Fill( area);              
cout<<" CA8PRUNE_NTRACKS               "<<endl;   CA8PRUNE_NTRACKS               ->Fill( ntracks);              
cout<<" CA8PRUNE_FLAVOUR               "<<endl;   CA8PRUNE_FLAVOUR               ->Fill( flavour);              
cout<<" CA8PRUNE_MASS_UNCORR           "<<endl;   CA8PRUNE_MASS_UNCORR           ->Fill( massuncorr);
cout<<" CA8PRUNE_MAXSUBJETBDISC        "<<endl;   CA8PRUNE_MAXSUBJETBDISC        ->Fill( maxbdisc);
cout<<" CA8PRUNE_MAXSUBJETBDISCFLAVOUR "<<endl;   CA8PRUNE_MAXSUBJETBDISCFLAVOUR ->Fill( maxbdiscpartonflavour);


// cout<<" CA8PRUNE_CH_MULT               "<<endl;   CA8PRUNE_CH_MULT                   ->Fill( ijet->chargedMultiplicity() );              
// cout<<" CA8PRUNE_NE_MULT               "<<endl;   CA8PRUNE_NE_MULT                   ->Fill( ijet->neutralMultiplicity() );              
// cout<<" CA8PRUNE_CHEF                  "<<endl;   CA8PRUNE_CHEF                      ->Fill( ijet->chargedHadronEnergyFraction() );              
// cout<<" CA8PRUNE_NHEF                  "<<endl;   CA8PRUNE_NHEF                      ->Fill( ijet->neutralHadronEnergyFraction() );              
// cout<<" CA8PRUNE_CEEF                  "<<endl;   CA8PRUNE_CEEF                      ->Fill( ijet->chargedEmEnergyFraction() );              
// cout<<" CA8PRUNE_NEEF                  "<<endl;   CA8PRUNE_NEEF                      ->Fill( ijet->neutralEmEnergyFraction() );              
// cout<<" CA8PRUNE_CMEF                  "<<endl;   CA8PRUNE_CMEF                      ->Fill( ijet->chargedMuEnergyFraction() );              
 

    }

    edm::Handle<std::vector<pat::Jet> > CMSTTJETS;
    iEvent.getByLabel( "patJetsCMSTopTagCHSPacked", CMSTTJETS );
    int count_CMSTTJETS=0;
    for ( std::vector<pat::Jet>::const_iterator jetBegin = CMSTTJETS->begin(), jetEnd = CMSTTJETS->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
    {
      if (count_CMSTTJETS >2 ) break;
      double pt            = ijet->pt();
      double phi           = ijet->phi();
      double eta           = ijet->eta();
      double rapidity      = ijet->rapidity();
      double mass          = ijet->mass();
      int    nconst        = ijet->numberOfDaughters();
      double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
      double area          = ijet->jetArea();
      double ntracks       = ijet->associatedTracks().size();
      double flavour       = abs(ijet->partonFlavour());

      reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
      double massuncorr = uncorrJet.mass();


      double minmass   = -1;
      double maxbdisc  = -20;
      double maxbdiscpartonflavour = -999;
      math::XYZTLorentzVector SumSubjets;

      if (nconst>=3){
        reco::Candidate const * c0 =  ijet->daughter(0);
        reco::Candidate const * c1 =  ijet->daughter(1);
        reco::Candidate const * c2 =  ijet->daughter(2);

        pat::Jet const * subjet0 = dynamic_cast<pat::Jet const *>(c0);
        pat::Jet const * subjet1 = dynamic_cast<pat::Jet const *>(c1);
        pat::Jet const * subjet2 = dynamic_cast<pat::Jet const *>(c2);

        math::XYZTLorentzVector Pair01(subjet0->p4()+subjet1->p4());
        math::XYZTLorentzVector Pair02(subjet0->p4()+subjet2->p4());
        math::XYZTLorentzVector Pair12(subjet1->p4()+subjet2->p4());

        double min3 = std::min(Pair01.mass(), Pair02.mass() );
        minmass = std::min(min3, Pair12.mass() );             
      }

      if (verbose) cout<<"CMS "<<count_CMSTTJETS<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<" minmass "<<minmass<<endl;
      
      for (int i = 0; i < nconst; i++ ) {
        reco::Candidate const * subjet =  ijet->daughter(i);
        pat::Jet const * patsubjet = dynamic_cast<pat::Jet const *>(subjet);        
        double sbdisc = patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags");
        double spartonflavour = abs(patsubjet->partonFlavour());
        if (sbdisc  > maxbdisc ) {maxbdisc = sbdisc; maxbdiscpartonflavour =spartonflavour; }
        SumSubjets+=patsubjet->p4();
        if (verbose) cout<<"  subjet pt "<<patsubjet->pt()<<"  mass "<<patsubjet->mass()<<"  bdisc "<<patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags") <<" area "<< patsubjet->jetArea() <<endl;
      }
      double groommass = SumSubjets.mass();
      if (verbose)cout<<" groommass"<<groommass<<" maxbdisc "<<maxbdisc<<endl;

      count_CMSTTJETS++;
      if (pt<400) continue;

      CMS_PT                    ->Fill( pt);              
      CMS_PHI                   ->Fill( phi);              
      CMS_ETA                   ->Fill( eta);              
      CMS_RAP                   ->Fill( rapidity);              
      CMS_MASS                  ->Fill( mass);              
      CMS_NCONST                ->Fill( nconst);              
      CMS_BDISC                 ->Fill( bdisc);              
      CMS_AREA                  ->Fill( area);              
      CMS_NTRACKS               ->Fill( ntracks);              
      CMS_FLAVOUR               ->Fill( flavour);              
      CMS_MASS_UNCORR           ->Fill( massuncorr);
      CMS_MAXSUBJETBDISC        ->Fill( maxbdisc);
      CMS_MAXSUBJETBDISCFLAVOUR ->Fill( maxbdiscpartonflavour);
      CMS_GROOMMASS             ->Fill( groommass);      
      CMS_MINMASS               ->Fill( minmass);      

    }


    edm::Handle<std::vector<pat::Jet> > CMSTTJETSFJ;
    iEvent.getByLabel( "patJetsCMSTopTagFJCHSPacked", CMSTTJETSFJ );
    int count_CMSTTJETSFJ=0;
    for ( std::vector<pat::Jet>::const_iterator jetBegin = CMSTTJETSFJ->begin(), jetEnd = CMSTTJETSFJ->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
    {
      if (count_CMSTTJETSFJ >2 ) break;
      double pt            = ijet->pt();
      double phi           = ijet->phi();
      double eta           = ijet->eta();
      double rapidity      = ijet->rapidity();
      double mass          = ijet->mass();
      int    nconst        = ijet->numberOfDaughters();
      double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
      double area          = ijet->jetArea();
      double ntracks       = ijet->associatedTracks().size();
      double flavour       = abs(ijet->partonFlavour());

      reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
      double massuncorr = uncorrJet.mass();


      double minmass   = -1;
      double maxbdisc  = -20;
      double maxbdiscpartonflavour = -999;
      math::XYZTLorentzVector SumSubjets;

      if (nconst>=3){
        reco::Candidate const * c0 =  ijet->daughter(0);
        reco::Candidate const * c1 =  ijet->daughter(1);
        reco::Candidate const * c2 =  ijet->daughter(2);

        pat::Jet const * subjet0 = dynamic_cast<pat::Jet const *>(c0);
        pat::Jet const * subjet1 = dynamic_cast<pat::Jet const *>(c1);
        pat::Jet const * subjet2 = dynamic_cast<pat::Jet const *>(c2);

        math::XYZTLorentzVector Pair01(subjet0->p4()+subjet1->p4());
        math::XYZTLorentzVector Pair02(subjet0->p4()+subjet2->p4());
        math::XYZTLorentzVector Pair12(subjet1->p4()+subjet2->p4());

        double min3 = std::min(Pair01.mass(), Pair02.mass() );
        minmass = std::min(min3, Pair12.mass() );             
      }

      if (verbose) cout<<"CMSFJ "<<count_CMSTTJETSFJ<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<" minmass "<<minmass<<endl;
      
      for (int i = 0; i < nconst; i++ ) {
        reco::Candidate const * subjet =  ijet->daughter(i);
        pat::Jet const * patsubjet = dynamic_cast<pat::Jet const *>(subjet);        
        double sbdisc = patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags");
        double spartonflavour = abs(patsubjet->partonFlavour());
        if (sbdisc  > maxbdisc ) {maxbdisc = sbdisc; maxbdiscpartonflavour =spartonflavour; }
        SumSubjets+=patsubjet->p4();
        if (verbose) cout<<"  subjet pt "<<patsubjet->pt()<<"  mass "<<patsubjet->mass()<<"  bdisc "<<patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags") <<" area "<< patsubjet->jetArea() <<endl;
      }
      double groommass = SumSubjets.mass();
      if (verbose)cout<<" groommass"<<groommass<<" maxbdisc "<<maxbdisc<<endl;

      count_CMSTTJETSFJ++;
      if (pt<400) continue;

      CMSFJ_PT                    ->Fill( pt);              
      CMSFJ_PHI                   ->Fill( phi);              
      CMSFJ_ETA                   ->Fill( eta);              
      CMSFJ_RAP                   ->Fill( rapidity);              
      CMSFJ_MASS                  ->Fill( mass);              
      CMSFJ_NCONST                ->Fill( nconst);              
      CMSFJ_BDISC                 ->Fill( bdisc);              
      CMSFJ_AREA                  ->Fill( area);              
      CMSFJ_NTRACKS               ->Fill( ntracks);              
      CMSFJ_FLAVOUR               ->Fill( flavour);              
      CMSFJ_MASS_UNCORR           ->Fill( massuncorr);
      CMSFJ_MAXSUBJETBDISC        ->Fill( maxbdisc);
      CMSFJ_MAXSUBJETBDISCFLAVOUR ->Fill( maxbdiscpartonflavour);
      CMSFJ_GROOMMASS             ->Fill( groommass);      
      CMSFJ_MINMASS               ->Fill( minmass);      

    }

    edm::Handle<std::vector<pat::Jet> > HEPTTJETS;
    iEvent.getByLabel( "patJetsHEPTopTagCHSPacked", HEPTTJETS );
    int count_HEPTTJETS=0;
    for ( std::vector<pat::Jet>::const_iterator jetBegin = HEPTTJETS->begin(), jetEnd = HEPTTJETS->end(), ijet = jetBegin; ijet != jetEnd; ++ijet ) 
    {
      if (count_HEPTTJETS >2 ) break;
      double pt            = ijet->pt();
      double phi           = ijet->phi();
      double eta           = ijet->eta();
      double rapidity      = ijet->rapidity();
      double mass          = ijet->mass();
      int    nconst        = ijet->numberOfDaughters();
      double bdisc         = ijet->bDiscriminator("combinedSecondaryVertexBJetTags");
      double area          = ijet->jetArea();
      double ntracks       = ijet->associatedTracks().size();
      double flavour       = abs(ijet->partonFlavour());

      reco::Candidate::LorentzVector uncorrJet = ijet->correctedP4(0);
      double massuncorr = uncorrJet.mass();

      if (verbose) cout<<"HEP "<<count_HEPTTJETS<<"  pt "<<pt<<" mass "<<mass<<" nconst "<<nconst<<" ntracks "<<ntracks<<" bdisc "<<bdisc<< " area "<<area<<endl;
      double maxbdisc=-20;
      double maxbdiscpartonflavour = -999;
      for (int i = 0; i < nconst; i++ ) {
        reco::Candidate const * subjet =  ijet->daughter(i);
        pat::Jet const * patsubjet = dynamic_cast<pat::Jet const *>(subjet);
        double sbdisc = patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags");
        double spartonflavour = abs(patsubjet->partonFlavour());
        if (sbdisc  > maxbdisc ) {maxbdisc = sbdisc; maxbdiscpartonflavour =spartonflavour; }
        if (verbose) cout<<"  subjet pt "<<patsubjet->pt()<<"  mass "<<patsubjet->mass()<<"  bdisc "<<patsubjet->bDiscriminator("combinedSecondaryVertexBJetTags") <<" area "<< patsubjet->jetArea() <<endl;
      }
      count_HEPTTJETS++;
      if (pt<200) continue;

      if (nconst>=3){
        reco::Candidate const * c0 =  ijet->daughter(0);
        reco::Candidate const * c1 =  ijet->daughter(1);
        reco::Candidate const * c2 =  ijet->daughter(2);

        pat::Jet const * subjet0 = dynamic_cast<pat::Jet const *>(c0);
        pat::Jet const * subjet1 = dynamic_cast<pat::Jet const *>(c1);
        pat::Jet const * subjet2 = dynamic_cast<pat::Jet const *>(c2);

        math::XYZTLorentzVector Trip012(subjet0->p4()+subjet1->p4()+subjet2->p4());
        math::XYZTLorentzVector Pair01( subjet0->p4()+subjet1->p4() );
        math::XYZTLorentzVector Pair02( subjet0->p4()+subjet2->p4() );
        math::XYZTLorentzVector Pair12( subjet1->p4()+subjet2->p4() );

        double M12M012 = -1;
        double AtanM02M01 = -1;
        if ( Trip012.mass()!=0 ) M12M012     =  Pair12.mass() / Trip012.mass();
        if ( Pair01.mass()!=0 )  AtanM02M01  = atan( Pair02.mass() / Pair01.mass() ) ; 
        if (verbose) cout<<" Trip012 "<<Trip012.mass()<<" M12M012 "<<M12M012<<" AtanM02M01 "<<AtanM02M01<<endl;


         // print ', ptda1 = {0:6.2f}, ptda2 = {1:6.2f}, ptda3 = {2:6.2f}'.format( jet.daughter(0).pt(), jet.daughter(1).pt(), jet.daughter(2).pt() )
         //    print ', massda1 = {0:6.2f}, massda2 = {1:6.2f}, massda3 = {2:6.2f}'.format( jet.daughter(0).mass(), jet.daughter(1).mass(), jet.daughter(2).mass() )
         //    print ', bdiscda0 = {0:6.3f}, bdiscda1 = {1:6.3f}, bdiscda2 = {2:6.3f}'.format( jet.daughter(0).bDiscriminator("combinedSecondaryVertexBJetTags"), jet.daughter(1).bDiscriminator("combinedSecondaryVertexBJetTags"), jet.daughter(2).bDiscriminator("combinedSecondaryVertexBJetTags") ),
         //    m123 = (jet.daughter(0).p4()+jet.daughter(1).p4()+jet.daughter(2).p4()).mass()
         //    m12 = (jet.daughter(0).p4()+jet.daughter(1).p4()).mass()
         //    m13 = (jet.daughter(0).p4()+jet.daughter(2).p4()).mass()
         //    m23 = (jet.daughter(1).p4()+jet.daughter(2).p4()).mass()       
         //    print ', m123 = {0:6.2f}, m12 = {1:6.2f}, m13 = {2:6.2f}, m23 = {3:6.2f}'.format( m123, m12, m13, m23 )

        HEP_PT                    ->Fill( pt);                                
        HEP_PHI                   ->Fill( phi);                               
        HEP_ETA                   ->Fill( eta);                               
        HEP_RAP                   ->Fill( rapidity);                           
        HEP_MASS                  ->Fill( mass);                                                   
        HEP_NCONST                ->Fill( nconst);                                                   
        HEP_BDISC                 ->Fill( bdisc);                                      
        HEP_AREA                  ->Fill( area);                                      
        HEP_NTRACKS               ->Fill( ntracks);                                      
        HEP_FLAVOUR               ->Fill( flavour);                                      
        HEP_MASS_UNCORR           ->Fill( massuncorr);                        
        HEP_MAXSUBJETBDISC        ->Fill( maxbdisc);                        
        HEP_MAXSUBJETBDISCFLAVOUR ->Fill( maxbdiscpartonflavour);                        
        HEP_M123                  ->Fill( Trip012.mass());                         
        HEP_M12                   ->Fill( Pair01.mass());                         
        HEP_M13                   ->Fill( Pair02.mass());                         
        HEP_M23                   ->Fill( Pair12.mass());                         
        HEP_M23M123               ->Fill( M12M012);                         
        HEP_ATANM13M12            ->Fill( AtanM02M01);                         
        HEP_ATANM13M12_M23M123    ->Fill( AtanM02M01,M12M012);                        




      }
      if (verbose) cout<<"last line of hep loop"<<endl;

    }
    if (verbose) cout<<"done analyze"<<endl;

}


// ------------ method called once each job just before starting event loop  ------------
void 
JetSubstructureTester::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetSubstructureTester::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetSubstructureTester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetSubstructureTester);
