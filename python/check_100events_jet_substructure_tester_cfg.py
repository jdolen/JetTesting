import FWCore.ParameterSet.Config as cms

###############################################
# SETUP
process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
   #eventsToSkip = cms.untracked.VEventRange('1:1:1368-1:1:1368'), 
   fileNames = cms.untracked.vstring(
      #'file:/uscmst1b_scratch/lpc1/lpcphys/jdolen/TestMiniAnalyzer706p1_1000events.root'
      #'file:/uscmst1b_scratch/lpc1/lpcphys/jdolen/TestJetsMiniAOD_706p1_big.root'
#      'file:/eos/uscms/store/user/jdolen/TestJetsAOD_706p1_verybig.root'
        'file:/eos/uscms/store/user/jdolen/check_100events_TestJetsMiniAOD_706p1.root'
		#'file:/uscmst1b_scratch/lpc1/lpcphys/jdolen/TestJetsAOD_706p1_big.root'
  )
)

###############################################
# HISTOGRAM MAKER
process.ana = cms.EDAnalyzer('JetSubstructureTester',
       jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    pfCands = cms.InputTag("packedPFCandidates"),
  packedGen = cms.InputTag("packedGenParticles"),   #packedGenParticles: all status == 1 particles (for gen jets)
  prunedGen = cms.InputTag("prunedGenParticles")    #prunedGenParticles: the interesting particles (for matching)
)

process.TFileService = cms.Service("TFileService",
      #fileName = cms.string("Jet_Substructure_Tester_MiniAOD_big.root"),
      fileName = cms.string("Jet_Substructure_Tester_MiniAOD_100events.root"),
      closeFileFast = cms.untracked.bool(True)
  )

###############################################

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
   process.ana
  )
#process.end = cms.EndPath(process.out)
