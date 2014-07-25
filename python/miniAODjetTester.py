import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
   eventsToSkip = cms.untracked.VEventRange('1:1:1368-1:1:1368'), 
   fileNames = cms.untracked.vstring(
    	'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_1.root',
    	'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_2.root',
    	'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_3.root',
    	'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_4.root',
    	'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_5.root',
    	'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_6.root',
    	'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_7.root'
	)
)



process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START70_V6::All'

process.ak4PFJets.src = 'packedPFCandidates'


process.demo = cms.EDAnalyzer('MiniAODjetTester.cc',
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    pfCands = cms.InputTag("packedPFCandidates"),
	packedGen = cms.InputTag("packedGenParticles"),
	prunedGen = cms.InputTag("prunedGenParticles")
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("CheckMA.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.demo)
