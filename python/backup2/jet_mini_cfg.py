## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")



process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoJets.Configuration.GenJetParticles_cff")
from TopQuarkAnalysis.TopPairBSM.filters_cff import applyFilters

process.source.fileNames = ['file:/eos/uscms/store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_1.root']
process.maxEvents.input = 10
process.out.outputCommands += [
    'keep *_*_*_*'
    ]
process.out.fileName = 'test_jet_mini.root'
