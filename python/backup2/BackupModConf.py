import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.options.allowUnscheduled = cms.untracked.bool(True)
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
   eventsToSkip = cms.untracked.VEventRange('1:1:1368-1:1:1368'), 
   fileNames = cms.untracked.vstring(
        #'root://cmsxrootd-site.fnal.gov//store/user/jstupak/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140703_030111/0000/miniAOD-prod_PAT_10.root'
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

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'

# GEN JETS
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

process.ak15GenJets = ak5GenJets.clone( rParam = cms.double(1.5),
                                           src = cms.InputTag("packedGenParticles"))

process.ak3GenJets = ak5GenJets.clone( rParam = cms.double(0.3),
                                           src = cms.InputTag("packedGenParticles"))

process.ak15GenJetsFiltered = ak5GenJets.clone(
  rParam = cms.double(1.5),
  src = cms.InputTag("packedGenParticles"),
  useFiltering = cms.bool(True),
  nFilt = cms.int32(3),
  rFilt = cms.double(0.3),
  writeCompound = cms.bool(True),
  jetCollInstanceName=cms.string("SubJets")
  )



process.ak4PFJets.src = 'packedPFCandidates'
process.ak15PFJets = process.ak4PFJets.clone(rParam = 1.5,  doAreaFastjet = True)

process.load('RecoJets.JetProducers.ak4PFJetsFiltered_cfi')
process.ak15PFJetsFiltered=process.ak4PFJetsFiltered.clone(rParam=1.5,  doAreaFastjet = True)

# As tracks are not stored in miniAOD, and b-tag fwk for CMSSW < 72X does not accept candidates
# we need to recreate tracks and pv for btagging in standard reco format:
process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

process.ak5JetTracksAssociatorAtVertexPF.jets = cms.InputTag("ak4PFJets")
process.ak5JetTracksAssociatorAtVertexPF.tracks = cms.InputTag("unpackedTracksAndVertices")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag("unpackedTracksAndVertices","secondary","")
process.combinedSecondaryVertex.trackMultiplicityMin = 1
process.ak15JetTracksAssociatorAtVertexPF=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak15PFJets'), coneSize = 1.5)

# for module in [process.patJetsAK15PF]:
#     module.addJetCharge = False
#     module.addBTagInfo = False #For some reason this has to be False
#     module.getJetMCFlavour = False
#     module.addAssociatedTracks = False



from PhysicsTools.PatAlgos.tools.jetTools import *
# from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
# from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

addJetCollection(
    process,
    labelName = 'AK15PF',
    jetSource = cms.InputTag('ak15PFJets'),
    algo = 'ak15',
    rParam = 1.5,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],

    )

process.patJetPartonMatchAK15PF.matched='prunedGenParticles'
process.patJetCorrFactorsAK15PF.primaryVertices = 'offlineSlimmedPrimaryVertices'
process.patJetGenJetMatchAK15PF.matched = 'ak15GenJets'#'slimmedGenJets'
process.patJetPartons.particles = "prunedGenParticles"

addJetCollection(
    process,
    labelName = 'AK15Filtered',
    jetSource = cms.InputTag('ak15PFJetsFiltered'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )

process.patJetPartonMatchAK15Filtered.matched='prunedGenParticles'
process.patJetCorrFactorsAK15Filtered.primaryVertices = 'offlineSlimmedPrimaryVertices'
process.patJetGenJetMatchAK15Filtered.matched = 'ak15GenJetsFiltered'#'slimmedGenJets'
process.patJetPartonMatchAK15Filtered.matched = "prunedGenParticles"

addJetCollection(
    process,
    labelName = 'AK15FilteredSubjets',
    jetSource = cms.InputTag('ak15PFJetsFiltered','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )

process.patJetPartonMatchAK15FilteredSubjets.matched='prunedGenParticles'
process.patJetCorrFactorsAK15FilteredSubjets.primaryVertices = 'offlineSlimmedPrimaryVertices'
process.patJetGenJetMatchAK15FilteredSubjets.matched = 'ak3GenJets'#slimmedGenJets'
#process.patJetGenJetMatchAK15PF.matched = "ak4GenJets"
process.patJetPartonMatchAK15FilteredSubjets.matched = "prunedGenParticles"


process.patJetsAK15FilteredPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsAK15Filtered" ),
    subjetSrc=cms.InputTag("patJetsAK15FilteredSubjets")
      )

# for ilabel in ['PatJetsCA8CMSTopTag',
#                'PatJetsCA8Pruned',
#                'PatJetsCA15HEPTopTag'] :
#     imerger = cms.EDProducer("BoostedJetMerger",
#                             jetSrc=cms.InputTag("good" + ilabel ),
#                             subjetSrc=cms.InputTag("selected" + ilabel + "Subjets")
#     )
#     setattr( process, 'good' + ilabel + 'Packed', imerger )



# patJetsAK15Filtered
# patJetsAK15FilteredSubjets





# #-------------------------------------
# ## N-subjettiness
# from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

# process.Njettiness = Njettiness.clone(
#     src = cms.InputTag("PFJetsCHS"),
#     cone = cms.double(options.jetRadius)
# )

# process.patJets.userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']




process.ana = cms.EDAnalyzer('MiniAnalyzer',
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    pfCands = cms.InputTag("packedPFCandidates"),
	packedGen = cms.InputTag("packedGenParticles"),   #packedGenParticles: all status == 1 particles (for gen jets)
	prunedGen = cms.InputTag("prunedGenParticles")    #prunedGenParticles: the interesting particles (for matching)
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("CheckMA.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('TestMod.root'),
                               outputCommands = cms.untracked.vstring([
                                'keep patJets_patJetsAK15PF_*_*',
                                'keep patJets_patJetsAK15Filtered_*_*',
                                'keep patJets_patJetsAK15FilteredSubjets_*_*',
                                'keep patJets_patJetsAK15FilteredPacked_*_*',
                                                                       ])
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
  process.ak3GenJets
  *process.ak15GenJets
  *process.ak15GenJetsFiltered
  *process.ak4PFJets
  *process.ak15PFJets
  *process.ak15PFJetsFiltered
  *process.patJetsAK15PF
  *process.patJetsAK15Filtered
  *process.patJetsAK15FilteredSubjets
  *process.patJetsAK15FilteredPacked
  *process.content
  *process.ana
  )
#process.end = cms.EndPath(process.out)
