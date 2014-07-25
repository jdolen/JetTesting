#! /usr/bin/env python
import math
import ROOT
import sys
from DataFormats.FWLite import Events, Handle
#files = ["patTuple_tlbsm_train_TT_Tune4C_13TeV-pythia8-tauola_Spring14dr_PU20bx25_POSTLS170_V5-v1_tlbsm_71x_v1.root"]
#files = ["patTuple_tlbsm_train_RSGluonToTT_M-3000_Tune4C_13TeV-pythia8_PU20bx25_POSTLS170_V5-v1_tlbsm_71x_v1_File1.root"]
#files = ["file:/eos/uscms/store/user/jdolen/CSA14/patTuple_tlbsm_train_RSGluonToTT_M-3000_Tune4C_13TeV-pythia8_Fall13dr_tsg_PU40bx25_POSTLS162_V2-v1_tlbsm_71x_v1_File1.root"]
#files = ["file:/eos/uscms/store/user/jdolen/CSA14/patTuple_tlbsm_train_RSGluonToTT_M-3000_Tune4C_13TeV-pythia8_Spring14dr_PU_S14_POSTLS170_V6-v1_tlbsm_71x_v1_File1.root"]
#files = ["patTuple_tlbsm_train_RSGluonToTT_M-3000_Tune4C_13TeV-pythia8_PU20bx25_POSTLS170_V5-v1_tlbsm_71x_v1.root"]
#files = ["file:/eos/uscms/store/user/jdolen/CSA14/patTuple_tlbsm_train_RSGluonToTT_M-3000_Tune4C_13TeV-pythia8_Spring14dr_PU_S14_POSTLS170_V6-v1_tlbsm_71x_v1_File1Big.root"]
#files = ["root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_1.root"]
#files = ["TestMod.root"]
files = ["/uscmst1b_scratch/lpc1/lpcphys/jdolen/TestMiniAnalyzer706p1_1000events.root"]
#files = ["patTuple_tlbsm_train_tlbsm_71x_v1.root"]
#printGen = True
events = Events (files)
handle0  = Handle ("std::vector<pat::Jet>")
handle1  = Handle ("std::vector<pat::Jet>")
handle2  = Handle ("std::vector<pat::Jet>")
handle2  = Handle ("std::vector<pat::Jet>")
# handle3  = Handle ("std::vector<pat::Jet>")
# handle4  = Handle ("std::vector<pat::Muon>")
# handle5  = Handle ("std::vector<pat::Electron>")

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
label0 = ("patJetsAK15PF")
label1 = ("patJetsAK15Filtered")
label2 = ("patJetsAK15FilteredSubjets")
label2 = ("patJetsAK15FilteredSubjets")
# label3 = ("goodPatJetsCA15HEPTopTagPacked")
# label4 = ("selectedPatMuons")
# label5 = ("selectedPatElectrons")


# #f = ROOT.TFile("outplots_RSGluonToTT_M-3000_PU40bx25_POSTLS162_tsg.root", "RECREATE")
# f = ROOT.TFile("outplots_RSGluonToTT_M-3000_PU40bx50_POSTLS170_4680events.root", "RECREATE")
# #f = ROOT.TFile("outplots_RSGluonToTT_M-3000_PU20bx25_POSTLS170_4797events.root", "RECREATE")
# f.cd()
# HTT_m123        = ROOT.TH1F ("HTT_m123",        "HTT_m123",      100, 0, 400)
# HTT_pt          = ROOT.TH1F ("HTT_pt",          "HTT_pt",        100, 0, 1000)
# HTT_pt_tagged   = ROOT.TH1F ("HTT_pt_tagged",   "HTT_pt_tagged", 100, 0, 1000)
# CMSTT_mass      = ROOT.TH1F ("CMSTT_mass",      "CMSTT_mass",    100, 0, 400)
# CMSTT_minmass   = ROOT.TH1F ("CMSTT_minmass",   "CMSTT_minmass", 100, 0, 400)
# CMSTT_pt        = ROOT.TH1F ("CMSTT_pt",        "CMSTT_pt",      100, 0, 1000)
# CMSTT_pt_tagged = ROOT.TH1F ("CMSTT_pt_tagged", "CMSTT_pt",      100, 0, 1000)

# loop over events
i = 0
for event in events:
    i = i + 1
    print  '--------- Processing Event ' + str(i)

    print '---- ' + label0
    event.getByLabel (label0, handle0)
    jets0 = handle0.product()
    ijet = 0
    for jet in jets0 :
        print ("Jet {0:4.0f}, pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, " +
               "ndaughters = {5:3.0f}, vtxmass = {6:6.2f}, area = {7:6.2f}, L1 = {8:6.2f}, L2 = {9:6.2f}, L3 = {10:6.2f}, " +
               "currLevel = {11:s}, bdic = {12:10.2f},").format(
            ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters(), jet.userFloat('secvtxMass'),
            jet.jetArea(), jet.jecFactor("L1FastJet"), jet.jecFactor("L2Relative"), jet.jecFactor("L3Absolute"), jet.currentJECLevel(), jet.bDiscriminator("combinedSecondaryVertexBJetTags")
            )
        ijet += 1


    print '---- ' + label1
    event.getByLabel (label1, handle1)
    jets1 = handle1.product()
    ijet = 0
    for jet in jets1 :
        print ("Jet {0:4.0f}, pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, " +
               "nda = {5:3.0f}, vtxmass = {6:6.2f}, area = {7:6.2f}, L1 = {8:6.2f}, L2 = {9:6.2f}, L3 = {10:6.2f}, " +
               "currLevel = {11:s}").format(
            ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters(), jet.userFloat('secvtxMass'),
            jet.jetArea(), jet.jecFactor("L1FastJet"), jet.jecFactor("L2Relative"), jet.jecFactor("L3Absolute"), jet.currentJECLevel()
            )
        if jet.numberOfDaughters() >= 3 :
            print ', ptda1 = {0:6.2f}, ptda2 = {1:6.2f}, ptda3 = {2:6.2f}'.format( jet.daughter(0).pt(), jet.daughter(1).pt(), jet.daughter(2).pt() )
            print ', massda1 = {0:6.2f}, massda2 = {1:6.2f}, massda3 = {2:6.2f}'.format( jet.daughter(0).mass(), jet.daughter(1).mass(), jet.daughter(2).mass() )
            print ', bdiscda0 = {0:6.3f}, bdiscda1 = {1:6.3f}, bdiscda2 = {2:6.3f}'.format( jet.daughter(0).bDiscriminator("combinedSecondaryVertexBJetTags"), jet.daughter(1).bDiscriminator("combinedSecondaryVertexBJetTags"), jet.daughter(2).bDiscriminator("combinedSecondaryVertexBJetTags") ),
            m123 = (jet.daughter(0).p4()+jet.daughter(1).p4()+jet.daughter(2).p4()).mass()
            m12 = (jet.daughter(0).p4()+jet.daughter(1).p4()).mass()
            m13 = (jet.daughter(0).p4()+jet.daughter(2).p4()).mass()
            m23 = (jet.daughter(1).p4()+jet.daughter(2).p4()).mass()       
            print ', m123 = {0:6.2f}, m12 = {1:6.2f}, m13 = {2:6.2f}, m23 = {3:6.2f}'.format( m123, m12, m13, m23 )
        ijet += 1


#         if printGen :
#             genPt = 0.
#             if jet.genJetFwdRef().isNonnull() and jet.genJetFwdRef().isAvailable() :
#                 genPt = jet.genJetFwdRef().pt()
#             else :
#                 genPt = -1.0
#             print (", gen pt = {0:6.2f}").format( genPt )
#         else :
#             print ''
#         ijet += 1
    
#     print '---- ' + label1
#     # use getByLabel, just like in cmsRun
#     event.getByLabel (label1, handle1)
#     # get the product
#     jets1 = handle1.product()

#     ijet = 0
#     for jet in jets1 :
#         print 'Jet {0:4.0f}, pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, nda = {5:3.0f}'.format(
#             ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters()
#             ),
#         if jet.numberOfDaughters() > 1 :
#             print ', ptda1 = {0:6.2f}, ptda1 = {1:6.2f}'.format( jet.daughter(0).pt(), jet.daughter(1).pt() )
#         else :
#             print ''
#         ijet += 1


#     print '---- ' + label2
#     # use getByLabel, just like in cmsRun
#     event.getByLabel (label2, handle2)
#     # get the product
#     jets2 = handle2.product()

#     ijet = 0
#     for jet in jets2 :
#         tagged = 0
#         CMSTT_pt.Fill( jet.pt() )
#         if (jet.pt() >350): 
#           CMSTT_mass.Fill( jet.mass() )
#           CMSTT_minmass.Fill( jet.tagInfo('CATop').properties().minMass)
#         if  jet.mass() > 140 and  jet.mass() < 250 and jet.tagInfo('CATop').properties().minMass >50:
#             tagged =1
#             CMSTT_pt_tagged.Fill( jet.pt() )
#         print 'Jet {0:4.0f}, pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, nda = {5:3.0f}, topmass = {6:6.2f}, minmass = {7:6.2f}, tagged = {8:1d}'.format(
#             ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters(), jet.tagInfo('CATop').properties().topMass, jet.tagInfo('CATop').properties().minMass, tagged
#             )
#         ijet += 1
   
        
#     print '---- ' + label3
#     # use getByLabel, just like in cmsRun
#     event.getByLabel (label3, handle3)
#     # get the product
#     jets3 = handle3.product()

#     ijet = 0
#     for jet in jets3 :
#         print 'Jet {0:4.0f}, pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, nda = {5:3.0f}'.format(
#             ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters()
#             ),
#         HTT_pt.Fill( jet.pt() )       
#         if jet.numberOfDaughters() > 2 :
#             print ', ptda0 = {0:6.2f}, ptda1 = {1:6.2f}, ptda2 = {2:6.2f}'.format( jet.daughter(0).pt(), jet.daughter(1).pt(), jet.daughter(2).pt() ),
#             print ', bdiscda0 = {0:6.3f}, bdiscda1 = {1:6.3f}, bdiscda2 = {2:6.3f}'.format( jet.daughter(0).bDiscriminator("combinedSecondaryVertexBJetTags"), jet.daughter(1).bDiscriminator("combinedSecondaryVertexBJetTags"), jet.daughter(2).bDiscriminator("combinedSecondaryVertexBJetTags") ),
#             #Tag jet
#             m123 = (jet.daughter(0).p4()+jet.daughter(1).p4()+jet.daughter(2).p4()).mass()
#             m12 = (jet.daughter(0).p4()+jet.daughter(1).p4()).mass()
#             m13 = (jet.daughter(0).p4()+jet.daughter(2).p4()).mass()
#             m23 = (jet.daughter(1).p4()+jet.daughter(2).p4()).mass()
#             HTT_m123.Fill(m123)
#             massWindowLower, massWindowUpper = 0.85, 1.14
#             cutCondition2 = cutCondition3= 0.35
#             mlow, mhigh = 140.0, 250.0
#             topmass, wmass = 172.3, 80.4
#             rmin=massWindowLower*wmass/topmass
#             rmax=massWindowUpper*wmass/topmass
#             keep = 0
#             if(math.atan(m13/m12)>0.2 and math.atan(m13/m12)<1.3 and m23/m123>rmin and m23/m123<rmax):
#                 keep = 1
#             cond2left=pow(rmin,2)*(1+pow((m13/m12),2))
#             cond2cent=1-pow(m23/m123,2)
#             cond2right=pow(rmax,2)*(1+pow(m13/m12,2))
#             if(cond2left<cond2cent and cond2cent<cond2right and  m23/m123>cutCondition2):
#                 keep = 1
#             cond3left=pow(rmin,2)*(1+pow((m12/m13),2))
#             cond3cent=1-pow(m23/123,2)
#             cond3right=pow(rmax,2)*(1+pow(m12/m13,2))
#             if(cond3left<cond3cent and cond3cent<cond3right and m23/m123>cutCondition3):
#                 keep = 1
#             if( m123 < mlow or m123 > mhigh):
#                 keep = 0
#             print ', m123 = {0:6.2f}, m12 = {1:6.2f}, m13 = {2:6.2f}, m23 = {3:6.2f}, tagged = {4:1d}'.format( m123, m12, m13, m23, keep )
            
#             if keep == 1:
#                 HTT_pt_tagged.Fill( jet.pt() )         

#         else :
#             print ''
#         ijet += 1


#     print '---- ' + label4
#     # use getByLabel, just like in cmsRun
#     event.getByLabel (label4, handle4)
#     # get the product
#     muons1 = handle4.product()

#     imuon = 0
#     for muon in muons1 :
#         if not muon.isGlobalMuon() :
#             continue
#         print 'Muon {0:4.0f}, pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, nda = {5:3.0f}, chi2/dof = {6:6.2f}'.format(
#             imuon, muon.pt(), muon.eta(), muon.phi(), muon.mass(), muon.numberOfDaughters(), muon.normChi2()
#             )
#         imuon += 1

#     print '---- ' + label5
#     # use getByLabel, just like in cmsRun
#     event.getByLabel (label5, handle5)
#     # get the product
#     electrons1 = handle5.product()

#     ielectron = 0
#     for electron in electrons1 :
#         print 'Electron {0:4.0f}, pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, nda = {5:3.0f}, eidTight = {6:6.2f}'.format(
#             ielectron, electron.pt(), electron.eta(), electron.phi(), electron.mass(), electron.numberOfDaughters(), electron.electronID('eidTight')
#             )
#         ielectron += 1 


# f.cd()
# HTT_m123        .Write()
# HTT_pt          .Write()
# HTT_pt_tagged   .Write()
# CMSTT_mass      .Write()
# CMSTT_minmass   .Write()
# CMSTT_pt        .Write()
# CMSTT_pt_tagged .Write()
# f.Close()
