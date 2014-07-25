from WMCore.Configuration import Configuration config = Configuration()

config.section_("General") 
config.General.requestName = 'TTbar'

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'ModConf.py' 
config.JobType.allowNonProductionCMSSW = False

config.section_("Data") 
config.Data.inputDataset = ' /RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/jdolen-Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1-814812ec83fce2f620905d2bb30e9100/USER' 
config.Data.dbsUrl = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/' 
config.Data.splitting = 'FileBased' 
config.Data.unitsPerJob = 3 
config.Data.publication = False 
config.Data.publishDbsUrl = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' 
config.Data.publishDataName = 'Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1' 
config.Data.outlfn = '/store/user/jdolen/' 
config.Data.ignoreLocality = False

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC' 
