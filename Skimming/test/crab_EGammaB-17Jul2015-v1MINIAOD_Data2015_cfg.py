from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.workArea = '/user/padeken/CMSSW/CMSSW_7_4_6_patch2/src/PxlSkimmer/Skimming/test'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'EGammaB-17Jul2015-v1MINIAOD_Data2015_v2'

config.section_("JobType")
config.JobType.psetName = '/user/padeken/CMSSW/CMSSW_7_4_6_patch2/src/PxlSkimmer/Skimming/test/configs/data_miniAOD_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['EGammaB-17Jul2015-v1MINIAOD_Data2015.pxlio']
config.JobType.pyCfgParams = ['name=EGammaB-17Jul2015-v1MINIAOD_Data2015', 'datasetpath=/EGamma/Run2015B-17Jul2015-v1/MINIAOD', 'globalTag=74X_dataRun2_Prompt_v1']

config.section_("Data")
config.Data.inputDataset = '/EGamma/Run2015B-17Jul2015-v1/MINIAOD'
config.Data.publishDBS = 'phys03'
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.publishDataName = '74X_dataRun2_Prompt_v1'
config.Data.splitting = 'FileBased'
config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = '/store/user/padeken/PxlSkim/CMSSW_7_4_v1.2/EGammaB-17Jul2015-v1MINIAOD_Data2015/'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
config.Data.runRange = '251244-251883'

config.section_("Site")
config.Site.storageSite = 'T2_DE_RWTH'

config.section_("User")
config.User.voGroup = 'dcms'
