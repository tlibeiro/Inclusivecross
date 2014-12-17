##template for data
import FWCore.ParameterSet.Config as cms
runOnMC=True
useData= not runOnMC
outputFile='ntuples_incxsec.root'
numEventsToRun=10

#process = cms.Process("JetAnalyzerProcess")
##################get ak7 pu chs################## 
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True,jetAlgo='AK7',
					runOnMC=runOnMC, postfix=postfix,
					jetCorrections=('AK5PF', []),
					pvCollection=cms.InputTag('goodOfflinePrimaryVertices')
          )
postfix2=postfix+"AK5"
usePF2PAT(process,runPF2PAT=True,jetAlgo='AK5',
					runOnMC=runOnMC, postfix=postfix2,
					jetCorrections=('AK5PF', []),
					pvCollection=cms.InputTag('goodOfflinePrimaryVertices')
          )
process.pfPileUpPFlow.checkClosestZVertex = False
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )
#########end PF2PAT----------------------------------

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if  useData :
#  process.GlobalTag.globaltag = cms.string( 'FT_53_V6C_AN3::All' )
  process.GlobalTag.globaltag = cms.string( 'GR_P_V43D::All' )
else :
  process.GlobalTag.globaltag = cms.string( 'STARTHI53_V28::All' )

if useData :
  fileNames = [
       # 'file:/eos/uscms/store/user/tlibeiro/Run2013APP.root'
#			'/store/data/Run2013A/PPJet/RECO/PromptReco-v1/000/211/831/00000/5AF48420-1478-E211-A19F-001D09F28D54.root',
       'file:/eos/uscms/store/user/tlibeiro/Run211792_lumi915_917.root'
  ]
else :
  fileNames = [
        'file:/eos/uscms/store/user/tlibeiro/QCDpt30_276TeV_Reco.root',
#       '/store/himc/HiWinter13/QCD_Pt_220_TuneZ2_2p76TeV_pythia6/GEN-SIM-RECO/STARTHI53_V28-v1/00000/0018EE61-BF7A-E311-8F4C-848F69FD29B8.root'
  ]
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(fileNames)
)
from JetAnalyzer.JetAnalyzer.jetanalyzer_cfi import *
#jetanalyzer.srcak7Jets    = cms.InputTag('ak7PFJets::RECO')
jetanalyzer.srcak7Jets    = cms.InputTag('pfJetsPFlow::PAT')
jetanalyzer.srcGenJets    = cms.InputTag('ak7GenJets')
jetanalyzer.srcPrimaryVertex = cms.InputTag('goodOfflinePrimaryVertices')
jetanalyzer.metLabel         = cms.InputTag('pfMet')
jetanalyzer.runOnMC          = cms.bool(runOnMC)
jetanalyzer.hlTriggersTags   = cms.InputTag('TriggerResults::HLT')
jetanalyzer.l1TriggersTags   = cms.InputTag("gtDigis::RECO")
jetanalyzer.l1TriggerList    = cms.vstring("L1_SingleJet16_BptxAND","L1_SingleJet36")
jetanalyzer.hlTriggerList    = cms.vstring("HLT_PAJet20_NoJetID_v1","HLT_PAJet40_NoJetID_v1",
																					 "HLT_PAJet60_NoJetID_v1","HLT_PAJet80_NoJetID_v1",
																					 "HLT_PAJet100_NoJetID_v1","HLT_PAJet120_NoJetID_v1")
if runOnMC:
#	jetanalyzer.JEC_GlobalTag    = cms.string('JEC_STARTHI53')
	jetanalyzer.JEC_GlobalTag    = cms.string('Winter14_V5_MC')
else:
#	jetanalyzer.JEC_GlobalTag    = cms.string('JEC_STARTHI53')
#	jetanalyzer.JEC_GlobalTag    = cms.string('Summer13_V4_DATA')
	jetanalyzer.JEC_GlobalTag    = cms.string('Winter14_V5_DATA')

jetanalyzerak5 = jetanalyzer.clone(
											saveHLTInfo=cms.untracked.bool(True),
											useak5Corr =cms.untracked.bool(True),
											srcak7Jets = cms.InputTag('pfJetsPFlowAK5::PAT'),
#											srcak7Jets = cms.InputTag('ak5PFJets::RECO'),
											srcGenJets = cms.InputTag('ak5GenJets')
)
process.jetanalyzer = jetanalyzer
process.jetanalyzerak5 = jetanalyzerak5
##HLT######---------------------------------
process.HLT =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring(
#							 						  'HLT_PAJet20_NoJetID_v1',
														'HLT_PAJet40_NoJetID_v1',
														'HLT_PAJet60_NoJetID_v1','HLT_PAJet80_NoJetID_v1',
														'HLT_PAJet100_NoJetID_v1','HLT_PAJet120_NoJetID_v1'   
													),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(True) # throw exception on unknown path names
)
###output TFile--------------------------------- 
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outputFile),
    closeFileFast = cms.untracked.bool(True)
)
#process.options.wantSummary = False
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(numEventsToRun) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.ERROR = cms.untracked.PSet(limit = cms.untracked.int32(5))

#stupid filters---------------------------------------------------
process.load('RecoMET.METFilters.metFilters_cff')
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
     applyfilter = cms.untracked.bool  (True),
     debugOn     = cms.untracked.bool  (False),
     numtrack    = cms.untracked.uint32(10),
     thresh      = cms.untracked.double(0.25)
   )
#stupid filters---------------------------------------------------
process.myseq = cms.Sequence(
									process.HLT* 
#									process.metFilters*process.scrapingVeto*
									process.goodOfflinePrimaryVertices*
                  getattr(process,"patPF2PATSequence"+postfix)*
                  getattr(process,"patPF2PATSequence"+postfix2)*
									process.jetanalyzer
								  *process.jetanalyzerak5
								)
if runOnMC:
	process.myseq.remove(process.HLT)
process.p = cms.Path(process.myseq)
#####------Writing to data stream---------------
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('dummyOutput.root'),
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(
                          SelectEvents = cms.vstring('p')
                                     )
)
process.e = cms.EndPath(process.out)

