import FWCore.ParameterSet.Config as cms
runOnMC=False
useData= not runOnMC
outputFile='ntuples_incxsec.root'
numEventsToRun=-1

###################get ak7 pu chs################## 
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True,jetAlgo='AK7',
					runOnMC=runOnMC, postfix=postfix,
					jetCorrections=('AK5PFchs', []),
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
#process = cms.Process("JetAnalyzerProcess")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(numEventsToRun) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.ERROR = cms.untracked.PSet(limit = cms.untracked.int32(5))

if  useData :
  process.GlobalTag.globaltag = cms.string( 'FT_53_V6C_AN3::All' )
else :
  process.GlobalTag.globaltag = cms.string( 'START53_V7G::All' )

if useData :
  fileNames = [
				'/store/data/Run2012C/JetHT/AOD/PromptReco-v2/000/198/969/58FEC838-16D1-E111-A403-001D09F2983F.root',
				'/store/data/Run2012C/JetHT/AOD/PromptReco-v2/000/198/969/4EAEA287-E3D0-E111-83F1-BCAEC5364C4C.root',
				'file:/uscmst1b_scratch/lpc1/3DayLifetime/tlibeiro/Run2012CJetMonAOD.root',
        '/store/data/Run2012C/JetMon/AOD/PromptReco-v2/000/198/955/CCA5E4B6-99CF-E111-A2E2-003048D2BCA2.root',
				'/store/data/Run2012C/JetMon/AOD/PromptReco-v2/000/198/941/7E513165-9BCF-E111-9AA4-5404A640A643.root'
  ]
else :
  fileNames = [
				'file:/uscmst1b_scratch/lpc1/3DayLifetime/tlibeiro/82B52E54-5E81-E111-A9A6-0018F3D096BA.root',
				'/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/0076C8E3-9AE1-E111-917C-003048D439AA.root',
  ]
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(fileNames)
)
from JetAnalyzer.JetAnalyzer.jetanalyzer_cfi import *
jetanalyzer.srcak7Jets    = cms.InputTag('pfJetsPFlow::PAT')
jetanalyzer.srcPrimaryVertex = cms.InputTag('goodOfflinePrimaryVertices')
jetanalyzer.metLabel         = cms.InputTag('pfMet')
jetanalyzer.runOnMC          = cms.bool(runOnMC)
jetanalyzer.hlTriggersTags   = cms.InputTag('TriggerResults::HLT')
jetanalyzer.l1TriggersTags   = cms.InputTag("gtDigis::RECO")
jetanalyzer.l1TriggerList    = cms.vstring("L1_SingleJet16", "L1_SingleJet36", "L1_SingleJet68", 
																					 "L1_SingleJet92", "L1_SingleJet128")
jetanalyzer.hlTriggerList    = cms.vstring("HLT_PFJet40","HLT_PFJet80","HLT_PFJet140", 
																						"HLT_PFJet200","HLT_PFJet260","HLT_PFJet320",
																						"HLT_PFJet400")
if runOnMC:
	jetanalyzer.JEC_GlobalTag    = cms.string('START53_V15')
else:
	jetanalyzer.JEC_GlobalTag    = cms.string('FT_53_V10_AN3')

process.jetanalyzer = jetanalyzer
##HLT######---------------------------------
process.HLT =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring(
														'HLT_PFJet40_v*','HLT_PFJet80_v*',   
														'HLT_PFJet140_v*','HLT_PFJet200_v*',   
														'HLT_PFJet260_v*','HLT_PFJet320_v*',   
													),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
)
###output TFile--------------------------------- 
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outputFile),
    closeFileFast = cms.untracked.bool(True)
)
process.options.wantSummary = False
process.p = cms.Path(
									process.HLT* 
									process.goodOfflinePrimaryVertices*
                  getattr(process,"patPF2PATSequence"+postfix)*
									process.jetanalyzer)
#####------Writing to data stream---------------
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('dummyOutput.root'),
    outputCommands = cms.untracked.vstring('drop *'),
    SelectEvents = cms.untracked.PSet(
                          SelectEvents = cms.vstring('p')
                                     )
)
process.e = cms.EndPath(process.out)

