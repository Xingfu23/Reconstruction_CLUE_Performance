import FWCore.ParameterSet.Config as cms
process = cms.Process("DPG")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet( limit = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        #modify the following string to file:/afs/cern.ch/user/w/who/your/path/to/step3.root
            'file:/afs/cern.ch/user/x/xisu/private/local/Reconstruction_WorkSpace/Ana_File/Low_Density_SinglePhoton.root')
        )

process.analyzer = cms.EDAnalyzer('clue_performance', 
    layerclusters = cms.untracked.InputTag('hgcalLayerClusters'),
    tracks = cms.untracked.InputTag('generalTracks'),
    photons = cms.untracked.InputTag('photonsFromMultiCl'),
    multis = cms.untracked.InputTag('hgcalMultiClusters'),
    pfcandsticl = cms.untracked.InputTag('pfTICL'), 
    pfcands = cms.untracked.InputTag('particleFlow'),
    tracksters = cms.untracked.InputTag('ticlTrackstersMerge')
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("hgcal_clue_performance_LD.root"))
process.p = cms.Path(process.analyzer)