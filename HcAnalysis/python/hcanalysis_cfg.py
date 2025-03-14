# Configuration file for running the Hc analysis

# imports
import os
import sys
import FWCore.ParameterSet.Config as cms

# initialize process
process = cms.Process("HcAnalysis")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

# set input file and number of events to process
# (note: use -1 to process all events in the input file)
inputfile = os.path.abspath(sys.argv[1])
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'file:{}'.format(inputfile)
        )
)

# define the processing steps and objects to use
process.analyzer = cms.EDAnalyzer('HcAnalysis',
  primaryVertices = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  prunedGenParticles = cms.untracked.InputTag("prunedGenParticles"),
  packedGenParticles = cms.untracked.InputTag("packedGenParticles"),
  packedPFCandidates = cms.untracked.InputTag('packedPFCandidates'),
  lostTracks = cms.untracked.InputTag('lostTracks')
)

# set output file
outputfile = sys.argv[2]
process.TFileService = cms.Service("TFileService",
  fileName = cms.string(outputfile),
)

process.p = cms.Path(process.analyzer)
