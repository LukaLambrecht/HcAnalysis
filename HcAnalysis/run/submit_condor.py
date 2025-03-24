import os
import sys
import six
import argparse

thisdir = os.path.abspath(os.path.dirname(__file__))
topdir = os.path.abspath(os.path.join(thisdir, '../'))
sys.path.append(topdir)

import tools.condortools as ct
from tools.datasettools import get_files


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samplelist', required=True, nargs='+')
    parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
    parser.add_argument('-m', '--maxfiles', default=-1, type=int,
      help='Number of files to process per entry in samplelist.')
    parser.add_argument('-n', '--nentries', default=-1, type=int,
      help='Number of entries to process per file.')
    parser.add_argument('--proxy', default=None)
    args = parser.parse_args()

    # set CMSSW config file
    config = os.path.abspath('../python/hcanalysis_cfg.py')

    # get CMSSW
    cmssw_version = os.getenv('CMSSW_BASE')
    print('Found following CMSSW base:')
    print(f'  - {cmssw_version}')
    if cmssw_version is None:
        msg = 'CMSSW version not set. Do cmsenv.'
        raise Exception(msg)

    # get proxy
    proxy = None
    if args.proxy is not None:
        proxy = os.path.abspath(args.proxy)
        print('Found following proxy:')
        print(f'  - {proxy}')

    # read samplelists
    datasets = []
    for samplelist in args.samplelist:
        with open(samplelist, 'r') as f:
            samples = [l.strip('\n') for l in f.readlines()]
            datasets += samples
    print(f'Read following datasets from provided samplelist ({len(datasets)}):')
    for d in datasets: print(f'  - {d}')

    # get input files
    inputfiles = get_files(datasets, maxfiles=args.maxfiles, verbose=True)

    # ask for confirmation
    print(f'Will submit {len(inputfiles)} jobs. Continue? (y/n)')
    go = six.moves.input()
    if go!='y': sys.exit()

    # loop over input files
    cmds = []
    for idx, f in enumerate(inputfiles):

        # make output file
        outputfile = os.path.join(args.outputdir, f'output_{idx+1}.root')
        
        # make the command
        cmd = f'cmsRun {config}'
        cmd += f' {f}'
        cmd += f' {args.nentries}'
        cmd += f' {outputfile}'
        cmds.append(cmd)

    # make output dir
    if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

    # submit the commands
    #for cmd in cmds: os.system(cmd)
    name = 'cjob_cmsRun'
    ct.submitCommandsAsCondorCluster(name, cmds,
        proxy=proxy,
        cmssw_version=cmssw_version,
        jobflavour='workday')
