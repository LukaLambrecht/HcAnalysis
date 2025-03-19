import os
import sys
import argparse

thisdir = os.path.abspath(os.path.dirname(__file__))
topdir = os.path.abspath(os.path.join(thisdir, '../'))
sys.path.append(topdir)

import tools.condortools as ct


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfiles', required=True, nargs='+')
    parser.add_argument('-o', '--outputdir', required=True, type=os.path.abspath)
    parser.add_argument('-c', '--config', required=True)
    parser.add_argument('--proxy', default=None)
    args = parser.parse_args()

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

    # loop over input files
    cmds = []
    for idx, f in enumerate(args.inputfiles):

        # make output file
        outputfile = os.path.join(args.outputdir, f'output_{idx+1}.root')
        
        # make the command
        cmd = f'cmsRun {args.config}'
        cmd += f' {f}'
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
