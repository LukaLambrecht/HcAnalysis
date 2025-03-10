# Ntuplizer for Hc analysis

### How to setup
Download and run the installation script as follows:

```
wget https://raw.githubusercontent.com/LukaLambrecht/HcAnalysis/refs/heads/main/setup.sh
bash setup.sh
```

### How to obtain a test file for quick testing and development
See instructions [here](https://github.com/LukaLambrecht/HcAnalysis/blob/main/HcAnalysis/testfiles/README.md).

### How to run
Go to `$CMSSW_BASE/src/HcAnalysis/HcAnalysis/python` and run the following command:
```
cmsRun hcanalysis_cfg.py <path to input file> <name of output file>
```
where both the input file and output file should be `.root` files.
The input file should be in MiniAOD format, while the output file contains a plain TTree.

Note: make sure to have done `cmsenv` in the CMSSW `src` directory containing the ntuplizer before running.

### How to make modifications
Modify the analyzers in the `HcAnalysis/HcAnalysis/src` directory.
If needed, also edit the corresponding headers in the `HcAnalysis/HcAnalysis/interface` directory.
Then recompile with `scramv1 b`.

Note: make sure to have done `cmsenv` in the CMSSW `src` directory containing the ntuplizer before recompiling.

### References
- [example analyzer](https://github.com/nikhsub/HZZ4lc/blob/master/CMSSW_12_6_0/src/hplusc/HcAnalyzer/plugins/HcAnalyzer.cc) by Nikilesh.
- Writing your own EDAnalyzer [tutorial](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule).
- MiniAOD analysis [documentation](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015).
