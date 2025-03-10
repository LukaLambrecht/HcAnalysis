# setup script for HcAnalysis package within CMSSW

# git settings
BRANCH=main
GIT_USER=LukaLambrecht
REPO_NAME=HcAnalysis

# other settings
RELEASE=CMSSW_14_0_4

# install CMSSW
echo "Installing $RELEASE..."
scram project CMSSW $RELEASE
cd $RELEASE/src
eval `scram runtime -sh`

# git clone this repository
echo "Cloning the $REPO_NAME repository..."
git cms-init
git clone https://github.com/$GIT_USER/$REPO_NAME
cd $CMSSW_BASE/src/$REPO_NAME
git checkout --track origin/$BRANCH

# compile and move into package
scramv1 b
echo "Setup finished"
