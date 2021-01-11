#!/bin/bash
# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
NON_DIST_REQS=$(cat ${IBM_CODE_PATH}/non_distributed_requirements.txt)
echo "
The following packages were not distributed in this container image because of their license:
${NON_DIST_REQS}
, as well as the packages installed using the following commands:
$(cat ${IBM_CODE_PATH}/non_distributed_requirements.sh)

By selecting '1) Install Python packages' you agree to trigger the download and installation of these packages now.
"

select answer in "Install Python packages" "Exit";
do
  case $answer in
    'Install Python packages' ) echo "Installing packages now..."; conda run -n base pip install ${NON_DIST_REQS}; /bin/bash ${IBM_CODE_PATH}/non_distributed_requirements.sh; exit;;
    * ) echo "Not installing packages. Exiting now..."; exit;;
  esac
done
