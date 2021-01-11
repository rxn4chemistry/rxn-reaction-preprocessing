#!/bin/bash
# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
set -Eeuo pipefail

current_user=$(whoami)
echo "Entering docker-entrypoint.sh as user:${current_user}..."

#app_user=rxn
#
## check if user is root
#if [[ "$(id -u)" = '0' ]]; then
#	# make sure we can write to stdout and stderr as $app_user
#	chown --dereference ${app_user} "/proc/$$/fd/1" "/proc/$$/fd/2" || :
#	exec gosu ${app_user} "$BASH_SOURCE" "$@"
#fi

echo "Running as user:${current_user}..."

# Install non-distributed dependencies
/bin/bash ${IBM_CODE_PATH}/download_non_distributed_requirements.sh

exec "$@"
