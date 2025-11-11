#!/usr/bin/env bash
# shellcheck disable=SC1090,SC1091

set -euo pipefail

if (( $# != 4 )); then
    echo "usage: $(basename "$0") <environment directory> <system> <config> <njobs>" 1>&2
    exit 1
fi
env_dir=$1
sys=$2
cfg=$3
njobs=$4

# Input validation
script_dir="$(dirname "$(readlink -f "${BASH_SOURCE:-$0}")")"
if [ ! -d "${script_dir}/systems/${sys}" ]; then
    echo "error: system directory '${sys}' does not exist"
    exit 1
fi
entry=$(jq ".configs[]|select(.name==\"${cfg}\")" "${script_dir}/systems/${sys}"/grid-config.json)
if [ -z "$entry" ]; then
  echo "error: config \"${cfg}\" does not exist for system '${sys}'"
  configs=$(jq -r ".configs[]|.name" "${script_dir}/systems/${sys}"/grid-config.json)
  echo "available configs:"
  for cfgname in ${configs[@]}; do
    echo "  ${cfgname}"
  done
  exit 1
fi

# Begin build
if [ -f "${script_dir}/systems/${sys}/files/shell-wrapper.sh" ]; then
    shellwrapper="${script_dir}/systems/${sys}/files/shell-wrapper.sh "
else
    shellwrapper=""
fi
${shellwrapper}./bootstrap-env.sh $1 $2
${shellwrapper}./build-grid.sh $1 $3 $4
${shellwrapper}./build-benchmark.sh $1 $3 $4

