#!/usr/bin/env bash
set -euo pipefail

# --- Parse args (pass-through to the executable) ---
# We accept the same flags you use locally and forward them to the binary.
ARGS=("$@")

echo "[INFO] Hostname: $(hostname)"
echo "[INFO] PWD: $(pwd)"
echo "[INFO] LS: "; ls -alh > /dev/null 2>&1

# --- Proxy checks ---
if [[ -n "${X509_USER_PROXY:-}" ]]; then
  echo "[INFO] X509_USER_PROXY is set to: ${X509_USER_PROXY}"
  if command -v voms-proxy-info > /dev/null 2>&1; then
    voms-proxy-info -all || echo "[WARN] Could not read proxy details" > /dev/null 2>&1
  else
    echo "[INFO] voms-proxy-info not found; skipping verbose proxy dump"
  fi
else
  echo "[WARN] X509_USER_PROXY is not set. If your job reads protected data (e.g. via XRootD), set and forward a proxy."
fi

# --- Make sure the binary is executable ---
chmod +x process_HiggsAnalysis_wCutsExtended_NoPairing_v0

# If you need CMSSW/LCG env, uncomment and adapt one of the blocks below.
# --- Example: CMSSW runtime (adjust to your CMSSW release) ---
# source /cvmfs/cms.cern.ch/cmsset_default.sh
# export SCRAM_ARCH=el9_amd64_gcc13   # example; change to your arch
# CMSSW_BASE_LOCAL="$PWD/CMSSW_14_1_0"
# if [ -d "$CMSSW_BASE_LOCAL" ]; then
#   echo "[INFO] Using local CMSSW at $CMSSW_BASE_LOCAL"
# else
#   echo "[WARN] CMSSW not bundled; relying on static binary or LCG"
# fi
# eval "$(scramv1 runtime -sh)" || true

# --- Example: LCG view (adjust version if needed) ---
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh || true
set -u

echo "[INFO] Starting job at $(date)"
echo "[INFO] Command: ./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 ${ARGS[*]}"

./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 "${ARGS[@]}"
EXIT_CODE=$?

echo "[INFO] Finished with exit code: ${EXIT_CODE} at $(date)"
exit ${EXIT_CODE}
