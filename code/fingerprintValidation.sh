#!/bin/bash
# Fingerprint the files the FE/BYM2 validation needs, to verify the laptop and the
# cluster (syvert1) copies match after rsync. Run from the repo root on BOTH
# machines, then compare the COMBINED/CODE/INPUTS hash lines (or `diff` the saved
# manifests). Works in Git Bash (Windows) and on Linux.
#
#   bash code/fingerprintValidation.sh                 # print manifest + hashes
#   bash code/fingerprintValidation.sh > fp_local.txt  # save to diff across hosts
#   # then on the cluster: bash code/fingerprintValidation.sh > fp_cluster.txt
#   #                      diff fp_local.txt fp_cluster.txt
#
# CODE files (code/*.R, *.cpp) are hashed with CR stripped so Windows CRLF vs Linux
# LF never causes a false mismatch. INPUT .RData are binary, hashed as-is. Compiled
# TMB artifacts (.o/.so/.dll) are intentionally EXCLUDED (legitimately differ per
# platform; rebuilt on the cluster by prepValidationFE.R).

cd "$(dirname "$0")/.." || exit 1

hasher() { if command -v sha256sum >/dev/null 2>&1; then sha256sum; else shasum -a 256; fi; }

# ---- file groups -------------------------------------------------------------
# CODE: validation + dependency sources (everything sourced via setup.R + the
# validation chain). Hash the whole code dir's R/cpp so nothing is missed.
codeFiles=$(ls code/*.R code/*.cpp 2>/dev/null | sort)
# INPUTS: ONLY the non-git data the validation actually needs (not all of
# savedOutput/global). This is exactly what setup.R load()s (19 files) + the 4
# extra files the validation chain reads (adm2Mat/admFinalMat/popMeanSDCal/
# predGridInputs) + the precomputed TMB inputs. Update this list if setup.R's or
# the validation chain's load() set changes (see: grep load code/setup.R).
inputFiles=$(cat <<'EOF'
savedOutput/global/NigeriaMapData.RData
savedOutput/global/covariates.RData
savedOutput/global/datMICS.RData
savedOutput/global/easpaNGA.RData
savedOutput/global/easpaNGAMICS.RData
savedOutput/global/easpaNGAed.RData
savedOutput/global/easpaNGAedMICS.RData
savedOutput/global/ed.RData
savedOutput/global/edMICS.RData
savedOutput/global/popMatNGA.RData
savedOutput/global/popMatNGAThresh.RData
savedOutput/global/popMatNGAed.RData
savedOutput/global/popMatNGAedThresh.RData
savedOutput/global/poppStratMICS.RData
savedOutput/global/poppStratMICSThresh.RData
savedOutput/global/poppaNGA.RData
savedOutput/global/poppsubNGA.RData
savedOutput/global/poppsubNGAThresh.RData
savedOutput/global/urbProps.RData
savedOutput/global/adm2Mat.RData
savedOutput/global/admFinalMat.RData
savedOutput/global/popMeanSDCal.RData
savedOutput/global/predGridInputs.RData
savedOutput/validation/fullMDMInputs_FE.RData
EOF
)

# ---- per-file lines ----------------------------------------------------------
# CODE (CR-stripped), INPUTS (raw). Format: <sha256>  <bytes>  <path>
codeLines=""
for f in $codeFiles; do
  [ -f "$f" ] || { codeLines+="MISSING  0  $f"$'\n'; continue; }
  h=$(tr -d '\r' < "$f" | hasher | awk '{print $1}')
  sz=$(tr -d '\r' < "$f" | wc -c | tr -d ' ')
  codeLines+="$h  $sz  $f"$'\n'
done
inputLines=""
for f in $inputFiles; do
  [ -f "$f" ] || { inputLines+="MISSING  0  $f"$'\n'; continue; }
  h=$(hasher < "$f" | awk '{print $1}')
  sz=$(wc -c < "$f" | tr -d ' ')
  inputLines+="$h  $sz  $f"$'\n'
done

# ---- output ------------------------------------------------------------------
echo "===== CODE (code/*.R, *.cpp; CR-stripped) ====="
printf "%s" "$codeLines"
echo "===== INPUTS (savedOutput .RData) ====="
printf "%s" "$inputLines"
echo "----- summary -----"
codeHash=$(printf "%s" "$codeLines"  | hasher | awk '{print $1}')
inpHash=$( printf "%s" "$inputLines" | hasher | awk '{print $1}')
allHash=$(printf "%s%s" "$codeLines" "$inputLines" | hasher | awk '{print $1}')
nCode=$(printf "%s" "$codeLines"  | grep -c '^[0-9a-f]')
nInp=$( printf "%s" "$inputLines" | grep -c '^[0-9a-f]')
nMiss=$(printf "%s%s" "$codeLines" "$inputLines" | grep -c '^MISSING')
echo "CODE    $codeHash  ($nCode files)"
echo "INPUTS  $inpHash  ($nInp files)"
echo "COMBINED $allHash"
[ "$nMiss" -gt 0 ] && echo "WARNING: $nMiss expected file(s) MISSING"
exit 0
