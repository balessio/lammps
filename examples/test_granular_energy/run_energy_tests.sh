#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LAMMPS_BIN="${LMP_BIN:-${LAMMPS:-${ROOT_DIR}/build/lmp}}"
if [[ ! -f "${LAMMPS_BIN}" || ! -x "${LAMMPS_BIN}" ]]; then
  LAMMPS_BIN="${ROOT_DIR}/build/lmp"
fi
OUT_BASE="${OUT_BASE:-${SCRIPT_DIR}/runs}"
PYTHON="${PYTHON:-python3.11}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/lammps_granular_energy_mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/lammps_granular_energy_xdg}"
export MPLBACKEND="${MPLBACKEND:-Agg}"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"
DUMP="${DUMP:-0}"
DUMP_EVERY="${DUMP_EVERY:-5000}"
LAMMPS_LOGS="${LAMMPS_LOGS:-0}"
QUICK="${QUICK:-0}"

if [[ "${QUICK}" == "1" ]]; then
  STEPS_TWO="${STEPS_TWO:-80000}"
  STEPS_SETTLE="${STEPS_SETTLE:-100000}"
else
  STEPS_TWO="${STEPS_TWO:-300000}"
  STEPS_SETTLE="${STEPS_SETTLE:-2000000}"
fi

AVEWINDOW="${AVEWINDOW:-100}"
DT_TWO="${DT_TWO:-1.0e-7}"
DT_SETTLE="${DT_SETTLE:-1.0e-7}"
PAIR_CUTOFF="${PAIR_CUTOFF:-0.006}"
WALL_TOP="${WALL_TOP:-0.080}"

if [[ -n "${CASES:-}" ]]; then
  read -r -a CASE_LIST <<< "${CASES}"
else
  CASE_LIST=(hooke hooke_history hertz_history mindlin mindlin_force mindlin_rescale_force)
fi

if [[ -n "${SCENARIOS:-}" ]]; then
  read -r -a SCENARIO_LIST <<< "${SCENARIOS}"
else
  SCENARIO_LIST=(pair wall)
fi

law_for_case() {
  case "$1" in
    elastic)
      echo "hooke 2.0e4 0.0 tangential linear_nohistory 0.0 0.0 damping velocity dissipative_heat 1.0 1.0 1.0"
      ;;
    hooke)
      echo "hooke 2.0e4 10.0 tangential linear_nohistory 0.0 0.0 damping velocity dissipative_heat 1.0 1.0 1.0 limit_damping"
      ;;
    hooke_history)
      echo "hooke 2.0e4 10.0 tangential linear_history 5.7e3 1.0 0.3 damping velocity dissipative_heat 1.0 1.0 1.0 limit_damping"
      ;;
    hertz_history)
      echo "hertz 2.0e5 5.0 tangential linear_history 5.7e3 0.5 0.3 damping velocity dissipative_heat 1.0 1.0 1.0 limit_damping"
      ;;
    mindlin)
      echo "hertz/material 5.0e5 0.2 0.3 tangential mindlin NULL 0.5 0.4 damping tsuji dissipative_heat 1.0 1.0 1.0 limit_damping"
      ;;
    mindlin_force)
      echo "hertz/material 5.0e5 0.2 0.3 tangential mindlin/force NULL 0.5 0.4 damping tsuji dissipative_heat 1.0 1.0 1.0 limit_damping"
      ;;
    mindlin_rescale_force)
      echo "hertz/material 5.0e5 0.2 0.3 tangential mindlin_rescale/force NULL 0.5 0.4 damping tsuji dissipative_heat 1.0 1.0 1.0 limit_damping"
      ;;
    dissipative)
      echo "hertz/material 5.0e5 0.2 0.3 tangential mindlin_rescale/force NULL 0.5 0.4 damping tsuji dissipative_heat 1.0 1.0 1.0 limit_damping"
      ;;
    *)
      echo "Unknown case: $1" >&2
      return 1
      ;;
  esac
}

write_case_files() {
  local scenario="$1"
  local case_name="$2"
  local outdir="$3"
  local law
  law="$(law_for_case "${case_name}")"
  mkdir -p "${outdir}"
  if [[ "${DUMP}" != "0" ]]; then
    mkdir -p "${outdir}/dumpfiles"
  fi
  printf 'pair_coeff * * %s cutoff %s\n' "${law}" "${PAIR_CUTOFF}" > "${outdir}/pair.inc"
  printf 'fix wall all wall/gran granular %s zplane 0.0 %s contacts\n' "${law}" "${WALL_TOP}" > "${outdir}/wall.inc"
  cat > "${outdir}/metadata.json" <<EOF
{
  "scenario": "${scenario}",
  "case": "${case_name}",
  "law": "${law}",
  "pair_cutoff": ${PAIR_CUTOFF},
  "wall_top": ${WALL_TOP},
  "dump": ${DUMP},
  "avewindow": ${AVEWINDOW}
}
EOF
}

run_case() {
  local scenario="$1"
  local case_name="$2"
  local outdir="${OUT_BASE}/${scenario}/${case_name}"
  write_case_files "${scenario}" "${case_name}" "${outdir}"

  echo "Running ${scenario}/${case_name}"
  local -a lmp_io_args
  if [[ "${LAMMPS_LOGS}" == "1" ]]; then
    lmp_io_args=(-log "${outdir}/log.lammps" -screen "${outdir}/screen.txt")
  else
    lmp_io_args=(-log none -screen none)
  fi

  if [[ "${scenario}" == "pair" ]]; then
    "${LAMMPS_BIN}" -in "${SCRIPT_DIR}/in.two_grains_periodic" \
      "${lmp_io_args[@]}" \
      -var outdir "${outdir}" -var pair_file "${outdir}/pair.inc" \
      -var runSteps "${STEPS_TWO}" -var avewindow "${AVEWINDOW}" \
      -var dump "${DUMP}" -var dumpEvery "${DUMP_EVERY}" -var dt "${DT_TWO}"
    "${PYTHON}" "${SCRIPT_DIR}/plot_energy.py" "${outdir}" --dt "${DT_TWO}" --avewindow "${AVEWINDOW}" --title "${scenario}/${case_name}"
  elif [[ "${scenario}" == "wall" ]]; then
    "${LAMMPS_BIN}" -in "${SCRIPT_DIR}/in.settling_wall" \
      "${lmp_io_args[@]}" \
      -var outdir "${outdir}" -var pair_file "${outdir}/pair.inc" -var wall_file "${outdir}/wall.inc" \
      -var runSteps "${STEPS_SETTLE}" -var avewindow "${AVEWINDOW}" \
      -var dump "${DUMP}" -var dumpEvery "${DUMP_EVERY}" -var dt "${DT_SETTLE}"
    "${PYTHON}" "${SCRIPT_DIR}/plot_energy.py" "${outdir}" --dt "${DT_SETTLE}" --avewindow "${AVEWINDOW}" --title "${scenario}/${case_name}"
  else
    echo "Unknown scenario: ${scenario}" >&2
    return 1
  fi
}

mkdir -p "${OUT_BASE}"

for scenario in "${SCENARIO_LIST[@]}"; do
  for case_name in "${CASE_LIST[@]}"; do
    run_case "${scenario}" "${case_name}"
  done
done

echo "Wrote results under ${OUT_BASE}"
