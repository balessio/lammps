#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LAMMPS_BIN="${LMP_BIN:-${LAMMPS:-${ROOT_DIR}/build/lmp}}"
if [[ ! -f "${LAMMPS_BIN}" || ! -x "${LAMMPS_BIN}" ]]; then
  LAMMPS_BIN="${ROOT_DIR}/build/lmp"
fi

OUT_BASE="${OUT_BASE:-${SCRIPT_DIR}/runs_shear}"
PYTHON="${PYTHON:-python3.11}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/lammps_granular_energy_mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/lammps_granular_energy_xdg}"
export MPLBACKEND="${MPLBACKEND:-Agg}"

mkdir -p "${OUT_BASE}" "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

DUMP="${DUMP:-0}"
DUMP_EVERY="${DUMP_EVERY:-10000}"
LAMMPS_LOGS="${LAMMPS_LOGS:-0}"
DT="${DT:-2.0e-7}"
PREP_STEPS="${PREP_STEPS:-${COMPACT_STEPS:-240000}}"
SHEAR_STEPS="${SHEAR_STEPS:-440000}"
AVEWINDOW="${AVEWINDOW:-200}"
THERMO_EVERY="${THERMO_EVERY:-10000}"
PAIR_CUTOFF="${PAIR_CUTOFF:-0.007}"
PRESSURE_2D="${PRESSURE_2D:-0.45}"
PRESSURE_3D="${PRESSURE_3D:-500.0}"
SHEAR_GAP_2D="${SHEAR_GAP_2D:-0.044}"
SHEAR_GAP_3D="${SHEAR_GAP_3D:-0.052}"
SERVO_GAIN_2D="${SERVO_GAIN_2D:-5.0}"
SERVO_GAIN_3D="${SERVO_GAIN_3D:-0.0050}"
SERVO_MAX="${SERVO_MAX:-1.0}"
REUSE_PACKING="${REUSE_PACKING:-0}"

LAW="hertz/material 5.0e5 0.2 0.3 tangential mindlin_rescale/force NULL 0.5 0.4 damping tsuji dissipative_heat 1.0 1.0 1.0 limit_damping"

run_case() {
  local scenario="$1"
  local prep_input="$2"
  local shear_input="$3"
  local pressure="$4"
  local shear_gap="$5"
  local servo_gain="$6"
  local outdir="${OUT_BASE}/${scenario}"
  local restart_file="${outdir}/packed.restart"
  local data_file="${outdir}/packed.data"

  mkdir -p "${outdir}"
  if [[ "${DUMP}" != "0" ]]; then
    mkdir -p "${outdir}/prep_dumpfiles" "${outdir}/dumpfiles"
  fi

  printf 'pair_coeff * * %s cutoff %s\n' "${LAW}" "${PAIR_CUTOFF}" > "${outdir}/pair.inc"
  cat > "${outdir}/metadata.json" <<EOF
{
  "scenario": "${scenario}",
  "law": "${LAW}",
  "dt": ${DT},
  "prep_steps": ${PREP_STEPS},
  "shear_steps": ${SHEAR_STEPS},
  "total_shear": 1.0,
  "time_unit": "s",
  "heat_unit": "J",
  "heat_rate_unit": "W",
  "pressure": ${pressure},
  "servo_gain": ${servo_gain},
  "servo_max": ${SERVO_MAX},
  "restart_file": "packed.restart",
  "data_file": "packed.data",
  "dump": ${DUMP},
  "avewindow": ${AVEWINDOW}
}
EOF

  if [[ "${REUSE_PACKING}" != "1" || ! -f "${restart_file}" ]]; then
    echo "Preparing ${scenario}"
    local -a prep_io_args
    if [[ "${LAMMPS_LOGS}" == "1" ]]; then
      prep_io_args=(-log "${outdir}/prep.log.lammps" -screen "${outdir}/prep.screen.txt")
    else
      prep_io_args=(-log none -screen none)
    fi
    "${LAMMPS_BIN}" -in "${SCRIPT_DIR}/${prep_input}" \
      "${prep_io_args[@]}" \
      -var outdir "${outdir}" -var pair_file "${outdir}/pair.inc" \
      -var dt "${DT}" -var pressure "${pressure}" -var prepSteps "${PREP_STEPS}" \
      -var avewindow "${AVEWINDOW}" -var thermoEvery "${THERMO_EVERY}" \
      -var dump "${DUMP}" -var dumpEvery "${DUMP_EVERY}" \
      -var servoGain "${servo_gain}" -var servoMax "${SERVO_MAX}" \
      -var restart_file "${restart_file}" -var data_file "${data_file}"
  else
    echo "Reusing ${scenario} packing"
  fi

  echo "Shearing ${scenario}"
  local -a shear_io_args
  if [[ "${LAMMPS_LOGS}" == "1" ]]; then
    shear_io_args=(-log "${outdir}/log.lammps" -screen "${outdir}/screen.txt")
  else
    shear_io_args=(-log none -screen none)
  fi
  "${LAMMPS_BIN}" -in "${SCRIPT_DIR}/${shear_input}" \
    "${shear_io_args[@]}" \
    -var outdir "${outdir}" -var pair_file "${outdir}/pair.inc" \
    -var dt "${DT}" -var pressure "${pressure}" -var shearGap "${shear_gap}" \
    -var restart_file "${restart_file}" \
    -var shearSteps "${SHEAR_STEPS}" -var avewindow "${AVEWINDOW}" \
    -var thermoEvery "${THERMO_EVERY}" -var dump "${DUMP}" -var dumpEvery "${DUMP_EVERY}" \
    -var servoGain "${servo_gain}" -var servoMax "${SERVO_MAX}"

  "${PYTHON}" "${SCRIPT_DIR}/plot_shear_heat.py" "${outdir}" \
    --dt "${DT}" --avewindow "${AVEWINDOW}" --compact-steps 0 \
    --title "${scenario}: pressure-loaded shear heat"
}

if [[ -n "${SHEAR_CASES:-}" ]]; then
  read -r -a SHEAR_CASE_LIST <<< "${SHEAR_CASES}"
else
  SHEAR_CASE_LIST=(shear)
fi

for shear_case in "${SHEAR_CASE_LIST[@]}"; do
  case "${shear_case}" in
    shear)
      run_case shear in.prepare_shear_cell_2d in.shear_cell_2d "${PRESSURE_2D}" "${SHEAR_GAP_2D}" "${SERVO_GAIN_2D}"
      ;;
    shear3d)
      run_case shear3d in.prepare_shear_cell_3d in.shear_cell_3d "${PRESSURE_3D}" "${SHEAR_GAP_3D}" "${SERVO_GAIN_3D}"
      ;;
    *)
      echo "Unknown shear case: ${shear_case}" >&2
      exit 1
      ;;
  esac
done

echo "Wrote shear heat results under ${OUT_BASE}"
