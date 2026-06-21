/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "gran_sub_mod_tangential.h"

#include "error.h"
#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_normal.h"
#include "granular_model.h"
#include "math_extra.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathExtra;

static constexpr double EPSILON = 1e-10;

/* ----------------------------------------------------------------------
   Default model
------------------------------------------------------------------------- */

GranSubModTangential::GranSubModTangential(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  allow_synchronization = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangential::init()
{
  if (gm->dissipative_heat)
    nsvector = 2;
}

/* ---------------------------------------------------------------------- */

double GranSubModTangential::calculate_heat()
{
  if (gm->heat_tang_damp != 0 || gm->heat_tang_fric != 0)
    error->one(FLERR, "Granular tangential model {} does not calculate a dissipative heat term", name);

  if (gm->calculate_svector) {
    gm->svector[index_svector] = 0.0;
    gm->svector[index_svector + 1] = 0.0;
  }

  return 0.0;
}

/* ---------------------------------------------------------------------- */

double GranSubModTangential::elastic_potential()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModTangentialNone::GranSubModTangentialNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  allow_synchronization = 1;
}

/* ----------------------------------------------------------------------
   Linear model with no history
------------------------------------------------------------------------- */

GranSubModTangentialLinearNoHistory::GranSubModTangentialLinearNoHistory(GranularModel *gm,
                                                                         LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 2;
  size_history = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::coeffs_to_local()
{
  k = 0.0;    // No tangential stiffness with no history
  xt = coeffs[0];
  mu = coeffs[1];

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal linear no history tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::calculate_forces()
{
  // classic pair gran/hooke (no history)
  damp = xt * gm->damping_model->get_damp_prefactor();

  double Fscrit = mu * gm->normal_model->get_fncrit();
  double fsmag = damp * gm->vrel;

  if (gm->vrel != 0.0)
    gm->Ft = MIN(Fscrit, fsmag) / gm->vrel;
  else
    gm->Ft = 0.0;

  scale3(-gm->Ft, gm->vtr, gm->fs);
  gm->StrainEnergyTang = 0.0;
}

/* ---------------------------------------------------------------------- */

double GranSubModTangentialLinearNoHistory::calculate_heat()
{
  double dq_damp, dq_friction;
  double Fscrit = mu * gm->normal_model->get_fncrit();

  if (Fscrit <= gm->Ft) {
    dq_damp = 0.0;
    dq_friction = gm->Ft * gm->vrel * gm->dt;
  }
  else {
    dq_damp = gm->Ft * gm->vrel * gm->dt;
    dq_friction = 0.0;
  }

  dq_damp *= gm->heat_tang_damp;
  dq_friction *= gm->heat_tang_fric;

  if (gm->calculate_svector) {
    gm->svector[index_svector] = dq_damp;
    gm->svector[index_svector + 1] = dq_friction;
  }

  return dq_damp + dq_friction;
}

/* ----------------------------------------------------------------------
   Linear model with history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistory::GranSubModTangentialLinearHistory(GranularModel *gm,
                                                                     LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  allow_synchronization = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || xt < 0.0 || mu < 0.0) error->all(FLERR, "Illegal linear tangential model");
}

/* ---------------------------------------------------------------------- */

/*void GranSubModTangentialLinearHistory::calculate_forces()
{
  // Note: this is the same as the base Mindlin calculation except k isn't scaled by contact radius
  double magfs, magfs_inv, rsht, shrmag, temp_array[3], vtr2[3];
  int frame_update = 0;

  double Ftangdamp, Ftangelasvec[3], Ftangelas_prevvec[3], StrEngPrev;

  damp = xt * gm->damping_model->get_damp_prefactor();

  double Fscrit = gm->normal_model->get_fncrit() * mu;
  double *history = &gm->history[history_index];
  
  scale3(-k, history, Ftangelas_prevvec);

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);

    if (frame_update) rotate_rescale_vec(history, gm->nx);
    scale3(-k, history, Ftangelas_prevvec);

    // update history, tangential force using velocities at half step
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    scale3(gm->dt, gm->vtr, temp_array);
    add3(history, temp_array, history);

    if(gm->synchronized_verlet == 1) {
      rsht = dot3(history, gm->nx_unrotated);
      frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);
      //Second projection to nx (t+\Delta t)
      if (frame_update) rotate_rescale_vec(history, gm->nx_unrotated);
    }
  }

  // tangential forces = history + tangential velocity damping
  scale3(-k, history, gm->fs);
  scale3(-k, history, Ftangelasvec);

  //Rotating vtr for damping term in nx direction
  if (frame_update && gm->synchronized_verlet == 1) {
    copy3(gm->vtr, vtr2);
    rotate_rescale_vec(vtr2, gm->nx_unrotated);
  } else {
    copy3(gm->vtr, vtr2);
  }
  scale3(damp, vtr2, temp_array);
  Ftangdamp = len3(temp_array);
  sub3(gm->fs, temp_array, gm->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  bool fric_on = magfs > Fscrit;

  if (fric_on) {
    shrmag = len3(history);

    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);
      scale3(-1.0 / k, history);
      scale3(Fscrit * magfs_inv, gm->fs);
      scale3(Fscrit * magfs_inv, Ftangelasvec);
      Ftangdamp = Fscrit * magfs_inv * Ftangdamp;
    } else {
      zero3(gm->fs);
      zero3(Ftangelasvec);
      Ftangdamp = 0.0;
    }
  }
  else gm->dq_friction_hold = 0.0;

  double slip[3];
  sub3(Ftangelasvec, Ftangelas_prevvec, temp_array);
  scale3(1.0 / k, temp_array);
  scale3(gm->dt, gm->vtr, slip);
  sub3(slip, temp_array, slip);
  add3(Ftangelas_prevvec, Ftangelasvec, temp_array);
  scale3(0.5, temp_array);
  gm->dq_friction_hold = -dot3(slip, temp_array);
  if (!fric_on) gm->dq_friction_hold = 0.0;

  gm->dq_damp_hold = Ftangdamp * gm->vrel * gm->dt;
  gm->StrainEnergyTang = 0.5 * dot3(Ftangelasvec, Ftangelasvec) / k;

}*/

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::calculate_forces()
{
  // Note: this is the same as the base Mindlin calculation except k isn't scaled by contact radius
  double magfs, magfs_inv, rsht, shrmag, temp_array[3], vtr2[3];
  int frame_update = 0;

  double F_tang_elas[3], F_tang_elas_prev[3], history_prev[3], slip[3];

  damp = xt * gm->damping_model->get_damp_prefactor();

  double Fscrit = gm->normal_model->get_fncrit() * mu;
  double *history = &gm->history[history_index];

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);

    if (frame_update) rotate_rescale_vec(history, gm->nx);

    copy3(history, history_prev);

    // update history, tangential force using velocities at half step
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    scale3(gm->dt, gm->vtr, temp_array);
    add3(history, temp_array, history);

    if(gm->synchronized_verlet == 1) {
      rsht = dot3(history, gm->nx_unrotated);
      frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);
      //Second projection to nx (t+\Delta t)
      if (frame_update) rotate_rescale_vec(history, gm->nx_unrotated);
    }
  } else {
    copy3(history, history_prev);
  }
  scale3(-k, history_prev, F_tang_elas_prev);

  // tangential forces = history + tangential velocity damping
  scale3(-k, history, gm->fs);

  //Rotating vtr for damping term in nx direction
  if (frame_update && gm->synchronized_verlet == 1) {
    copy3(gm->vtr, vtr2);
    rotate_rescale_vec(vtr2, gm->nx_unrotated);
  } else {
    copy3(gm->vtr, vtr2);
  }
  scale3(damp,vtr2, temp_array);

  sub3(gm->fs, temp_array, gm->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  bool fric_on = magfs > Fscrit;
  gm->dq_friction_hold = 0.0;
  double damp_scale = 1.0;
  if (fric_on) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      damp_scale = Fscrit * magfs_inv;
      scale3(damp_scale, history);
      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
      damp_scale = 0.0;
    }
  }

  scale3(-k, history, F_tang_elas);

  double tangential_step[3];
  scale3(gm->dt, vtr2, tangential_step);

  double dq_friction = 0.0;
  add3(F_tang_elas, F_tang_elas_prev, temp_array);
  scale3(0.5, temp_array, temp_array);

  copy3(tangential_step, slip);
  sub3(slip, history, slip);
  add3(slip, history_prev, slip);

  if (fric_on) dq_friction = MAX(0.0, -dot3(slip, temp_array));
  gm->dq_friction_hold = dq_friction;
  gm->dq_damp_hold = damp_scale * damp * dot3(vtr2, vtr2) * gm->dt;

  gm->StrainEnergyTang = 0.5*dot3(F_tang_elas, F_tang_elas)/k;

}

/* ---------------------------------------------------------------------- */

double GranSubModTangentialLinearHistory::calculate_heat()
{
  
  double dq_friction = gm->dq_friction_hold; 
  double dq_damp = gm->dq_damp_hold;

  dq_damp *= gm->heat_tang_damp;
  dq_friction *= gm->heat_tang_fric;

  if (gm->calculate_svector) {
    gm->svector[index_svector] = dq_damp;
    gm->svector[index_svector + 1] = dq_friction;
  }

  return dq_damp + dq_friction;
}

/* ---------------------------------------------------------------------- */

double GranSubModTangentialLinearHistory::elastic_potential()
{
  double *history = &gm->history[history_index];
  return 0.5 * k * dot3(history, history);
}

/* ----------------------------------------------------------------------
   Linear model with history from pair gran/hooke/history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistoryClassic::GranSubModTangentialLinearHistoryClassic(
    GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangentialLinearHistory(gm, lmp)
{
  contact_radius_flag = 1;    // Sets gran/hooke/history behavior
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistoryClassic::calculate_forces()
{
  double magfs, magfs_inv, rsht, shrmag;
  double temp_array[3];

  double Ftangdamp, Ftangelasvec[3], Ftangelas_prevvec[3];

  damp = xt * gm->damping_model->get_damp_prefactor();

  double Fscrit = gm->normal_model->get_fncrit() * mu;
  double *history = &gm->history[history_index];

  scale3(-k, history, Ftangelas_prevvec);

  // update history
  if (gm->history_update) {
    scale3(gm->dt, gm->vtr, temp_array);
    add3(history, temp_array, history);
  }

  shrmag = len3(history);

  // rotate shear displacements
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    scale3(rsht, gm->nx, temp_array);
    sub3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  if (contact_radius_flag) {
    scale3(-k * gm->contact_radius, history, gm->fs);
    scale3(-k * gm->contact_radius, history, Ftangelasvec);
  }
  else {
    scale3(-k, history, gm->fs);
    scale3(-k, history, Ftangelasvec);
  }
  scale3(damp, gm->vtr, temp_array);
  Ftangdamp = len3(temp_array);
  sub3(gm->fs, temp_array, gm->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);
      scale3(-1.0 / k, history);
      scale3(Fscrit * magfs_inv, gm->fs);
      scale3(Fscrit * magfs_inv, Ftangelasvec);
    } else {
      zero3(gm->fs);
      zero3(Ftangelasvec);
      Ftangdamp = 0.0;
    }
  }
  Ftangdamp = len3(gm->vtr) * damp;
  gm->dq_damp_hold = damp * dot3(gm->vtr, gm->vtr) * gm->dt; 

  gm->dq_friction_hold = 0.0; // to implement

  gm->StrainEnergyTang = 0.5 * dot3(Ftangelasvec, Ftangelasvec) / k;

}

/* ---------------------------------------------------------------------- */

double GranSubModTangentialLinearHistoryClassic::calculate_heat()
{
  double dq_damp = gm->dq_damp_hold;
  double dq_friction = gm->dq_friction_hold;

  dq_damp *= gm->heat_tang_damp;
  dq_friction *= gm->heat_tang_fric;

  if (gm->calculate_svector) {
    gm->svector[index_svector] = dq_damp;
    gm->svector[index_svector + 1] = dq_friction;
  }

  return dq_damp + dq_friction;
}

/* ----------------------------------------------------------------------
   Mindlin from pair gran/hertz/history
------------------------------------------------------------------------- */

GranSubModTangentialMindlinClassic::GranSubModTangentialMindlinClassic(GranularModel *gm,
                                                                       LAMMPS *lmp) :
    GranSubModTangentialLinearHistoryClassic(gm, lmp)
{
  contact_radius_flag = 1;    // Sets gran/hertz/history behavior
}

/* ----------------------------------------------------------------------
   Mindlin model
------------------------------------------------------------------------- */

GranSubModTangentialMindlin::GranSubModTangentialMindlin(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  mindlin_force = 0;
  mindlin_rescale = 0;
  contact_radius_flag = 1;
  allow_synchronization = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k == -1) {
    if (!gm->normal_model->get_material_properties())
      error->all(FLERR,
                 "Must either specify tangential stiffness or material properties for normal model "
                 "for the Mindlin tangential style");

    double Emod = gm->normal_model->get_emod();
    double poiss = gm->normal_model->get_poiss();

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu < 0.0) error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1)
    coeffs[0] = -1;
  else
    coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag;
  double temp_array[3], vtr2[3];
  int frame_update = 0;

  double Ftangdamp, Ftangelasvec[3], Ftangelas_prevvec[3];

  damp = xt * gm->damping_model->get_damp_prefactor();

  double *history = &gm->history[history_index];
  double Fscrit = gm->normal_model->get_fncrit() * mu;

  k_scaled = k * gm->contact_radius;

  zero3(Ftangelas_prevvec);
  if (!mindlin_force)
    scale3(-k_scaled, history, Ftangelas_prevvec);
  else
    add3(Ftangelas_prevvec, history, Ftangelas_prevvec);

  // on unloading, rescale the shear displacements/force
  if (mindlin_rescale)
    if (gm->contact_radius < history[3]) scale3(gm->contact_radius / history[3], history);

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) rotate_rescale_vec(history, gm->nx);

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * gm->dt, gm->vtr, temp_array);
    } else {
      scale3(gm->dt, gm->vtr, temp_array);
    }
    add3(history, temp_array, history);

    if (mindlin_rescale) history[3] = gm->contact_radius;

    if (gm->synchronized_verlet == 1) {
      // second projection to full step normal
      rsht = dot3(history, gm->nx_unrotated);
      if (mindlin_force) {
        frame_update = fabs(rsht) > (EPSILON * Fscrit);
      } else {
        frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
      }
      if (frame_update) rotate_rescale_vec(history, gm->nx_unrotated);
    }
  }

  // tangential forces = history + tangential velocity damping
  // Rotating vtr for damping term in nx direction
  if (frame_update && gm->synchronized_verlet) {
    copy3(gm->vtr, vtr2);
    rotate_rescale_vec(vtr2, gm->nx_unrotated);
  } else {
    copy3(gm->vtr, vtr2);
  }
  scale3(-damp, vtr2, gm->fs);
  Ftangdamp = len3(gm->fs);

  zero3(Ftangelasvec);
  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(gm->fs, temp_array, gm->fs);
    scale3(-k_scaled, history, Ftangelasvec);
  } else {
    add3(gm->fs, history, gm->fs);
    add3(Ftangelasvec, history, Ftangelasvec);
  }

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      const double clip_scale = Fscrit * magfs_inv;
      scale3(clip_scale, history);
      scale3(Fscrit * magfs_inv, gm->fs);
      scale3(Fscrit * magfs_inv, Ftangelasvec);
      Ftangdamp *= clip_scale;

    } else {
      zero3(gm->fs);

      zero3(Ftangelasvec);
      Ftangdamp = 0.0;
    
    }
  }

  gm->dq_damp_hold = Ftangdamp * len3(gm->vtr) * gm->dt;

  gm->dq_friction_hold = 0.0; // to implement

  gm->StrainEnergyTang = 0.5 * dot3(Ftangelasvec, Ftangelasvec) / k_scaled;
}

/* ---------------------------------------------------------------------- */

double GranSubModTangentialMindlin::calculate_heat()
{
  double dq_damp = gm->dq_damp_hold;
  double dq_friction = gm->dq_friction_hold;

  dq_damp *= gm->heat_tang_damp;
  dq_friction *= gm->heat_tang_fric;

  if (gm->calculate_svector) {
    gm->svector[index_svector] = dq_damp;
    gm->svector[index_svector + 1] = dq_friction;
  }

  return dq_damp + dq_friction;
}

/* ----------------------------------------------------------------------
   Mindlin force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinForce::GranSubModTangentialMindlinForce(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  mindlin_force = 1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescale::GranSubModTangentialMindlinRescale(GranularModel *gm,
                                                                       LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescaleForce::GranSubModTangentialMindlinRescaleForce(GranularModel *gm,
                                                                                 LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_force = 1;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}
