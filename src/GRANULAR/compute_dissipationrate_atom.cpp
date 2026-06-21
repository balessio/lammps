/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_dissipationrate_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDissipationrateAtom::ComputeDissipationrateAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), dissipationrate(nullptr)
{
  if (narg < 0) error->all(FLERR, "Illegal compute dissipationrate/atom command");

  peratom_flag = 1;
  size_peratom_cols = 3;
  comm_reverse = 3;

  nmax = 0;

  // error checks

  if (!atom->radius_flag) error->all(FLERR, "Compute dissipationrate/atom requires atom attribute radius");
}

/* ---------------------------------------------------------------------- */

ComputeDissipationrateAtom::~ComputeDissipationrateAtom()
{
  memory->destroy(dissipationrate);
}

/* ---------------------------------------------------------------------- */

void ComputeDissipationrateAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute dissipationrate/atom requires a pair style be defined");

  if (force->pair->single_extra < 15)
    error->all(FLERR, "Compute dissipationrate/atom requires pair style granular with dissipative_heat");

  if (modify->get_compute_by_style("dissipationrate/atom").size() > 1 && comm->me == 0)
    error->warning(FLERR, "More than one compute dissipationrate/atom");

  // need an occasional neighbor list

  neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeDissipationrateAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeDissipationrateAtom::compute_peratom()
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  double radi, fpair;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *type = atom->type;

  invoked_peratom = update->ntimestep;

  // grow dissipationrate array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(dissipationrate);
    nmax = atom->nmax;
    memory->create(dissipationrate, nmax, 3, "dissipationrate/atom:dissipationrate");
    array_atom = dissipationrate;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  Pair *pair = force->pair;

  // compute number of dissipationrates for each atom in group
  // dissipationrate if distance <= sum of radii
  // tally for both I and J

  double **x = atom->x;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  bool update_i_flag, update_j_flag;

  for (i = 0; i < nall; i++)
    for (j = 0; j < 3; j++) dissipationrate[i][j] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    itype = type[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      // Only tally for atoms in compute group (groupbit)
      update_i_flag = (mask[i] & groupbit);
      update_j_flag = (mask[j] & groupbit);
      if (!update_i_flag && !update_j_flag) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      // Do not pre-filter by radsum here. Granular history contacts can expose
      // deleted-contact heat through pair->single() just after separation.
      pair->single(i, j, itype, jtype, rsq, 1.0, 1.0, fpair);

      if (update_i_flag) {
        dissipationrate[i][0] += 0.5 * force->pair->svector[12] / update->dt;
        dissipationrate[i][1] += 0.5 * force->pair->svector[13] / update->dt;
        dissipationrate[i][2] += 0.5 * force->pair->svector[14] / update->dt;
      }

      if (update_j_flag) {
        dissipationrate[j][0] += 0.5 * force->pair->svector[12] / update->dt;
        dissipationrate[j][1] += 0.5 * force->pair->svector[13] / update->dt;
        dissipationrate[j][2] += 0.5 * force->pair->svector[14] / update->dt;
      }
    }
  }

  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm(this);

}

/* ---------------------------------------------------------------------- */

int ComputeDissipationrateAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dissipationrate[i][0];
    buf[m++] = dissipationrate[i][1];
    buf[m++] = dissipationrate[i][2];
  }
  return m;

}

/* ---------------------------------------------------------------------- */

void ComputeDissipationrateAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    dissipationrate[j][0] += buf[m++];
    dissipationrate[j][1] += buf[m++];
    dissipationrate[j][2] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDissipationrateAtom::memory_usage()
{
  double bytes = (double) nmax * 3 * sizeof(double);
  return bytes;
}
