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
#include "pair_granular.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDissipationrateAtom::ComputeDissipationrateAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), group2(nullptr), dissipationrate(nullptr)
{
  if (narg < 0) error->all(FLERR, "Illegal compute dissipationrate/atom command");

  jgroup = group->find("all");
  jgroupbit = group->bitmask[jgroup];

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_reverse = 1;

  nmax = 0;

  // error checks

  if (!atom->radius_flag) error->all(FLERR, "Compute dissipationrate/atom requires atom attribute radius");
}

/* ---------------------------------------------------------------------- */

ComputeDissipationrateAtom::~ComputeDissipationrateAtom()
{
  memory->destroy(dissipationrate);
  delete[] group2;
}

/* ---------------------------------------------------------------------- */

void ComputeDissipationrateAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute dissipationrate/atom requires a pair style be defined");

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
  int i, j, ii, jj, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;

  invoked_peratom = update->ntimestep;

  // grow dissipationrate array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(dissipationrate);
    nmax = atom->nmax;
    memory->create(dissipationrate, nmax, "dissipationrate/atom:dissipationrate");
    vector_atom = dissipationrate;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute number of dissipationrates for each atom in group
  // dissipationrate if distance <= sum of radii
  // tally for both I and J

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  bool update_i_flag, update_j_flag;

  for (i = 0; i < nall; i++) dissipationrate[i] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // Only proceed if i is either part of the compute group or will contribute to dissipationrates
    if (!(mask[i] & groupbit) && !(mask[i] & jgroupbit)) continue;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // Only tally for atoms in compute group (groupbit) if neighbor is in group2 (jgroupbit)
      update_i_flag = (mask[i] & groupbit) && (mask[j] & jgroupbit);
      update_j_flag = (mask[j] & groupbit) && (mask[i] & jgroupbit);
      if (!update_i_flag && !update_j_flag) continue;

      //if (update_i_flag) dissipationrate[i] += force->pair;
      //if (update_j_flag) dissipationrate[j] += 1.0;
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
  for (i = first; i < last; i++) buf[m++] = dissipationrate[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeDissipationrateAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    dissipationrate[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDissipationrateAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
