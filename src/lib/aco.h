/*

   Ant Colony Optimization search for AutoDock Vina (Vina-ACO).

   A MAX-MIN Ant System in the spirit of PLANTS (Korb, Stuetzle, Exner 2006/2007):
   every degree of freedom (position x/y/z, orientation, torsions) owns a discretized
   pheromone trail; ants sample each DOF independently from its trail, the constructed
   pose is refined with Vina's BFGS local search (Lamarckian), and the best refined
   solutions deposit pheromone. Trails are bounded by [tau_min, tau_max].

   Orientation is parameterized with Shoemake's (u1, u2, u3) in [0,1)^3, which maps a
   uniform distribution onto the Haar-uniform distribution on SO(3), so uniform trails
   do not bias the initial orientations.

   Licensed under the Apache License, Version 2.0.

*/

#ifndef VINA_ACO_H
#define VINA_ACO_H

#include <string>

#include "incrementable.h"
#include "model.h"

struct ant_colony {
	unsigned num_ants;          // ants per iteration
	unsigned max_iterations;    // hard cap on iterations
	unsigned patience;          // iterations without improvement before a restart / stop
	unsigned max_restarts;      // pheromone re-initializations allowed after stagnation
	unsigned max_evals;         // 0 = unlimited
	unsigned local_steps;       // BFGS steps for every ant
	unsigned elite_local_steps; // BFGS steps for the iteration-best ant
	fl evaporation;             // rho
	fl tau_ratio;               // tau_min = tau_max / (tau_ratio * num_bins)
	unsigned num_elite;         // number of ranked ants that deposit pheromone
	fl position_bin;            // Angstrom
	unsigned angle_bins;        // bins for torsions and orientation parameters
	vec hunt_cap;
	fl min_rmsd;
	sz num_saved_mins;
	// variant: 0 = MMAS with discretized per-DOF pheromone trails (PLANTS-like),
	//          1 = ACO_R (Socha & Dorigo 2008): pheromone is a ranked solution archive,
	//              ants sample Gaussian kernels around a guide solution.
	int variant;
	unsigned archive_size;      // ACO_R: k
	fl acor_q;                  // ACO_R: locality of guide selection (smaller = greedier)
	fl acor_xi;                 // ACO_R: kernel width factor (like 1 - evaporation)
	fl acor_explore;            // ACO_R: expected number of DOFs per ant drawn from a broad distribution (tau_min analogue)
	bool ls_halving;            // ACO_R: successive halving of the ants' local search budget
	unsigned consensus;         // ACO_R: stop once this many later restart phases re-find the best pose (0 = off)
	fl rank_lambda;             // final ranking: energy + rank_lambda * -log((support+1)/(phases+1)); 0 = energy only
	fl budget_scale;            // default eval budget = budget_scale * (50 + movable atoms + 10 * DOF); 0 = none
	bool skip_ls;               // ACO_R: QuickVina 2-style first-order consistency check before each ant's local search

	// defaults = the "balanced" preset
	ant_colony() : num_ants(10), max_iterations(100000), patience(15), max_restarts(1000), max_evals(0),
		local_steps(20), elite_local_steps(30), evaporation(0.15), tau_ratio(4), num_elite(3),
		position_bin(0.375), angle_bins(72), hunt_cap(10, 10, 10), min_rmsd(0.5), num_saved_mins(20),
		variant(1), archive_size(10), acor_q(0.1), acor_xi(0.85),
		acor_explore(1.0), ls_halving(false), consensus(2), rank_lambda(0), budget_scale(600), skip_ls(false) {}

	// named parameter sets: "fast", "balanced" (default), "accurate"; returns false for unknown names
	static bool preset(const std::string& name, ant_colony& out);

	// returns an empty string if the parameters are usable, else a description of the problem
	std::string validate() const;

	// out is sorted. If phases is given, the best pose of every restart phase (ACO_R) is appended to it
	// (one representative per phase, used for basin-population ranking).
	void operator()(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
	                const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,
	                output_container* phases = NULL) const;
private:
	void run_mmas(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
	              const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	void run_acor(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
	              const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,
	              output_container* phases) const;
};

struct parallel_aco {
	ant_colony aco;
	sz num_tasks;   // independent colonies (exhaustiveness)
	sz num_threads;
	bool display_progress;
	parallel_aco() : num_tasks(1), num_threads(1), display_progress(true) {}
	void operator()(const model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
	                const vec& corner1, const vec& corner2, rng& generator, std::function<void(double)>* progress_callback,
	                output_container* phases = NULL) const;
};

#endif
