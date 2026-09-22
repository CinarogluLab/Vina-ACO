/*

   Ant Colony Optimization search for AutoDock Vina (Vina-ACO).
   See aco.h for an overview.

   Licensed under the Apache License, Version 2.0.

*/

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "aco.h"
#include "coords.h"
#include "quasi_newton.h"
#include "parallel.h"
#include "parallel_progress.h"

namespace {

const fl two_pi = 2 * pi;

// One pheromone trail over a discretized interval [lo, hi).
struct trail {
	fl lo, width;
	bool periodic;
	flv tau;
	flv cumulative;

	trail(fl lo_, fl hi_, sz bins, bool periodic_) : lo(lo_), periodic(periodic_), tau(bins, 0), cumulative(bins, 0) {
		width = (hi_ - lo_) / bins;
	}
	sz size() const { return tau.size(); }

	sz bin_of(fl v) const {
		const int n = int(tau.size());
		int b = int(std::floor((v - lo) / width));
		if(periodic) {
			b %= n;
			if(b < 0) b += n;
		}
		else
			b = std::max(0, std::min(n - 1, b));
		return sz(b);
	}
	fl sample(rng& generator) {
		fl acc = 0;
		VINA_FOR_IN(i, tau) {
			acc += tau[i];
			cumulative[i] = acc;
		}
		const fl r = random_fl(0, acc, generator);
		const sz b = std::min(sz(std::lower_bound(cumulative.begin(), cumulative.end(), r) - cumulative.begin()), tau.size() - 1);
		return lo + (b + random_fl(0, 1, generator)) * width;
	}
	// deposit on the bin of v, with half the amount spilled on each neighbour (smoothing)
	void deposit(fl v, fl amount) {
		const int n = int(tau.size());
		const int b = int(bin_of(v));
		tau[b] += amount;
		for(int d = -1; d <= 1; d += 2) {
			int nb = b + d;
			if(periodic) nb = (nb + n) % n;
			else if(nb < 0 || nb >= n) continue;
			tau[nb] += 0.5 * amount;
		}
	}
	void evaporate(fl rho) {
		VINA_FOR_IN(i, tau)
			tau[i] *= (1 - rho);
	}
	void clamp(fl tmin, fl tmax) {
		VINA_FOR_IN(i, tau)
			tau[i] = std::max(tmin, std::min(tmax, tau[i]));
	}
	void reset(fl value) {
		std::fill(tau.begin(), tau.end(), value);
	}
};

// Shoemake: uniform (u1, u2, u3) in [0,1)^3 -> Haar-uniform unit quaternion
qt shoemake_to_quaternion(fl u1, fl u2, fl u3) {
	u1 = std::max(fl(0), std::min(fl(1), u1));
	const fl s1 = std::sqrt(1 - u1), s2 = std::sqrt(u1);
	const fl t1 = two_pi * u2, t2 = two_pi * u3;
	return qt(s2 * std::cos(t2), s1 * std::sin(t1), s1 * std::cos(t1), s2 * std::sin(t2));
}

void quaternion_to_shoemake(const qt& q_in, fl& u1, fl& u2, fl& u3) {
	qt q = q_in;
	quaternion_normalize_approx(q, 1e-3);
	const fl w = q.R_component_1(), x = q.R_component_2(), y = q.R_component_3(), z = q.R_component_4();
	u1 = z * z + w * w;
	u2 = std::atan2(x, y) / two_pi; if(u2 < 0) u2 += 1;
	u3 = std::atan2(z, w) / two_pi; if(u3 < 0) u3 += 1;
}

struct dof_codec {
	std::vector<trail> trails;

	dof_codec(const conf_size& s, const vec& corner1, const vec& corner2, fl position_bin, unsigned angle_bins) {
		VINA_FOR_IN(i, s.ligands) {
			VINA_FOR(d, 3) {
				const fl len = corner2[d] - corner1[d];
				const sz bins = std::max(sz(4), std::min(sz(1024), sz(std::ceil(len / position_bin))));
				trails.push_back(trail(corner1[d], corner2[d], bins, false));
			}
			trails.push_back(trail(0, 1, angle_bins, false)); // u1
			trails.push_back(trail(0, 1, angle_bins, true));  // u2
			trails.push_back(trail(0, 1, angle_bins, true));  // u3
			VINA_FOR(t, s.ligands[i])
				trails.push_back(trail(-pi, pi, angle_bins, true));
		}
		VINA_FOR_IN(i, s.flex)
			VINA_FOR(t, s.flex[i])
				trails.push_back(trail(-pi, pi, angle_bins, true));
	}

	void construct(conf& c, rng& generator) {
		sz k = 0;
		VINA_FOR_IN(i, c.ligands) {
			ligand_conf& lig = c.ligands[i];
			VINA_FOR(d, 3)
				lig.rigid.position[d] = trails[k++].sample(generator);
			const fl u1 = trails[k++].sample(generator);
			const fl u2 = trails[k++].sample(generator);
			const fl u3 = trails[k++].sample(generator);
			lig.rigid.orientation = shoemake_to_quaternion(u1, u2, u3);
			VINA_FOR_IN(t, lig.torsions)
				lig.torsions[t] = trails[k++].sample(generator);
		}
		VINA_FOR_IN(i, c.flex)
			VINA_FOR_IN(t, c.flex[i].torsions)
				c.flex[i].torsions[t] = trails[k++].sample(generator);
	}

	void deposit(const conf& c, fl amount) {
		sz k = 0;
		VINA_FOR_IN(i, c.ligands) {
			const ligand_conf& lig = c.ligands[i];
			VINA_FOR(d, 3)
				trails[k++].deposit(lig.rigid.position[d], amount);
			fl u1, u2, u3;
			quaternion_to_shoemake(lig.rigid.orientation, u1, u2, u3);
			trails[k++].deposit(u1, amount);
			trails[k++].deposit(u2, amount);
			trails[k++].deposit(u3, amount);
			VINA_FOR_IN(t, lig.torsions)
				trails[k++].deposit(normalized_angle(lig.torsions[t]), amount);
		}
		VINA_FOR_IN(i, c.flex)
			VINA_FOR_IN(t, c.flex[i].torsions)
				trails[k++].deposit(normalized_angle(c.flex[i].torsions[t]), amount);
	}

	void evaporate_deposit_clamp(fl rho, const std::vector<const conf*>& depositors, const flv& amounts, fl tau_max, fl tau_ratio) {
		VINA_FOR_IN(i, trails)
			trails[i].evaporate(rho);
		VINA_FOR_IN(i, depositors)
			deposit(*depositors[i], amounts[i]);
		VINA_FOR_IN(i, trails)
			trails[i].clamp(tau_max / (tau_ratio * trails[i].size()), tau_max);
	}

	void reset(fl value) {
		VINA_FOR_IN(i, trails)
			trails[i].reset(value);
	}
};

struct by_energy {
	bool operator()(const output_type* a, const output_type* b) const { return a->e < b->e; }
};

} // namespace

bool ant_colony::preset(const std::string& name, ant_colony& out) {
	out = ant_colony(); // balanced
	if(name == "balanced")
		return true;
	if(name == "fast") {
		out.budget_scale = 300;
		out.consensus = 1;
		return true;
	}
	if(name == "accurate") {
		out.budget_scale = 3000;
		out.consensus = 2;
		return true;
	}
	return false;
}

std::string ant_colony::validate() const {
	if(num_ants < 1) return "number of ants must be >= 1";
	if(max_iterations < 1) return "number of iterations must be >= 1";
	if(patience < 1) return "patience must be >= 1";
	if(num_elite < 1) return "number of elite ants must be >= 1";
	if(!(std::isfinite(evaporation) && evaporation > 0 && evaporation < 1)) return "evaporation must be in (0, 1)";
	if(!(std::isfinite(tau_ratio) && tau_ratio > 0)) return "tau_ratio must be finite and > 0";
	if(!(std::isfinite(position_bin) && position_bin > 0)) return "position bin width must be finite and > 0";
	if(angle_bins < 2 || angle_bins > 4096) return "angle bins must be in [2, 4096]";
	if(!(std::isfinite(min_rmsd) && min_rmsd >= 0)) return "min_rmsd must be finite and >= 0";
	if(num_saved_mins < 1) return "number of saved poses must be >= 1";
	if(variant != 0 && variant != 1) return "variant must be 0 (mmas) or 1 (acor)";
	if(archive_size < 2) return "archive size must be >= 2";
	if(!(std::isfinite(budget_scale) && budget_scale >= 0)) return "budget scale must be finite and >= 0";
	if(!(std::isfinite(acor_explore) && acor_explore >= 0)) return "explore must be finite and >= 0";
	if(!(std::isfinite(rank_lambda) && rank_lambda >= 0)) return "rank lambda must be finite and >= 0";
	if(!(std::isfinite(acor_q) && acor_q > 0)) return "q must be finite and > 0";
	if(!(std::isfinite(acor_xi) && acor_xi > 0)) return "xi must be finite and > 0";
	return "";
}

void ant_colony::operator()(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                            const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,
                            output_container* phases) const {
	const std::string problem = validate();
	if(!problem.empty())
		throw std::invalid_argument("ACO: " + problem);
	if(variant == 0)
		run_mmas(m, out, p, ig, corner1, corner2, increment_me, generator);
	else
		run_acor(m, out, p, ig, corner1, corner2, increment_me, generator, phases);
}

void ant_colony::run_mmas(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                          const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
	const vec authentic_v(1000, 1000, 1000);
	const conf_size s = m.get_size();
	change g(s);
	dof_codec codec(s, corner1, corner2, position_bin, angle_bins);

	quasi_newton ant_ls;   ant_ls.max_steps = local_steps;
	quasi_newton elite_ls; elite_ls.max_steps = elite_local_steps;

	const sz n_elite = std::max(sz(1), std::min(sz(num_elite), sz(num_ants)));
	flv amounts; // global best gets 1, iteration ranks get 1/(r+1)
	fl deposit_total = 1;
	VINA_FOR(r, n_elite) deposit_total += 1.0 / (r + 1);
	const fl tau_max = deposit_total / evaporation;
	codec.reset(tau_max);

	boost::ptr_vector<output_type> ants;
	VINA_FOR(a, num_ants)
		ants.push_back(new output_type(conf(s), max_fl));
	std::vector<output_type*> ranked(num_ants);

	output_type global_best(conf(s), max_fl);
	bool have_global_best = false;
	unsigned since_improvement = 0;
	unsigned restarts = 0;
	int evalcount = 0;
	unsigned iterations_done = 0, best_iteration = 0;
	static const bool debug = std::getenv("VINA_ACO_DEBUG") != NULL;

	VINA_U_FOR(iteration, max_iterations) {
		iterations_done = iteration + 1;
		if(increment_me)
			++(*increment_me);
		if(max_evals > 0 && evalcount >= int(max_evals) && !out.empty())
			break;

		// construction + short Lamarckian local search
		sz built = 0;
		VINA_FOR(a, num_ants) {
			if(max_evals > 0 && evalcount >= int(max_evals) && built > 0)
				break;
			output_type& ant = ants[a];
			codec.construct(ant.c, generator);
			ant_ls(m, p, ig, ant, g, hunt_cap, evalcount);
			ranked[built++] = &ant;
		}
		std::sort(ranked.begin(), ranked.begin() + built, by_energy());

		// deeper refinement of the iteration elite, with authentic (uncapped) energies
		bool improved = false;
		const sz n_refine = std::min(n_elite, built);
		VINA_FOR(r, n_refine) {
			output_type& e = *ranked[r];
			elite_ls(m, p, ig, e, g, authentic_v, evalcount);
			e.coords = m.get_heavy_atom_movable_coords();
			add_to_output_container(out, e, min_rmsd, num_saved_mins);
			if(e.e < global_best.e - 1e-3) {
				if(!out.empty() && e.e <= out.front().e) best_iteration = iteration;
				global_best = e;
				have_global_best = true;
				improved = true;
			}
		}
		std::sort(ranked.begin(), ranked.begin() + n_refine, by_energy());

		// pheromone update (MMAS)
		std::vector<const conf*> depositors;
		amounts.clear();
		if(have_global_best) {
			depositors.push_back(&global_best.c);
			amounts.push_back(1);
		}
		VINA_FOR(r, n_refine) {
			depositors.push_back(&ranked[r]->c);
			amounts.push_back(1.0 / (r + 1));
		}
		codec.evaporate_deposit_clamp(evaporation, depositors, amounts, tau_max, tau_ratio);

		// stagnation handling
		if(improved)
			since_improvement = 0;
		else if(++since_improvement >= patience) {
			if(restarts >= max_restarts)
				break;
			++restarts;
			since_improvement = 0;
			codec.reset(tau_max);
			global_best.e = max_fl; // the output container keeps the old optimum
			have_global_best = false;
		}
	}
	if(debug)
		std::cerr << "ACO colony: iterations " << iterations_done << " evals " << evalcount << " restarts " << restarts
		          << " best " << out.front().e << " found_at_iter " << best_iteration << "\n";
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e);
}

namespace {

fl angle_distance(fl a, fl b) { return std::abs(normalized_angle(a - b)); }

fl orientation_distance(const qt& a, const qt& b) {
	qt qa = a, qb = b;
	quaternion_normalize_approx(qa, 1e-3);
	quaternion_normalize_approx(qb, 1e-3);
	fl dot = std::abs(qa.R_component_1() * qb.R_component_1() + qa.R_component_2() * qb.R_component_2() +
	                  qa.R_component_3() * qb.R_component_3() + qa.R_component_4() * qb.R_component_4());
	dot = std::min(fl(1), dot);
	return 2 * std::acos(dot); // rotation angle in [0, pi]
}

// Gaussian-kernel sampling around the guide solution archive[l] (ACO_R).
// With probability p_explore a DOF ignores the archive and is drawn from a broad distribution
// (uniform torsion, uniform random rotation, position within 2 A of the guide): the continuous
// counterpart of the MMAS lower pheromone bound tau_min, which keeps every value reachable.
void acor_sample(conf& c, const std::vector<const output_type*>& archive, sz l, fl xi, fl p_explore,
                 const vec& corner1, const vec& corner2, rng& generator) {
	const conf& guide = archive[l]->c;
	const sz k = archive.size();
	const fl denom = fl(std::max(sz(1), k - 1));
	c = guide;
	VINA_FOR_IN(i, c.ligands) {
		ligand_conf& lig = c.ligands[i];
		const ligand_conf& gl = guide.ligands[i];
		const bool explore_position = random_fl(0, 1, generator) < p_explore;
		const vec jump = explore_position ? 2.0 * random_inside_sphere(generator) : zero_vec;
		VINA_FOR(d, 3) {
			fl sum = 0;
			VINA_FOR(e, k) sum += std::abs(archive[e]->c.ligands[i].rigid.position[d] - gl.rigid.position[d]);
			const fl sigma = xi * sum / denom;
			fl x = gl.rigid.position[d] + jump[d] + (sigma > 0 ? random_normal(0, sigma, generator) : 0);
			lig.rigid.position[d] = std::max(corner1[d], std::min(corner2[d], x));
		}
		if(random_fl(0, 1, generator) < p_explore)
			lig.rigid.orientation = random_orientation(generator);
		else {
			fl sum = 0;
			VINA_FOR(e, k) sum += orientation_distance(archive[e]->c.ligands[i].rigid.orientation, gl.rigid.orientation);
			const fl sigma = xi * sum / denom / std::sqrt(fl(3)); // per-component sigma of the rotation vector
			if(sigma > 0) {
				vec rot(random_normal(0, sigma, generator), random_normal(0, sigma, generator), random_normal(0, sigma, generator));
				quaternion_increment(lig.rigid.orientation, rot);
			}
		}
		VINA_FOR_IN(t, lig.torsions) {
			if(random_fl(0, 1, generator) < p_explore) {
				lig.torsions[t] = random_fl(-pi, pi, generator);
				continue;
			}
			fl sum = 0;
			VINA_FOR(e, k) sum += angle_distance(archive[e]->c.ligands[i].torsions[t], gl.torsions[t]);
			const fl sigma = xi * sum / denom;
			lig.torsions[t] = normalized_angle(gl.torsions[t] + (sigma > 0 ? random_normal(0, sigma, generator) : 0));
		}
	}
	VINA_FOR_IN(i, c.flex)
		VINA_FOR_IN(t, c.flex[i].torsions) {
			if(random_fl(0, 1, generator) < p_explore) {
				c.flex[i].torsions[t] = random_fl(-pi, pi, generator);
				continue;
			}
			fl sum = 0;
			VINA_FOR(e, k) sum += angle_distance(archive[e]->c.flex[i].torsions[t], guide.flex[i].torsions[t]);
			const fl sigma = xi * sum / denom;
			c.flex[i].torsions[t] = normalized_angle(guide.flex[i].torsions[t] + (sigma > 0 ? random_normal(0, sigma, generator) : 0));
		}
}

// QuickVina 2 start filter (Alhossary et al., Bioinformatics 2015, 31, 2214). The colony remembers the
// last 10n points visited by the ants' local searches (n = number of DOFs), with energy and gradient sign.
// A new ant is minimized only if one of its 2n nearest remembered points q passes, in every DOF i where
// both gradients are nonzero with the same sign, the test sign(g_i) * sign((f(p) - f(q)) (p_i - q_i)) <= 0,
// i.e. the data allow a stationary point between p and q. Otherwise the ant most likely lies on a slope
// the colony has already descended, and its BFGS is skipped after the single starting evaluation.
// Unlike QuickVina, torsion differences are taken modulo 2 pi and there is no limit on the number of DOFs.
struct visited_history : public bfgs_observer<conf, change> {
	sz n, capacity, next;
	std::vector<flv> xs;
	flv fs;
	std::vector<std::vector<signed char> > ds;
	std::vector<bool> periodic;
	bool allow_skip, skipped;
	unsigned n_checked, n_skipped;
	flv x_tmp, dist;
	std::vector<sz> order;
	std::vector<signed char> d_tmp;

	explicit visited_history(const conf_size& s) : n(0), capacity(0), next(0), allow_skip(true), skipped(false), n_checked(0), n_skipped(0) {
		VINA_FOR_IN(i, s.ligands) {
			VINA_FOR(d, 6) periodic.push_back(false);
			VINA_FOR(t, s.ligands[i]) periodic.push_back(true);
			n += 6 + s.ligands[i];
		}
		VINA_FOR_IN(i, s.flex)
			VINA_FOR(t, s.flex[i]) { periodic.push_back(true); ++n; }
		capacity = 10 * n;
	}
	void flatten(const conf& c, const change& g) {
		x_tmp.clear();
		VINA_FOR_IN(i, c.ligands) {
			const ligand_conf& lig = c.ligands[i];
			VINA_FOR(d, 3) x_tmp.push_back(lig.rigid.position[d]);
			const vec r = quaternion_to_angle(lig.rigid.orientation);
			VINA_FOR(d, 3) x_tmp.push_back(r[d]);
			VINA_FOR_IN(t, lig.torsions) x_tmp.push_back(lig.torsions[t]);
		}
		VINA_FOR_IN(i, c.flex)
			VINA_FOR_IN(t, c.flex[i].torsions) x_tmp.push_back(c.flex[i].torsions[t]);
		d_tmp.resize(n);
		VINA_FOR(i, n) {
			const fl gi = g(i);
			d_tmp[i] = gi > 0 ? 1 : (gi < 0 ? -1 : 0);
		}
	}
	fl diff(sz i, fl a, fl b) const { return periodic[i] ? normalized_angle(a - b) : a - b; }
	void add(fl f) {
		if(xs.size() < capacity) {
			xs.push_back(x_tmp); fs.push_back(f); ds.push_back(d_tmp);
		}
		else {
			xs[next] = x_tmp; fs[next] = f; ds[next] = d_tmp;
			next = (next + 1) % capacity;
		}
	}
	bool consistent_with(sz e, fl f) const { // true = a stationary point may lie between the two points
		const fl df = f - fs[e];
		VINA_FOR(i, n) {
			const signed char a = d_tmp[i], b = ds[e][i];
			if(a == 0 || b == 0 || a != b) continue;
			const fl dx = diff(i, x_tmp[i], xs[e][i]);
			// Equation 2 is satisfied when (f(p)-f(q))(p_i-q_i) is zero.
			if(df == 0 || dx == 0) continue;
			const bool f_higher = df > 0;
			const bool x_higher = dx > 0;
			if(a > 0 ? !(x_higher ^ f_higher) : (x_higher ^ f_higher))
				return false;
		}
		return true;
	}
	bool interesting(fl f) {
		if(xs.size() < capacity) return true;
		dist.resize(xs.size());
		order.resize(xs.size());
		VINA_FOR_IN(e, xs) {
			fl d2 = 0;
			VINA_FOR(i, n) d2 += sqr(diff(i, x_tmp[i], xs[e][i]));
			dist[e] = d2;
			order[e] = e;
		}
		const sz m = std::min(order.size(), 2 * n);
		std::partial_sort(order.begin(), order.begin() + m, order.end(), [&](sz a, sz b) { return dist[a] < dist[b]; });
		VINA_FOR(j, m)
			if(consistent_with(order[j], f))
				return true;
		return false;
	}
	bool start(const conf& x, fl f, const change& g) {
		flatten(x, g);
		skipped = false;
		if(allow_skip) {
			++n_checked;
			if(!interesting(f)) {
				skipped = true;
				++n_skipped;
				return false;
			}
		}
		add(f);
		return true;
	}
	void step(const conf& x, fl f, const change& g) {
		flatten(x, g);
		add(f);
	}
};

} // namespace

void ant_colony::run_acor(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                          const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator,
                          output_container* phases) const {
	const vec authentic_v(1000, 1000, 1000);
	const conf_size s = m.get_size();
	change g(s);
	quasi_newton ant_ls;   ant_ls.max_steps = local_steps;
	quasi_newton elite_ls; elite_ls.max_steps = elite_local_steps;
	const sz k = std::max(sz(2), sz(archive_size));
	const sz n_elite = std::max(sz(1), std::min(sz(num_elite), sz(num_ants)));
	static const bool debug = std::getenv("VINA_ACO_DEBUG") != NULL;
	int evalcount = 0;
	// expected number of DOFs explored per ant = acor_explore (position/orientation count as one each)
	const fl p_explore = std::min(fl(1), acor_explore / fl(std::max(sz(1), s.num_degrees_of_freedom() - 4)));
	const auto budget_left = [&]() { return max_evals == 0 || evalcount < int(max_evals); };

	// Short local search of an ant under the capped ("hunt") potential, then a rescore under the
	// authentic potential, so that every archive energy is measured with the same objective.
	// With skip_ls, an ant whose local search is vetoed by the start filter gets e = max_fl and is discarded;
	// it costs one evaluation. Random (re)initialization ants are never vetoed but feed the history.
	visited_history history(s);
	const auto ant_search = [&](output_type& o) {
		if(skip_ls) {
			ant_ls(m, p, ig, o, g, hunt_cap, evalcount, &history);
			if(history.skipped) {
				o.e = max_fl;
				return false;
			}
		}
		else
			ant_ls(m, p, ig, o, g, hunt_cap, evalcount);
		o.e = m.eval_deriv(p, ig, authentic_v, g); // model is already set to o.c by quasi_newton
		++evalcount;
		return true;
	};
	// successive halving: every ant gets half of the local steps, the better half gets the rest
	quasi_newton half_ls; half_ls.max_steps = std::max(1u, local_steps / 2);
	const auto ant_search_half = [&](output_type& o) {
		half_ls(m, p, ig, o, g, hunt_cap, evalcount);
		o.e = m.eval_deriv(p, ig, authentic_v, g);
		++evalcount;
	};
	// Deeper local search under the authentic potential; the result is a candidate output pose.
	const auto elite_search = [&](output_type& o) {
		elite_ls(m, p, ig, o, g, authentic_v, evalcount);
		o.coords = m.get_heavy_atom_movable_coords();
		add_to_output_container(out, o, min_rmsd, num_saved_mins);
	};

	fl init_best = max_fl; // best energy after the first initialization, for the degenerate-case check
	boost::ptr_vector<output_type> pool; // owns archive members
	std::vector<output_type*> archive;
	flv cumulative;
	fl best_e = max_fl; // best of the current phase
	unsigned since_improvement = 0, restarts = 0, iterations_done = 0, best_iteration = 0;
	// consensus stopping: best pose of each finished phase vs the best pose of all earlier phases
	output_type phase_best(conf(s), max_fl), overall_best(conf(s), max_fl);
	unsigned confirmations = 0;
	const fl consensus_energy_tol = 1.0, consensus_rmsd_tol = 2.0;
	// returns true when the search may stop because the optimum was re-found enough times
	const auto finish_phase = [&]() {
		if(phase_best.e >= max_fl) return false;
		if(phases) phases->push_back(new output_type(phase_best));
		if(overall_best.e < max_fl && phase_best.e > overall_best.e - consensus_energy_tol &&
		   phase_best.e < overall_best.e + consensus_energy_tol &&
		   rmsd_upper_bound(phase_best.coords, overall_best.coords) < consensus_rmsd_tol)
			++confirmations;
		else if(phase_best.e < overall_best.e)
			confirmations = 0;
		if(phase_best.e < overall_best.e)
			overall_best = phase_best;
		phase_best.e = max_fl;
		return consensus > 0 && confirmations >= consensus;
	};

	// rank weights of the archive (Gaussian in the rank, width q*k)
	const auto update_weights = [&]() {
		cumulative.resize(archive.size());
		fl acc = 0;
		VINA_FOR_IN(l, archive) {
			acc += std::exp(-fl(l * l) / (2 * sqr(acor_q * k)));
			cumulative[l] = acc;
		}
	};
	// fills the archive with up to k random solutions after a short local search;
	// the best of them is refined and offered to the output container
	const auto initialize = [&]() {
		archive.clear();
		VINA_FOR(a, k) {
			if(a > 0 && !budget_left())
				break;
			output_type* o = new output_type(conf(s), max_fl);
			o->c.randomize(corner1, corner2, generator);
			history.allow_skip = false;
			ant_search(*o);
			history.allow_skip = true;
			pool.push_back(o);
			archive.push_back(o);
		}
		std::sort(archive.begin(), archive.end(), by_energy());
		output_type best(*archive.front());
		elite_search(best);
		if(init_best >= max_fl) init_best = out.empty() ? best.e : out.front().e;
		update_weights();
	};
	initialize();

	std::vector<const output_type*> guides;
	boost::ptr_vector<output_type> ants;
	VINA_FOR(a, num_ants)
		ants.push_back(new output_type(conf(s), max_fl));
	std::vector<output_type*> ranked(num_ants);

	VINA_U_FOR(iteration, max_iterations) {
		if(!budget_left())
			break;
		iterations_done = iteration + 1;
		if(increment_me)
			++(*increment_me);

		guides.assign(archive.begin(), archive.end());
		sz built = 0;
		VINA_FOR(a, num_ants) {
			if(built > 0 && !budget_left())
				break;
			const fl r = random_fl(0, cumulative.back(), generator);
			const sz l = std::min(sz(std::lower_bound(cumulative.begin(), cumulative.end(), r) - cumulative.begin()), guides.size() - 1);
			output_type& ant = ants[a];
			acor_sample(ant.c, guides, l, acor_xi, p_explore, corner1, corner2, generator);
			if(ls_halving) ant_search_half(ant);
			else if(!ant_search(ant)) continue;
			ranked[built++] = &ant;
		}
		std::sort(ranked.begin(), ranked.begin() + built, by_energy());
		if(ls_halving) {
			VINA_FOR(a, (built + 1) / 2)
				ant_search_half(*ranked[a]);
			std::sort(ranked.begin(), ranked.begin() + built, by_energy());
		}

		// deeper (Lamarckian) refinement of the iteration elite -> output poses
		bool improved = false;
		VINA_FOR(r, std::min(n_elite, built)) {
			output_type& e = *ranked[r];
			elite_search(e);
			if(e.e < phase_best.e)
				phase_best = e;
			if(e.e < best_e - 1e-3) {
				best_e = e.e;
				best_iteration = iteration;
				improved = true;
			}
		}

		// archive update: keep the best k of archive + ants (all energies authentic)
		VINA_FOR(a, built) {
			if(archive.size() < k || ranked[a]->e < archive.back()->e) {
				output_type* o = new output_type(*ranked[a]);
				pool.push_back(o);
				if(archive.size() < k) archive.push_back(o);
				else archive.back() = o;
				std::sort(archive.begin(), archive.end(), by_energy());
			}
		}
		update_weights();

		// drop pool members no longer referenced, to keep memory bounded
		if(pool.size() > 4 * k) {
			boost::ptr_vector<output_type> kept;
			VINA_FOR_IN(l, archive) {
				output_type* o = new output_type(*archive[l]);
				kept.push_back(o);
				archive[l] = o;
			}
			pool.swap(kept);
		}

		if(improved)
			since_improvement = 0;
		else if(++since_improvement >= patience) {
			if(finish_phase() || restarts >= max_restarts || !budget_left())
				break;
			++restarts;
			since_improvement = 0;
			best_e = max_fl;
			initialize();
		}
	}
	if(phases && phase_best.e < max_fl)
		phases->push_back(new output_type(phase_best));
	if(debug)
		std::cerr << "ACO_R colony: iterations " << iterations_done << " evals " << evalcount << " restarts " << restarts
		          << " confirmations " << confirmations
		          << " skipped_ls " << history.n_skipped << "/" << history.n_checked
		          << " init_best " << init_best
		          << " best " << out.front().e << " last_improvement_iter " << best_iteration << "\n";
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e);
}

namespace {

struct aco_task {
	model m;
	output_container out;
	output_container phases;
	rng generator;
	aco_task(const model& m_, int seed) : m(m_), generator(static_cast<rng::result_type>(seed)) {}
};
typedef boost::ptr_vector<aco_task> aco_task_container;

struct aco_aux {
	const ant_colony* aco;
	const precalculate_byatom* p;
	const igrid* ig;
	const vec* corner1;
	const vec* corner2;
	parallel_progress* pg;
	aco_aux(const ant_colony* aco_, const precalculate_byatom* p_, const igrid* ig_, const vec* corner1_, const vec* corner2_, parallel_progress* pg_)
		: aco(aco_), p(p_), ig(ig_), corner1(corner1_), corner2(corner2_), pg(pg_) {}
	void operator()(aco_task& t) const {
		(*aco)(t.m, t.out, *p, *ig, *corner1, *corner2, pg, t.generator, &t.phases);
	}
};

} // namespace

void parallel_aco::operator()(const model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                              const vec& corner1, const vec& corner2, rng& generator, std::function<void(double)>* progress_callback,
                              output_container* phases) const {
	parallel_progress pp(progress_callback);
	aco_aux aux(&aco, &p, &ig, &corner1, &corner2, (display_progress ? (&pp) : NULL));
	aco_task_container tasks;
	VINA_FOR(i, num_tasks)
		tasks.push_back(new aco_task(m, random_int(0, 1000000, generator)));
	if(display_progress)
		pp.init(num_tasks * aco.max_iterations);
	parallel_iter<aco_aux, aco_task_container, aco_task, true> iter(&aux, num_threads);
	iter.run(tasks);
	VINA_FOR_IN(i, tasks)
		VINA_FOR_IN(j, tasks[i].out)
			add_to_output_container(out, tasks[i].out[j], aco.min_rmsd, aco.num_saved_mins);
	out.sort();
	if(phases)
		VINA_FOR_IN(i, tasks)
			VINA_FOR_IN(j, tasks[i].phases)
				phases->push_back(new output_type(tasks[i].phases[j]));
}
