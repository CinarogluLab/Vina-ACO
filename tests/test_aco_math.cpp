// Independent invariants for the actual internal pheromone/orientation implementation.
// Includes implementation solely to access private helpers; do not also link aco.o.
#include <cassert>
#include <iostream>
#include "../src/lib/aco.cpp"

int main() {
    rng gen(143);
    trail periodic(-pi, pi, 72, true);
    periodic.reset(1);
    assert(periodic.bin_of(-pi) == periodic.bin_of(pi));
    assert(periodic.bin_of(-pi-0.001) == 71);
    periodic.deposit(-pi, 4);
    assert(periodic.tau[0] == 5 && periodic.tau[71] == 3 && periodic.tau[1] == 3);
    periodic.evaporate(0.5);
    periodic.clamp(0.75, 2);
    for(auto t: periodic.tau) assert(t >= 0.75 && t <= 2);
    trail position(-10, 10, 80, false);
    position.reset(1);
    for(int i=0;i<10000;++i) {
        const fl x=position.sample(gen);
        assert(x>=-10 && x<=10);
    }
    position.reset(0.001);
    position.deposit(0,100);
    int near=0;
    for(int i=0;i<10000;++i) if(std::abs(position.sample(gen))<0.5) ++near;
    assert(near>9900); // Reinforcement must actually affect construction.
    double means[4]={0,0,0,0}, squares[4]={0,0,0,0};
    for(int i=0;i<20000;++i) {
        qt q=shoemake_to_quaternion(random_fl(0,1,gen),random_fl(0,1,gen),random_fl(0,1,gen));
        double u1,u2,u3;quaternion_to_shoemake(q,u1,u2,u3);
        qt restored=shoemake_to_quaternion(u1,u2,u3);
        double qv[]={q.R_component_1(),q.R_component_2(),q.R_component_3(),q.R_component_4()};
        double rv[]={restored.R_component_1(),restored.R_component_2(),restored.R_component_3(),restored.R_component_4()};
        double norm=0,dot=0;
        for(int j=0;j<4;++j){norm+=qv[j]*qv[j];dot+=qv[j]*rv[j];means[j]+=qv[j];squares[j]+=qv[j]*qv[j];}
        assert(std::abs(norm-1)<1e-10);
        assert(std::abs(std::abs(dot)-1)<1e-8); // same rotation modulo q/-q
    }
    for(int j=0;j<4;++j){assert(std::abs(means[j]/20000)<0.025);assert(std::abs(squares[j]/20000-0.25)<0.015);}
    assert(angle_distance(pi-0.01,-pi+0.01)<0.021);
    qt identity(1,0,0,0), opposite(-1,0,0,0);
    assert(orientation_distance(identity,opposite)<1e-8);
    conf_size size; size.ligands.push_back(2);
    output_type guide(conf(size),-1);
    guide.c.ligands[0].rigid.position=vec(1,2,3);
    guide.c.ligands[0].rigid.orientation=identity;
    guide.c.ligands[0].torsions[0]=pi-0.001;
    guide.c.ligands[0].torsions[1]=-pi+0.001;
    std::vector<const output_type*> archive{&guide,&guide};
    conf sampled(size);
    acor_sample(sampled,archive,0,0.85,0.0,vec(-10,-10,-10),vec(10,10,10),gen);
    for(int d=0;d<3;++d) assert(sampled.ligands[0].rigid.position[d]==guide.c.ligands[0].rigid.position[d]);
    assert(orientation_distance(sampled.ligands[0].rigid.orientation,identity)<1e-8);
    for(int t=0;t<2;++t) assert(angle_distance(sampled.ligands[0].torsions[t],guide.c.ligands[0].torsions[t])<1e-8);
    ant_colony params;
    assert(params.validate().empty());
    params.num_ants=0; assert(!params.validate().empty());
    // the initial orientations must be Haar-uniform on SO(3), i.e. the rotation angle theta has
    // density (1 - cos theta)/pi on [0, pi]. One-sample Kolmogorov-Smirnov against the analytic CDF
    // F(theta) = (theta - sin theta)/pi, for Vina's random_orientation and for the Shoemake map used by MMAS.
    {
        const int N = 20000;
        const double crit = 1.63 / std::sqrt(double(N)); // ~1% two-sided critical value
        for(int which = 0; which < 2; ++which) {
            std::vector<double> ang;
            ang.reserve(N);
            for(int i = 0; i < N; ++i) {
                qt q = (which == 0) ? random_orientation(gen)
                                    : shoemake_to_quaternion(random_fl(0, 1, gen), random_fl(0, 1, gen), random_fl(0, 1, gen));
                ang.push_back(orientation_distance(q, qt_identity)); // rotation angle in [0, pi]
            }
            std::sort(ang.begin(), ang.end());
            double d = 0;
            for(int i = 0; i < N; ++i) {
                const double F = (ang[i] - std::sin(ang[i])) / pi;
                d = std::max(d, std::max(std::abs(F - double(i) / N), std::abs(double(i + 1) / N - F)));
            }
            assert(d < crit); // both samplers follow the Haar law
        }
    }
    // the SO(3) kernel width must scale with the archive spread, as eq 2 states. With an archive
    // whose members sit at a known rotation angle from the guide, the sampled perturbation scales with xi.
    {
        const conf_size s2{{0}, {}};
        output_type guide(conf(s2), 0), other(conf(s2), 0);
        guide.c.ligands[0].rigid.orientation = qt_identity;
        const double spread = 0.4;
        other.c.ligands[0].rigid.orientation = angle_to_quaternion(vec(spread, 0, 0));
        std::vector<const output_type*> archive{&guide, &other};
        double prev = 0;
        for(double xi : {0.25, 0.5, 1.0}) {
            double sum = 0;
            const int M = 4000;
            for(int i = 0; i < M; ++i) {
                conf c(s2);
                acor_sample(c, archive, 0, xi, 0.0, vec(-10, -10, -10), vec(10, 10, 10), gen);
                sum += orientation_distance(c.ligands[0].rigid.orientation, qt_identity);
            }
            const double mean_angle = sum / M;
            assert(mean_angle > prev); // wider kernel -> larger perturbations
            // expected per-component sigma is xi * spread / sqrt(3); the mean rotation angle tracks it
            assert(mean_angle > 0.25 * xi * spread && mean_angle < 4.0 * xi * spread);
            prev = mean_angle;
        }
    }
    // QuickVina 2 start filter: 1 ligand, 2 torsions -> n = 8 variables, history capacity 80
    {
        conf_size s2; s2.ligands.push_back(2);
        const auto point = [&](fl pos, fl tors, fl grad_rigid, fl grad_tors, conf& c, change& g) {
            c = conf(s2); g = change(s2);
            c.ligands[0].rigid.position = vec(pos, pos, pos);
            c.ligands[0].rigid.orientation = qt_identity;
            for(int t=0;t<2;++t) { c.ligands[0].torsions[t] = tors; g.ligands[0].torsions[t] = grad_tors; }
            // orientation gradient 0: those variables are ignored (identical orientations would count as "not further along")
            for(int d=0;d<3;++d) { g.ligands[0].rigid.position[d] = grad_rigid; g.ligands[0].rigid.orientation[d] = 0; }
        };
        // history of 80 copies of q; returns whether p is skipped
        const auto skipped = [&](fl q_pos, fl q_tors, fl q_f, fl q_gr, fl q_gt, fl p_pos, fl p_tors, fl p_f, fl p_gr, fl p_gt) {
            visited_history h(s2);
            assert(h.n == 8 && h.capacity == 80);
            conf c; change g(s2);
            point(q_pos, q_tors, q_gr, q_gt, c, g);
            assert(h.start(c, q_f, g) && !h.skipped); // not full yet: always minimized
            for(int i=1;i<80;++i) h.step(c, q_f, g);
            point(p_pos, p_tors, p_gr, p_gt, c, g);
            const bool go = h.start(c, p_f, g);
            assert(go == !h.skipped);
            return h.skipped;
        };
        // same positive slope, p further along it and higher: monotone, no stationary point between -> skip (Eq. 1 and 2 fail)
        assert(skipped(0, 0, 0, 1, 1,   1, 1, 1, 1, 1));
        // same positive slope but p further along and LOWER: a stationary point must lie between -> minimize (Eq. 2)
        assert(!skipped(0, 0, 0, 1, 1,   1, 1, -1, 1, 1));
        // opposite gradient signs in every variable -> minimize (Eq. 1)
        assert(!skipped(0, 0, 0, 1, 1,   1, 1, 1, -1, -1));
        // torsion displacement across the +-pi seam: 3.0 -> -3.0 is +0.28 rad, so p is further along the slope and higher -> skip
        // (without modulo-2pi differences it would look like -6 rad and p would be minimized)
        assert(skipped(0, 3.0, 0, 0, 1,   0, -3.0, 1, 0, 1));
        // Equation 2 has a zero product when either displacement or energy difference is zero; both must pass.
        assert(!skipped(0, 0, 0, 1, 1,   0, 0, 1, 1, 1));
        assert(!skipped(0, 0, 0, 1, 1,   1, 1, 0, 1, 1));
    }
    std::cout << "PASS: trail bounds, wrapping, evaporation, reinforcement, quaternion roundtrip/uniformity, SO3 sign equivalence, collapsed ACOR kernel, API validation, QuickVina 2 start filter, Haar uniformity (KS), SO(3) kernel scaling\n";
}
