#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "ControlAllocation.h"

namespace {

using Rows = std::vector<std::vector<double>>;

struct InstanceInput {
    Rows B;                 // 6 x m, rows [Mx My Mz Fx Fy Fz]
    Rows v_sp;              // N x 6
    std::vector<double> umin;
    std::vector<double> umax;
    std::vector<int> active_rows;
};

struct InstanceResult {
    Rows u;                 // N x m
    Rows v_achieved;        // N x 6
    double allocator_s{0.0};
    double restore_s{0.0};
    int fail_count{0};
    int fallback_count{0};
};

struct MethodResult {
    Rows u;                 // N x total actuator count
    Rows v_achieved;        // N x 6
    Rows residual;          // v_sp - v_achieved
    double total_s{0.0};
    double allocator_s{0.0};
    double restore_s{0.0};
    int fail_count{0};
    int fallback_count{0};
};

std::string path_join(const std::string &dir, const std::string &name)
{
    return dir.empty() || dir.back() == '/' ? dir + name : dir + "/" + name;
}

Rows read_csv(const std::string &path)
{
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("cannot open " + path);
    }

    Rows rows;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            if (!value.empty()) {
                row.push_back(std::stod(value));
            }
        }

        if (!row.empty()) {
            rows.push_back(row);
        }
    }

    return rows;
}

std::vector<double> read_row(const std::string &path)
{
    Rows rows = read_csv(path);
    if (rows.size() != 1) {
        throw std::runtime_error("expected one-row vector: " + path);
    }
    return rows[0];
}

void write_csv(const std::string &path, const Rows &rows)
{
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("cannot write " + path);
    }

    file << std::setprecision(17);
    for (const auto &row : rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i] << (i + 1 < row.size() ? "," : "\n");
        }
    }
}

int ncol(const Rows &A)
{
    return A.empty() ? 0 : static_cast<int>(A[0].size());
}

std::vector<int> active_rows(const Rows &B)
{
    std::vector<int> rows;

    for (int r = 0; r < static_cast<int>(B.size()); ++r) {
        for (double value : B[r]) {
            if (std::fabs(value) > 1.0e-9) {
                rows.push_back(r);
                break;
            }
        }
    }

    return rows;
}

Rows active_B(const InstanceInput &inst)
{
    Rows B;
    for (int r : inst.active_rows) {
        B.push_back(inst.B[r]);
    }
    return B;
}

std::vector<double> active_v(const InstanceInput &inst, int sample)
{
    std::vector<double> v;
    for (int r : inst.active_rows) {
        v.push_back(inst.v_sp[sample][r]);
    }
    return v;
}

std::vector<double> B_times_u(const Rows &B, const std::vector<double> &u)
{
    std::vector<double> v(B.size(), 0.0);

    for (size_t r = 0; r < B.size(); ++r) {
        for (size_t c = 0; c < B[r].size(); ++c) {
            v[r] += B[r][c] * u[c];
        }
    }

    return v;
}

std::vector<double> clamp_u(std::vector<double> u,
                            const std::vector<double> &umin,
                            const std::vector<double> &umax)
{
    for (size_t i = 0; i < u.size(); ++i) {
        u[i] = std::min(std::max(u[i], umin[i]), umax[i]);
    }
    return u;
}

std::vector<double> pinv_alloc(const Rows &B,
                               const std::vector<double> &v,
                               const std::vector<double> &umin,
                               const std::vector<double> &umax)
{
    const int k = static_cast<int>(B.size());
    const int m = ncol(B);
    Eigen::MatrixXf B_eig(k, m);
    Eigen::VectorXf v_eig(k);

    for (int r = 0; r < k; ++r) {
        v_eig(r) = static_cast<float>(v[r]);
        for (int c = 0; c < m; ++c) {
            B_eig(r, c) = static_cast<float>(B[r][c]);
        }
    }

    Eigen::JacobiSVD<Eigen::MatrixXf> svd(B_eig, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXf u_eig = svd.solve(v_eig);
    std::vector<double> u(m, 0.0);

    for (int i = 0; i < m; ++i) {
        u[i] = static_cast<double>(u_eig(i));
    }

    return clamp_u(u, umin, umax);
}

template<int K, int M>
InstanceResult run_lpca_instance(const InstanceInput &inst,
                                 const std::string &method,
                                 bool use_restoring)
{
    Rows B_eff = active_B(inst);
    const int N = static_cast<int>(inst.v_sp.size());
    InstanceResult result;
    result.u.assign(N, std::vector<double>(M, 0.0));
    result.v_achieved.assign(N, std::vector<double>(6, 0.0));

    float B_array[K][M];
    float umin_array[M];
    float umax_array[M];

    for (int r = 0; r < K; ++r) {
        for (int c = 0; c < M; ++c) {
            B_array[r][c] = static_cast<float>(B_eff[r][c]);
        }
    }

    for (int c = 0; c < M; ++c) {
        umin_array[c] = static_cast<float>(inst.umin[c]);
        umax_array[c] = static_cast<float>(inst.umax[c]);
    }

    Aircraft<K, M> aircraft(B_array, umax_array, umin_array);
    DP_LP_ControlAllocator<K, M> allocator(aircraft);

    for (int sample = 0; sample < N; ++sample) {
        float v[K];
        float raw[M];
        float u[M];
        int err = 0;
        float rho = 0.0f;
        std::vector<double> v_eff = active_v(inst, sample);

        for (int r = 0; r < K; ++r) {
            v[r] = static_cast<float>(v_eff[r]);
        }

        const auto alloc_start = std::chrono::high_resolution_clock::now();
        if (method == "pca_dir") {
            allocator.DP_LPCA(v, raw, err, rho);
        } else {
            allocator.DPscaled_LPCA(v, raw, err, rho);
        }
        const auto alloc_end = std::chrono::high_resolution_clock::now();

        if (use_restoring) {
            const auto restore_start = std::chrono::high_resolution_clock::now();
            allocator.restoring(raw, u);
            const auto restore_end = std::chrono::high_resolution_clock::now();
            result.restore_s += std::chrono::duration<double>(restore_end - restore_start).count();
        } else {
            for (int i = 0; i < M; ++i) {
                u[i] = raw[i];
            }
        }

        result.allocator_s += std::chrono::duration<double>(alloc_end - alloc_start).count();

        for (int i = 0; i < M; ++i) {
            result.u[sample][i] = std::min(std::max(static_cast<double>(u[i]), inst.umin[i]), inst.umax[i]);
        }

        result.v_achieved[sample] = B_times_u(inst.B, result.u[sample]);
        result.fail_count += err != 0;
    }

    return result;
}

InstanceResult run_pinv_instance(const InstanceInput &inst)
{
    const int N = static_cast<int>(inst.v_sp.size());
    const int m = ncol(inst.B);
    Rows B_eff = active_B(inst);
    InstanceResult result;
    result.u.assign(N, std::vector<double>(m, 0.0));
    result.v_achieved.assign(N, std::vector<double>(6, 0.0));
    result.fallback_count = N;

    for (int sample = 0; sample < N; ++sample) {
        std::vector<double> v_eff = active_v(inst, sample);
        const auto start = std::chrono::high_resolution_clock::now();
        result.u[sample] = pinv_alloc(B_eff, v_eff, inst.umin, inst.umax);
        const auto finish = std::chrono::high_resolution_clock::now();
        result.allocator_s += std::chrono::duration<double>(finish - start).count();
        result.v_achieved[sample] = B_times_u(inst.B, result.u[sample]);
    }

    return result;
}

InstanceResult run_instance(const InstanceInput &inst,
                            const std::string &method,
                            bool use_restoring)
{
    const int k = static_cast<int>(inst.active_rows.size());
    const int m = ncol(inst.B);

    if (k <= 1 || m <= k) {
        return run_pinv_instance(inst);
    }

    if (k == 3 && m == 4) return run_lpca_instance<3, 4>(inst, method, use_restoring);
    if (k == 3 && m == 6) return run_lpca_instance<3, 6>(inst, method, use_restoring);
    if (k == 3 && m == 8) return run_lpca_instance<3, 8>(inst, method, use_restoring);
    if (k == 4 && m == 5) return run_lpca_instance<4, 5>(inst, method, use_restoring);
    if (k == 4 && m == 7) return run_lpca_instance<4, 7>(inst, method, use_restoring);

    return run_pinv_instance(inst);
}

std::vector<InstanceInput> read_inputs(const std::string &input_dir)
{
    Rows meta = read_csv(path_join(input_dir, "case_meta.csv"));
    const int num_instances = static_cast<int>(std::lround(meta[0][2]));
    std::vector<InstanceInput> instances;

    for (int i = 0; i < num_instances; ++i) {
        const std::string prefix = path_join(input_dir, "inst_" + std::to_string(i));
        InstanceInput inst;
        inst.B = read_csv(prefix + "_B.csv");
        inst.v_sp = read_csv(prefix + "_Y.csv");
        inst.umin = read_row(prefix + "_umin.csv");
        inst.umax = read_row(prefix + "_umax.csv");
        inst.active_rows = active_rows(inst.B);
        instances.push_back(inst);
    }

    return instances;
}

Rows combined_v_sp(const std::vector<InstanceInput> &instances)
{
    const int N = static_cast<int>(instances[0].v_sp.size());
    Rows v_sp(N, std::vector<double>(6, 0.0));

    for (const auto &inst : instances) {
        for (int sample = 0; sample < N; ++sample) {
            for (int axis = 0; axis < 6; ++axis) {
                v_sp[sample][axis] += inst.v_sp[sample][axis];
            }
        }
    }

    return v_sp;
}

MethodResult run_method(const std::vector<InstanceInput> &instances,
                        const std::string &method,
                        bool use_restoring)
{
    const int N = static_cast<int>(instances[0].v_sp.size());
    int total_u_dim = 0;
    for (const auto &inst : instances) {
        total_u_dim += ncol(inst.B);
    }

    MethodResult result;
    result.u.assign(N, std::vector<double>(total_u_dim, 0.0));
    result.v_achieved.assign(N, std::vector<double>(6, 0.0));
    const Rows v_sp = combined_v_sp(instances);
    int col0 = 0;

    const auto total_start = std::chrono::high_resolution_clock::now();
    for (const auto &inst : instances) {
        InstanceResult inst_result = run_instance(inst, method, use_restoring);

        for (int sample = 0; sample < N; ++sample) {
            for (int c = 0; c < ncol(inst.B); ++c) {
                result.u[sample][col0 + c] = inst_result.u[sample][c];
            }
            for (int axis = 0; axis < 6; ++axis) {
                result.v_achieved[sample][axis] += inst_result.v_achieved[sample][axis];
            }
        }

        col0 += ncol(inst.B);
        result.allocator_s += inst_result.allocator_s;
        result.restore_s += inst_result.restore_s;
        result.fail_count += inst_result.fail_count;
        result.fallback_count += inst_result.fallback_count;
    }
    const auto total_end = std::chrono::high_resolution_clock::now();

    result.total_s = std::chrono::duration<double>(total_end - total_start).count();
    result.residual = v_sp;
    for (int sample = 0; sample < N; ++sample) {
        for (int axis = 0; axis < 6; ++axis) {
            result.residual[sample][axis] -= result.v_achieved[sample][axis];
        }
    }

    return result;
}

void write_method_result(const std::string &output_dir,
                         const std::string &method,
                         const MethodResult &result)
{
    write_csv(path_join(output_dir, method + "_u.csv"), result.u);
    write_csv(path_join(output_dir, method + "_v_achieved.csv"), result.v_achieved);
    write_csv(path_join(output_dir, method + "_residual.csv"), result.residual);
    write_csv(path_join(output_dir, method + "_timing.csv"),
              {{result.total_s, result.allocator_s, result.restore_s,
                static_cast<double>(result.fail_count),
                static_cast<double>(result.fallback_count)}});
}

} // namespace

int main(int argc, char **argv)
{
    try {
        if (argc != 3) {
            std::cerr << "usage: ca_offline_benchmark <input_dir> <output_dir>\n";
            return 2;
        }

        const std::string input_dir = argv[1];
        const std::string output_dir = argv[2];
        mkdir(output_dir.c_str(), 0755);

        Rows meta = read_csv(path_join(input_dir, "case_meta.csv"));
        const bool use_restoring = meta[0].size() >= 4 && std::lround(meta[0][3]) != 0;
        const std::vector<InstanceInput> instances = read_inputs(input_dir);

        write_method_result(output_dir, "pca_dir", run_method(instances, "pca_dir", use_restoring));
        write_method_result(output_dir, "pca_dpscaled", run_method(instances, "pca_dpscaled", use_restoring));

        std::cout << "Wrote C++ allocation CSVs to " << output_dir << "\n";
        return 0;

    } catch (const std::exception &e) {
        std::cerr << "ca_offline_benchmark failed: " << e.what() << "\n";
        return 1;
    }
}
