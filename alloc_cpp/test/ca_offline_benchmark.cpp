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
#include <sys/types.h>
#include <vector>

#include "ControlAllocation.h"

namespace {

using Rows = std::vector<std::vector<double>>;

struct InstanceInput {
    Rows B;
    Rows Y;
    std::vector<double> umin;
    std::vector<double> umax;
    std::vector<int> active_rows;
};

struct SeriesResult {
    Rows u;
    Rows raw;
    Rows y_achieved;
    Rows errout;
    Rows rho;
    double allocator_s{0.0};
    double restore_s{0.0};
    int fail_count{0};
    int fallback_count{0};
};

std::string join_path(const std::string &a, const std::string &b)
{
    if (a.empty()) {
        return b;
    }

    if (a.back() == '/') {
        return a + b;
    }

    return a + "/" + b;
}

void ensure_dir(const std::string &path)
{
    mkdir(path.c_str(), 0755);
}

bool file_exists(const std::string &path)
{
    std::ifstream file(path);
    return file.good();
}

Rows read_csv(const std::string &path)
{
    std::ifstream file(path);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open CSV: " + path);
    }

    Rows rows;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream line_stream(line);
        std::string value;
        std::vector<double> row;

        while (std::getline(line_stream, value, ',')) {
            if (!value.empty()) {
                row.push_back(std::stod(value));
            }
        }

        if (!row.empty()) {
            rows.push_back(row);
        }
    }

    if (rows.empty()) {
        throw std::runtime_error("Empty CSV: " + path);
    }

    return rows;
}

std::vector<double> read_row_vector(const std::string &path)
{
    Rows rows = read_csv(path);

    if (rows.size() != 1) {
        throw std::runtime_error("Expected one-row vector: " + path);
    }

    return rows[0];
}

void write_csv(const std::string &path, const Rows &rows)
{
    std::ofstream file(path);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot write CSV: " + path);
    }

    file << std::setprecision(17);

    for (const auto &row : rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i] << (i + 1 < row.size() ? "," : "\n");
        }
    }
}

void write_lines(const std::string &path, const std::vector<std::string> &lines)
{
    std::ofstream file(path);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot write text file: " + path);
    }

    for (const auto &line : lines) {
        file << line << "\n";
    }
}

int num_cols(const Rows &rows)
{
    return rows.empty() ? 0 : static_cast<int>(rows[0].size());
}

std::vector<int> active_rows_from_B(const Rows &B)
{
    std::vector<int> rows;

    for (int r = 0; r < static_cast<int>(B.size()); ++r) {
        bool active = false;

        for (double value : B[r]) {
            if (std::fabs(value) > 1.0e-9) {
                active = true;
                break;
            }
        }

        if (active) {
            rows.push_back(r);
        }
    }

    return rows;
}

Rows extract_active_B(const InstanceInput &inst)
{
    Rows B_eff;
    B_eff.reserve(inst.active_rows.size());

    for (int row : inst.active_rows) {
        B_eff.push_back(inst.B[row]);
    }

    return B_eff;
}

std::vector<double> extract_active_y(const InstanceInput &inst, int sample)
{
    std::vector<double> y;
    y.reserve(inst.active_rows.size());

    for (int row : inst.active_rows) {
        y.push_back(inst.Y[sample][row]);
    }

    return y;
}

std::vector<double> mat_vec(const Rows &A, const std::vector<double> &x)
{
    std::vector<double> y(A.size(), 0.0);

    for (size_t r = 0; r < A.size(); ++r) {
        for (size_t c = 0; c < A[r].size(); ++c) {
            y[r] += A[r][c] * x[c];
        }
    }

    return y;
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
                               const std::vector<double> &y,
                               const std::vector<double> &umin,
                               const std::vector<double> &umax)
{
    const int k = static_cast<int>(B.size());
    const int m = num_cols(B);

    Eigen::MatrixXf B_eig(k, m);
    Eigen::VectorXf y_eig(k);

    for (int r = 0; r < k; ++r) {
        y_eig(r) = static_cast<float>(y[r]);

        for (int c = 0; c < m; ++c) {
            B_eig(r, c) = static_cast<float>(B[r][c]);
        }
    }

    Eigen::JacobiSVD<Eigen::MatrixXf> svd(B_eig, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXf u_eig = svd.solve(y_eig);

    std::vector<double> u(m, 0.0);

    for (int i = 0; i < m; ++i) {
        u[i] = static_cast<double>(u_eig(i));
    }

    return clamp_u(u, umin, umax);
}

template<int K, int M>
SeriesResult run_lpca_template(const InstanceInput &inst,
                               const std::string &method,
                               bool use_restoring)
{
    Rows B_eff = extract_active_B(inst);
    const int N = static_cast<int>(inst.Y.size());
    SeriesResult result;
    result.u.assign(N, std::vector<double>(M, 0.0));
    result.raw.assign(N, std::vector<double>(M, 0.0));
    result.y_achieved.assign(N, std::vector<double>(6, 0.0));
    result.errout.assign(N, std::vector<double>(1, 0.0));
    result.rho.assign(N, std::vector<double>(1, 0.0));

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
        float y[K];
        float raw[M];
        float u[M];
        int err = 0;
        float rho = 0.0f;
        std::vector<double> y_eff = extract_active_y(inst, sample);

        for (int r = 0; r < K; ++r) {
            y[r] = static_cast<float>(y_eff[r]);
        }

        const auto alloc_start = std::chrono::high_resolution_clock::now();

        if (method == "pca_dir") {
            allocator.DP_LPCA(y, raw, err, rho);
        } else if (method == "pca_dpscaled") {
            allocator.DPscaled_LPCA(y, raw, err, rho);
        } else {
            throw std::runtime_error("Unsupported C++ method: " + method);
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

        std::vector<double> u_row(M, 0.0);

        for (int i = 0; i < M; ++i) {
            result.raw[sample][i] = raw[i];
            u_row[i] = std::min(std::max(static_cast<double>(u[i]), inst.umin[i]), inst.umax[i]);
        }

        for (int i = 0; i < M; ++i) {
            result.u[sample][i] = u_row[i];
        }

        result.y_achieved[sample] = mat_vec(inst.B, u_row);
        result.errout[sample][0] = static_cast<double>(err);
        result.rho[sample][0] = static_cast<double>(rho);

        if (err != 0) {
            result.fail_count += 1;
        }
    }

    return result;
}

SeriesResult run_pinv_series(const InstanceInput &inst)
{
    const int N = static_cast<int>(inst.Y.size());
    const int m = num_cols(inst.B);
    SeriesResult result;
    result.u.assign(N, std::vector<double>(m, 0.0));
    result.raw.assign(N, std::vector<double>(m, 0.0));
    result.y_achieved.assign(N, std::vector<double>(6, 0.0));
    result.errout.assign(N, std::vector<double>(1, 0.0));
    result.rho.assign(N, std::vector<double>(1, 0.0));
    result.fallback_count = N;

    Rows B_eff = extract_active_B(inst);

    for (int sample = 0; sample < N; ++sample) {
        std::vector<double> y_eff = extract_active_y(inst, sample);
        const auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> u = pinv_alloc(B_eff, y_eff, inst.umin, inst.umax);
        const auto finish = std::chrono::high_resolution_clock::now();
        result.allocator_s += std::chrono::duration<double>(finish - start).count();

        result.u[sample] = u;
        result.raw[sample] = u;
        result.y_achieved[sample] = mat_vec(inst.B, u);
    }

    return result;
}

SeriesResult run_instance_series(const InstanceInput &inst,
                                 const std::string &method,
                                 bool use_restoring)
{
    const int k = static_cast<int>(inst.active_rows.size());
    const int m = num_cols(inst.B);

    if (k <= 1 || m <= k) {
        return run_pinv_series(inst);
    }

    if (k == 3 && m == 4) {
        return run_lpca_template<3, 4>(inst, method, use_restoring);
    }

    if (k == 3 && m == 6) {
        return run_lpca_template<3, 6>(inst, method, use_restoring);
    }

    if (k == 3 && m == 8) {
        return run_lpca_template<3, 8>(inst, method, use_restoring);
    }

    if (k == 4 && m == 5) {
        return run_lpca_template<4, 5>(inst, method, use_restoring);
    }

    if (k == 4 && m == 7) {
        return run_lpca_template<4, 7>(inst, method, use_restoring);
    }

    SeriesResult result = run_pinv_series(inst);
    result.fallback_count = static_cast<int>(inst.Y.size());
    return result;
}

std::vector<InstanceInput> read_instances(const std::string &input_dir)
{
    Rows meta = read_csv(join_path(input_dir, "case_meta.csv"));

    if (meta[0].size() < 3) {
        throw std::runtime_error("case_meta.csv must contain [N,total_u_dim,num_instances,use_restoring]");
    }

    const int num_instances = static_cast<int>(std::lround(meta[0][2]));
    std::vector<InstanceInput> instances;
    instances.reserve(num_instances);

    for (int i = 0; i < num_instances; ++i) {
        const std::string prefix = join_path(input_dir, "inst_" + std::to_string(i));
        InstanceInput inst;
        inst.B = read_csv(prefix + "_B.csv");
        inst.Y = read_csv(prefix + "_Y.csv");
        inst.umin = read_row_vector(prefix + "_umin.csv");
        inst.umax = read_row_vector(prefix + "_umax.csv");
        inst.active_rows = active_rows_from_B(inst.B);

        if (inst.B.size() != 6) {
            throw std::runtime_error(prefix + "_B.csv must be 6 x m");
        }

        if (num_cols(inst.B) != static_cast<int>(inst.umin.size())
            || num_cols(inst.B) != static_cast<int>(inst.umax.size())) {
            throw std::runtime_error(prefix + " B/limit dimension mismatch");
        }

        instances.push_back(inst);
    }

    return instances;
}

Rows sum_command(const std::vector<InstanceInput> &instances)
{
    const int N = static_cast<int>(instances[0].Y.size());
    Rows y(N, std::vector<double>(6, 0.0));

    for (const auto &inst : instances) {
        for (int sample = 0; sample < N; ++sample) {
            for (int axis = 0; axis < 6; ++axis) {
                y[sample][axis] += inst.Y[sample][axis];
            }
        }
    }

    return y;
}

Rows subtract_rows(const Rows &a, const Rows &b)
{
    Rows out = a;

    for (size_t r = 0; r < out.size(); ++r) {
        for (size_t c = 0; c < out[r].size(); ++c) {
            out[r][c] -= b[r][c];
        }
    }

    return out;
}

void append_instance_output(Rows &u_all,
                            Rows &y_achieved,
                            int col0,
                            const SeriesResult &inst_result)
{
    for (size_t sample = 0; sample < u_all.size(); ++sample) {
        for (size_t c = 0; c < inst_result.u[sample].size(); ++c) {
            u_all[sample][col0 + static_cast<int>(c)] = inst_result.u[sample][c];
        }

        for (int axis = 0; axis < 6; ++axis) {
            y_achieved[sample][axis] += inst_result.y_achieved[sample][axis];
        }
    }
}

void run_method(const std::vector<InstanceInput> &instances,
                const std::string &method,
                bool use_restoring,
                const std::string &output_dir)
{
    const int N = static_cast<int>(instances[0].Y.size());
    int total_u_dim = 0;

    for (const auto &inst : instances) {
        total_u_dim += num_cols(inst.B);
    }

    Rows u_all(N, std::vector<double>(total_u_dim, 0.0));
    Rows y_achieved(N, std::vector<double>(6, 0.0));
    Rows y_command = sum_command(instances);

    double allocator_s = 0.0;
    double restore_s = 0.0;
    int fail_count = 0;
    int fallback_count = 0;
    int col0 = 0;

    const auto total_start = std::chrono::high_resolution_clock::now();

    for (const auto &inst : instances) {
        SeriesResult inst_result = run_instance_series(inst, method, use_restoring);
        append_instance_output(u_all, y_achieved, col0, inst_result);
        col0 += num_cols(inst.B);
        allocator_s += inst_result.allocator_s;
        restore_s += inst_result.restore_s;
        fail_count += inst_result.fail_count;
        fallback_count += inst_result.fallback_count;
    }

    const auto total_end = std::chrono::high_resolution_clock::now();
    const double total_s = std::chrono::duration<double>(total_end - total_start).count();
    Rows residual = subtract_rows(y_command, y_achieved);
    Rows timing = {{total_s, allocator_s, restore_s,
                    static_cast<double>(fail_count), static_cast<double>(fallback_count)}};

    write_csv(join_path(output_dir, method + "_u.csv"), u_all);
    write_csv(join_path(output_dir, method + "_y_achieved.csv"), y_achieved);
    write_csv(join_path(output_dir, method + "_residual.csv"), residual);
    write_csv(join_path(output_dir, method + "_timing.csv"), timing);

    // Per-instance diagnostics.  These keep the primary LP failure output
    // visible instead of hiding it behind a fallback.
    for (size_t i = 0; i < instances.size(); ++i) {
        SeriesResult inst_result = run_instance_series(instances[i], method, use_restoring);
        const std::string prefix = method + "_inst_" + std::to_string(i);
        write_csv(join_path(output_dir, prefix + "_raw.csv"), inst_result.raw);
        write_csv(join_path(output_dir, prefix + "_err.csv"), inst_result.errout);
        write_csv(join_path(output_dir, prefix + "_rho.csv"), inst_result.rho);
    }
}

} // namespace

int main(int argc, char **argv)
{
    try {
        if (argc < 3) {
            std::cerr << "usage: ca_offline_benchmark <input_dir> <output_dir>\n";
            return 2;
        }

        const std::string input_dir = argv[1];
        const std::string output_dir = argv[2];
        Rows meta = read_csv(join_path(input_dir, "case_meta.csv"));
        const bool use_restoring = meta[0].size() >= 4 && std::lround(meta[0][3]) != 0;

        ensure_dir(output_dir);
        std::vector<InstanceInput> instances = read_instances(input_dir);
        const std::vector<std::string> methods = {"pca_dir", "pca_dpscaled"};

        for (const auto &method : methods) {
            run_method(instances, method, use_restoring, output_dir);
        }

        write_lines(join_path(output_dir, "methods.txt"), methods);
        std::cout << "Wrote C++ offline allocation result to " << output_dir << "\n";
        return 0;

    } catch (const std::exception &e) {
        std::cerr << "ca_offline_benchmark failed: " << e.what() << "\n";
        return 1;
    }
}
