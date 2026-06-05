#include <cmath>
#include <cfloat>
#include <cstdlib>
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

using Rows = std::vector<std::vector<float>>;

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

Rows read_csv(const std::string &path, int expected_cols)
{
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open CSV: " + path);
    }

    Rows rows;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream line_stream(line);
        std::string value;
        std::vector<float> row;
        while (std::getline(line_stream, value, ',')) {
            row.push_back(std::stof(value));
        }

        if (expected_cols > 0 && static_cast<int>(row.size()) != expected_cols) {
            throw std::runtime_error("Bad column count in " + path);
        }
        rows.push_back(row);
    }

    if (rows.empty()) {
        throw std::runtime_error("Empty CSV: " + path);
    }
    return rows;
}

std::vector<int> read_index_row(const std::string &path, int expected_cols)
{
    Rows rows = read_csv(path, expected_cols);
    if (rows.size() != 1) {
        throw std::runtime_error("Index CSV must have exactly one row: " + path);
    }

    std::vector<int> out(rows[0].size());
    for (size_t i = 0; i < rows[0].size(); ++i) {
        out[i] = static_cast<int>(std::lround(rows[0][i]));
    }
    return out;
}

template<int R, int C>
void rows_to_array(const Rows &rows, float (&out)[R][C], const std::string &name)
{
    if (static_cast<int>(rows.size()) != R) {
        throw std::runtime_error("Bad row count for " + name);
    }
    for (int r = 0; r < R; ++r) {
        if (static_cast<int>(rows[r].size()) != C) {
            throw std::runtime_error("Bad column count for " + name);
        }
        for (int c = 0; c < C; ++c) {
            out[r][c] = rows[r][c];
        }
    }
}

template<int N>
void row_to_vector(const Rows &rows, float (&out)[N], const std::string &name)
{
    if (rows.size() != 1 || static_cast<int>(rows[0].size()) != N) {
        throw std::runtime_error("Bad vector shape for " + name);
    }
    for (int i = 0; i < N; ++i) {
        out[i] = rows[0][i];
    }
}

template<int N>
void copy_row(const std::vector<float> &row, float (&out)[N], const std::string &name)
{
    if (static_cast<int>(row.size()) != N) {
        throw std::runtime_error("Bad row width for " + name);
    }
    for (int i = 0; i < N; ++i) {
        out[i] = row[i];
    }
}

template<int N>
void append_array_row(Rows &rows, const float (&values)[N])
{
    std::vector<float> row(N);
    for (int i = 0; i < N; ++i) {
        row[i] = values[i];
    }
    rows.push_back(row);
}

void append_int_row(Rows &rows, int value)
{
    rows.push_back(std::vector<float>{static_cast<float>(value)});
}

void append_float_row(Rows &rows, float value)
{
    rows.push_back(std::vector<float>{value});
}

void write_csv(const std::string &path, const Rows &rows)
{
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open output CSV: " + path);
    }
    file << std::setprecision(17);
    for (const auto &row : rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i] << (i + 1 < row.size() ? "," : "\n");
        }
    }
}

template<int ControlSize, int EffectorSize>
struct MethodOutput {
    Rows u;
    Rows raw;
    Rows errout;
    Rows rho;
};

template<int ControlSize, int EffectorSize>
void run_dp_dir(DP_LP_ControlAllocator<ControlSize, EffectorSize> &allocator,
                const Rows &inputs,
                MethodOutput<ControlSize, EffectorSize> &out)
{
    for (const auto &row : inputs) {
        float input[ControlSize];
        float raw[EffectorSize];
        float u[EffectorSize];
        int err = 0;
        float rho = 0.0f;

        copy_row(row, input, "DP_LPCA input");
        allocator.DP_LPCA(input, raw, err, rho);
        allocator.restoring(raw, u);

        append_array_row(out.raw, raw);
        append_array_row(out.u, u);
        append_int_row(out.errout, err);
        append_float_row(out.rho, rho);
    }
}

template<int ControlSize, int EffectorSize>
void run_dp_scaled(DP_LP_ControlAllocator<ControlSize, EffectorSize> &allocator,
                   const Rows &inputs,
                   MethodOutput<ControlSize, EffectorSize> &out)
{
    for (const auto &row : inputs) {
        float input[ControlSize];
        float raw[EffectorSize];
        float u[EffectorSize];
        int err = 0;
        float rho = 0.0f;

        copy_row(row, input, "DPscaled_LPCA input");
        allocator.DPscaled_LPCA(input, raw, err, rho);
        allocator.restoring(raw, u);

        append_array_row(out.raw, raw);
        append_array_row(out.u, u);
        append_int_row(out.errout, err);
        append_float_row(out.rho, rho);
    }
}

Rows embed_split_output(const MethodOutput<3, 6> &servo_output,
                        const Rows &u_force,
                        const std::vector<int> &servo_cols_zero_based)
{
    if (u_force.size() != servo_output.u.size()) {
        throw std::runtime_error("u_force and split output length mismatch");
    }
    if (servo_cols_zero_based.size() != 6) {
        throw std::runtime_error("servo_cols_zero_based must have 6 entries");
    }

    Rows full;
    full.reserve(u_force.size());
    for (size_t sample = 0; sample < u_force.size(); ++sample) {
        if (u_force[sample].size() != 7 || servo_output.u[sample].size() != 6) {
            throw std::runtime_error("Bad split row width");
        }
        std::vector<float> row = u_force[sample];
        for (int i = 0; i < 6; ++i) {
            const int col = servo_cols_zero_based[i];
            if (col < 0 || col >= 7) {
                throw std::runtime_error("servo column index out of range");
            }
            row[col] += servo_output.u[sample][i];
        }
        full.push_back(row);
    }
    return full;
}

void write_method(const std::string &output_dir,
                  const std::string &name,
                  const Rows &u,
                  const Rows &raw,
                  const Rows &errout,
                  const Rows &rho)
{
    write_csv(join_path(output_dir, "output_cpp_shc09_log_" + name + ".csv"), u);
    write_csv(join_path(output_dir, "output_cpp_shc09_log_" + name + "_raw.csv"), raw);
    write_csv(join_path(output_dir, "output_cpp_shc09_log_" + name + "_errout.csv"), errout);
    write_csv(join_path(output_dir, "output_cpp_shc09_log_" + name + "_rho.csv"), rho);
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const std::string project_root = argc > 1 ? argv[1] : "/Users/mch/Proj/control_allocation";
        const std::string input_dir = join_path(project_root, "results/cpp_inputs");
        const std::string output_root = join_path(project_root, "results");
        const std::string output_dir = join_path(project_root, "results/cpp_outputs");
        ensure_dir(output_root);
        ensure_dir(output_dir);

        const Rows bpar_rows = read_csv(join_path(input_dir, "shc09_log_bpar.csv"), 7);
        const Rows vpar_rows = read_csv(join_path(input_dir, "shc09_log_v_par.csv"), 4);
        const Rows umin_rows = read_csv(join_path(input_dir, "shc09_log_umin.csv"), 7);
        const Rows umax_rows = read_csv(join_path(input_dir, "shc09_log_umax.csv"), 7);

        const Rows btorque_rows = read_csv(join_path(input_dir, "shc09_log_b_torque_servo.csv"), 6);
        const Rows vtorque_rows = read_csv(join_path(input_dir, "shc09_log_v_torque_left.csv"), 3);
        const Rows umin_servo_rows = read_csv(join_path(input_dir, "shc09_log_umin_servo.csv"), 6);
        const Rows umax_servo_rows = read_csv(join_path(input_dir, "shc09_log_umax_servo.csv"), 6);
        const Rows u_force_rows = read_csv(join_path(input_dir, "shc09_log_u_force.csv"), 7);
        const std::vector<int> servo_cols_zero_based =
            read_index_row(join_path(input_dir, "shc09_log_servo_cols_zero_based.csv"), 6);

        if (vpar_rows.size() != vtorque_rows.size() || vpar_rows.size() != u_force_rows.size()) {
            throw std::runtime_error("Input sample counts differ");
        }

        float Bpar[4][7];
        float BtorqueServo[3][6];
        float umin[7];
        float umax[7];
        float uminServo[6];
        float umaxServo[6];
        rows_to_array(bpar_rows, Bpar, "Bpar");
        rows_to_array(btorque_rows, BtorqueServo, "BtorqueServo");
        row_to_vector(umin_rows, umin, "umin");
        row_to_vector(umax_rows, umax, "umax");
        row_to_vector(umin_servo_rows, uminServo, "uminServo");
        row_to_vector(umax_servo_rows, umaxServo, "umaxServo");

        Aircraft<4, 7> aircraftBpar(Bpar, umax, umin);
        Aircraft<3, 6> aircraftSplit(BtorqueServo, umaxServo, uminServo);
        DP_LP_ControlAllocator<4, 7> allocatorBpar(aircraftBpar);
        DP_LP_ControlAllocator<3, 6> allocatorSplit(aircraftSplit);

        MethodOutput<4, 7> bparDir;
        MethodOutput<4, 7> bparScaled;
        MethodOutput<3, 6> splitDirServo;
        MethodOutput<3, 6> splitScaledServo;

        run_dp_dir(allocatorBpar, vpar_rows, bparDir);
        run_dp_scaled(allocatorBpar, vpar_rows, bparScaled);
        run_dp_dir(allocatorSplit, vtorque_rows, splitDirServo);
        run_dp_scaled(allocatorSplit, vtorque_rows, splitScaledServo);

        write_method(output_dir, "pca_dir_bpar", bparDir.u, bparDir.raw, bparDir.errout, bparDir.rho);
        write_method(output_dir, "pca_dpscaled_bpar", bparScaled.u, bparScaled.raw, bparScaled.errout, bparScaled.rho);
        write_method(output_dir, "split_pca_dir",
                     embed_split_output(splitDirServo, u_force_rows, servo_cols_zero_based),
                     splitDirServo.raw, splitDirServo.errout, splitDirServo.rho);
        write_method(output_dir, "split_pca_dpscaled",
                     embed_split_output(splitScaledServo, u_force_rows, servo_cols_zero_based),
                     splitScaledServo.raw, splitScaledServo.errout, splitScaledServo.rho);

        std::cout << "SHC09 log replay samples: " << vpar_rows.size() << std::endl;
        std::cout << "Output dir: " << output_dir << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "shc09_log_replay failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
