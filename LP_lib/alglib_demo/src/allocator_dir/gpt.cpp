#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> Vector;
typedef vector<Vector> Matrix;
// Certainly! Below is a simple implementation of the simplex algorithm in C++ to solve a standard linear programming (LP) problem. This implementation assumes the LP problem is given in standard form:

// Minimize: c^T * x
// Subject to: Ax = b, x >= 0
class Simplex {
private:
    Matrix tableau;
    int m, n; // Number of constraints and number of variables

public:
    Simplex(const Matrix &A, const Vector &b, const Vector &c) {
        m = A.size();
        n = A[0].size() + 1;

        tableau = Matrix(m + 1, Vector(n + m + 1, 0));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n - 1; ++j) {
                tableau[i][j] = A[i][j];
            }
            tableau[i][n - 1] = 1;
            tableau[i][n + i] = 1;
            tableau[i][n + m] = b[i];
        }

        for (int j = 0; j < n - 1; ++j) {
            tableau[m][j] = c[j] * -1;
        }
    }

    void pivot(int row, int column) {
        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= n + m; ++j) {
                if (i != row && j != column) {
                    tableau[i][j] -= tableau[row][j] * tableau[i][column] / tableau[row][column];
                }
            }
        }
        for (int i = 0; i <= m; ++i) {
            if (i != row) tableau[i][column] = 0;
        }
        for (int j = 0; j <= n + m; ++j) {
            if (j != column) tableau[row][j] /= tableau[row][column];
        }
        tableau[row][column] = 1;
    }

    bool simplex(int &col, int &row) {
        col = -1;
        row = -1;
        double min_val = 0;
        for (int j = 0; j < n + m; ++j) {
            if (tableau[m][j] < min_val) {
                min_val = tableau[m][j];
                col = j;
            }
        }
        if (col < 0) return false;

        double min_ratio = INFINITY;
        for (int i = 0; i < m; ++i) {
            if (tableau[i][col] > 0) {
                double ratio = tableau[i][n + m] / tableau[i][col];
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    row = i;
                }
            }
        }
        if (row < 0) return false;

        pivot(row, col);
        return true;
    }

    void solve() {
        while (true) {
            int col, row;
            if (!simplex(col, row)) break;
        }
    }

    double getOptimalValue() {
        return -tableau[m][n + m];
    }

    Vector getOptimalSolution() {
        Vector solution(n - 1);
        for (int i = 0; i < n - 1; ++i) {
            bool is_basic = false;
            int basic_index;
            for (int j = 0; j < m; ++j) {
                if (abs(tableau[j][i] - 1) < 1e-6 && !is_basic) {
                    is_basic = true;
                    basic_index = j;
                } else if (abs(tableau[j][i]) > 1e-6) {
                    is_basic = false;
                    break;
                }
            }
            if (is_basic) {
                solution[i] = tableau[basic_index][n + m];
            } else {
                solution[i] = 0;
            }
        }
        return solution;
    }
};

int main() {
    Matrix A = {{1, 1, 3}, {2, 2, 5}, {4, 1, 2}};
    Vector b = {30, 24, 36};
    Vector c = {3, 1, 2};

    Simplex simplex(A, b, c);
    cout << "Optimal Value: " << simplex.getOptimalValue() << endl;
    auto start = std::chrono::high_resolution_clock::now();
    simplex.solve();
    Vector solution = simplex.getOptimalSolution();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Execution time: " << elapsed.count() << "s\n";
    cout << "Optimal Solution: ";
    for (double val : solution) {
        cout << val << " ";
    }
    cout << endl;

    return 0;
}
