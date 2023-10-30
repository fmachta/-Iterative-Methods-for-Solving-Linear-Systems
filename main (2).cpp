#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

bool is_diagonally_dominant(const vector<vector<double>>& A) {
    for (int i = 0; i < A.size(); i++) {
        double sum = 0;
        for (int j = 0; j < A[i].size(); j++) {
            if (i != j)
                sum += abs(A[i][j]);
        }
        if (abs(A[i][i]) <= sum)
            return false;
    }
    return true;
}

vector<double> jacobi(const vector<vector<double>>& A, const vector<double>& b, vector<double> x0, double error, int max_iter) {
    int n = A.size();
    vector<double> x(n, 0.0);

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            x[i] = b[i];
            for (int j = 0; j < n; j++) {
                if (i != j)
                    x[i] -= A[i][j] * x0[j];
            }
            x[i] /= A[i][i];
        }

        double norm = 0;
        for (int i = 0; i < n; i++) {
            norm += (x[i] - x0[i]) * (x[i] - x0[i]);
        }

        if (sqrt(norm) < error) {
            return x;
        }

        x0 = x;
    }

    return x;
}

vector<double> gauss_seidel(const vector<vector<double>>& A, const vector<double>& b, vector<double> x0, double error, int max_iter) {
    int n = A.size();
    vector<double> x(n, 0.0);

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            x[i] = b[i];
            for (int j = 0; j < n; j++) {
                if (i != j)
                    x[i] -= A[i][j] * ((j < i) ? x[j] : x0[j]);
            }
            x[i] /= A[i][i];
        }

        double norm = 0;
        for (int i = 0; i < n; i++) {
            norm += (x[i] - x0[i]) * (x[i] - x0[i]);
        }

        if (sqrt(norm) < error) {
            return x;
        }

        x0 = x;
    }

    return x;
}

int main() {
    int n;
    cout << "Enter the number of equations (n <= 10): ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    cout << "Enter coefficients through:\n1. Command Line\n2. File\nChoose (1/2): ";
    int choice;
    cin >> choice;

    if (choice == 1) {
        cout << "Enter the augmented coefficient matrix (including b values):\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cin >> A[i][j];
            }
            cin >> b[i];
        }
    }
    else {
        string filename;
        cout << "Enter the filename: ";
        cin >> filename;

        ifstream file(filename);
        for (int i = 0; i < n && file; i++) {
            for (int j = 0; j < n && file; j++) {
                file >> A[i][j];
            }
            file >> b[i];
        }
    }

    if (!is_diagonally_dominant(A)) {
        cout << "Matrix is not strictly diagonally dominant. Methods may not converge.\n";
    }

    double error;
    cout << "Enter the desired stopping error: ";
    cin >> error;

    vector<double> x0(n, 0.0);
    cout << "Enter the starting solution:\n";
    for (int i = 0; i < n; i++) {
        cin >> x0[i];
    }

    // Jacobi method
    vector<double> x_jacobi = jacobi(A, b, x0, error, 50);

    // Gauss-Seidel method
    vector<double> x_gauss_seidel = gauss_seidel(A, b, x0, error, 50);

    cout << "Solution using Jacobi: [";
    for (int i = 0; i < n; i++) {
        cout << fixed << setprecision(4) << x_jacobi[i] << " ";
    }
    cout << "]T\n";

    cout << "Solution using Gauss-Seidel: [";
    for (int i = 0; i < n; i++) {
        cout << fixed << setprecision(4) << x_gauss_seidel[i] << " ";
    }
    cout << "]T\n";

    return 0;
}