//FATHAN IRFANI
//21120122120009
//METODE NUMERIK
//KELAS D

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Mencetak matriks
void printMatrix(vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
}

// Mencetak vektor
void printVector(vector<double>& vec) {
    for (double val : vec) {
        cout << val << "\t";
    }
    cout << endl;
}

// Alokasi matriks
vector<vector<double>> allocateMatrix(int rows, int cols) {
    return vector<vector<double>>(rows, vector<double>(cols, 0));
}

// Matriks balikan
vector<double> inverseMatrixMethod(vector<vector<double>>& A, vector<double>& b) {
    // Cek apakah matriks memiliki balikan
    if (A.size() != A[0].size()) {
        cerr << "Matriks A bukan matriks persegi, tidak memiliki balikan." << endl;
        return {};
    }

    int n = A.size();

    // Hitung determinan
    double det = 1.0;
    for (int i = 0; i < n; ++i) {
        det *= A[i][i];
    }

    if (det == 0) {
        cerr << "Determinan nol, matriks tidak memiliki balikan." << endl;
        return {};
    }

    // Hitung matriks balikan
    vector<vector<double>> A_inv = allocateMatrix(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Matriks balikan adalah adjoin(A) / det(A)
            // Di sini diasumsikan det(A) tidak nol
            A_inv[j][i] = pow(-1, i + j) * A[i][j] / det;
        }
    }

    // Perkalian matriks balikan dengan vektor B
    vector<double> x(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            x[i] += A_inv[i][j] * b[j];
        }
    }

    return x;
}

// Dekomposisi LU Gauss
vector<double> luGaussMethod(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();

    // Deklarasi matriks LU
    vector<vector<double>> L = allocateMatrix(n, n);
    vector<vector<double>> U = allocateMatrix(n, n);

    // Inisialisasi matriks L dan U
    for (int i = 0; i < n; ++i) {
        L[i][i] = 1; // Diagonal utama matriks L diisi dengan 1
    }

    // Proses dekomposisi LU
    for (int k = 0; k < n; ++k) {
        U[k][k] = A[k][k];
        for (int i = k + 1; i < n; ++i) {
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }
        for (int i = k + 1; i < n; ++i) {
            for (int j = k + 1; j < n; ++j) {
                A[i][j] = A[i][j] - L[i][k] * U[k][j];
            }
        }
    }

    // Proses substitusi maju
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // Proses substitusi mundur
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

// Dekomposisi Crout
vector<double> croutMethod(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();

    // Deklarasi matriks L dan U
    vector<vector<double>> L = allocateMatrix(n, n);
    vector<vector<double>> U = allocateMatrix(n, n);

    // Proses dekomposisi Crout
    for (int i = 0; i < n; ++i) {
        // Hitung matriks U
        for (int k = i; k < n; ++k) {
            double sum = 0;
            for (int j = 0; j < i; ++j) {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = A[i][k] - sum;
        }

        // Hitung matriks L
        for (int k = i; k < n; ++k) {
            if (i == k) {
                L[i][i] = 1; // Diagonal utama matriks L diisi dengan 1
            } else {
                double sum = 0;
                for (int j = 0; j < i; ++j) {
                    sum += L[k][j] * U[j][i];
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }

    // Proses substitusi maju
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // Proses substitusi mundur
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

    int main() {
    int n;
    cout <<"================="<<endl;
    cout <<"FATHAN IRFANI"<<endl;
    cout <<"21120122120009"<<endl;
    cout <<"METODE NUMERIK D"<<endl;
    cout <<"================="<<endl<<endl;
    cout << "Masukkan ukuran matriks (n x n): ";
    cin >> n;

    // Input matriks A
    cout << endl << "Masukkan matriks A (" << n << " x " << n << "):" << endl;
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        cout << "Baris ke-" << i + 1 << ": ";
        for (int j = 0; j < n; ++j) {
            cin >> A[i][j];
        }
    }

    // Input vektor B
    cout << endl<< "Masukkan vektor B (" << n << "):" << endl;
    vector<double> b(n);
    cout << "Elemen vektor B: ";
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }

    // Pilih metode
    int method;
    cout << endl<<"==========================================================="<<endl;
    cout << "Pilih metode untuk mencari solusi sistem persamaan linear:" << endl;
    cout << "1. Metode matriks balikan" << endl;
    cout << "2. Metode dekomposisi LU Gauss" << endl;
    cout << "3. Metode dekomposisi Crout" << endl;
    cout << "==========================================================="<<endl;
    cout << "Masukkan pilihan (1/2/3): ";
    cin >> method;


    vector<double> x;
    switch (method) {
        case 1:
            x = inverseMatrixMethod(A, b);
            break;
        case 2:
            x = luGaussMethod(A, b);
            break;
        case 3:
            x = croutMethod(A, b);
            break;
        default:
            cerr << "Pilihan tidak valid." << endl;
            return 1;
    }

    // Output solusi
    cout <<endl<< "Solusi sistem persamaan linear:" << endl;
    printVector(x);

    return 0;
}
