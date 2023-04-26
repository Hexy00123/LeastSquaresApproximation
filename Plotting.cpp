#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

class Matrix {
public:
    int rows;
    int cols;
    vector<vector<double>> data;

    // Constructor that initializes the matrix with zeros
    Matrix(int rows, int cols) : rows(rows), cols(cols), data(rows, vector<double>(cols, 0.0)) {}

    // Constructor that initializes the matrix with given data
    Matrix(vector<vector<double>> data) : rows(data.size()), cols(data[0].size()), data(data) {}

    // Overloading the + operator to add two matrices
    Matrix operator+(const Matrix &other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("Matrices must have the same dimensions to be added");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    // Overloading the - operator to subtract two matrices
    Matrix operator-(const Matrix &other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("Matrices must have the same dimensions to be subtracted");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    // Overloading the = operator to assign one matrix to another
    Matrix &operator=(const Matrix &other) {
        if (this != &other) {
            rows = other.rows;
            cols = other.cols;
            data = other.data;
        }
        return *this;
    }

    // Overloading the * operator to multiply two matrices
    Matrix operator*(const Matrix &other) const {
        if (cols != other.rows) {
            throw invalid_argument(
                    "The number of columns in the first matrix must match the number of rows in the second matrix to be multiplied");
        }
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    // Overloading the [] operator to access elements of the matrix
    vector<double> &operator[](int index) {
        if (index < 0 || index >= rows) {
            throw out_of_range("Index out of range");
        }
        return data[index];
    }

    // Overloading the >> operator to input matrix data from a stream
    friend istream &operator>>(istream &in, Matrix &matrix) {
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.cols; j++) {
                in >> matrix.data[i][j];
            }
        }
        return in;
    }

    // Overloading the << operator to output matrix data to a stream
    friend ostream &operator<<(ostream &out, const Matrix &matrix) {
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.cols; j++) {
                if (abs(matrix.data[i][j]) <= 1e-10)
                    out << "0.0000" << " ";
                else
                    out << matrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }

    // Transposes the matrix
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int size) : Matrix(size, size) {}

    SquareMatrix(const Matrix &matrix) : Matrix(matrix.rows, matrix.cols) {
        if (matrix.rows != matrix.cols) {
            throw invalid_argument("Matrix must be square to cast to SquareMatrix");
        }
        data = matrix.data;
    }

    // Returns the transpose of the square matrix
    SquareMatrix transpose() const {
        SquareMatrix result(rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[i][j] = data[j][i];
            }
        }
        return result;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int size) : SquareMatrix(size) {
        for (int i = 0; i < size; i++) {
            data[i][i] = 1.0;
        }
    }

    IdentityMatrix(const SquareMatrix &matrix) : SquareMatrix(matrix.rows) {
        if (matrix.rows != matrix.cols) {
            throw runtime_error("Cannot create identity matrix from non-square matrix");
        }
        for (int i = 0; i < rows; i++) {
            data[i][i] = 1.0;
        }
    }
};

class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int size, int row1, int row2) : SquareMatrix(size) {
        if (row1 >= size || row2 >= size || row1 < 0 || row2 < 0) {
            throw runtime_error("Invalid row indices");
        }
        for (int i = 0; i < size; i++) {
            data[i][i] = 1.0;
        }
        (*this)[row1][row1] = 0;
        (*this)[row2][row2] = 0;
        (*this)[row2][row1] = 1;
        (*this)[row1][row2] = 1;
    }

    PermutationMatrix(const SquareMatrix &matrix, int row1, int row2) : SquareMatrix(matrix.rows) {
        if (matrix.rows != matrix.cols) {
            throw runtime_error("Cannot create permutation matrix from non-square matrix");
        }
        if (row1 >= rows || row2 >= rows || row1 < 0 || row2 < 0) {
            throw runtime_error("Invalid row indices");
        }
        for (int i = 0; i < rows; i++) {
            data[i][i] = 1.0;
        }
        swap(data[row1], data[row2]);
    }
};

class EliminationMatrix : public SquareMatrix {
public:
    EliminationMatrix(const SquareMatrix &matrix, int row, int col) : SquareMatrix(matrix.rows) {
        if (matrix.rows != matrix.cols) {
            throw runtime_error("Cannot create elimination matrix from non-square matrix");
        }
        if (row >= rows || col >= cols || row < 0 || col < 0) {
            throw runtime_error("Invalid row or column index");
        }
        for (int i = 0; i < rows; i++) {
            data[i][i] = 1.0;
        }
        double factor = matrix.data[row][col] / matrix.data[col][col];
        data[row][col] = -factor;
    }
};

class ColumnVector : public Matrix {
public:
    // Constructor that initializes the column vector with zeros
    ColumnVector(int rows) : Matrix(rows, 1) {}

    // Constructor that initializes the column vector with given data
    ColumnVector(vector<double> data) : Matrix(data.size(), 1) {
        for (int i = 0; i < data.size(); i++) {
            this->data[i][0] = data[i];
        }
    }

    // Constructor that casts a Matrix to a ColumnVector
    ColumnVector(const Matrix &matrix) : Matrix(matrix.rows, 1) {
        if (matrix.cols != 1) {
            throw invalid_argument("Matrix must have only one column to be cast to a ColumnVector");
        }
        for (int i = 0; i < rows; i++) {
            this->data[i][0] = matrix.data[i][0];
        }
    }

    // Overloading the [] operator to access elements of the column vector
    double &operator[](int index) {
        if (index < 0 || index >= rows) {
            throw out_of_range("Index out of range");
        }
        return data[index][0];
    }
};

int findMaxPivot(SquareMatrix &m, int currentRow) {
    int maxRow = currentRow;
    for (int j = maxRow; j < m.rows; ++j) {
        if (abs(m[maxRow][currentRow]) < abs(m[j][currentRow])) {
            maxRow = j;
        }
    }
    return maxRow;
}

void eliminateSystem(SquareMatrix &m, SquareMatrix &inverse, ColumnVector &b) {
    int n = m.rows;

    // direct way
    for (int i = 0; i < n; ++i) {
        // permutation
        int maxRow = findMaxPivot(m, i);
        if (maxRow != i) {
            PermutationMatrix p(m, i, maxRow);
            m = p * m;
            inverse = p * inverse;
            b = p * b;
        }

        //elimination
        for (int row = i + 1; row < n; ++row) {
            EliminationMatrix e(m, row, i);
            if (abs(e[row][i]) > 1e-10) {
                m = e * m;
                inverse = e * inverse;
                b = e * b;
            }
        }
    }
    // way back
    for (int i = n - 1; i >= 0; --i) {
        for (int row = i - 1; row >= 0; --row) {
            EliminationMatrix e(m, row, i);
            if (abs(e[row][i]) > 1e-10) {
                m = e * m;
                inverse = e * inverse;
                b = e * b;
            }
        }
    }

    // diagonal normalization
    IdentityMatrix dn(n);
    for (int i = 0; i < n; ++i) {
        if (abs(m[i][i]) > 1e-10) {
            dn[i][i] = 1.0 / m[i][i];
        }
    }
    m = dn * m;
    inverse = dn * inverse;
    b = dn * b;
}

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

int main() {
#ifdef WIN32
    FILE *pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif

    cout << fixed;
    cout << setprecision(4);

    int datasetLength;
    cin >> datasetLength;

    Matrix data = Matrix(datasetLength, 2);
    cin >> data;

    int degree;
    cin >> degree;

    Matrix A = Matrix(datasetLength, degree + 1);
    ColumnVector b = ColumnVector(datasetLength);

    for (int i = 0; i < datasetLength; ++i) {
        for (int d = 0; d <= degree; ++d)
            A[i][d] = pow(data[i][0], d);

        b[i] = data[i][1];
    }
    cout << "A:\n";
    cout << A;

    SquareMatrix AtA = A.transpose() * A;
    cout << "A_T*A:\n";
    cout << AtA;

    ColumnVector bT = A.transpose() * b;
    SquareMatrix AtAinv = IdentityMatrix(AtA);
    eliminateSystem(AtA, AtAinv, bT);
    cout << "(A_T*A)^-1:\n";
    cout << AtAinv;

    b = A.transpose() * b;
    cout << "A_T*b:\n";
    cout << b;

    cout << "x~:\n";
    cout << bT;

    fprintf(pipe, "plot [-20 : 20] [-20 : 20] %lf*x**3 %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n",
            bT[3], bT[2], bT[1], bT[0]);

    for (double i = 0; i < datasetLength; ++i) {
        fprintf(pipe, "%f\t%f\n", data[i][0], data[i][1]);
    }

    fprintf(pipe, "e\n");
    fflush(pipe);

#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
    return 0;
}
