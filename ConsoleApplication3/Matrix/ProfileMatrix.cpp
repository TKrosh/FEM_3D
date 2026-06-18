// ProfileMatrix.cpp
#include "ProfileMatrix.h"

void ProfileMatrix::ReservedMemory_SpareStorage(std::vector<std::vector<int>> rowCols)
{
    std::cout << "we already preperae memory? how are you even get here?" << std::endl;
}

ProfileMatrix::ProfileMatrix(int n,
    std::vector<int> indexStart,
    std::map<int, double> dirichletValues,
    std::vector<int> freeIndex)
{
    DirichletValues = dirichletValues;
    FreeIndex = freeIndex;
    size = n;
    di.resize(n, 0.0);
    rightPart.resize(n, 0.0);
    ia.resize(n + 1, 0);

    for (int i = 1; i <= n; ++i)
    {
        ia[i] = ia[i - 1] + (i - 1 - indexStart[i - 1]);
    }
    al.resize(ia[n], 0.0);
    au.resize(ia[n], 0.0);
}

void ProfileMatrix::multiplicationByVector(std::vector<double>* v)
{
    std::vector<double> tmp = *v;
    for (int i = 0; i < size; ++i)
    {
        double sum = di[i] * tmp[i];

        // нижний треугольник 
        int start = ia[i];
        int end = ia[i + 1];
        int lineStart = i - (end - start);
        for (int pos = start; pos < end; ++pos)
        {
            int j = lineStart + (pos - start);
            sum += al[pos] * tmp[j];
        }

        // верхний треугольник
        for (int j = i + 1; j < size; ++j)
        {
            int js = ia[j];
            int je = ia[j + 1];
            int jLineStart = j - (je - js);
            if (i >= jLineStart)
            {
                int offset = i - jLineStart;
                sum += au[js + offset] * tmp[j];
            }
        }
        (*v)[i] = sum;
    }
}
void ProfileMatrix::multiplicationByVector_T(std::vector<double>* v)
{
    std::vector<double> tmp = *v;
    for (int i = 0; i < size; i++)
    {
        double sum = di[i] * tmp[i];

        // верхний треугольник 
        int start = ia[i];
        int end = ia[i + 1];
        int lineStart = i - (end - start);
        for (int pos = start; pos < end; pos++)
        {
            int j = lineStart + (pos - start);
            sum += au[pos] * tmp[j];
        }

        // нижний треугольник
        for (int j = i + 1; j < size; ++j)
        {
            int js = ia[j];
            int je = ia[j + 1];
            int jLineStart = j - (je - js);
            if (i >= jLineStart)
            {
                int offset = i - jLineStart;
                sum += al[js + offset] * tmp[j];
            }
        }
        (*v)[i] = sum;
    }
}

void ProfileMatrix::matrixClear()
{
    for (int i = 0; i < size; ++i)
    {
        di[i] = 0.0;
        rightPart[i] = 0.0;
    }
    for (int i = 0; i < ia[size]; ++i)
    {
        al[i] = 0.0;
        au[i] = 0.0;
    }
}

void ProfileMatrix::test()
{
    std::cout << "Global Matrix (lower triangle + diag):\n";
    for (int i = 0; i < size; ++i)
    {
        int count = ia[i + 1] - ia[i];
        int lineStart = i - count;
        for (int j = 0; j < i; ++j)
        {
            if (j >= lineStart)
                std::cout << al[ia[i] + (j - lineStart)] << " ";
            else
                std::cout << 0 << " ";
        }
        std::cout << di[i] << "\n";
    }

    std::cout << "\nGlobal Matrix (upper triangle):\n";
    for (int i = 0; i < size; ++i)
    {
        for (int j = i + 1; j < size; ++j)
        {
            int count = ia[j + 1] - ia[j];
            int lineStart = j - count;
            if (i >= lineStart)
                std::cout << au[ia[j] + (i - lineStart)] << " ";
            else
                std::cout << 0;
        }
        std::cout << "\n";
    }

    std::cout << "\nRight part:\n";
    for (int i = 0; i < size; ++i)
        std::cout << rightPart[i] << " ";
    std::cout << "\n";
}

void ProfileMatrix::AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix)
{
    for (int i = 0; i < elem.n; ++i)
    {
        int global_i = elem.CoordsIndexes[i];
        int free_i = FreeIndex[global_i];
        if (free_i != -1) di[free_i] += LocalMatrix[i][i];

        for (int j = 0; j < i; ++j)
        {
            int global_j = elem.CoordsIndexes[j];
            int free_j = FreeIndex[global_j];
            double val = LocalMatrix[i][j];

            if (free_i != -1 && free_j != -1)
            {
                if (free_i > free_j)
                {
                    int count = ia[free_i + 1] - ia[free_i];
                    int startIdx = free_i - count;
                    int offset = free_j - startIdx;
                    al[ia[free_i] + offset] += val;
                }
                else
                {
                    int count = ia[free_j + 1] - ia[free_j];
                    int startIdx = free_j - count;
                    int offset = free_i - startIdx;
                    al[ia[free_j] + offset] += val;
                }
            }
            else if (free_i != -1 && free_j == -1)
            {
                rightPart[free_i] -= val * DirichletValues[global_j];
            }
            else if (free_i == -1 && free_j != -1)
            {
                rightPart[free_j] -= val * DirichletValues[global_i];
            }
        }
    }
}

void ProfileMatrix::AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v)
{
    for (int i = 0; i < elem.n; ++i)
    {
        int global_i = elem.CoordsIndexes[i];
        int free_i = FreeIndex[global_i];
        if (free_i != -1) rightPart[free_i] += v[i];
    }
}

void ProfileMatrix::AddLocalMatrix_Garmonic(TreeLinearLagrange& elem, double** LocalMatrix)
{
    for (int i = 0; i < elem.n; ++i)
    {
        int global_i = elem.CoordsIndexes[i];
        int s_orig_i = 2 * global_i;
        int c_orig_i = s_orig_i + 1;
        int s_free_i = FreeIndex[s_orig_i];
        int c_free_i = FreeIndex[c_orig_i];

        for (int j = 0; j < elem.n; ++j)
        {
            int global_j = elem.CoordsIndexes[j];
            int s_orig_j = 2 * global_j;
            int c_orig_j = s_orig_j + 1;
            int s_free_j = FreeIndex[s_orig_j];
            int c_free_j = FreeIndex[c_orig_j];

            // Вспомогательная лямбда для добавления в профиль
            auto add_to_profile = [&](int row, int col, double val) {
                if (row == col) {
                    di[row] += val;
                }
                else if (row > col) {
                    // нижний треугольник (al)
                    int count = ia[row + 1] - ia[row];
                    int startIdx = row - count;
                    int offset = col - startIdx;
                    if (offset >= 0 && offset < count) {
                        al[ia[row] + offset] += val;
                    }
                    // else: связь не предусмотрена профилем – игнорируем (не должно случаться)
                }
                else { // row < col
                    // верхний треугольник (au), хранимый по столбцам
                    int count = ia[col + 1] - ia[col];
                    int startIdx = col - count;
                    int offset = row - startIdx;
                    if (offset >= 0 && offset < count) {
                        au[ia[col] + offset] += val;
                    }
                }
                };

            // Блок ss (p_ij)
            double p_ij = LocalMatrix[2 * i][2 * j];
            if (s_free_i != -1 && s_free_j != -1)
                add_to_profile(s_free_i, s_free_j, p_ij);
            else if (s_free_i != -1 && s_free_j == -1)
                rightPart[s_free_i] -= p_ij * DirichletValues[s_orig_j];
            //else if (s_free_i == -1 && s_free_j != -1)
            //    rightPart[s_free_j] -= p_ij * DirichletValues[s_orig_i];

            // Блок sc (-c_ij)
            double mc_ij = LocalMatrix[2 * i][2 * j + 1];
            if (s_free_i != -1 && c_free_j != -1)
                add_to_profile(s_free_i, c_free_j, mc_ij);
            else if (s_free_i != -1 && c_free_j == -1)
                rightPart[s_free_i] -= mc_ij * DirichletValues[c_orig_j];
            //else if (s_free_i == -1 && c_free_j != -1)
            //    rightPart[c_free_j] -= mc_ij * DirichletValues[s_orig_i];

            // Блок cs (+c_ij)
            double pc_ij = LocalMatrix[2 * i + 1][2 * j];
            if (c_free_i != -1 && s_free_j != -1)
                add_to_profile(c_free_i, s_free_j, pc_ij);
            else if (c_free_i != -1 && s_free_j == -1)
                rightPart[c_free_i] -= pc_ij * DirichletValues[s_orig_j];
            //else if (c_free_i == -1 && s_free_j != -1)
            //    rightPart[s_free_j] -= pc_ij * DirichletValues[c_orig_i];

            // Блок cc (p_ij)
            double p_ij2 = LocalMatrix[2 * i + 1][2 * j + 1];
            if (c_free_i != -1 && c_free_j != -1)
                add_to_profile(c_free_i, c_free_j, p_ij2);
            else if (c_free_i != -1 && c_free_j == -1)
                rightPart[c_free_i] -= p_ij2 * DirichletValues[c_orig_j];
            //else if (c_free_i == -1 && c_free_j != -1)
            //    rightPart[c_free_j] -= p_ij2 * DirichletValues[c_orig_i];
        }
    }
}

void ProfileMatrix::AddToRightPart_Garmonic(TreeLinearLagrange& elem, std::vector<double> v)
{
    for (int i = 0; i < elem.n; ++i)
    {
        int global_i = elem.CoordsIndexes[i];
        int s_free_i = FreeIndex[2 * global_i];
        if (s_free_i != -1) rightPart[s_free_i] += v[2 * i];

        int c_free_i = FreeIndex[2 * global_i + 1];
        if (c_free_i != -1) rightPart[c_free_i] += v[2 * i + 1];
    }
}

void ProfileMatrix::GetDiag(std::vector<double>* d)
{
    for (int i = 0; i < size; ++i)
        (*d)[i] = di[i];
}

void ProfileMatrix::LU_decompose()
{
    for (int i = 0; i < size; i++)
    {
        int i0 = ia[i];
        int i1 = ia[i + 1];
        int len = i1 - i0;
        int firstCol = i - len;
        double sumD = 0;

        int j = firstCol;
        for (int ij = i0; ij < i1; ij++, j++)
        {
            double sumL = 0;
            double sumU = 0;

            int j0 = ia[j];
            int j1 = ia[j + 1];
            int ki = i0;
            int kj = j0;
            int kol_umn_i = ij - i0;
            int kol_umn_j = j1 - j0;
            int kol_umn_r = kol_umn_i - kol_umn_j;

            if (kol_umn_r > 0)
                ki += kol_umn_r;
            else
                kj -= kol_umn_r;

            for (; ki < ij; ki++, kj++)
            {
                sumL += al[ki] * au[kj];
                sumU += au[ki] * al[kj];
            }

            al[ij] = (al[ij] - sumL) / di[j];
            au[ij] = au[ij] - sumU;
            sumD += al[ij] * au[ij];
        }
        di[i] -= sumD;
    }
}

void ProfileMatrix::Solve_L(std::vector<double>& y) {
    for (int i = 0; i < size; ++i) {
        double sum = 0.0;
        int start = i - (ia[i + 1] - ia[i]);
        for (int p = ia[i]; p < ia[i + 1]; ++p) {
            int j = start + (p - ia[i]);
            sum += al[p] * y[j];
        }
        y[i] = rightPart[i] - sum;
    }
}

void ProfileMatrix::Solve_U(std::vector<double>& x) {
    for (int i = size - 1; i >= 0; --i) {
        x[i] /= di[i];

        int start = i - (ia[i + 1] - ia[i]);
        for (int p = ia[i]; p < ia[i + 1]; ++p) {
            int j = start + (p - ia[i]);
            x[j] -= au[p] * x[i];
        }
    }
}