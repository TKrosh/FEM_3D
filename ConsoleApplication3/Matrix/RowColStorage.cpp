#include "RowColStorage.h"

RowColStorage::RowColStorage(int n,
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
}

void RowColStorage::ReservedMemory_SpareStorage(std::vector<std::vector<int>> rowCols)
{

    // Строим ia
    ia[0] = 0;
    for (int i = 0; i < size; i++)
        ia[i + 1] = ia[i] + rowCols[i].size();

    // Строим ja 
    ja.resize(ia[size], 0);
    for (int i = 0; i < size; i++)
        std::copy(rowCols[i].begin(), rowCols[i].end(), ja.begin() + ia[i]);

    al.resize(ia[size]);
    au.resize(ia[size]);
}

void RowColStorage::multiplicationByVector(std::vector<double>* v)
{
    std::vector<double> tmp = *v;
    for (int i = 0; i < size; i++)
    {
        double sum = 0.0;
        for (int k = ia[i]; k < ia[i + 1]; k++)
        {
            int j = ja[k]; // is a index in vector x that we need to multiply on
            if (j < i)
                sum += al[k] * tmp[j];
            else if (j > i)
                sum += au[k] * tmp[j];
        }
        (*v)[i] = tmp[i] * di[i] + sum;
    }
}

void RowColStorage::multiplicationByVector_T(std::vector<double>* v)
{
    std::vector<double> tmp = *v;
    std::vector<double> res(size, 0.0);
    for (int i = 0; i < size; ++i) {
        res[i] += di[i] * tmp[i];
        for (int k = ia[i]; k < ia[i + 1]; ++k) {
            int j = ja[k];
            if (j < i)
                res[j] += al[k] * tmp[i];
            else if (j > i)
                res[j] += au[k] * tmp[i];    
        }
    }
    *v = res;
}

void RowColStorage::matrixClear()
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

void RowColStorage::test()
{
    std::cout << "Global Matrix (size " << size << "):\n";
    for (int i = 0; i < size; ++i)
    {
        int row_start = ia[i];
        int row_end = ia[i + 1];
        int curr = row_start;  // текущая позиция в строке i

        for (int j = 0; j < size; ++j)
        {
            if (j == i)
            {
                std::cout << di[i] << " ";
            }
            else
            {
                // Ищем столбец j в строке i.
                // ja в строке упорядочен по возрастанию – двигаем указатель curr.
                while (curr < row_end && ja[curr] < j)
                    ++curr;

                if (curr < row_end && ja[curr] == j)
                {
                    double val = (j < i) ? al[curr] : au[curr];
                    std::cout << val << " ";
                    ++curr;  // готовимся к следующему j
                }
                else
                {
                    std::cout << "0 ";
                }
            }
        }
        std::cout << "\n";
    }

    std::cout << "\nRight part:\n";
    for (int i = 0; i < size; ++i)
        std::cout << rightPart[i] << " ";
    std::cout << "\n";
}

void RowColStorage::AddLocalMatrix_Garmonic(TreeLinearLagrange& elem, double** LocalMatrix)
{
    int n = elem.n; // число узлов в элементе
    for (int iloc = 0; iloc < n; ++iloc) {
        int glob_i = elem.CoordsIndexes[iloc];
        // Две степени свободы для узла iloc
        int s_free_i = FreeIndex[2 * glob_i];
        int c_free_i = FreeIndex[2 * glob_i + 1];

        for (int jloc = 0; jloc < n; ++jloc) {
            int glob_j = elem.CoordsIndexes[jloc];
            int s_free_j = FreeIndex[2 * glob_j];
            int c_free_j = FreeIndex[2 * glob_j + 1];

            // Блок 2x2 локальной матрицы
            double a_ss = LocalMatrix[2 * iloc][2 * jloc];
            double a_sc = LocalMatrix[2 * iloc][2 * jloc + 1];
            double a_cs = LocalMatrix[2 * iloc + 1][2 * jloc];
            double a_cc = LocalMatrix[2 * iloc + 1][2 * jloc + 1];

            // Вспомогательная лямбда для обработки одной пары (i_free, j_free, val)
            auto add_pair = [&](int i_free, int j_free, double val, int comp_j) {
                if (i_free != -1 && j_free != -1) {
                    if (i_free == j_free) {
                        di[i_free] += val;
                    }
                    else {
                        // Ищем j_free в строке i_free
                        int beg = ia[i_free];
                        int end = ia[i_free + 1];
                        auto it = std::lower_bound(ja.begin() + beg, ja.begin() + end, j_free);
                        if (it != ja.begin() + end && *it == j_free) {
                            int k = it - ja.begin();
                            if (j_free < i_free)
                                al[k] += val;
                            else
                                au[k] += val;
                        }
                    }
                }
                else if (i_free != -1 && j_free == -1) {
                    double u_j = 0.0;
                    if (comp_j == 0) // s-компонента
                        u_j = DirichletValues[2 * glob_j];
                    else             // c-компонента
                        u_j = DirichletValues[2 * glob_j + 1];
                    rightPart[i_free] -= val * u_j;
                }
                else if (i_free == -1 && j_free != -1) {
                }
                };
            {
                int i_free = s_free_i, j_free = s_free_j;
                double val = a_ss;
                if (i_free != -1 && j_free != -1) {
                    if (i_free == j_free) di[i_free] += val;
                    else {
                        int beg = ia[i_free], end = ia[i_free + 1];
                        auto it = std::lower_bound(ja.begin() + beg, ja.begin() + end, j_free);
                        if (it != ja.begin() + end && *it == j_free) {
                            int k = it - ja.begin();
                            if (j_free < i_free) al[k] += val;
                            else au[k] += val;
                        }
                    }
                }
                else if (i_free != -1 && j_free == -1) {
                    double u = DirichletValues[2 * glob_j];
                    rightPart[i_free] -= val * u;
                }
                //else if (i_free == -1 && j_free != -1) {
                //    double u = DirichletValues[2 * glob_i];
                //    rightPart[j_free] -= val * u;
                //}
            }
            // (s_i, c_j)
            {
                int i_free = s_free_i, j_free = c_free_j;
                double val = a_sc;
                if (i_free != -1 && j_free != -1) {
                    if (i_free == j_free) di[i_free] += val;
                    else {
                        int beg = ia[i_free], end = ia[i_free + 1];
                        auto it = std::lower_bound(ja.begin() + beg, ja.begin() + end, j_free);
                        if (it != ja.begin() + end && *it == j_free) {
                            int k = it - ja.begin();
                            if (j_free < i_free) al[k] += val;
                            else au[k] += val;
                        }
                    }
                }
                else if (i_free != -1 && j_free == -1) {
                    double u = DirichletValues[2 * glob_j + 1];
                    rightPart[i_free] -= val * u;
                }
                //else if (i_free == -1 && j_free != -1) {
                //    double u = DirichletValues[2 * glob_i];
                //    rightPart[j_free] -= val * u;
                //}
            }
            // (c_i, s_j)
            {
                int i_free = c_free_i, j_free = s_free_j;
                double val = a_cs;
                if (i_free != -1 && j_free != -1) {
                    if (i_free == j_free) di[i_free] += val;
                    else {
                        int beg = ia[i_free], end = ia[i_free + 1];
                        auto it = std::lower_bound(ja.begin() + beg, ja.begin() + end, j_free);
                        if (it != ja.begin() + end && *it == j_free) {
                            int k = it - ja.begin();
                            if (j_free < i_free) al[k] += val;
                            else au[k] += val;
                        }
                    }
                }
                else if (i_free != -1 && j_free == -1) {
                    double u = DirichletValues[2 * glob_j];
                    rightPart[i_free] -= val * u;
                }
                //else if (i_free == -1 && j_free != -1) {
                //    double u = DirichletValues[2 * glob_i + 1];
                //    rightPart[j_free] -= val * u;
                //}
            }
            // (c_i, c_j)
            {
                int i_free = c_free_i, j_free = c_free_j;
                double val = a_cc;
                if (i_free != -1 && j_free != -1) {
                    if (i_free == j_free) di[i_free] += val;
                    else {
                        int beg = ia[i_free], end = ia[i_free + 1];
                        auto it = std::lower_bound(ja.begin() + beg, ja.begin() + end, j_free);
                        if (it != ja.begin() + end && *it == j_free) {
                            int k = it - ja.begin();
                            if (j_free < i_free) al[k] += val;
                            else au[k] += val;
                        }
                    }
                }
                else if (i_free != -1 && j_free == -1) {
                    double u = DirichletValues[2 * glob_j + 1];
                    rightPart[i_free] -= val * u;
                }
                //else if (i_free == -1 && j_free != -1) {
                //    double u = DirichletValues[2 * glob_i + 1];
                //    rightPart[j_free] -= val * u;
                //}
            }
        }
    }
}

void RowColStorage::AddToRightPart_Garmonic(TreeLinearLagrange& elem, std::vector<double> v)
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

void RowColStorage::GetDiag(std::vector<double>* d)
{
    for (int i = 0; i < size; ++i)
        (*d)[i] = di[i];
}

void RowColStorage::LU_preconditioning()
{
    L.resize(ia[size], 0.0);
    U.resize(ia[size], 0.0);
    D.resize(size, 0.0);
    for (int i = 0; i < size; i++)
    {
        double sumD = 0;
        int i0 = ia[i];
        int i1 = ia[i + 1];
        for (int ik = i0; ik < i1; ik++)
        {
            double sumL = 0, sumU = 0;
            int ij = ja[ik]; // is a column that we need to multiply on
            int j0 = ia[ij];
            int j1 = ia[ij + 1];

            int i_start = i0;
            int j_start = j0;

            int amount_mult_in_col = j1 - j0;
            int amount_mult_in_line = ik - i0;
            int razn_mult = amount_mult_in_col - amount_mult_in_line;

            if (razn_mult > 0)
            {
                j_start += razn_mult;
            }
            else
            {
                i_start -= razn_mult;
            }

            for (; i_start < ik; i_start++, j_start++)
            {
                //cout << i_start << "(" << ja[i_start] << ")" << " " << j_start << "(" << ja[j_start] << ")" << endl;
                if (ja[i_start] == ja[j_start])
                {
                    sumL += L[i_start] * U[j_start];
                    sumU += L[j_start] * U[i_start];
                }
            }
            U[ik] = au[ik] - sumU;
            L[ik] = (al[ik] - sumL) / di[ij];
            sumD += U[ik] * L[ik];
        }
        D[i] = di[i] - sumD;
    }
}

void RowColStorage::Solve_L(std::vector<double>& y)
{
    for (int i = 0; i < size; i++)
    {
        double sum = 0.0;
        for (int k = ia[i]; k < ia[i + 1]; ++k) 
        {
            int j = ja[k];
            if (j < i)
                sum += L[k] * y[j];
        }
        y[i] = (y[i] - sum);
    }
}

void RowColStorage::Solve_U(std::vector<double>& x)
{
    for (int i = size - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int k = ia[i]; k < ia[i + 1]; ++k) {
            int j = ja[k];
            if (j > i)
                sum += U[k] * x[j];
        }
        x[i] = (x[i] - sum) / D[i];
    }
}

void RowColStorage::Solve_L_Transope(std::vector<double>& y)
{
    //we change L to U, becouse matrix A was transpose
    for (int i = 0; i < size; i++)
    {
        int i0 = ia[i];
        int i1 = ia[i + 1];
        for (int j = i0; j < i1; j++)
        {
            int ij = ja[j];
            y[i] -= U[j] * y[ij];
        }
        y[i] /= D[i];
    }
}

void RowColStorage::Solve_U_Transope(std::vector<double>& x)
{
    //we change U to L, becouse matrix A was transpose
    for (int i = size - 1; i >= 0; i--)
    {
        int i0 = ia[i];
        int i1 = ia[i + 1];
        for (int j = i0; j < i1; j++)
        {
            int ij = ja[j];
            x[ij] -= L[j] * x[i];
        }
    }
}

