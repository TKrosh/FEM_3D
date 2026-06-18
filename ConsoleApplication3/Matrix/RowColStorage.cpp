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
    for (int i = 0; i < size; i++)
    {
        double sum = 0.0;
        for (int k = ia[i]; k < ia[i + 1]; k++)
        {
            int j = ja[k]; // is a index in vector x that we need to multiply on
            if (j < i)
                sum += au[k] * tmp[j];
            else if (j > i)
                sum += al[k] * tmp[j];
        }
        (*v)[i] = tmp[i] * di[i] + sum;
    }
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


#include <algorithm>  // для std::lower_bound

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


