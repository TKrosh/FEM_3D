#pragma once
#include <vector>
#include "../Matrix/GlobalMatrix.h"

template<typename Storage>
class MatrixSolver
{
public:
	int MaxIters, Size;
	double Eps, fNorm;
	GlobalMatrix<Storage> Matrix;

	std::vector<double> _r;
	std::vector<double> _Ax;
	std::vector<double> _z;
	std::vector<double> _p;
	std::vector<double> _x;

	MatrixSolver(int maxIters,
		double eps,
		GlobalMatrix<Storage> matrix) :
		MaxIters(maxIters),
		Eps(eps),
		Matrix(matrix),
		Size(matrix.storage.size)
	{
		_r.resize(Size, 0.0);
		_Ax.resize(Size, 0.0);
		_z.resize(Size, 0.0);
		_p.resize(Size, 0.0);
		_x.resize(Size, 0.0);
	};

	MatrixSolver() : MaxIters(10000), Eps(1e-11), Size(0) {};

	std::vector<double>LOC()
	{
		fNorm = scalarMult(Matrix.storage.rightPart, Matrix.storage.rightPart);
		//r0
		CalcResidual(); //(f - Ax0)

		//z0
		for (int i = 0; i < Size; i++) _z[i] = _r[i];

		//p0
		for (int i = 0; i < Size; i++) _p[i] = _z[i];
		Matrix.multiplicationByVector(&_p);

		double mistake = sqrt(scalarMult(_r, _r) / fNorm);
		for (int k = 0; k < MaxIters && mistake > Eps; k++)
		{
			double mistake_r = scalarMult(_r, _r);
			double scalc_p = scalarMult(_p, _p);

			if (scalc_p < Eps) break;
			double alfa_k = scalarMult(_p, _r) / scalc_p;

			for (int i = 0; i < Size; i++) 
			{
				_x[i] += alfa_k * _z[i];
				_r[i] -= alfa_k * _p[i];
			}

			// A U-1 r
			for (int i = 0; i < Size; i++) _Ax[i] = _r[i];
			Matrix.multiplicationByVector(&_Ax);
			double betta = -(scalarMult(_p, _Ax) / scalc_p);

			for (int i = 0; i < Size; i++) 
			{
				_z[i] = _r[i] + betta * _z[i];
				_p[i] = _Ax[i] + betta * _p[i];
			}
			mistake = mistake_r - alfa_k * alfa_k * scalc_p;

			if (alfa_k < 1e-15) break;
		}
		//std::cout << std::uppercase << std::scientific << std::setprecision(15) << mistake << std::endl;
		return _x;
	}

	std::vector<double>LOC_D()
	{
		std::vector du(Size, 0.0);
		Matrix.GetDiag(&du);
		fNorm = scalarMult(Matrix.storage.rightPart, Matrix.storage.rightPart);
		//r0
		CalcResidual(); //(f - Ax0)
		for (int i = 0; i < Size; i++) _r[i] = _r[i] / sqrt(du[i]);
		//z0
		for (int i = 0; i < Size; i++) _z[i] = _r[i] / sqrt(du[i]);
		//p0
		for (int i = 0; i < Size; i++) _p[i] = _z[i];
		Matrix.multiplicationByVector(&_p);
		for (int i = 0; i < Size; i++) _p[i] /= sqrt(du[i]);

		double mistake = sqrt(scalarMult(_r, _r) / fNorm);

		for (int k = 0; k < MaxIters && mistake > Eps; k++)
		{
			double mistake_r = scalarMult(_r, _r);
			double scalc_p = scalarMult(_p, _p);

			if (scalc_p < Eps) break;
			double alfa_k = scalarMult(_p, _r) / scalc_p;

			for (int i = 0; i < Size; i++)
			{
				_x[i] += alfa_k * _z[i];
				_r[i] -= alfa_k * _p[i];
			}

			// L-1 A U-1 r
			for (int i = 0; i < Size; i++) _Ax[i] = _r[i] / sqrt(du[i]);
			Matrix.multiplicationByVector(&_Ax);
			for (int i = 0; i < Size; i++) _Ax[i] = _Ax[i] / sqrt(du[i]);
			double betta = -(scalarMult(_p, _Ax) / scalc_p);

			for (int i = 0; i < Size; i++) 
			{
				_z[i] = _r[i] / sqrt(du[i]) + betta * _z[i];
				_p[i] = _Ax[i] + betta * _p[i];
			}
			mistake = mistake_r - alfa_k * alfa_k * scalc_p;

			if (alfa_k < 1e-15) break;
		}
		return _x;
	}

	std::vector<double>CGM_D()
	{
		std::vector du(Size, 0.0);
		Matrix.GetDiag(&du);
		fNorm = scalarMult(Matrix.storage.rightPart, Matrix.storage.rightPart);

		//r0
		CalcResidual(); //(f - Ax0)
		for (int i = 0; i < Size; i++) _r[i] = _r[i] / du[i];
		Matrix.multiplicationByVector(&_r);
		for (int i = 0; i < Size; i++) _r[i] = _r[i] / sqrt(du[i]);

		for (int i = 0; i < Size; i++) _z[i] = _r[i];
		for (int i = 0; i < Size; i++) _x[i] *= sqrt(du[i]);
		double mistake = sqrt(scalarMult(_r, _r) / fNorm);

		for (int k = 0; k < MaxIters && mistake > Eps; k++)
		{
			double mistake_r = scalarMult(_r, _r);

			for (int i = 0; i < Size; i++) _Ax[i] = _z[i] / sqrt(du[i]);
			Matrix.multiplicationByVector(&_Ax);
			for (int i = 0; i < Size; i++) _Ax[i] /= du[i];
			Matrix.multiplicationByVector(&_Ax);
			for (int i = 0; i < Size; i++) _Ax[i] /= sqrt(du[i]);

			double alfa_k = mistake_r / scalarMult(_Ax, _z);
			if (alfa_k < 1e-15) break;

			for (int i = 0; i < Size; i++)
			{
				_x[i] += alfa_k * _z[i];
				_r[i] -= alfa_k * _Ax[i];
			}

			double betta = scalarMult(_r, _r) / mistake_r;
			for (int i = 0; i < Size; i++)
			{
				_z[i] = _r[i] + betta * _z[i];
			}

			mistake = sqrt(scalarMult(_r, _r) / fNorm);
		}

		for (int i = 0; i < Size; i++) _x[i] /= sqrt(du[i]);
		return _x;
	}

	//f - Ax
	void CalcResidual()
	{
		std::vector<double> r;

		Matrix.multiplicationByVector(&_r);
		for (int i = 0; i < Size; i++)
		{
			_r[i] = Matrix.storage.rightPart[i] - _r[i];
		};
	}

	//(x, x)
	double scalarMult(std::vector<double> x, std::vector<double> y)
	{
		int size = x.size();
		double res = 0.0;
		for (int i = 0; i < size; i++)
			res += x[i] * y[i];
		return res;
	};
};
