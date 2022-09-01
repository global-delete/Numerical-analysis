using System;
using System.Collections.Generic;
using System.Text;

namespace ChislMet
{
    class LinAlg
    {
        public static Vector LU_NT(Matrix A, Vector B)
        {
            if (A.GetSizeM() != A.GetSizeN() || A.GetSizeN() != B.GetSize()) return null;
            for (int i = 0; i < A.GetSizeN(); i++)
            {
                if (A[i, i] == 0) return null;
                for (int j = i + 1; j < A.GetSizeN(); j++)
                {
                    if (A[i, j] != 0) return null;
                }
            }
            Vector x = new Vector(A.GetSizeN());
            x[0] = B[0] / A[0, 0];
            for (int i = 1; i < A.GetSizeN(); i++)
            {
                double s = 0;
                for (int k = 0; k < i; k++) s += A[i, k] * x[k];
                x[i] = (B[i] - s) / A[i, i];

            }
            return x;

        }


        public static Vector LU_UT(Matrix A, Vector B)
        {
            int An = A.GetSizeN(), Am = A.GetSizeM(), BSize = B.GetSize();
            if (An != Am || An != BSize)
            {
                return null;
            }
            for (int i = An - 1; i >= 0; i--)
            {
                if (A[i, i] == 0)
                    return null;
                for (int j = i - 1; j >= 0; j--)
                {
                    if (A[i, j] != 0) return null;
                }

            }

            Vector X = new Vector(An);
            for (int i = 0; i < An; i++)
            {
                int m = An - i - 1;
                double s = 0;
                for (int j = m; j < An; j++)
                {
                    s += A[m, j] * X[j];
                }
                X[m] = (B[m] - s) / A[m, m];
            }
            X[An - 1] = B[An - 1] / A[An - 1, An - 1];
            return X;
        }

        public static Vector Gauss(Matrix A, Vector B)
        {
            int An = A.GetSizeN(), Am = A.GetSizeM(), BSize = B.GetSize();
            if (An != Am || An != BSize)
            {
                return null;
            }
            Vector X = new Vector(An);
            int n = 0;
            // прямой ход
            for (int j = 0; j < Am; j++)
            {
                double s = 0;

                for (int i = j; i < An; i++)
                {
                    Vector vec = A.GetColumnVector(j);
                    if (Math.Abs(vec[i]) > s)
                    {
                        s = Math.Abs(vec[i]);
                        n = i;
                    }
                }
                A.SwapRows(j, n);
                B.SwapRows(j, n);

                for (int i = j + 1; i < Am; i++) // minus, tmn
                {
                    double d = A[i, j] / A[j, j];
                    for (int k = 0; k < An; k++)
                    {
                        A[i, k] -= (d * A[j, k]);
                    }
                    B[i] -= (d * B[i]);
                }

            }
            //обратный ход
            for (int i = An - 1; i >= 0; i--)
                if (A[i, i] == 0) return null;

            for (int i = 0; i < An; i++)
            {
                int m = An - i - 1;
                double s = 0;
                for (int j = m; j < An; j++)
                {
                    s += A[m, j] * X[j];
                }
                X[m] = (B[m] - s) / A[m, m];
            }
            X[An - 1] = B[An - 1] / A[An - 1, An - 1];
            return X;
        }
        // Решение СЛУ с трех диагональной матрицей методом прогонки

        // e - верхняя диагональ 
        // d - главная(средняя) диагональ 
        // c - нижняя диагональ
        // b - правая часть (столбец)

        public Vector Progonka(Vector e, Vector d, Vector c, Vector b)
        {
            int n = e.size;
            Vector x = new Vector(n);
            Vector u = new Vector(n);
            Vector v = new Vector(n);

            u[1] = -e[0] / d[0];//альфа1
            v[1] = b[0] / d[0];//бета1

            for (int i = 1; i < n - 1; i++)
            {
                u[i + 1] = -e[i] / (d[i] + c[i] * u[i]);
                double ddd = (d[i] + c[i] * u[i]);

                v[i + 1] = (-c[i] * v[i] + b[i]) / ddd;
            }
            double fff = (d[n - 1] + c[n - 1] * u[n - 1]);
            x[n - 1] = (-c[n - 1] * v[n - 1] + b[n - 1]) / fff;

            for (int i = n - 2; i >= 0; i--)
                x[i] = u[i + 1] * x[i + 1] + v[i + 1];

            return x;
        }

        //метод Грамма-Шмидта
        public Vector GramSchmidt(Matrix mx, Vector vc, int size)
        {
            if (size < 1)
            {
                return null;
            }

            //вспомогательные матрицы
            Matrix mu = new Matrix(new double[size, size]);
            Matrix mv = new Matrix(new double[size, size]);

            //вспомогательные векторы
            Vector temp = new Vector(size);
            Vector t = new Vector(size);

            //вектор результатов
            Vector result = new Vector(size);
            double y;

            //обнуление вспомогательных матриц
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    mu[i, j] = 0;
                    mv[i, j] = 0;
                }

            //перенос нулевой строки данной матрицы в нулевую вспомогательную матрицу
            for (int i = 0; i < size; i++)
                mu[0, i] = mx[0, i];

            //добавление во вторую вспомогательную матрицу нормализованного
            //вектора(нулевую строки) из первой вспомогательной матрицы
            mv.SetRowVector(mu.GetRowVector(0).Normalize(), 0);

            //проведение ортогонализации нулевого элемента в векторе a/|b|
            t[0] = vc[0] / mu.GetRowVector(0).Norma1();

            int count = 1;
            do
            {
                double temp_h;
                temp_h = 0;

                //зануление вспомогателльного вектора
                for (int i = 0; i < size; i++)
                    temp[i] = 0;

                for (int j = 0; j < count; j++)
                {
                    y = mx.GetRowVector(count) * mv.GetRowVector(j);

                    for (int i = 0; i < size; i++)
                        temp[i] += y * mv[j, i];

                    temp_h += y * t[j];
                }

                //вычисление матрицы mu
                for (int i = 0; i < size; i++)
                    mu[count, i] = mx[count, i] - temp[i];

                //ортогонализация остальных элементов вектора(деление вектора на длину (a / |a|) )
                t[count] = (vc[count] - temp_h) / mu.GetRowVector(count).Norma1();
                //создание ортонормированного базиса(деление вектора на двойную длину (a / ||a||) )
                mv.SetRowVector(mu.GetRowVector(count).Normalize().Normalize(), count);

                count++;

            } while (count < size);

            //ортонормированный базис умножаем на нормированный вектор
            result = Matrix.Transpose(mv) * t;

            return result;
        }

        // Метод последовательных приближений
        public Vector Iterations(Matrix matrix, Vector v, double epsilon)
        {


            var _matrix = DominantDiag(matrix, v);

            var _work = BuildMatrix(_matrix);
            var _matrixSize = _matrix.GetRows();



            var solution = new Vector(new double[_matrixSize]);

            var current_sol = new Vector(new double[_matrixSize]);

            for (int i = 0; i < _matrixSize; i++)
                current_sol[i] = _work[i, _matrixSize];

            var prev_iteration = new Vector(new double[_matrixSize]);


            while (IterationStop(current_sol, prev_iteration, epsilon))
            {
                prev_iteration = current_sol;
                current_sol = new Vector(new double[_matrixSize]);

                for (int i = 0; i < _matrixSize; i++)
                {
                    for (int j = 0; j < _matrixSize; j++)
                    {
                        current_sol[i] += _work[i, j] * prev_iteration[j];
                    }

                    current_sol[i] += _work[i, _matrixSize];
                }
            }

            solution = current_sol;

            return solution;
        }



        //построение матрицы C и вектора d Ax=B =>  x= Bx+g  g=D^(-1)*B
        public static Matrix BuildMatrix(Matrix matrix)
        {
            Matrix _work = new Matrix(matrix);

            var _matrixSize = matrix.GetRows();

            for (int i = 0; i < _matrixSize; i++)
            {
                int a_ii = (int)_work[i, i];

                if (a_ii != 0)
                {
                    _work[i, i] = 0;

                    _work[i, _matrixSize] /= a_ii;
                    for (int j = 0; j < _matrixSize; j++)
                    {
                        if (i != j)
                        {
                            _work[i, j] /= -a_ii;
                        }
                    }
                }
                else
                {
                    _work[i, i] = 1;

                    for (int j = 0; j < _matrixSize; j++)
                    {
                        if (i != j)
                        {
                            _work[i, j] *= -1;
                        }
                    }
                }
            }

            return _work;
        }
        //построение доминантной диагонали в матрице значений
        public static Matrix DominantDiag(Matrix matrix, Vector B)
        {
            Matrix items = new Matrix(matrix);
            Vector b = new Vector(B.data);
            double max;
            int indmax;

            for (int i = 0; i < b.size; i++)
            {
                max = Math.Abs(items[i, i]);
                indmax = i;
                for (int j = i; j < b.size; j++)
                    if (Math.Abs(items[j, i]) > max)
                    {
                        max = Math.Abs(items[j, i]);
                        indmax = j;
                    }

                if (i != indmax)
                {
                    items.SwapRows(i, indmax);
                    double temp = b[i];
                    b[i] = b[indmax];
                    b[indmax] = temp;
                }
            }

            return AddVectorMatrix(items, b);
        }
        //добавление вектора в новый столбец матрицы
        public static Matrix AddVectorMatrix(Matrix m, Vector v)
        {
            double[,] q = new double[m.GetRows(), m.GetColumns() + 1];
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetColumns() + 1; j++)
                {
                    if (j < m.GetColumns())
                        q[i, j] = m[i, j];
                    else
                        q[i, j] = v[i];
                }
            }

            return new Matrix(q);
        }
        //нормализация вектора(максимальный элемент)
        private static double MaxAbs(Vector vc)
        {
            double normal = 0;
            double abs_vc_i = 0;
            for (int i = 0; i < vc.size; i++)
            {
                abs_vc_i = Math.Abs(vc[i]);

                if (normal < abs_vc_i)
                    normal = abs_vc_i;
            }

            return normal;
        }

        //критерий точности для остановки
        public static bool IterationStop(Vector a, Vector b, double epsilon)
        {
            Vector c = new Vector(a.data);


            c = c - b;

            if (MaxAbs(c) <= epsilon) return false;
            else return true;
        }



    }
}