using System;
using System.Collections.Generic;
using System.Text;

namespace ChislMet
{
    class Matrix
    {
        protected int sizeN, sizeM;
        protected double[,] data;
        public Matrix(int sizeN, int sizeM)
        {
            this.sizeN = sizeN;
            this.sizeM = sizeM;
            data = new double[sizeN, sizeM];
        }

        public Matrix(double[,] matrix)
        {
            sizeN = matrix.GetLength(0);
            sizeM = matrix.GetLength(1);
            data = new double[sizeN, sizeM];
            for (int i = 0; i < sizeN; i++)
            {
                for (int j = 0; j < sizeM; j++)
                {
                    data[i, j] = matrix[i, j];
                }

            }

        }
        public Matrix Transp()
        {
            Matrix rez = new Matrix(this.sizeM, this.sizeN);
            for (int i = 0; i < this.sizeM; i++)
            {
                for (int j = 0; j < this.sizeN; j++)
                {
                    rez[i, j] = this[j, i];
                }
            }

            return rez;
        }
        public Matrix(Matrix m)
        {
            this.sizeN = m.sizeN;
            this.sizeM = m.sizeM;
            data = new double[sizeN, sizeM];
            for (int i = 0; i < sizeN; i++)
                for (int j = 0; j < sizeM; j++)
                    data[i, j] = m.data[i, j];
        }
        public int GetRows()
        {
            return sizeN;
        }
        public int GetColumns()
        {
            return sizeM;
        }

        public bool SetElement(double el, int indexI, int indexJ)
        {
            if (indexI < 0 || indexI >= sizeN || indexJ < 0 || indexJ >= sizeM) return false;
            data[indexI, indexJ] = el;
            return true;
        }
        public double GetElement(int indexI, int indexJ)
        {
            if (indexI < 0 || indexI >= sizeN || indexJ < 0 || indexJ >= sizeM) return default(double);
            return data[indexI, indexJ];
        }

        public Matrix Copy()
        {
            Matrix rez = new Matrix(data);
            return rez;
        }
        public Vector GetRowVector(int indexN)
        {
            Vector rezRow = new Vector(sizeM);
            for (int i = 0; i < sizeM; i++)
                rezRow[i] = data[indexN, i];
            return rezRow;

        }
        //задать вектор-строку матрицы
        public void SetRowVector(Vector vc, int place)
        {
            for (int j = 0; j < sizeM; j++)
                data[place, j] = vc.GetElement(j);
        }

        public Vector GetColumnVector(int indexM)
        {
            Vector rezCol = new Vector(sizeN);
            for (int i = 0; i < sizeN; i++)
                rezCol[i] = data[i, indexM];
            return rezCol;

        }
        public static Matrix operator *(double c, Matrix m)
        {
            Matrix rez = new Matrix(m.sizeN, m.sizeM);

            for (int i = 0; i < m.sizeN; i++)
                for (int j = 0; j < m.sizeM; j++)
                    rez.data[i, j] = m.data[i, j] * c;
            return rez;

        }
        public static Matrix operator *(Matrix m, double c)
        {
            Matrix rez = new Matrix(m.sizeN, m.sizeM);

            for (int i = 0; i < m.sizeN; i++)
                for (int j = 0; j < m.sizeM; j++)
                    rez.data[i, j] = m.data[i, j] * c;
            return rez;

        }

        public static Matrix operator +(Matrix u, Matrix v)
        {
            if (u.sizeN != v.sizeN || u.sizeM != v.sizeM) return null;
            Matrix rez = new Matrix(v.sizeN, v.sizeM);
            for (int i = 0; i < v.sizeN; i++)
                for (int j = 0; j < v.sizeM; j++)
                    rez.data[i, j] = u.data[i, j] + v.data[i, j];
            return rez;

        }
        public static Matrix operator -(Matrix u, Matrix v)
        {
            if (u.sizeN != v.sizeN || u.sizeM != v.sizeM) return null;
            Matrix rez = new Matrix(v.sizeN, v.sizeM);
            for (int i = 0; i < v.sizeN; i++)
                for (int j = 0; j < v.sizeM; j++)
                    rez.data[i, j] = u.data[i, j] - v.data[i, j];
            return rez;

        }




        public Matrix UMinus()
        {
            Matrix rez = new Matrix(sizeN, sizeM);
            for (int i = 0; i < sizeN; i++)
                for (int j = 0; j < sizeM; j++)
                    rez.data[i, j] = -data[i, j];
            return rez;
        }



        public override string ToString()
        {
            for (int i = 0; i < data.GetLength(0); i++)
            {
                for (int j = 0; j < data.GetLength(1); j++)
                {
                    Console.Write(data[i, j] + "\t");
                }
                Console.WriteLine("");
            }
            Console.WriteLine("\n");

            return "";
        }
        public int GetSizeM()
        {
            return sizeM;
        }
        public int GetSizeN()
        {
            return sizeN;
        }

        public void SwapRows(int f, int s)
        {
            Vector v = GetRowVector(f);
            for (int i = 0; i < sizeM; i++)
            {
                this[f, i] = this[s, i];
            }
            for (int i = 0; i < sizeM; i++)
            {
                this[s, i] = v[i];
            }


        }

        public double this[int indexN, int indexM]
        {
            get
            {
                if (indexN < 0 || indexN >= sizeN || indexM < 0 || indexM >= sizeM) return default(double);
                return data[indexN, indexM];
            }
            set
            {
                if (indexN >= 0 && indexN < sizeN && indexM >= 0 && indexM < sizeM)
                    data[indexN, indexM] = value;
            }
        }
        //умножение матрицы на матрицу
        public static Matrix operator *(Matrix u, Matrix v)
        {
            if (u.sizeM != v.sizeN) return null;
            Matrix rez = new Matrix(v.sizeN, v.sizeM);
            for (int i = 0; i < u.sizeN; i++)
                for (int j = 0; j < v.sizeM; j++)
                    for (int k = 0; k < v.sizeN; k++)
                        rez[i, j] += u[i, k] * v[k, j];

            return rez;
        }
        //умножение матрицы на вектор
        public static Vector operator *(Matrix u, Vector v)
        {
            if (u.sizeM != v.GetSize()) return null;
            Vector rez = new Vector(u.sizeN);
            for (int i = 0; i < u.sizeN; i++)
            {
                rez[i] = 0;
                for (int j = 0; j < u.sizeM; j++)
                    rez[i] += u.data[i, j] * v[j];
            }
            return rez;
        }

        //транспонирование
        public static Matrix Transpose(Matrix mx)
        {
            double[,] t = new double[mx.sizeM, mx.sizeN];
            for (int i = 0; i < mx.sizeN; i++)
                for (int j = 0; j < mx.sizeM; j++)
                    t[j, i] = mx[i, j];
            return new Matrix(t);
        }
    }
}
