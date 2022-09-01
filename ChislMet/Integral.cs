using System;
using System.Collections.Generic;
using System.Text;

namespace ChislMet
{
    class Integral
    {
        //функция для интеграла
        public delegate double xFunc(double x);
        //функция для двойного интеграла
        public delegate double yxFunc(double x, double y);



        //интегрирование методом прямоугольников
        public double IntegralRectangle(double down, double up, double eps, xFunc f)
        {
            int n = 1;
            double s = 0;
            double prev = 0, h = 0;
            do
            {
                prev = s;
                n *= 2;
                h = (up - down) / n;
                for (int i = 0; i < n; i++)
                {
                    s += f(down + h * (i + 0.5));
                }
                s *= h;
            } while (Math.Abs(s - prev) > eps);
            return s;
        }


        //интегрирование методом трапеций
        public double IntegralTrap(double down, double up, double eps, xFunc f)
        {
            double x = down;
            int n = 1;
            double h = (up - down) / n;
            double prev, s = 0, res = 0;
            double sn = (f(down) + f(up)) / 2.0;
            do
            {
                prev = res;
                n *= 2;
                h = (up - down) / n;
                for (int i = 1; i < n; i = i + 2)
                {
                    s += f(down + i * h);
                }
                res = h * (sn + s);
            } while (Math.Abs(res - prev) > eps);
            return res;
        }
        //интегрирование методом Симпсона
        public double IntegralSimpson(double down, double up, double eps, xFunc f)
        {
            //prev-предыдущий вычисленый интеграл
            //cur-новый интеграл, с большим N.
            double prev = eps + 1, cur = 0;

            for (int n = 2; (n <= 4) || (Math.Abs(cur - prev) > eps); n *= 2)
            {
                double h, sumChet = 0, sumNech = 0, sum = 0;
                //Шаг интегрирования вложили в метод
                h = (up - down) / (2 * n);
                for (int i = 1; i <= 2 * n - 1; i += 2)
                {
                    //Значения с нечётными индексами
                    sumNech += f(down + h * i);
                    //Значения с чётными индексами
                    sumChet += f(down + h * (i + 1));
                }
                sum = f(down) + 4 * sumNech + 2 * sumChet - f(up);
                prev = cur;
                cur = (h / 3) * sum;
            }
            return cur;
        }


        //интегрирование методом Симпсона для двойного интеграла
        public double IntegralDoubleSimpson(double down, double up, double yDown, double yUp, int N, yxFunc func)
        {
            //матрица для проверки
            // Matrix testMat = new Matrix(2 * N, 2 * N);

            double h = (up - down) / (2 * N); //Шаг по х
            double k = (yUp - yDown) / (2 * N); //Шаг по у

            double[] str = new double[2 * N];

            str[0] = 1;

            bool flag = false;

            for (int i = 1; i < str.Length - 1; i++)
            {
                if (flag)
                    str[i] = 2;
                else
                    str[i] = 4;

                flag = !flag;
            }

            str[str.Length - 1] = 1;

            double sum = 0;

            for (int i = 0; i < 2 * N; i++)
            {
                for (int j = 0; j < 2 * N; j++)
                {
                    var l = str[i];

                    if (j > 0 && j < (2 * N) - 1)
                    {
                        if ((j & 1) == 1)
                            l *= 4;
                        else
                            l *= 2;
                    }

                    //testMat[i, j] = l;

                    sum += l * func(down + h * i, yDown + k * j);
                }
            }

            //Console.WriteLine(testMat);

            return (h * k * sum) / 9;
        }


    }
}
