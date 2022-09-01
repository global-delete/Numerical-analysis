using System;
using System.Collections.Generic;
using System.Text;

namespace ChislMet
{
    public delegate double Func(double x); //делегат созданный, чтобы в метод можно было передавать функцию 
    class NonLinRoot
    {
        //корни уравнения методом половинного деления
        public double Half(double left, double right, double eps, Func func)
        {
            // left - левая граница   right-правая граница   eps- заданная точность ɛ
            double delta = right - left; //середина отрезка
            double curDelta = delta;
            double Fmin = func(left); //Найдем значение функции в левой границе
            double Fmax = func(right); //Найдем значение функции в правой границе
            int k = 2;

            while (curDelta > eps)  // пока не будет достигнута заданная точность ɛ
            {
                double x = (left + right) / 2; //Разобьем отрезок пополам
                double Fx = func(x); //Найдем значение функции в точке х
                k++;
                /*
                 Проверим условие Fmin*F(x) < 0. Если условие выполнено, то корень расположен на отрезке [left,х].
                 В этом случае необходимо точку right переместить в точку х. Если условие не выполнено, то
                 корень расположен на отрезке [х,right]. В этом случае необходимо точку left переместить в точку х.
                 */
                if (Fmin * Fx < 0)
                {
                    right = x;
                    Fmax = Fx;
                }
                else
                {
                    left = x;
                    Fmin = Fx;
                }
                curDelta = (right - left);
            }

            return (left + right) / 2;
        }

        //корни уравнения методом Ньютона
        public double Newton(double left, double right, double eps, Func func)
        {
            // left - левая граница   right-правая граница   eps-заданная точность ɛ   

            double cur = left + (right - left) / 2;
            bool flag = false;
            double nextCur = 0; // следующее приближение
            double oldDelta = 0;
            double delta;
            do
            {
                if (flag)
                    cur = nextCur;
                double fcur = func(cur);
                double der = (func(cur + eps / 2) - fcur) / (eps / 2); // линейная зависимость в окрестности точки
                nextCur = cur - func(cur) / der; // формула метода простой итерации
                delta = Math.Abs(nextCur - cur); // величина поправки

                if ((delta <= oldDelta) || !flag)
                    oldDelta = delta;
                else return double.NaN;
                flag = true;

            } while (delta > eps); // пока не будет достигнута заданная точность ɛ

            return nextCur;
        }



        //корни уравнения методом последовательных приближений
        public double Priblij(double left, double right, double eps, Func func)
        {
            // left - левая граница   right-правая граница   eps-заданная точность ɛ   

            double x1 = left + (right - left) / 2; //xk+1
            double x0;
            double oDelta, nDelta = 0; //δхк, δхк+1
            int iterations = 0;
            

            do
            {

                oDelta = nDelta;
                x0 = x1;

                //перобразование f(x) к φ(x) 
                x1 = Transform(func, eps, x0);

                iterations++;
              
               
                nDelta = x1 - x0; // вычисляем оценку

                if (!(Math.Abs(nDelta) >= eps) && (iterations < 2 || Math.Abs(nDelta) < Math.Abs(oDelta)))
                {
                    return double.NaN;
                }

            } while ((Math.Abs(nDelta) >= eps) && (iterations < 2 || Math.Abs(nDelta) < Math.Abs(oDelta)));

            return x1;

        }

        //перобразование f(x) к φ(x) 
        public double Transform(Func func, double eps, double x)
        {
            double M = -(func(x + eps) - func(x - eps)) / (2 * eps);
            return x + func(x) / M;
        }
    }
}
