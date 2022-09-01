using System;

namespace ChislMet
{
    class Program
    {
        static void Main(string[] args)
        {
            Vector v1 = new Vector(3);
            v1[0] = 2.5; v1[1] = -0.3;
            Vector v3 = new Vector(3);
            v3[0] = 1; v3[1] = 2; v3[2] = 3;
            double[] dd = new double[2] { 0.7, 1.25 };
            Vector v2 = new Vector(dd);
            Matrix m1 = new Matrix(3, 3);
            m1[0, 0] = 9.0; m1[0, 1] = 3.0; m1[0, 2] = 4.0;
            m1[1, 0] = 0.0; m1[1, 1] = 6.0; m1[1, 2] = 1.0;
            m1[2, 0] = 0.0; m1[2, 1] = 7.0; m1[2, 2] = 8.0;
            Console.WriteLine(m1);
            //   Vector v4 = LinAlg.LU_UT(m1,v3);
            Console.WriteLine(v3);

            LinAlg.Gauss(m1,v3);
            Console.WriteLine(v3);

            //  Console.WriteLine(v4);

            LinAlg LA = new LinAlg();
            //Прогонка
            double[,] LAdm8 = new double[3, 3] { {1, 1, 0},
                                                 {3,-1, 2},
                                                 {0, 1,-1},
                                             };
            Matrix LAm8 = new Matrix(LAdm8);

            double[] e1 = new double[3] { 1, 2, 0 };
            double[] d1 = new double[3] { 1, -1, -1 };
            double[] c1 = new double[3] { 0, 3, 1 };
            double[] b1 = new double[3] { 3, 1, 2 };
            Vector d = new Vector(d1);
            Vector e = new Vector(e1);
            Vector c = new Vector(c1);
            Vector b = new Vector(b1);
            Console.WriteLine(LA.Progonka(e, d, c, b));

        
         
            Integral integ = new Integral();
            double xDown = 5;
            double xUp = 10;
            double yDown = 5;
            double yUp = 10;

            double f(double x)
            {
                return x * x;
            }

            Console.WriteLine("\nМетод прямоугольников");
            Console.WriteLine(integ.IntegralRectangle(xDown, xUp, 0.001, f));

            Console.WriteLine("\nМетод трапеций");
            Console.WriteLine(integ.IntegralTrap(xDown, xUp, 0.0001, f));

            Console.WriteLine("\nМетод Cимпсона");
            Console.WriteLine(integ.IntegralSimpson(xDown, xUp, 0.00001, f));

            double xy(double x, double y)
            {
                return x * x + y * y;
            }

            Console.WriteLine("\nМетод Cимпсона двойного интеграла\n");

            int N = 5;

            Console.WriteLine(integ.IntegralDoubleSimpson(xDown, xUp, yDown, yUp, N, xy));


            //поиск корней уравнения
            Console.WriteLine("\nКорни ур-ий:");
            NonLinRoot kr = new NonLinRoot();

            double l = -0.5, r = 5.0, epsi = 0.01;

            Console.WriteLine(); ;

            Console.Write("Метод половинного деления: ");
            Console.Write(kr.Half(l, r, epsi, x => (x * x - 9)));

            Console.WriteLine();

            Console.Write("Метод Ньютона: ");
            Console.Write(kr.Newton(l, r, epsi, x => (x * x - 9)));

            Console.WriteLine();

            Console.Write("Метод Последовательного приближения: ");
            Console.Write(kr.Priblij(l, r, epsi, x => (x * x - 9)));

            Console.WriteLine();



            Console.ReadKey();
        }
    }
}










/*    Vector v1 = new Vector(2);
            v1[0] = 2.5; v1[1] = -0.3;
            double[] dd = new double[2] { 0.7, 1.25 };
            Vector v2 = new Vector(dd);
            double sc = v1 * v2;
            Console.WriteLine(sc);
            Console.WriteLine(v1);
            Console.WriteLine(v2);
            Vector v3 = v1.Add(v2);
            Console.WriteLine(v3);

            Console.WriteLine("{0}+{1}={2}", v1, v2, v3, v1 + v2);
            Vector v4 = v3 * 5.0;
            Console.WriteLine(v4);
*/