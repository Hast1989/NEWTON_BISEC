#include <iostream>
#include<math.h>
#include<vector>
#include <fstream>
#include<cmath>
#include <iomanip>
#include <fstream>
int iter=0;
double PI = 3.141592653589793;
double delta = 0.00001;
double Function(double x)
{
    return std::sqrt(40 - x * x) * std::sin(x * PI) * std::cos(std::sqrt(40 - x * x) * PI) - x * std::sin(std::sqrt(40 - x * x) * PI) * std::cos(x * PI);
    double r1, r2;
    r1 = std::pow((sqrt(3) * x * x * x - 2 * x + 5)/(7+sqrt(7)), (log(3) / log(10)));
    r2 = std::asin((x * x + x + sqrt(3)) / (2 * x - 2));
    //return r1 + r2;
    return std::tan(x) - x;
    return sin(x * 10.);
    return(x + 0.1) * (x + 0.22) * (x + 0.55) * (x + 0.75)*(x+0.7);
    return 35 * x*x*x - 67 * x*x - 3 * x + 3;
    return sqrt(x + 1) - 1;
}
double Derivative(double x)
{
    return (Function(x + delta) - Function(x)) / delta;
}
std::vector<double> local(double a, double b, int n)
{
    std::vector<double> points;
    double x1, x2;
    double h = (b - a) / n;
    for (int i = 0; i < n; i++)
    {
        x1 = a + i * h;
        x2 = a + (i + 1) * h;
        //std::cout << x1 << ' ' << x2 << std::endl;
       // std::cout << Function(x1) << ' ' << Function(x2) << std::endl;
        if (Function(x1) * Function(x2) < 0)
        {
            
            points.push_back(x1);
            points.push_back(x2);
        }
    }
    return points;
}
double bisec(double a, double b,double eps)
{
    double a1, b1,h;
    a1 = a;
    b1 = b;
    iter = 0;
    while (b1-a1 > eps)
    {
        iter++;
        h = (a1 + b1) / 2;
        if (Function(a1) * Function(h) < 0)
        {
            b1 = h;
        }
        else
        {
            a1 = h;
        }
        if (fabs(Function(a1)) < eps)
        {
            return a1;
        }
        if (fabs(Function(b1)) < eps)
        {
            return b1;
        }
        //std::cout << a1 << ' ' << b1 << std::endl;
       //std::cout << Function(a1) << ' ' << Function(b1) << std::endl;

    }
    return (a1 + b1) / 2;
}
double newton(double a, double b, double eps)
{
    double x0, xk;
    x0 = (Function(a) * b - Function(b) * a) / (Function(a) - Function(b));
   //x0 = 8;
    xk = x0;
    iter = 0;
    while ( fabs(Function(xk))> fabs(Derivative(xk))*eps)
    {
        iter++;
        xk = x0 - Function(x0)/ Derivative(x0);
        if ((xk > a) && (xk < b))
        {
            x0 = xk;
        }
        else
        {
            if (Function(x0) * Function(b) < 0)
            {
                xk = (Function(x0) * b - Function(b) * x0) / (Function(x0) - Function(b));
            }
            else
            {
                xk = (Function(a) * x0 - Function(x0) * a) / (Function(a) - Function(x0));
            }
            x0 = xk;
        }
        //std::cout << xk << ' ' << Function(xk) << std::endl;
        
    }
    
    return xk;
}
void test1()
{
    std::cout<< std::setprecision(15);
    std::vector<double> points;
    double a, b, f;
    double eps = 0.00000000000001;
    int n;
    a = -std::sqrt(40);
    b = 0;
    n = 1;
    points = local(a, b, 100);
    for (int i = 0; i < points.size(); i++)
    {
        f = bisec(points[i], points[i + 1], eps);
        //std::cout << points[i] << ' '<< points[i+1] << std::endl;
        std::cout << n << ' ' << "Bisec(eps=" << eps << ") " <<  f<< ' ' << "iter " << iter<<' '<<Function(f)<< std::endl;
        i++;
        n++;
    }
    n = 1;
    for (int i = 0; i < points.size(); i++)
    {
        //std::cout << points[i] << ' '<< points[i+1] << std::endl;
        f = newton(points[i], points[i + 1], eps);
        std::cout << n << ' ' << "Newton(eps=" << eps << ") " << f << ' ' << "iter " << iter << ' ' << Function(f) << std::endl;
        i++;
        n++;
    }
}
double Function1(double x, double y)
{
    return y * y - 0.3 * x * (x - 4) * (x - 4);
}
double Function2(double x, double y)
{
    return y-9/(x*x+5);
}
double AD11(double x, double y)
{
    return -0.9 * x * x + 4.8 * x - 4.8;
}
double AD12(double x, double y)
{
    return 2*y;
}
double AD21(double x, double y)
{
    return (18*x)/((5+x*x)* (5 + x * x));
}
double AD22(double x, double y)
{
    return 1;
}
double D11(double x, double y)
{
    return (Function1(x+delta,y)- Function1(x, y))/delta;
}
double D12(double x, double y)
{
    return (Function1(x, y + delta) - Function1(x, y)) / delta;
}
double D21(double x, double y)
{
    return (Function2(x + delta, y) - Function2(x, y)) / delta;
}
double D22(double x, double y)
{
    return (Function2(x , y + delta) - Function2(x, y)) / delta;
}
std::vector<double> setka(double L1, double L2, int n)
{
    std::vector<double> points;
    double h1, h2;
    h1 = 2 * L1 / n;
    h2 = 2 * L2 / n;
    for (int i = 0; i < n+1; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            points.push_back(-L1 + i * h1);
            points.push_back(-L2 + j * h2);
        }
    }
    return points;
}
void newtonAD(double x,double y,double eps,double* result)
{
    double df11,f1,f2,df12,df21,df22,xk, yk, xkl, ykl,det;
    xkl = x;
    ykl = y;
    xk = x + 1;
    yk = y + 1;
    iter = 0;
    while (sqrt((xkl - xk) * (xkl - xk) + (ykl - yk) * (ykl - yk)) > eps)
    {
        xkl = xk;
        ykl = yk;
        iter++;
        df11 = AD11(xkl, ykl);
        df12 = AD12(xkl, ykl);
        df21 = AD21(xkl, ykl);
        df22 = AD22(xkl, ykl);
        f1 = Function1(xkl, ykl);
        f2 = Function2(xkl, ykl);
       // std::cout << df11 << ' ' << df12 <<' '<<df21<<' '<<df22<< std::endl;
        det = df11 * df22 - df21 * df12;
        //std::cout << det<< std::endl;
        xk = xkl - (df22 * f1 - df12 * f2) / det;
        yk = ykl + (df21 * f1 - df11 * f2) / det;
        //std::cout << xk << ' ' << yk << std::endl;
        //std::cout << iter << std::endl;
        if (iter == 30)
        {
            break;
        }
    }
    result[0] = xk;
    result[1] = yk;
}
void newtonD(double x, double y, double eps, double* result)
{
    double df11, f1, f2, df12, df21, df22, xk, yk, xkl, ykl, det;
    xkl = x;
    ykl = y;
    xk = x + 1;
    yk = y + 1;
    iter = 0;
    while (sqrt((xkl - xk) * (xkl - xk) + (ykl - yk) * (ykl - yk)) > eps)
    {
        xkl = xk;
        ykl = yk;
        iter++;
        df11 = D11(xkl, ykl);
        df12 = D12(xkl, ykl);
        df21 = D21(xkl, ykl);
        df22 = D22(xkl, ykl);
        f1 = Function1(xkl, ykl);
        f2 = Function2(xkl, ykl);
        //std::cout << df11 << ' ' << df12 << ' ' << df21 << ' ' << df22 << std::endl;
        det = df11 * df22 - df21 * df12;
        //std::cout << det << std::endl;
        xk = xkl - (df22 * f1 - df12 * f2) / det;
        yk = ykl + (df21 * f1 - df11 * f2) / det;
        //std::cout << xk << ' ' << yk << std::endl;
        //std::cout << iter << std::endl;
        if (iter == 30)
        {
            break;
        }
    }
    result[0] = xk;
    result[1] = yk;
}
void newtonAD1Z(double x, double y, double eps, double* result)
{
    double df11, f1, f2, df12, df21, df22, xk, yk, xkl, ykl, det, omega, zx, zy;
    xkl = x;
    ykl = y;
    xk = x + 1;
    yk = y + 1;
    iter = 0;
    omega = 0.95;
    while (sqrt((xkl - xk) * (xkl - xk) + (ykl - yk) * (ykl - yk)) > eps)
    {
        zx = xk - xkl;
        zy = yk - ykl;
        xkl = xk;
        ykl = yk;
        iter++;
        df11 = AD11(xkl, ykl);
        df12 = AD12(xkl, ykl);
        df21 = AD21(xkl, ykl);
        df22 = AD22(xkl, ykl);
        f1 = -Function1(xkl, ykl);
        f2 = -Function2(xkl, ykl);
        // std::cout << df11 << ' ' << df12 <<' '<<df21<<' '<<df22<< std::endl;
        //det = df11 * df22 - df21 * df12;
        //std::cout << det<< std::endl;
        zx = (1 - omega) * zx - omega * (df12 * zy) / df11 + omega * f1 / df11;
        zy = -omega * df21 * zx / df22 + (1 - omega) * zy + omega * f2 / df22;
        xk = zx + xkl;
        yk = zy + ykl;
        //std::cout << xk << ' ' << yk << std::endl;
        //std::cout << iter << std::endl;
        if (iter == 60)
        {
            break;
        }
    }
    result[0] = xk;
    result[1] = yk;
}
void newtonD1Z(double x, double y, double eps, double* result)
{
    double df11, f1, f2, df12, df21, df22, xk, yk, xkl, ykl, det,omega, zx, zy;
    xkl = x;
    ykl = y;
    xk = x + 1;
    yk = y + 1;
    iter = 0;
    omega = 0.95;
    while (sqrt((xkl - xk) * (xkl - xk) + (ykl - yk) * (ykl - yk)) > eps)
    {
        zx = xk - xkl;
        zy = yk - ykl;
        xkl = xk;
        ykl = yk;
        iter++;
        df11 = D11(xkl, ykl);
        df12 = D12(xkl, ykl);
        df21 = D21(xkl, ykl);
        df22 = D22(xkl, ykl);
        f1 = -Function1(xkl, ykl);
        f2 = -Function2(xkl, ykl);
        //std::cout << df11 << ' ' << df12 << ' ' << df21 << ' ' << df22 << std::endl;
        //std::cout << det << std::endl;
        zx = (1 - omega) * zx - omega * (df12 * zy) / df11 + omega * f1 / df11;
        zy = -omega * df21 * zx / df22 + (1 - omega) * zy + omega * f2 / df22;
        xk = zx + xkl;
        yk = zy + ykl;
        //std::cout << xk << ' ' << yk << std::endl;
        //std::cout << iter << std::endl;
        if (iter == 60)
        {
            break;
        }
    }
    result[0] = xk;
    result[1] = yk;
}
void test2()
{
    std::ofstream ans1,ans2;
    ans1.open("newtonD200.txt");
    ans2.open("newtonAD200.txt");
    double L1, L2,f1,f2;
    double eps = 0.000001;
    double result[3][2] = { {0.83022,1.581927},{4.332158,0.378667},{3.487259,0.524446} };
    int n,indr;
    L1 = 10;
    L2 = 10;
    n = 200;
    std::cout << std::setprecision(15);
    ans1 << std::setprecision(15);
    ans2 << std::setprecision(15);
    std::vector<double> p;
    double* res;
    res = new double[2];
    p = setka(L1, L2, n);
    for (int i=0; i < p.size(); i++)
    {
        newtonD(p[i], p[i+1], eps, res);
        f1 = Function1(res[0], res[1]);
        f2 = Function2(res[0], res[1]);
        indr = 0;
        if (abs(result[0][0] - res[0]) + abs(result[0][1] - res[1]) < 0.0001)
        {
            indr = 1;
        }
        if (abs(result[1][0] - res[0]) + abs(result[1][1] - res[1]) < 0.0001)
        {
            indr = 2;
        }
        if (abs(result[2][0] - res[0]) + abs(result[2][1] - res[1]) < 0.0001)
        {
            indr = 3;
        }
       //std::cout << "newtonD " << p[i] << ' ' << p[i + 1] << ' ' << "result " <<  res[0] << ' ' << res[1]<<' ' << "iter " << iter << std::endl;
        ans1 << p[i] << ' ' << p[i + 1] << ' ' << iter <<' '<<sqrt(f1*f1+f2*f2)<<' ' << res[0] << ' ' << res[1] <<' ' <<indr<< std::endl;
        newtonAD(p[i], p[i + 1], eps, res);
        f1 = Function1(res[0], res[1]);
        f2 = Function2(res[0], res[1]);
        indr = 0;
        if (abs(result[0][0] - res[0]) + abs(result[0][1] - res[1]) < 0.0001)
        {
            indr = 1;
        }
        if (abs(result[1][0] - res[0]) + abs(result[1][1] - res[1]) < 0.0001)
        {
            indr = 2;
        }
        if (abs(result[2][0] - res[0]) + abs(result[2][1] - res[1]) < 0.0001)
        {
            indr = 3;
        }
        //std::cout << "newtonAD " << p[i] << ' ' << p[i + 1] << ' ' << "result " << res[0] << ' ' << res[1] << ' ' << "iter " << iter << std::endl;
        ans2 << p[i] << ' ' << p[i + 1] << ' ' << iter << ' ' << sqrt(f1 * f1 + f2 * f2) << ' ' << res[0] << ' ' << res[1] << ' ' << indr << std::endl;
        i++;
       // std::cout << std::endl;
       // std::cout << std::endl;
    }
    
}
void test3()
{
    std::ofstream ans1, ans2;
    ans1.open("newtonD2001Z.txt");
    ans2.open("newtonAD2001Z.txt");
    double L1, L2, f1, f2;
    double eps = 0.000001;
    double result[3][2] = { {0.83022,1.581927},{4.332158,0.378667},{3.487259,0.524446} };
    int n, indr;
    L1 = 10;
    L2 = 10;
    n = 200;
    std::cout << std::setprecision(15);
    ans1 << std::setprecision(15);
    ans2 << std::setprecision(15);
    std::vector<double> p;
    double* res;
    res = new double[2];
    p = setka(L1, L2, n);
    for (int i = 0; i < p.size(); i++)
    {
        newtonD1Z(p[i], p[i + 1], eps, res);
        f1 = Function1(res[0], res[1]);
        f2 = Function2(res[0], res[1]);
        indr = 0;
        if (abs(result[0][0] - res[0]) + abs(result[0][1] - res[1]) < 0.001)
        {
            indr = 1;
        }
        if (abs(result[1][0] - res[0]) + abs(result[1][1] - res[1]) < 0.0001)
        {
            indr = 2;
        }
        if (abs(result[2][0] - res[0]) + abs(result[2][1] - res[1]) < 0.0001)
        {
            indr = 3;
        }
        //std::cout << "newtonD " << p[i] << ' ' << p[i + 1] << ' ' << "result " <<  res[0] << ' ' << res[1]<<' ' << "iter " << iter << std::endl;
        ans1 << p[i] << ' ' << p[i + 1] << ' ' << iter << ' ' << sqrt(f1 * f1 + f2 * f2) << ' ' << res[0] << ' ' << res[1] << ' ' << indr << std::endl;
        newtonAD1Z(p[i], p[i + 1], eps, res);
        f1 = Function1(res[0], res[1]);
        f2 = Function2(res[0], res[1]);
        indr = 0;
        if (abs(result[0][0] - res[0]) + abs(result[0][1] - res[1]) < 0.0001)
        {
            indr = 1;
        }
        if (abs(result[1][0] - res[0]) + abs(result[1][1] - res[1]) < 0.0001)
        {
            indr = 2;
        }
        if (abs(result[2][0] - res[0]) + abs(result[2][1] - res[1]) < 0.001)
        {
            indr = 3;
        }
        //std::cout << "newtonAD " << p[i] << ' ' << p[i + 1] << ' ' << "result " << res[0] << ' ' << res[1] << ' ' << "iter " << iter << std::endl;
        ans2 << p[i] << ' ' << p[i + 1] << ' ' << iter << ' ' << sqrt(f1 * f1 + f2 * f2) << ' ' << res[0] << ' ' << res[1] << ' ' << indr << std::endl;
        i++;
        //std::cout << std::endl;
       // std::cout << std::endl;
    }

}
int main()
{
    //test3();
    //test2();
    test1();
    std::cout << "Hello World!\n";

}

