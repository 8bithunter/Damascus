using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using System.Net;
using System;
using UnityEngine.UI;

public class Funcs : MonoBehaviour
{
    public static double dx = 0.00000001;
    public static int integrationResolution = 1000;
    public Transform starti;

    static Complex integrationStart;

    public static Matrix matrix;
    public static Physics physics;

    public void Start()
    {
        matrix = GetComponent<Matrix>();
        physics = GetComponent<Physics>();
    }
    private void Update()
    {
        integrationStart = new Complex(starti.position.x, starti.position.y);
    }

    public static Complex Function(Complex unscaledz)
    {
        Complex z = unscaledz * Scaler.scale;
        Complex output;

        if (physics != null && physics.physicMode)
        {
            output = physics.CalculateWaveFunction(z);
        }
        else if (matrix != null && matrix.matrixMode)
        {
            output = matrix.Transform(z);
        }
        else
        { 
            //Change the right side of the following assignment to your desired function using "z" as your variable
            //eg. Complex output = Complex.Sin(z) + Complex.Pow(z, 3) + z;
            output = Zeta(z);
        }

        return output / Scaler.scale;
    }

    /*public static Complex Derivative(Complex z)
    {
        Complex fz = Function(z);
        Complex fz_plus_h = Function(z + dx);
        Complex fz_minus_h = Function(z - dx);

        double df_dx = (fz_plus_h.Real - fz_minus_h.Real) / (2 * dx);
        double df_dy = (fz_plus_h.Imaginary - fz_minus_h.Imaginary) / (2 * dx);

        return new Complex(df_dx, df_dy) / Scaler.scale;
    }*/

    public static Complex Derivative(Complex z)
    {
        Complex fz_plus_h = Function(z + new Complex(dx, 0));
        Complex fz_minus_h = Function(z - new Complex(dx, 0));

        return ((fz_plus_h - fz_minus_h) / (2 * dx)) / Scaler.scale;
    }


    public static Complex RiemannSum(Complex endPoint)
    {
        Complex antiResult = Complex.Zero;

        Complex deltaZ = (endPoint - integrationStart) / integrationResolution;

        for (int i = 0; i < integrationResolution; i++)
        {
            Complex z = integrationStart + deltaZ * i;
            Complex f_z = Function(z);
            antiResult += f_z * deltaZ;
        }
        return (antiResult * Scaler.scale);
    }

    public static Complex SimpsonsRule(Complex end)
    {
        Complex h = (end - integrationStart) / integrationResolution;
        Complex result = Function(integrationStart) + Function(end);

        for (int i = 1; i < integrationResolution; i += 2)
        {
            result += 4 * Function(integrationStart + i * h);
        }

        for (int i = 2; i < integrationResolution - 1; i += 2)
        {
            result += 2 * Function(integrationStart + i * h);
        }

        return (result * h / 3.0) * Scaler.scale;
    }

    public static Complex Mandelbrot(Complex z)
    {
        Complex output = Complex.Zero;
        for (int i = 0; i < 100; i++)
        {
            output = Complex.Pow(output, 2) + z;
        }
        return output;
    }
    public static Complex Zeta(Complex s)
    {
        const double eps = 1e-12;
        int N = 400;

        // ---------- helpers ----------
        Complex ZetaDirichlet(Complex z)
        {
            Complex sum = Complex.Zero;
            for (int n = 1; n <= N; n++)
                sum += Complex.Pow(n, -z);

            // crude Euler–Maclaurin tail correction
            Complex nN = new Complex(N, 0);
            sum += Complex.Pow(nN, 1 - z) / (z - 1);
            sum += 0.5 * Complex.Pow(nN, -z);
            return sum;
        }

        Complex ZetaEta(Complex z)
        {
            Complex sum = Complex.Zero;
            for (int n = 1; n <= N; n++)
            {
                double sign = (n & 1) == 1 ? 1.0 : -1.0;
                sum += sign * Complex.Pow(n, -z);
            }
            return sum;
        }

        // ---------- special values ----------
        if (s == Complex.Zero) return new Complex(-0.5, 0);
        if (s == Complex.One) return new Complex(double.PositiveInfinity, 0);

        // ---------- Re(s) > 1 ----------
        if (s.Real > 1 + eps)
            return ZetaDirichlet(s);

        // ---------- 0 < Re(s) ≤ 1 ----------
        if (s.Real > eps)
        {
            Complex eta = ZetaEta(s);
            Complex denom = Complex.One - Complex.Pow(2, 1 - s);

            if (Complex.Abs(denom) < 1e-10)
                throw new ArithmeticException("ζ(s) unstable near pole/zero of continuation");

            return eta / denom;
        }

        // ---------- Re(s) ≤ 0 : functional equation ----------
        Complex factor =
            Complex.Pow(2, s) *
            Complex.Pow(Math.PI, s - 1) *
            Complex.Sin(Math.PI * s / 2) *
            Gamma(1 - s);

        Complex t = Complex.One - s;

        // t now has Re(t) ≥ 1
        return factor * Zeta(t);
    }


    static int g = 7;
    static double[] p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

    public static Complex Gamma(Complex z)
    {
        if (z.Real < 0.5)
        {
            return Math.PI / (Complex.Sin(Math.PI * z) * Gamma(1 - z));
        }
        else
        {
            z -= 1;
            Complex x = p[0];
            for (var i = 1; i < g + 2; i++)
            {
                x += p[i] / (z + i);
            }
            Complex t = z + g + 0.5;
            return Complex.Sqrt(2 * Math.PI) * (Complex.Pow(t, z + 0.5)) * Complex.Exp(-t) * x;
        }
    }

    public static Complex LambertW(Complex z, int maxIterations = 100, double tolerance = 1e-10)
    {
        if (z == Complex.Zero) return Complex.Zero;

        Complex w = z.Magnitude < 1 ? z : Complex.Log(z);

        for (int i = 0; i < maxIterations; i++)
        {
            Complex ew = Complex.Exp(w);
            Complex wew = w * ew;
            Complex deltaW = (wew - z) / (ew * (w + 1) - (w + 2) * (wew - z) / (2 * w + 2));
            w -= deltaW;

            if (Complex.Abs(deltaW) < tolerance)
                break;
        }

        return w;
    }

    public static Complex CreateSymmetry(Complex z, Complex symmetryNumber)
    {
        return (Complex.Pow(z, symmetryNumber)) / Complex.Abs(Complex.Pow(z, symmetryNumber - 1));
    }
}

