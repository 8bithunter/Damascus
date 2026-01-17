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

    // ---- physics wavepacket parameters and diagnostic settings ----
    public static double phys_m = 1.0;               // mass
    public static double phys_sigma0 = 0.25;        // initial width (position uncertainty at t=0)
    public static double phys_x0 = 0.0, phys_y0 = 0.0; // initial center at t=0
    public static double phys_kx0 = 1.0, phys_ky0 = 1.0; // mean wavevector (momentum)
    public static int uncertaintyGrid = 81;          // grid resolution used for diagnostics (odd recommended)
    public static double uncertaintyRangeFactor = 6.0; // range = factor * sigma(t)
    public static double uncertaintyInterval = 0.5;  // seconds between diagnostic logs
    private static double lastUncertaintyTime = -1e9;


    public void Start()
    {
        matrix = GetComponent<Matrix>();
        physics = GetComponent<Physics>();
    }
    private void Update()
    {
        integrationStart = new Complex(starti.position.x, starti.position.y);

        // run diagnostics periodically while physics mode is active
        if (physics != null && physics.physicMode)
        {
            if (Time.time - lastUncertaintyTime >= uncertaintyInterval)
            {
                ComputeAndLogUncertainties(physics.time);
                lastUncertaintyTime = Time.time;
            }
        }
    }

    public static Complex Function(Complex unscaledz)
    {
        Complex z = unscaledz * Scaler.scale;
        Complex output;
        if (physics != null && physics.physicMode)
        {
            // Moving 2D Gaussian wavepacket: known position at t=0, speed uncertain.
            double t = physics.time;

            // use the shared phys_* parameters above
            double m = phys_m;
            double sigma0 = phys_sigma0;
            double x0 = phys_x0, y0 = phys_y0;
            double kx0 = phys_kx0, ky0 = phys_ky0;

            // Derived
            double omega = (kx0 * kx0 + ky0 * ky0) / (2.0 * m); // dispersion relation
            double vx = kx0 / m;
            double vy = ky0 / m;

            // Translate center by group velocity
            double cx = x0 + vx * t;
            double cy = y0 + vy * t;

            double dxr = z.Real - cx;
            double dyr = z.Imaginary - cy;

            // Time-dependent complex width: a = sigma0^2 + i * t / (2 m)
            Complex a = new Complex(sigma0 * sigma0, t / (2.0 * m));

            // Exponent (2D separable Gaussian centered at cx,cy)
            Complex exponent = -(dxr * dxr + dyr * dyr) / (4.0 * a);

            // Proper prefactor for separable 2D Gaussian (visualization; we'll renormalize in diagnostics)
            // For separable 1D factors, prefactor_2D = 1/(2*pi*a)
            Complex prefactor = 1.0 / (2.0 * Math.PI * a);

            // Phase including spatial plane-wave and temporal -omega t
            Complex planePhase = Complex.Exp(Complex.ImaginaryOne * (kx0 * (z.Real - cx) + ky0 * (z.Imaginary - cy) - omega * t));

            Complex psi = prefactor * Complex.Exp(exponent) * planePhase;

            // optional phase-chirp: unit-modulus quadratic phase in x (preserves |ψ|^2 -> σ_x unchanged,
            // but increases momentum spread σ_p). alpha can be time-dependent via physics.chirpRate.
            if (physics.phaseChirp)
            {
                double alpha = physics.chirpAlpha + physics.chirpRate * physics.time; // choose sign/magnitude as needed
                double dxCenter = z.Real - cx; // use coordinate relative to packet center
                // quadratic chirp only in x-direction; you can add y term similarly if desired
                Complex chirp = Complex.Exp(Complex.ImaginaryOne * alpha * dxCenter * dxCenter);
                psi *= chirp;
            }

            if (physics.modulusSquared)
            {
                double prob = Complex.Abs(psi);
                prob = prob * prob;             // |ψ|^2
                output = new Complex(prob, 0.0);
            }
            else
            {
                output = psi;
            }
        }
        else if (matrix != null && matrix.matrixMode)
        {
            output = matrix.Transform(z);
        }
        else

            //Change the right side of the following assignment to your desired function using "z" as your variable
            //eg. Complex output = Complex.Sin(z) + Complex.Pow(z, 3) + z;
            output = Mandelbrot(z);

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

    public static Complex Zeta(Complex z)
    {
        if (z.Real >= 0.5)
        {
            Complex output = Complex.Zero;
            for (int i = 1; i < 1000; i++)
            {
                output = output + 1 / Complex.Pow(i, z);
            }
            return output;
        }
        else
        {
            Complex output = Complex.Pow(2, z) * Complex.Pow(Math.PI, z - 1);
            return (output * Complex.Sin((Math.PI * z) / 2) * Gamma(1 - z) * Zeta(1 - z));
        }
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

    public static void ComputeAndLogUncertainties(double t)
    {
        int N = Math.Max(5, uncertaintyGrid | 1); // ensure odd and >=5 for 4th-order stencil
        double m = phys_m;
        double sigma0 = phys_sigma0;
        double x0 = phys_x0, y0 = phys_y0;
        double kx0 = phys_kx0, ky0 = phys_ky0;

        // analytic instantaneous width (1D) for guidance (ħ = 1)
        double sigma_t = sigma0 * Math.Sqrt(1.0 + Math.Pow(t / (2.0 * m * sigma0 * sigma0), 2.0));

        double cx = x0 + (kx0 / m) * t;
        double cy = y0 + (ky0 / m) * t;

        double range = uncertaintyRangeFactor * Math.Max(sigma_t, sigma0);
        double dxGrid = (2.0 * range) / (N - 1);

        // allocate psi grid
        Complex[,] psi = new Complex[N, N];
        double[,] prob = new double[N, N];

        double sumP = 0.0;
        for (int i = 0; i < N; i++)
        {
            double x = cx - range + i * dxGrid;
            for (int j = 0; j < N; j++)
            {
                double y = cy - range + j * dxGrid;
                Complex p = PhysicalPsi(new Complex(x, y), t, false); // raw psi (no Scaler division)
                psi[i, j] = p;
                double P = Complex.Abs(p);
                P = P * P;
                prob[i, j] = P;
                sumP += P;
            }
        }

        double areaElement = dxGrid * dxGrid;
        double norm = sumP * areaElement;
        if (norm <= 0) { Debug.Log("ComputeAndLogUncertainties: normalization zero."); return; }

        // compute position moments
        double sumX = 0.0, sumX2 = 0.0;
        for (int i = 0; i < N; i++)
        {
            double x = cx - range + i * dxGrid;
            for (int j = 0; j < N; j++)
            {
                double P = prob[i, j];
                sumX += x * P;
                sumX2 += x * x * P;
            }
        }
        double meanX = (sumX * areaElement) / norm;
        double meanX2 = (sumX2 * areaElement) / norm;
        double sigmaX = Math.Sqrt(Math.Max(0.0, meanX2 - meanX * meanX));

        // compute momentum moments along x using 4th-order finite differences (more accurate & stable)
        Complex integral_px = Complex.Zero;
        Complex integral_px2 = Complex.Zero;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // high-order derivative in x (4th-order central where possible)
                Complex dpsi_dx;
                if (i >= 2 && i <= N - 3)
                {
                    dpsi_dx = (-psi[i + 2, j] + 8.0 * psi[i + 1, j] - 8.0 * psi[i - 1, j] + psi[i - 2, j]) / (12.0 * dxGrid);
                }
                else if (i == 1)
                {
                    // fallback to 3-point forward
                    dpsi_dx = (-3.0 * psi[i, j] + 4.0 * psi[i + 1, j] - psi[i + 2, j]) / (2.0 * dxGrid);
                }
                else if (i == N - 2)
                {
                    dpsi_dx = (3.0 * psi[i, j] - 4.0 * psi[i - 1, j] + psi[i - 2, j]) / (2.0 * dxGrid);
                }
                else // i == 0 or i == N-1
                {
                    // one-sided 2-point
                    if (i == 0) dpsi_dx = (psi[i + 1, j] - psi[i, j]) / dxGrid;
                    else dpsi_dx = (psi[i, j] - psi[i - 1, j]) / dxGrid;
                }

                // second derivative using 4th-order central where possible
                Complex d2psi_dx2;
                if (i >= 2 && i <= N - 3)
                {
                    d2psi_dx2 = (-psi[i + 2, j] + 16.0 * psi[i + 1, j] - 30.0 * psi[i, j] + 16.0 * psi[i - 1, j] - psi[i - 2, j]) / (12.0 * dxGrid * dxGrid);
                }
                else if (i == 1)
                {
                    d2psi_dx2 = (psi[i + 2, j] - 2.0 * psi[i + 1, j] + psi[i, j]) / (dxGrid * dxGrid);
                }
                else if (i == N - 2)
                {
                    d2psi_dx2 = (psi[i, j] - 2.0 * psi[i - 1, j] + psi[i - 2, j]) / (dxGrid * dxGrid);
                }
                else
                {
                    d2psi_dx2 = (psi[Math.Min(i + 1, N - 1), j] - 2.0 * psi[i, j] + psi[Math.Max(i - 1, 0), j]) / (dxGrid * dxGrid);
                }

                Complex conjPsi = Complex.Conjugate(psi[i, j]);

                integral_px += conjPsi * dpsi_dx;
                integral_px2 += conjPsi * d2psi_dx2;
            }
        }

        Complex pxComplex = -Complex.ImaginaryOne * integral_px * areaElement / norm;
        double meanPx = pxComplex.Real;

        Complex px2Complex = -1.0 * integral_px2 * areaElement / norm;
        double meanPx2 = px2Complex.Real;

        double sigmaPxSq = meanPx2 - meanPx * meanPx;
        double sigmaPx = sigmaPxSq > 0 ? Math.Sqrt(sigmaPxSq) : 0.0;

        // analytic momentum width for initial Gaussian (ħ = 1): sigma_p_analytic = 1/(2 sigma0)
        double sigmaPAnalytic = 1.0 / (2.0 * sigma0);

        Debug.Log(string.Format(
            "Quantum diagnostics t={0:F3}s: ⟨x⟩={1:F6}, σ_x={2:F6}, ⟨p_x⟩={3:F6}, σ_p_x(numeric)={4:F6}, σ_p_x(analytic)={5:F6}, σ_xσ_p(numeric)={6:F6}",
            t, meanX, sigmaX, meanPx, sigmaPx, sigmaPAnalytic, sigmaX * sigmaPx));
    }

    private static Complex PhysicalPsi(Complex z, double t, bool modulusOnly)
    {
        double m = phys_m;
        double sigma0 = phys_sigma0;
        double x0 = phys_x0, y0 = phys_y0;
        double kx0 = phys_kx0, ky0 = phys_ky0;

        double vx = kx0 / m;
        double vy = ky0 / m;

        double cx = x0 + vx * t;
        double cy = y0 + vy * t;

        double dxr = z.Real - cx;
        double dyr = z.Imaginary - cy;

        Complex a = new Complex(sigma0 * sigma0, t / (2.0 * m));
        Complex exponent = -(dxr * dxr + dyr * dyr) / (4.0 * a);
        Complex prefactor = 1.0 / (2.0 * Math.PI * a);
        double omega = (kx0 * kx0 + ky0 * ky0) / (2.0 * m);
        Complex planePhase = Complex.Exp(Complex.ImaginaryOne * (kx0 * (z.Real - cx) + ky0 * (z.Imaginary - cy) - omega * t));

        Complex psi = prefactor * Complex.Exp(exponent) * planePhase;
        if (modulusOnly)
        {
            double P = Complex.Abs(psi);
            P = P * P;
            return new Complex(P, 0.0);
        }
        return psi;
    }
}

