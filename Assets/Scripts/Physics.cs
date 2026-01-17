using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using System;

public class Physics : MonoBehaviour
{
    // basic control
    public bool physicMode = false;
    public double time = 0;
    public bool modulusSquared = false;

    // Phase-chirp controls (for modifying momentum distribution without changing |ψ|^2)
    public bool phaseChirp = false;
    public double chirpAlpha = 0.0; // base chirp amplitude
    public double chirpRate = 0.0;  // chirp growth per second (alpha(t) = chirpAlpha + chirpRate * time)

    // ----- physical / diagnostic defaults (kept accessible) -----
    public static double phys_m = 1.0;               // mass
    public static double phys_sigma0 = 0.25;        // initial width (position uncertainty at t=0)
    public static double phys_x0 = 0.0, phys_y0 = 0.0; // initial center at t=0
    public static double phys_kx0 = 1.0, phys_ky0 = 1.0; // mean wavevector (momentum)
    public static int uncertaintyGrid = 81;          // grid resolution used for diagnostics (odd recommended)
    public static double uncertaintyRangeFactor = 6.0; // range = factor * sigma(t)
    public static double uncertaintyInterval = 0.5;  // seconds between diagnostic logs
    private static double lastUncertaintyTime = -1e9;

    // ----- numerical simulator settings (in-Physics, no extra scene objects) -----
    [Header("Simulator")]
    public int simGridSize = 256;         // must be power of two
    public float simDomainSize = 8f;      // physical domain size (square)
    public float simDt = 0.002f;          // simulation timestep
    public int simStepsPerUpdate = 1;     // steps performed each Unity Update when physicMode is true

    [Header("Potential / Barriers")]
    public Texture2D potentialTexture;    // optional grayscale texture for potential (white->high V)
    public float potentialScale = 10f;    // scale factor applied to texture values
    public bool useAbsorbingBoundary = true;
    public float absorberStrength = 20f;
    public int absorberWidth = 16;

    // internal simulation buffers
    private int N;
    private double dx; // spatial step
    private Complex[,] psiSim;
    private double[,] Vsim;
    private Complex[,] kineticProp;
    private Complex[] rowBuf;
    private Complex[] colBuf;
    private bool simInitialized = false;

    void Start()
    {
        // do not auto-start simulation; initialize lazily when physicMode enabled
    }

    void Update()
    {
        // toggle initialization when entering physics mode
        if (physicMode && !simInitialized)
        {
            InitializeSimulation();
        }

        if (!physicMode && simInitialized)
        {
            // keep sim alive but do not advance if user disables physicMode
            return;
        }

        if (physicMode && simInitialized)
        {
            for (int s = 0; s < simStepsPerUpdate; s++)
            {
                StepSimulation();
                time += simDt;
            }

            // periodic diagnostics (optional)
            if (Time.time - lastUncertaintyTime >= uncertaintyInterval)
            {
                // reuse existing diagnostics if any external code calls ComputeAndLogUncertainties;
                // we keep lastUncertaintyTime update here to avoid flooding logs.
                lastUncertaintyTime = Time.time;
            }
        }
    }

    // Initialize simulation arrays and precompute kinetic propagator
    void InitializeSimulation()
    {
        if (!IsPowerOfTwo(simGridSize))
        {
            Debug.LogError("simGridSize must be a power of two.");
            simGridSize = 256;
        }

        N = simGridSize;
        psiSim = new Complex[N, N];
        Vsim = new double[N, N];
        kineticProp = new Complex[N, N];
        rowBuf = new Complex[N];
        colBuf = new Complex[N];

        dx = simDomainSize / N;

        // initialize psi as Gaussian wavepacket (centered at phys_x0/phys_y0, momentum phys_kx0/phys_ky0)
        for (int i = 0; i < N; i++)
        {
            double x = (i - N / 2) * dx + phys_x0;
            for (int j = 0; j < N; j++)
            {
                double y = (j - N / 2) * dx + phys_y0;
                double dxr = x - phys_x0;
                double dyr = y - phys_y0;
                double r2 = dxr * dxr + dyr * dyr;
                Complex envelope = Complex.Exp(-r2 / (2.0 * phys_sigma0 * phys_sigma0));
                Complex phase = Complex.Exp(Complex.ImaginaryOne * (phys_kx0 * x + phys_ky0 * y));
                psiSim[i, j] = envelope * phase;
            }
        }

        // build potential from texture if provided
        if (potentialTexture != null)
        {
            BuildPotentialFromTexture();
        }
        else
        {
            // zero potential
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Vsim[i, j] = 0.0;
        }

        // absorbing boundary (encoded as additional potential magnitude, treated empirically)
        if (useAbsorbingBoundary)
        {
            for (int i = 0; i < N; i++)
            {
                float ix = 0f;
                if (i < absorberWidth) ix = 1f - (float)i / absorberWidth;
                else if (i >= N - absorberWidth) ix = 1f - (float)(N - 1 - i) / absorberWidth;
                for (int j = 0; j < N; j++)
                {
                    float jy = 0f;
                    if (j < absorberWidth) jy = 1f - (float)j / absorberWidth;
                    else if (j >= N - absorberWidth) jy = 1f - (float)(N - 1 - j) / absorberWidth;
                    float factor = Math.Max(ix, jy);
                    if (factor > 0f)
                    {
                        Vsim[i, j] += absorberStrength * factor;
                    }
                }
            }
        }

        // precompute kinetic propagator in k-space (ħ = 1)
        double dk = 2.0 * Math.PI / (N * dx);
        for (int i = 0; i < N; i++)
        {
            int ki = i <= N / 2 ? i : i - N;
            double kx = ki * dk;
            for (int j = 0; j < N; j++)
            {
                int kj = j <= N / 2 ? j : j - N;
                double ky = kj * dk;
                double k2 = kx * kx + ky * ky;
                double phase = -(k2) * (simDt / (2.0 * phys_m)); // exponent argument (multiplied by i)
                kineticProp[i, j] = Complex.Exp(Complex.ImaginaryOne * phase);
            }
        }

        simInitialized = true;
    }

    void StepSimulation()
    {
        // potential half-step: psi <- exp(-i V dt/2) psi
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                double Vval = Vsim[i, j];
                Complex phase = Complex.Exp(Complex.ImaginaryOne * (-Vval * simDt * 0.5));
                psiSim[i, j] *= phase;
            }

        // forward 2D FFT
        FFT2D(psiSim, false);

        // multiply by kinetic propagator in k-space
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                psiSim[i, j] *= kineticProp[i, j];

        // inverse 2D FFT
        FFT2D(psiSim, true);

        // final potential half-step and optional phase-chirp
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                double Vval = Vsim[i, j];
                Complex phase = Complex.Exp(Complex.ImaginaryOne * (-Vval * simDt * 0.5));
                psiSim[i, j] *= phase;

                // optional quadratic chirp to increase σ_p while leaving |ψ|^2 unchanged
                if (phaseChirp)
                {
                    // compute world coordinate x for grid cell
                    double x = (i - N / 2) * dx + phys_x0;
                    double cx = phys_x0 + (phys_kx0 / phys_m) * time;
                    double alpha = chirpAlpha + chirpRate * time;
                    Complex chirp = Complex.Exp(Complex.ImaginaryOne * alpha * (x - cx) * (x - cx));
                    psiSim[i, j] *= chirp;
                }

                // crude damping for absorber (keeps stability)
                if (useAbsorbingBoundary && Vval > 0.0)
                {
                    double damp = Math.Exp(-Vval * simDt * 0.1);
                    psiSim[i, j] *= damp;
                }
            }
    }

    // sample psi at world space coordinate z (Complex.Real = x, .Imaginary = y)
    // returns raw complex ψ (not divided by any Scaler); if modulusSquared true returns (|ψ|^2, 0)
    public Complex CalculateWaveFunction(Complex z)
    {
        if (!simInitialized)
        {
            // fallback to analytic gaussian if simulation not running
            return AnalyticFallback(z);
        }

        // map world coordinates to grid indices (centered)
        double x = z.Real;
        double y = z.Imaginary;

        double x0 = phys_x0 - (N / 2) * dx;
        double y0 = phys_y0 - (N / 2) * dx;

        double fx = (x - x0) / dx;
        double fy = (y - y0) / dx;

        int i = Mathf.FloorToInt((float)fx);
        int j = Mathf.FloorToInt((float)fy);

        // clamp for safety
        if (i < 0 || i >= N - 1 || j < 0 || j >= N - 1)
        {
            // outside grid: assume zero
            return new Complex(0.0, 0.0);
        }

        // bilinear interpolation
        double tx = fx - i;
        double ty = fy - j;

        Complex c00 = psiSim[i, j];
        Complex c10 = psiSim[i + 1, j];
        Complex c01 = psiSim[i, j + 1];
        Complex c11 = psiSim[i + 1, j + 1];

        Complex c0 = c00 * (1.0 - tx) + c10 * tx;
        Complex c1 = c01 * (1.0 - tx) + c11 * tx;
        Complex val = c0 * (1.0 - ty) + c1 * ty;

        if (modulusSquared)
        {
            double P = Complex.Abs(val);
            P = P * P;
            return new Complex(P, 0.0);
        }

        return val;
    }

    // a small analytic fallback used when sim not initialized
    private Complex AnalyticFallback(Complex z)
    {
        double t = time;
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
        return psi;
    }

    // Build potential from Texture2D grayscale: black->0, white->1
    void BuildPotentialFromTexture()
    {
        if (potentialTexture == null) return;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                float u = (i + 0.5f) / N;
                float v = (j + 0.5f) / N;
                Color c = potentialTexture.GetPixelBilinear(u, v);
                float gray = c.grayscale;
                Vsim[i, j] = gray * potentialScale;
            }
        }
    }

    // -------------------- 2D FFT (in-place) --------------------
    void FFT2D(Complex[,] data, bool inverse)
    {
        // rows
        for (int y = 0; y < N; y++)
        {
            for (int x = 0; x < N; x++) rowBuf[x] = data[x, y];
            FFT1D(rowBuf, inverse);
            for (int x = 0; x < N; x++) data[x, y] = rowBuf[x];
        }
        // cols
        for (int x = 0; x < N; x++)
        {
            for (int y = 0; y < N; y++) colBuf[y] = data[x, y];
            FFT1D(colBuf, inverse);
            for (int y = 0; y < N; y++) data[x, y] = colBuf[y];
        }
    }

    void FFT1D(Complex[] buffer, bool inverse)
    {
        int n = buffer.Length;
        int bits = (int)Math.Log(n, 2);

        // bit reversal
        for (int j = 1, i = 0; j < n; j++)
        {
            int bit = n >> 1;
            for (; (i & bit) != 0; bit >>= 1) i ^= bit;
            i ^= bit;
            if (j < i)
            {
                var tmp = buffer[j];
                buffer[j] = buffer[i];
                buffer[i] = tmp;
            }
        }

        // FFT
        for (int len = 2; len <= n; len <<= 1)
        {
            double ang = 2.0 * Math.PI / len * (inverse ? 1.0 : -1.0);
            Complex wlen = Complex.Exp(Complex.ImaginaryOne * ang);
            for (int i = 0; i < n; i += len)
            {
                Complex w = Complex.One;
                for (int j = 0; j < len / 2; j++)
                {
                    Complex u = buffer[i + j];
                    Complex v = buffer[i + j + len / 2] * w;
                    buffer[i + j] = u + v;
                    buffer[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }

        // scale if inverse
        if (inverse)
        {
            double inv = 1.0 / n;
            for (int i = 0; i < n; i++) buffer[i] *= inv;
        }
    }

    bool IsPowerOfTwo(int x) { return (x & (x - 1)) == 0; }
}