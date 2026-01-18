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
    public static double phys_x0 = 0.0, phys_y0 = 0.0; // initial center at t=0 (world units)
    public static double phys_kx0 = 1.0, phys_ky0 = 1.0; // mean wavevector (momentum)
    public static int uncertaintyGrid = 81;          // grid resolution used for diagnostics (odd recommended)
    public static double uncertaintyRangeFactor = 6.0; // range = factor * sigma(t)
    public static double uncertaintyInterval = 0.5;  // seconds between diagnostic logs
    private static double lastUncertaintyTime = -1e9;

    // ----- numerical simulator settings (in-Physics, no extra scene objects) -----
    [Header("Simulator")]
    public int simGridSize = 256;         // must be power of two
    public float simDomainSize = 8f;      // physical domain size (square) in world units (unscaled)
    public float simDt = 0.002f;          // simulation timestep
    public int simStepsPerUpdate = 1;     // steps performed each Unity Update when physicMode is true

    // If true, normalize ψ after each time step (overrides physical loss by absorbers)
    public bool autoNormalize = false;

    [Header("Potential / Barriers")]
    public Texture2D potentialTexture;    // optional grayscale texture for potential (white->high V)
    public float potentialScale = 10f;    // scale factor applied to texture values
    public bool useAbsorbingBoundary = true;
    public float absorberStrength = 20f;
    public int absorberWidth = 16;

    // Lasso tool settings (draw polygon with left mouse, release to apply)
    [Header("Lasso Tool")]
    public bool enableLasso = true;                   // enable/disable lasso input
    public float lassoMinSegment = 0.05f;             // minimum world distance between recorded points
    public float lassoPotential = 1000.0f;            // potential applied inside polygon
    public bool lassoReplaceExisting = true;          // replace or add to existing potential
    public KeyCode lassoToggleKey = KeyCode.L;        // optional key to toggle lasso mode

    // internal lasso state
    private List<UnityEngine.Vector2> lassoPoints = new List<UnityEngine.Vector2>();
    private bool isDrawingLasso = false;

    // internal simulation buffers
    private int N;
    private double dx; // spatial step (scaled)
    private double originX; // scaled world x coordinate of grid cell [0,0]
    private double originY; // scaled world y coordinate of grid cell [0,0]
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
        // toggle lasso on/off
        if (Input.GetKeyDown(lassoToggleKey))
        {
            enableLasso = !enableLasso;
            Debug.Log("Lasso " + (enableLasso ? "enabled" : "disabled"));
        }

        // reset time and reinitialize wavefunction when pressing R
        if (Input.GetKeyDown(KeyCode.R))
        {
            ResetSimulationTime();
        }

        // observe / collapse on T
        if (Input.GetKeyDown(KeyCode.T))
        {
            ObserveCollapseToPoint();
        }

        // lasso input handling (independent of sim state)
        if (enableLasso)
        {
            HandleLassoInput();
        }

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

                // optional auto-normalize each step (not physical if absorbers are meant to remove probability)
                if (autoNormalize)
                {
                    NormalizePsi();
                }
            }

            // periodic diagnostics (optional)
            if (Time.time - lastUncertaintyTime >= uncertaintyInterval)
            {
                lastUncertaintyTime = Time.time;
            }
        }
    }

    // Observe: collapse the wavefunction to a narrow Gaussian centered on a sampled cell
    // instead of a delta. This recreates the "starting particle" shape at the measured location
    // and samples a post-measurement mean momentum from the momentum distribution.
    void ObserveCollapseToPoint()
    {
        if (!simInitialized)
        {
            InitializeSimulation();
        }
        if (psiSim == null) return;

        // compute weights = |psi|^2
        double total = 0.0;
        double[,] weights = new double[N, N];
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                double mag = Complex.Abs(psiSim[i, j]);
                double w = mag * mag;
                weights[i, j] = w;
                total += w;
            }
        }

        if (total <= 0.0)
        {
            Debug.LogWarning("ObserveCollapse: total probability zero, cannot sample.");
            return;
        }

        // sample a cell by cumulative distribution
        double r = UnityEngine.Random.value;
        double accum = 0.0;
        int selI = 0, selJ = 0;
        bool found = false;
        for (int i = 0; i < N && !found; i++)
        {
            for (int j = 0; j < N; j++)
            {
                accum += weights[i, j] / total;
                if (r <= accum)
                {
                    selI = i; selJ = j;
                    found = true;
                    break;
                }
            }
        }
        if (!found)
        {
            selI = N - 1;
            selJ = N - 1;
        }

        // compute chosen center in scaled world coordinates
        double scale = (double)Scaler.scale;
        double centerX_scaled = originX + selI * dx;
        double centerY_scaled = originY + selJ * dx;

        // Build a narrow Gaussian (same shape as the starting particle) centered at the sampled cell.
        // Use the same convention as InitializeSimulation: envelope uses phys_sigma0 (in world units).
        // For scaled coordinates we convert sigma to scaled units.
        double sigma_scaled = phys_sigma0 * scale;
        double twoSigma2 = 2.0 * sigma_scaled * sigma_scaled;

        // Sample post-measurement mean momentum (kx,ky) from Gaussian momentum distribution
        // analytic sigma_p = 1/(2*sigma) (ħ=1). Use phys_kx0,phys_ky0 as mean.
        double momentumStd = 1.0 / (2.0 * phys_sigma0);
        double sampledKx = SampleNormal(phys_kx0, momentumStd);
        double sampledKy = SampleNormal(phys_ky0, momentumStd);

        // Build new psiSim as normalized Gaussian with plane-wave phase from sampled momentum
        // We'll compute unnormalized values and then normalize with NormalizePsi()
        for (int i = 0; i < N; i++)
        {
            double x = originX + i * dx;
            for (int j = 0; j < N; j++)
            {
                double y = originY + j * dx;
                double dxr = x - centerX_scaled;
                double dyr = y - centerY_scaled;
                double r2 = dxr * dxr + dyr * dyr;

                // envelope (in scaled coords)
                double envelope = Math.Exp(-r2 / twoSigma2);

                // phase: use unscaled momentum*position as in initialization (undo scale in argument)
                double phaseArg = sampledKx * (x / scale) + sampledKy * (y / scale);
                Complex phase = Complex.Exp(Complex.ImaginaryOne * phaseArg);

                psiSim[i, j] = envelope * phase;
            }
        }

        // Normalize so total probability = 1
        NormalizePsi();

        // Log chosen location (convert scaled world back to unscaled world units for convenience)
        double worldX_unscaled = centerX_scaled / scale;
        double worldY_unscaled = centerY_scaled / scale;

        Debug.LogFormat("ObserveCollapse: collapsed to Gaussian at cell [{0},{1}] -> world {2:F3},{3:F3} with sampled k=({4:F3},{5:F3})",
            selI, selJ, worldX_unscaled, worldY_unscaled, sampledKx, sampledKy);
    }

    // Box-Muller sampler for normal distribution (mean, std)
    double SampleNormal(double mean, double std)
    {
        // Use UnityEngine.Random.value for uniform RNG
        double u1 = Math.Max(1e-12, UnityEngine.Random.value);
        double u2 = Math.Max(1e-12, UnityEngine.Random.value);
        double z0 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        return mean + z0 * std;
    }

    // Reset simulation time to zero and reset psi to initial analytic packet
    void ResetSimulationTime()
    {
        time = 0.0;
        if (simInitialized)
        {
            // reinitialize psiSim (keep potentials)
            double scale = (double)Scaler.scale;
            for (int i = 0; i < N; i++)
            {
                double x = originX + i * dx;
                for (int j = 0; j < N; j++)
                {
                    double y = originY + j * dx;
                    double dxr = x - (phys_x0 * scale);
                    double dyr = y - (phys_y0 * scale);
                    double r2 = dxr * dxr + dyr * dyr;
                    Complex envelope = Complex.Exp(-r2 / (2.0 * phys_sigma0 * phys_sigma0));
                    Complex phase = Complex.Exp(Complex.ImaginaryOne * (phys_kx0 * (x / scale) + phys_ky0 * (y / scale)));
                    psiSim[i, j] = envelope * phase;
                }
            }

            // normalize initial packet to unit total probability
            NormalizePsi();

            Debug.Log("Physics simulation time reset and wavefunction reinitialized.");
        }
        else
        {
            Debug.Log("Physics time reset to zero (simulation not initialized).");
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

        // take current Scaler.scale into account so simulator coordinates match Funcs/visuals
        double scale = (double)Scaler.scale;

        // dx is scaled spatial step (simDomainSize is in unscaled world units)
        dx = (simDomainSize * scale) / N;

        // origin is scaled world coordinate of grid[0,0]
        originX = (phys_x0 * scale) - (N / 2.0) * dx;
        originY = (phys_y0 * scale) - (N / 2.0) * dx;

        // initialize psi as Gaussian wavepacket (centered at phys_x0/phys_y0, momentum phys_kx0/phys_ky0)
        for (int i = 0; i < N; i++)
        {
            double x = originX + i * dx;
            for (int j = 0; j < N; j++)
            {
                double y = originY + j * dx;
                double dxr = x - (phys_x0 * scale);
                double dyr = y - (phys_y0 * scale);
                double r2 = dxr * dxr + dyr * dyr;
                Complex envelope = Complex.Exp(-r2 / (2.0 * phys_sigma0 * phys_sigma0));
                Complex phase = Complex.Exp(Complex.ImaginaryOne * (phys_kx0 * (x / scale) + phys_ky0 * (y / scale)));
                psiSim[i, j] = envelope * phase;
            }
        }

        // normalize initial packet
        NormalizePsi();

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
        double scale = (double)Scaler.scale;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                double Vval = Vsim[i, j];
                Complex phase = Complex.Exp(Complex.ImaginaryOne * (-Vval * simDt * 0.5));
                psiSim[i, j] *= phase;

                // optional quadratic chirp to increase σ_p while leaving |ψ|^2 unchanged
                if (phaseChirp)
                {
                    double x = originX + i * dx;
                    double cx = (phys_x0 * scale) + (phys_kx0 / phys_m) * time * scale;
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

    // Normalize psiSim so ∫ |ψ|^2 dA = 1 over the simulation domain.
    // Note: with absorbing boundaries this will counteract intended loss if autoNormalize is enabled.
    void NormalizePsi()
    {
        if (psiSim == null) return;
        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                double mag = Complex.Abs(psiSim[i, j]);
                sum += mag * mag;
            }
        }
        double area = dx * dx;
        double norm = sum * area;
        if (norm <= 0.0) return;
        double factor = 1.0 / Math.Sqrt(norm);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                psiSim[i, j] *= factor;
            }
        }
    }

    // Lasso input: start on left mouse down, collect points while held, on release fill polygon into Vsim
    void HandleLassoInput()
    {
        // only allow drawing when left mouse is used
        if (Input.GetMouseButtonDown(0))
        {
            // start drawing
            isDrawingLasso = true;
            lassoPoints.Clear();
            AddLassoPointFromMouse();
        }
        else if (Input.GetMouseButton(0) && isDrawingLasso)
        {
            AddLassoPointFromMouse();
        }
        else if (Input.GetMouseButtonUp(0) && isDrawingLasso)
        {
            isDrawingLasso = false;
            // finalize polygon
            if (lassoPoints.Count >= 3)
            {
                if (!simInitialized)
                {
                    InitializeSimulation();
                }
                ApplyLassoToPotential();
            }
            lassoPoints.Clear();
        }

        // simple visual feedback using Debug.DrawLine (editor/Game view with Gizmos enabled)
        if (isDrawingLasso && lassoPoints.Count >= 2)
        {
            for (int i = 0; i < lassoPoints.Count - 1; i++)
            {
                Debug.DrawLine(lassoPoints[i], lassoPoints[i + 1], Color.yellow, 0f, false);
            }
            // draw current mouse to last point
            UnityEngine.Vector3 world = GetMouseWorldPoint();
            Debug.DrawLine(lassoPoints[lassoPoints.Count - 1], world, Color.yellow, 0f, false);
        }
    }

    // add current mouse world point to lasso if it is far enough from last
    void AddLassoPointFromMouse()
    {
        UnityEngine.Vector2 world = GetMouseWorldPoint();
        if (lassoPoints.Count == 0)
        {
            lassoPoints.Add(world);
            return;
        }
        if (UnityEngine.Vector2.Distance(lassoPoints[lassoPoints.Count - 1], world) >= lassoMinSegment)
        {
            lassoPoints.Add(world);
        }
    }

    // convert screen mouse position to world (z = 0 plane)
    UnityEngine.Vector2 GetMouseWorldPoint()
    {
        Camera cam = Camera.main;
        if (cam == null) cam = Camera.current;
        UnityEngine.Vector3 mp = Input.mousePosition;
        mp.z = cam != null ? -cam.transform.position.z : 0f; // put on z=0 plane
        if (cam != null)
        {
            UnityEngine.Vector3 wp = cam.ScreenToWorldPoint(mp);
            return new UnityEngine.Vector2(wp.x, wp.y);
        }
        return new UnityEngine.Vector2(mp.x, mp.y);
    }

    // Apply polygon to Vsim: set or add potential inside polygon
    void ApplyLassoToPotential()
    {
        if (lassoPoints == null || lassoPoints.Count < 3) return;

        // scale polygon points to simulator (scaled) coordinates so they match Vsim grid
        double scale = (double)Scaler.scale;
        List<UnityEngine.Vector2> polyScaled = new List<UnityEngine.Vector2>(lassoPoints.Count);
        foreach (var p in lassoPoints)
        {
            polyScaled.Add(p * (float)scale);
        }

        // If lassoPotential is negative we treat polygon as a sink:
        // - inside polygon: set (or add) the negative potential (sink)
        // - outside polygon: "push up" potentials to at least +abs(lassoPotential)
        bool isSink = lassoPotential < 0.0;
        double absVal = Math.Abs(lassoPotential);

        // for each grid cell, test cell center against polygon (all in scaled coords)
        for (int i = 0; i < N; i++)
        {
            double worldX = originX + i * dx;
            for (int j = 0; j < N; j++)
            {
                double worldY = originY + j * dx;
                UnityEngine.Vector2 p = new UnityEngine.Vector2((float)worldX, (float)worldY);
                bool inside = PointInPolygon(p, polyScaled);

                if (!isSink)
                {
                    // normal behavior: set/add potential inside polygon
                    if (inside)
                    {
                        if (lassoReplaceExisting) Vsim[i, j] = lassoPotential;
                        else Vsim[i, j] += lassoPotential;
                    }
                }
                else
                {
                    // sink behavior
                    if (inside)
                    {
                        // set or add the negative sink potential
                        if (lassoReplaceExisting) Vsim[i, j] = lassoPotential;
                        else Vsim[i, j] += lassoPotential;
                    }
                    else
                    {
                        // push outside cells up to at least +absVal
                        if (Vsim[i, j] < absVal)
                        {
                            Vsim[i, j] = absVal;
                        }
                    }
                }
            }
        }

        Debug.Log("Lasso applied: polygon " + (isSink ? "sink" : "barrier") + " (value=" + lassoPotential + ")");
    }

    // standard ray-casting point-in-polygon
    bool PointInPolygon(UnityEngine.Vector2 point, List<UnityEngine.Vector2> polygon)
    {
        bool inside = false;
        int count = polygon.Count;
        for (int i = 0, j = count - 1; i < count; j = i++)
        {
            UnityEngine.Vector2 pi = polygon[i];
            UnityEngine.Vector2 pj = polygon[j];
            bool intersect = ((pi.y > point.y) != (pj.y > point.y)) &&
                             (point.x < (pj.x - pi.x) * (point.y - pi.y) / (pj.y - pi.y + 1e-12f) + pi.x);
            if (intersect) inside = !inside;
        }
        return inside;
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

        // take current scale into account so lookup matches where user draws
        double scale = (double)Scaler.scale;

        // map world coordinates to grid indices (centered) in scaled space
        double xScaled = z.Real * scale;
        double yScaled = z.Imaginary * scale;

        double fx = (xScaled - originX) / dx;
        double fy = (yScaled - originY) / dx;

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