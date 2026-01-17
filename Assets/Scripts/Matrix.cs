using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;

public class Matrix : MonoBehaviour
{
    public bool matrixMode = false;
    //public bool eigenValueV = false;

    public double[] î = new double[2] {1,
                                       0};

    public double[] ĵ = new double[2] {0,
                                       1};

    public Complex Transform(Complex value)
    {
        return new Complex(value.Real * î[0] + value.Imaginary * ĵ[0], value.Real * î[1] + value.Imaginary * ĵ[1]);
    }
}
