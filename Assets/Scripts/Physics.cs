using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Physics : MonoBehaviour
{
    public bool physicMode = false;
    public double time = 0;
    public bool modulusSquared = false;

    // Phase-chirp controls (for modifying momentum distribution without changing |ψ|^2)
    public bool phaseChirp = false;
    public double chirpAlpha = 0.0; // base chirp amplitude
    public double chirpRate = 0.0;  // chirp growth per second (alpha(t) = chirpAlpha + chirpRate * time)

    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
