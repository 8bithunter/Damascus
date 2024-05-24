using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Setter : MonoBehaviour
{
    public Transform neightx;
    public Transform nsevenx;
    public Transform nsixx;
    public Transform nfivex;
    public Transform nfourx;
    public Transform nthreex;
    public Transform ntwox;
    public Transform nonex;
    public Transform ox;
    public Transform onex;
    public Transform twox;
    public Transform threex;
    public Transform fourx;
    public Transform fivex;
    public Transform sixx;
    public Transform sevenx;
    public Transform eightx;

    public Transform nfivey;
    public Transform nfoury;
    public Transform nthreey;
    public Transform ntwoy;
    public Transform noney;
    public Transform oy;
    public Transform oney;
    public Transform twoy;
    public Transform threey;
    public Transform foury;
    public Transform fivey;

    void Start()
    {
        nonex.position = new Vector3(-1, 0, ox.transform.position.z);
        ntwox.position = new Vector3(-2, 0, ox.transform.position.z);
        nthreex.position = new Vector3(-3, 0, ox.transform.position.z);
        nfourx.position = new Vector3(-4, 0, ox.transform.position.z);
        nfivex.position = new Vector3(-5, 0, ox.transform.position.z);
        nsixx.position = new Vector3(-6, 0, ox.transform.position.z);
        nsevenx.position = new Vector3(-7, 0, ox.transform.position.z);
        neightx.position = new Vector3(-8, 0, ox.transform.position.z);
        ox.position = new Vector3(0, 0, ox.transform.position.z);
        onex.position = new Vector3(1, 0, ox.transform.position.z);
        twox.position = new Vector3(2, 0, ox.transform.position.z);
        threex.position = new Vector3(3, 0, ox.transform.position.z);
        fourx.position = new Vector3(4, 0, ox.transform.position.z);
        fivex.position = new Vector3(5, 0, ox.transform.position.z);
        sixx.position = new Vector3(6, 0, ox.transform.position.z);
        sevenx.position = new Vector3(7, 0, ox.transform.position.z);
        eightx.position = new Vector3(8, 0, ox.transform.position.z);

        noney.position = new Vector3(0, -1, ox.transform.position.z);
        ntwoy.position = new Vector3(0, -2, ox.transform.position.z);
        nthreey.position = new Vector3(0, -3, ox.transform.position.z);
        nfoury.position = new Vector3(0, -4, ox.transform.position.z);
        nfivey.position = new Vector3(0, -5, ox.transform.position.z);
        oy.position = new Vector3(0, 0, ox.transform.position.z);
        oney.position = new Vector3(0, 1, ox.transform.position.z);
        twoy.position = new Vector3(0, 2, ox.transform.position.z);
        threey.position = new Vector3(0, 3, ox.transform.position.z);
        foury.position = new Vector3(0, 4, ox.transform.position.z);
        fivey.position = new Vector3(0, 5, ox.transform.position.z);
    }
}
