using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using System;

public class Scaler : MonoBehaviour
{
    public static double scale = 1;
    private double minScale = 0.1;
    private double maxScale = 10;

    public RectTransform Unitcircle1;
    public RectTransform Unitcircle2;

    public GameObject input;

    private float orginalscale1;
    private float orginalscale2;

    private HeatMap heatmap;
    private Matrix matrix;

    public TMP_Text posr;
    public TMP_Text negr;
    public TMP_Text posi;
    public TMP_Text negi;

    void Start()
    {
        orginalscale1 = Unitcircle1.localScale.x;
        orginalscale2 = Unitcircle2.localScale.x;

        heatmap = GetComponent<HeatMap>();
        matrix = GetComponent<Matrix>();
    }

    void Update()
    {
        float scrollWheelInput = Input.GetAxis("Mouse ScrollWheel");

        BoundedScale(scale * (1 + scrollWheelInput * 0.3));

        if ((Input.GetKeyDown(KeyCode.Equals) && (Input.GetKey(KeyCode.LeftControl) || Input.GetKey(KeyCode.RightControl))))
        {
            BoundedScale(scale * 0.5);
        }
        else if ((Input.GetKeyDown(KeyCode.Minus) && (Input.GetKey(KeyCode.LeftControl) || Input.GetKey(KeyCode.RightControl))))
        {
            BoundedScale(scale * 2);
        }

        Unitcircle2.localScale = new Vector3(orginalscale2 / (float)scale + (float)(orginalscale2 - orginalscale1 + 11), orginalscale2 / (float)scale + (float)(orginalscale2 - orginalscale1 + 11), orginalscale2 / (float)scale + (float)(orginalscale1 - orginalscale2 + 11));
        Unitcircle1.localScale = new Vector3(orginalscale2 / (float)scale + 11, orginalscale2 / (float)scale + 11, orginalscale2 / (float)scale + 11);

        if (!matrix.matrixMode)
        {
            posr.text = Math.Round(3 * scale, 1).ToString("0.0");
            negr.text = (-Math.Round(3 * scale, 1)).ToString("0.0");
            posi.text = Math.Round(3 * scale, 1).ToString("0.0") + "i";
            negi.text = (-Math.Round(3 * scale, 1)).ToString("0.0") + "i";
        }
        else
        {
            posr.text = Math.Round(3 * scale, 1).ToString("0.0") + " î";
            negr.text = (-Math.Round(3 * scale, 1)).ToString("0.0") + " î";
            posi.text = Math.Round(3 * scale, 1).ToString("0.0") + " ĵ";
            negi.text = (-Math.Round(3 * scale, 1)).ToString("0.0") + " ĵ";
        }
    }

    public void BoundedScale(double scaleinput)
    {
        if (scaleinput >= minScale && scaleinput <= maxScale)
        {
            scale = scaleinput;
        }
        else if (scaleinput < minScale)
        {
            scale = minScale;
        }
        else if (scaleinput > maxScale)
        {
            scale = maxScale;
        }
    }
}
