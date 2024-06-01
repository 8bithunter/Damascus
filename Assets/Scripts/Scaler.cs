using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Scaler : MonoBehaviour
{
    public static double scale = 1;

    public RectTransform Unitcircle1;
    public RectTransform Unitcircle2;

    public GameObject input;

    private float orginalscale1;
    private float orginalscale2;

    private HeatMap heatmap;

    void Start()
    {
        orginalscale1 = Unitcircle1.localScale.x;
        orginalscale2 = Unitcircle2.localScale.x;

        heatmap = GetComponent<HeatMap>();
    }

    void Update()
    {
        float scrollWheelInput = Input.GetAxis("Mouse ScrollWheel");

        BoundedScale(scale * (1 + scrollWheelInput * 0.3));

        if ((Input.GetKeyDown(KeyCode.Equals) && (Input.GetKey(KeyCode.LeftControl) || Input.GetKey(KeyCode.RightControl))))
        {
            BoundedScale(scale * 0.5);
            heatmap.CreateHeatMap();
        }
        else if ((Input.GetKeyDown(KeyCode.Minus) && (Input.GetKey(KeyCode.LeftControl) || Input.GetKey(KeyCode.RightControl))))
        {
            BoundedScale(scale * 2);
            heatmap.CreateHeatMap();
        }

        Unitcircle2.localScale = new Vector3(orginalscale2 / (float)scale + (float)(orginalscale2 - orginalscale1 + 11), orginalscale2 / (float)scale + (float)(orginalscale2 - orginalscale1 + 11), orginalscale2 / (float)scale + (float)(orginalscale1 - orginalscale2 + 11));
        Unitcircle1.localScale = new Vector3(orginalscale2 / (float)scale + 11, orginalscale2 / (float)scale + 11, orginalscale2 / (float)scale + 11);
    }

    public void BoundedScale(double scaleinput)
    {
        if (scaleinput >= 0.25 && scaleinput <= 4)
        {
            scale = scaleinput;
        }
    }
}
