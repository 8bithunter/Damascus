using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Setter : MonoBehaviour
{
    private int screenLength = 9;
    private int screenHeight = 5;

    private float lineThickness = 0.02f;
    public GameObject squarePrefab;
   
    private void Start()
    {
        for (int i = -screenLength; i <= screenLength; i++)
        {
            GameObject square = Instantiate(squarePrefab, new Vector3(i, 0, 0), Quaternion.identity);
            square.transform.localScale = new Vector3(lineThickness, 20, 1);
            if (i == 0) square.GetComponent<SpriteRenderer>().color = Color.white;

            square.name = "x_" + i.ToString();
        }

        for (int i = -screenHeight; i <= screenHeight; i++)
        {
            GameObject square = Instantiate(squarePrefab, new Vector3(0, i, 0), Quaternion.identity);
            square.transform.localScale = new Vector3(20, lineThickness, 1);
            if (i == 0) square.GetComponent<SpriteRenderer>().color = Color.white;

            square.name = "y_" + i.ToString();
        }
    }
}
