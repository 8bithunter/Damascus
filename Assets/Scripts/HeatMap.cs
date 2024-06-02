using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Numerics;
using Vector3 = UnityEngine.Vector3;
using System;
using Unity.VisualScripting;
using UnityEngine.UIElements;

public class HeatMap : MonoBehaviour
{
    public bool doHeatMap = true;
    private bool heatMapGone = false;

    public Sprite squareSprite; 
    public float opacity = 1f; 
    public float squareSize = 1f; 
    private int screenLength = 9;
    private int screenHeight = 5;

    public double calibrationDistance = 4;
    double lowReal = double.MaxValue;
    double HighReal = double.MinValue;
    double lowImag = double.MaxValue;
    double HighImag = double.MinValue;
    double HighMag = double.MinValue;

    GameObject[,] spriteObjects;

    public int zeroColorPower = 25;

    public Color realColor = Color.red;
    public Color imagColor = Color.green;
    public Color zeroColor = Color.blue;

    void Start()
    {
        spriteObjects = new GameObject[(int)Math.Round(screenLength / squareSize) * 2 + 1, (int)Math.Round(screenHeight / squareSize) * 2 + 1];
        CreateHeatMap();
        UpdateHeatMap();
    }
    private long frameCount = 0;
    private void Update()
    {

        if (doHeatMap && heatMapGone)
        {
            heatMapGone = false;
        }

        if(!doHeatMap && !heatMapGone)
        {
            EraseHeatMap();
            heatMapGone = true;
        }

        if (doHeatMap)
        {
            frameCount++;

            if (frameCount % 5 == 0)
                UpdateHeatMap();
        }

        if (doHeatMap && Input.GetKeyDown(KeyCode.R))
        {
            DistroyHeatMap();
            CreateHeatMap();
            UpdateHeatMap();
        }
    }
    public void CreateHeatMap()
    {
        for (int x = -(int)Math.Round(screenLength / squareSize); x <= (int)Math.Round(screenLength / squareSize); x++)
        {
            for (int y = -(int)Math.Round(screenHeight / squareSize); y <= (int)Math.Round(screenHeight / squareSize); y++)
            {
                GameObject spriteObject = new GameObject("Tile_" + x + "_" + y);

                SpriteRenderer spriteRenderer = spriteObject.AddComponent<SpriteRenderer>();

                spriteRenderer.sprite = squareSprite;

                Vector3 position = new Vector3(x * squareSize, y * squareSize, 0);

                spriteObject.transform.position = position;

                spriteObject.transform.localScale = new Vector3(squareSize, squareSize, 1f);

                Complex functionOutput = Funcs.function(new Complex(position.x, position.y));

                if (position.x > -calibrationDistance && position.x < calibrationDistance && position.y > -calibrationDistance && position.y < calibrationDistance)
                {
                    if (functionOutput.Real < lowReal) lowReal = functionOutput.Real;
                    if (functionOutput.Real > HighReal) HighReal = functionOutput.Real;
                    if (functionOutput.Imaginary < lowImag) lowImag = functionOutput.Imaginary;
                    if (functionOutput.Imaginary > HighImag) HighImag = functionOutput.Imaginary;
                    if (Complex.Abs(functionOutput) > HighMag) HighMag = Complex.Abs(functionOutput);
                }

                spriteObjects[x + (int)Math.Round(screenLength / squareSize), y + (int)Math.Round(screenHeight / squareSize)] = spriteObject;
            }
        }
    }

    public void UpdateHeatMap()
    {
        if (doHeatMap)
        {
            for (int x = 0; x < (int)Math.Round(screenLength / squareSize) * 2 + 1; x++)
            {
                for (int y = 0; y < (int)Math.Round(screenHeight / squareSize) * 2 + 1; y++)
                {
                    GameObject spriteObject = spriteObjects[x, y];

                    SpriteRenderer spriteRenderer = spriteObject.GetComponent<SpriteRenderer>();

                    Complex functionOutput = Funcs.function(new Complex(spriteObject.transform.position.x, spriteObject.transform.position.y));

                    Color color = ComplexToColor(functionOutput);

                    spriteRenderer.color = color;
                }
            }
        }
    }

    public void EraseHeatMap()
    {
        for (int x = 0; x < (int)Math.Round(screenLength / squareSize) * 2 + 1; x++)
        {
            for (int y = 0; y < (int)Math.Round(screenHeight / squareSize) * 2 + 1; y++)
            {
                 GameObject spriteObject = spriteObjects[x, y];

                 SpriteRenderer spriteRenderer = spriteObject.GetComponent<SpriteRenderer>();

                 spriteRenderer.color = new Color (0, 0, 0, 0);
            }
        }
    }

    public void DistroyHeatMap()
    {
        for (int x = 0; x < (int)Math.Round(screenLength / squareSize) * 2 + 1; x++)
        {
            for (int y = 0; y < (int)Math.Round(screenHeight / squareSize) * 2 + 1; y++)
            {
                GameObject spriteObject = spriteObjects[x, y];

                if (spriteObject != null)
                {
                    Destroy(spriteObject);
                }
            }
        }

        lowReal = double.MaxValue;
        HighReal = double.MinValue;
        lowImag = double.MaxValue;
        HighImag = double.MinValue;
        HighMag = double.MinValue;
    }

    Color ComplexToColor(Complex value)
    { 
        float magnitude = (float) Complex.Abs(value);
        float closeZero = Mathf.InverseLerp(0, (float)HighMag, magnitude);
        if (magnitude > HighMag) closeZero = 1;

        float realScale = Mathf.InverseLerp((float)lowReal, (float)HighReal, (float)value.Real);
        if (value.Real > HighReal) realScale = 1;
        if (value.Real < lowReal) realScale = 0;

        float imagScale = Mathf.InverseLerp((float)lowImag, (float)HighImag, (float)value.Imaginary);
        if (value.Imaginary > HighImag) imagScale = 1;
        if (value.Imaginary < lowImag) imagScale = 0;

        float actualCloseZero = Mathf.Pow((1 - closeZero), zeroColorPower);

        float currentOpacity = opacity;

        if (actualCloseZero > opacity)
        {
            currentOpacity = actualCloseZero;
        }
        realScale -= actualCloseZero / 2;
        imagScale -= actualCloseZero / 2;

        Color color1 = ScaleColor(realColor, realScale);
        Color color2 = ScaleColor(imagColor, imagScale);
        Color color3 = ScaleColor(zeroColor, actualCloseZero);

        return AverageColor(color1, color2, color3, currentOpacity);
    }

    Color ScaleColor(Color color, float scaleValue)
    {
        scaleValue = Math.Max(0, Math.Min(1, scaleValue));

        float newR = (color.r * scaleValue);
        float newG = (color.g * scaleValue);
        float newB = (color.b * scaleValue);

        return new Color(newR, newG, newB);
    }

    public static Color AverageColor(Color color1, Color color2, Color color3, float opacity)
    {
        float averageR = (color1.r + color2.r + color3.r) / 3;
        float averageG = (color1.g + color2.g + color3.g) / 3;
        float averageB = (color1.b + color2.b + color3.b) / 3;

        return new Color(averageR, averageG, averageB, opacity);
    }
}
