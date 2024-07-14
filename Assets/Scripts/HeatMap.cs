using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Numerics;
using Vector3 = UnityEngine.Vector3;
using System;
using Unity.VisualScripting;
using UnityEngine.UIElements;
using System.Drawing;
using Color = UnityEngine.Color;

public class HeatMap : MonoBehaviour
{
    public bool doHeatMap = true;
    public bool derivativeHeatMap = false;
    public bool integralHeatMap = false;

    private bool heatMapGone = false;

    public Sprite squareSprite; 
    public float opacity = 1f; 
    public float squareSize = 1f;
    private float tempSquareSize;
    private uint screenLength = 9;
    private uint screenHeight = 5;

    public double calibrationDistance = 4;
    private double lowReal = double.MaxValue;
    private double HighReal = double.MinValue;
    private double lowImag = double.MaxValue;
    private double HighImag = double.MinValue;
    private double HighMag = double.MinValue;

    private double previousScale;

    public int zeroColorPower = 25;

    public bool HSVZeroWhite = false;
    public double offsetDegrees;

    public bool DoRGB = false;
    public Color realColor = Color.red;
    public Color imagColor = Color.green;
    public Color zeroColor = Color.blue;

    private GameObject[,] spriteObjects;

    public static Matrix matrix;

    void Start()
    {
        matrix = GetComponent<Matrix>();

        spriteObjects = new GameObject[(int)Math.Round(screenLength / squareSize) * 2 + 1, (int)Math.Round(screenHeight / squareSize) * 2 + 1];
        CreateHeatMap();
        UpdateHeatMap();

        previousScale = Scaler.scale;
    }

    private void Update()
    {
        if (squareSize <= 0) squareSize = 0.1f;

        if (doHeatMap && heatMapGone)
        {
            heatMapGone = false;
        }

        if(!doHeatMap && !heatMapGone)
        {
            EraseHeatMap();
            heatMapGone = true;
        }

        if (doHeatMap && previousScale != Scaler.scale)
        {
            UpdateHeatMap(); 
        }

        previousScale = Scaler.scale;

        if (doHeatMap && Input.GetKeyDown(KeyCode.R))
        {
            DistroyHeatMap();
            CreateHeatMap();
            UpdateHeatMap();
        }
    }

    public Complex HeatMapFunction(Complex value)
    {
        if (derivativeHeatMap) return Funcs.Derivative(value);
        else if (integralHeatMap) return Funcs.SimpsonsRule(value);
        else return Funcs.Function(value);    
    }
    public void CreateHeatMap()
    {
        tempSquareSize = squareSize;

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

                Complex functionOutput = HeatMapFunction(new Complex(position.x, position.y));

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
            for (int x = 0; x < (int)Math.Round(screenLength / tempSquareSize) * 2 + 1; x++)
            {
                for (int y = 0; y < (int)Math.Round(screenHeight / tempSquareSize) * 2 + 1; y++)
                {
                   
                    GameObject spriteObject = spriteObjects[x, y];

                    SpriteRenderer spriteRenderer = spriteObject.GetComponent<SpriteRenderer>();

                    Complex functionOutput = HeatMapFunction(new Complex(spriteObject.transform.position.x, spriteObject.transform.position.y));

                    Color color = ComplexToColor(functionOutput);

                    /*if (matrix.matrixMode && matrix.eigenValueV)
                    {
                        if (Complex.Abs(new Complex(spriteObject.transform.position.x, spriteObject.transform.position.y).Phase - functionOutput.Phase) < 0.1)
                        {
                            color = Color.gray;
                        }
                    }*/

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
        for (int x = 0; x < (int)Math.Round(screenLength / tempSquareSize) * 2 + 1; x++)
        {
            for (int y = 0; y < (int)Math.Round(screenHeight / tempSquareSize) * 2 + 1; y++)
            {
                GameObject spriteObject = spriteObjects[x, y];

                if (spriteObject != null)
                {
                    Destroy(spriteObject);
                }
            }
        }

        tempSquareSize = squareSize;

        spriteObjects = new GameObject[(int)Math.Round(screenLength / squareSize) * 2 + 1, (int)Math.Round(screenHeight / squareSize) * 2 + 1];

        lowReal = double.MaxValue;
        HighReal = double.MinValue;
        lowImag = double.MaxValue;
        HighImag = double.MinValue;
        HighMag = double.MinValue;
    }

    Color ComplexToColor(Complex value)
    { 
        return DoRGB ? ComplexToRGB(value) : ComplexToHSV(value);
    }

    Color ComplexToRGB(Complex value)
    {
        float magnitude = (float)Complex.Abs(value);
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
            if (currentOpacity > 0.9) currentOpacity = 0.9f;
        }
        realScale -= actualCloseZero / 2;
        imagScale -= actualCloseZero / 2;

        Color color1 = ScaleColor(realColor, realScale);
        Color color2 = ScaleColor(imagColor, imagScale);
        Color color3 = ScaleColor(zeroColor, actualCloseZero);

        return AverageColor(color1, color2, color3, currentOpacity);
    }

    Color ComplexToHSV(Complex value)
    {
        double angleRadians = Math.Atan2(value.Imaginary, value.Real);

        double angleDegrees = angleRadians * (180 / Math.PI);

        angleDegrees = (angleDegrees + 360 + offsetDegrees) % 360;

        float hue = (float)(angleDegrees / 360.0);

        float magnitude = (float)Complex.Abs(value);
        float s = Mathf.InverseLerp(0, (float)HighMag, magnitude);
        if (magnitude > HighMag) s = 1;
        float actualS = Mathf.Sqrt(s);
        float v;
        if (!HSVZeroWhite)
        {
            v = 1 - Mathf.Pow((1 - s), zeroColorPower);
            if (v > 0.5) v = 0.5f;
        }
        else
        {
            v = Mathf.Pow((1 - s), zeroColorPower);
            if (v < 0.5) v = 0.5f;
        }

        Color color = Color.HSVToRGB(hue, actualS, v);
        color.a = opacity;
        return color;
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
