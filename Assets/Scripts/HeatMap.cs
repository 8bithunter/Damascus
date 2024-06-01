using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Numerics;
using Vector3 = UnityEngine.Vector3;
using System;
using Unity.VisualScripting;

public class HeatMap : MonoBehaviour
{
    public bool doHeatMap = true;

    public Sprite squareSprite; 
    public float opacity = 0.3f; 
    public float squareSize = 1f; 
    private int screenLength = 9;
    private int screenHeight = 5;

    double lowReal = double.MaxValue;
    double HighReal = double.MinValue;
    double lowImag = double.MaxValue;
    double HighImag = double.MinValue;
    double HighMag = double.MinValue;

    void Start()
    {
        CreateHeatMap();
    }
    public void CreateHeatMap()
    {
        if (doHeatMap)
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

                    Complex functionOutput = Funcs.function(new Complex(x, y));

                    //if (x > -4 / squareSize && x < 4 / squareSize && y > -4 / squareSize && y < 4 / squareSize)

                    if (functionOutput.Real < lowReal) lowReal = functionOutput.Real;
                    if (functionOutput.Real > HighReal) HighReal = functionOutput.Real;
                    if (functionOutput.Imaginary < lowImag) lowImag = functionOutput.Imaginary;
                    if (functionOutput.Imaginary > HighImag) HighImag = functionOutput.Imaginary;
                    if (Complex.Abs(functionOutput) > HighMag) HighMag = Complex.Abs(functionOutput);

                    Complex value = functionOutput;

                    Color color = ComplexToColor(value);

                    spriteRenderer.color = color;
                }
            }
        }
    }
    Color ComplexToColor(Complex value)
    { 
        float magnitude = (float) Complex.Abs(value);

        float preLerpOpacity = Mathf.InverseLerp(0, (float)HighMag, magnitude);
        //opacity = Mathf.Lerp(0f, 0.2f, preLerpOpacity);

        float normalizedReal = Mathf.InverseLerp((float)lowReal, (float)HighReal, (float)value.Real);
        if (value.Real > HighReal) normalizedReal = 1;
        if (value.Real < lowReal) normalizedReal = 0;
        float normalizedImaginary = Mathf.InverseLerp((float)lowImag, (float)HighImag, (float)value.Imaginary);
        if (value.Imaginary > HighImag) normalizedImaginary = 1;
        if (value.Imaginary < lowImag) normalizedImaginary = 0;

        float red = Mathf.Lerp(0f, 1f, normalizedReal);
        float green = 0;
        float blue = Mathf.Lerp(0f, 1f, normalizedImaginary);

        return new Color(red, green, blue, opacity);
    }
}
