using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using System.Numerics;
using UnityEngine;
using System.Net;
using TMPro;
using UnityEditor.UI;

public class Outputter : MonoBehaviour
{
    public int threshold = 1000;
    public bool doderivative = true;
    public bool dointegral = true;

    public Transform input;
    public Transform output;
    public Transform outputd;
    public Transform outputad;
    public Transform starti;

    public TMP_Text inputvalue;
    public TMP_Text outputvalue;
    public TMP_Text dvalue;
    public TMP_Text ivalue;
    public TMP_Text startivalue;

    void Update()
    {

        Complex complexNumber = new Complex(input.position.x, input.position.y);
        Complex result = Funcs.function(complexNumber);
        Complex derivative = new Complex(threshold - 1, threshold - 1);
        Complex antiderivative = new Complex(threshold - 1, threshold - 1);
        Complex startinumber = new Complex(starti.position.x, starti.position.y);

        if (doderivative)
        {
            derivative = Funcs.Derivative(complexNumber);
        }

        if (dointegral)
        {
            antiderivative = Funcs.SimpsonsRule(complexNumber);
            starti.gameObject.SetActive(true);
        }
        else
        {
            starti.gameObject.SetActive(false);
        }

        if ((result.Real < threshold && result.Imaginary < threshold) && (result.Real > -threshold && result.Imaginary > -threshold))
        {
            output.position = new UnityEngine.Vector3((float)result.Real, (float)result.Imaginary, output.position.z);
        }
        else
        {
            output.position = new UnityEngine.Vector3((float)threshold - 1, (float)threshold - 1, output.position.z);
        } 

        if ((derivative.Real < threshold && derivative.Imaginary < threshold) && (derivative.Real > -threshold && derivative.Imaginary > -threshold))
        {
            outputd.position = new UnityEngine.Vector3((float)derivative.Real, (float)derivative.Imaginary, output.position.z);
        }
        else
        {
            outputd.position = new UnityEngine.Vector3((float)threshold - 1, (float)threshold - 1, output.position.z);
        }

        if ((antiderivative.Real < threshold && antiderivative.Imaginary < threshold) && (antiderivative.Real > -threshold && antiderivative.Imaginary > -threshold))
        {
            outputad.position = new UnityEngine.Vector3((float)antiderivative.Real, (float)antiderivative.Imaginary, output.position.z);
        }
        else
        {
            outputad.position = new UnityEngine.Vector3((float)threshold - 1, (float)threshold - 1, output.position.z);
        }

        inputvalue.text = FormatComplexNumber(complexNumber, true);
        outputvalue.text = FormatComplexNumber(result, true);
        dvalue.text = FormatComplexNumber(derivative, doderivative);
        ivalue.text = FormatComplexNumber(antiderivative, dointegral);
        startivalue.text = FormatComplexNumber(startinumber, dointegral);
    }

    public string FormatComplexNumber(Complex complexNumber, bool dothing)
    {
        string realPart = Math.Abs(complexNumber.Real * Scaler.scale).ToString("0.000");
        string imagPart = Math.Abs(complexNumber.Imaginary * Scaler.scale).ToString("0.000");

        string formattedNumber;

        if (!dothing)
        {
            return "N/A";
        }
        else if ((Math.Abs(complexNumber.Real) > threshold || Math.Abs(complexNumber.Imaginary) > threshold) || (double.IsNaN(complexNumber.Real) || double.IsNaN(complexNumber.Imaginary)))
        {
            return "Too Big!";
        }
        else
        { 
            if (complexNumber.Real >= 0 && complexNumber.Imaginary >= 0)
            {
                formattedNumber = $"{realPart} + {imagPart}i";
            }
            else if (complexNumber.Real >= 0 && complexNumber.Imaginary < 0)
            {
                formattedNumber = $"{realPart} - {imagPart}i";
            }
            else if (complexNumber.Real < 0 && complexNumber.Imaginary >= 0)
            {
                formattedNumber = $"-{realPart} + {imagPart}i";
            }
            else
            {
                formattedNumber = $"-{realPart} - {imagPart}i";
            }
            return formattedNumber;
        }
    }
}
