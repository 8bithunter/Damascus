using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.Numerics;
using Vector3 = UnityEngine.Vector3;

public class CursorFollow : MonoBehaviour
{
    public GameObject objectToMove;
    public GameObject output;
    public GameObject starti;
    bool isoutput = true;
    public float moveSpeed = 5f;
    public float arrowSpeed = 0.001f;

    private bool isRecordingX = false;
    private bool isRecordingY = false;

    private string recordedXValue = "";
    private string recordedYValue = "";

    void Update()
    {
        if (Input.GetKeyDown(KeyCode.I))
        {
            objectToMove = isoutput ? starti : output;
            isoutput = !isoutput;
        }

        if (Input.GetKey(KeyCode.Return))
        {
            output.transform.position = new Vector3(0, 0, 0);
            starti.transform.position = new Vector3(0, 0, 0);
            Scaler.scale = 1;
        }

        if (Input.GetMouseButton(0))
        {
            Vector3 cursorScreenPosition = Input.mousePosition;
            Vector3 cursorWorldPosition = Camera.main.ScreenToWorldPoint(cursorScreenPosition);
            cursorWorldPosition.z = 0f;

            objectToMove.transform.position = Vector3.Lerp(objectToMove.transform.position, cursorWorldPosition, moveSpeed * Time.deltaTime);
        }

        if (Input.GetKey(KeyCode.LeftControl) || Input.GetKey(KeyCode.RightControl))
        {
            if (Input.GetKey(KeyCode.RightArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.right * arrowSpeed * 4f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.right * arrowSpeed);
                }
            }
            if (Input.GetKey(KeyCode.LeftArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.left * arrowSpeed * 4f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.left * arrowSpeed);
                }
            }
            if (Input.GetKey(KeyCode.UpArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.up * arrowSpeed * 4f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.up * arrowSpeed);
                }
            }
            if (Input.GetKey(KeyCode.DownArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.down * arrowSpeed * 4f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.down * arrowSpeed);
                }
            }
        }
        else
        {
            if (Input.GetKeyDown(KeyCode.RightArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.right * arrowSpeed * 10f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.right * arrowSpeed);
                }
            }
            if (Input.GetKeyDown(KeyCode.LeftArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.left * arrowSpeed * 10f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.left * arrowSpeed);
                }
            }
            if (Input.GetKeyDown(KeyCode.UpArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.up * arrowSpeed * 10f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.up * arrowSpeed);
                }
            }
            if (Input.GetKeyDown(KeyCode.DownArrow))
            {
                if (Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift))
                {
                    objectToMove.transform.Translate(Vector3.down * arrowSpeed * 10f);
                }
                else
                {
                    objectToMove.transform.Translate(Vector3.down * arrowSpeed);
                }
            }
        }

        if (Input.GetKey(KeyCode.LeftAlt) || Input.GetKey(KeyCode.RightAlt))
        {
            Complex UnitCircle = Complex.FromPolarCoordinates(1, Math.Atan2(objectToMove.transform.position.y, objectToMove.transform.position.x)) / Scaler.scale;

            Vector3 UnitCircleVector = new Vector3((float)UnitCircle.Real, (float)UnitCircle.Imaginary, objectToMove.transform.position.z);

            objectToMove.transform.position = UnitCircleVector;
        }

        if (Input.GetKeyDown(KeyCode.X))
        {
            isRecordingX = true;
            recordedXValue = "";
        }

        if (Input.GetKeyUp(KeyCode.X))
        {
            isRecordingX = false;
            float xCoordinate;
            if (float.TryParse(recordedXValue, out xCoordinate))
            {
                objectToMove.transform.position = new Vector3(xCoordinate / (float)Scaler.scale, objectToMove.transform.position.y, objectToMove.transform.position.z);
            }
        }

        if (Input.GetKeyDown(KeyCode.Y))
        {
            isRecordingY = true;
            recordedYValue = ""; 
        }

        if (Input.GetKeyUp(KeyCode.Y))
        {
            isRecordingY = false;
            float yCoordinate;
            if (float.TryParse(recordedYValue, out yCoordinate))
            {
                objectToMove.transform.position = new Vector3(objectToMove.transform.position.x, yCoordinate / (float)Scaler.scale, objectToMove.transform.position.z);
            }
        }

        if (isRecordingX)
        {
            RecordNumericInput(ref recordedXValue);
        }

        if (isRecordingY)
        {
            RecordNumericInput(ref recordedYValue);
        }
    }
    void RecordNumericInput(ref string recordedValue)
    {
        if (Input.anyKeyDown)
        {
            foreach (KeyCode keyCode in System.Enum.GetValues(typeof(KeyCode)))
            {
                if ((keyCode >= KeyCode.Alpha0 && keyCode <= KeyCode.Alpha9) ||
                    keyCode == KeyCode.Period || keyCode == KeyCode.KeypadPeriod || keyCode == KeyCode.Minus)
                {
                    if (Input.GetKeyDown(keyCode))
                    {
                        if (keyCode == KeyCode.Period || keyCode == KeyCode.KeypadPeriod)
                        {
                            if (!recordedValue.Contains("."))
                                recordedValue += ".";
                        }
                        else if (keyCode == KeyCode.Minus)
                        {
                            if (recordedValue == "")
                                recordedValue += "-";
                        }
                        else
                        {
                            recordedValue += (keyCode - KeyCode.Alpha0);
                        }
                    }
                }
            }
        }
    }
}
