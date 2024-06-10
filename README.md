![DamascusLogoRound](https://github.com/8bithunter/Damascus/assets/148516367/8bd5e2e5-c19f-4519-9ccd-7e57b1778626)

<div align="center">
  <p><strong>Damascus is a Unity-based software for graphing and visualizing complex functions, derivatives, and integrals.</strong></p>
</div>

# Installation
Step 1: Download and install the latest version of Unity Hub and the Unity editor.

Step 2: Download the latest version of Damascus, and unzip the file.

Step 3: At the top of Unity Hub, click "Add" and select the unzipped file.

Step 4: That's it! Just run the program, and you are good to go!

# Features
- All functions are supported! Use premade functions in "System.Numerics.Complex", or write your own!
- Easily input any desired value, and receive the output, derivative, and integral (from the starting value) at that point.
- Along with direct inputs, Damascus allows you to generate maps of your chosen function, enabling you to visualize the shape of functions and identify zeroes.
- Two by two matrices are now supporting, letting you see the effect of transformations on any vector.

# Usage
### Running the Program
Once the program is open in the Unity editor, select the "Damascus" scene in the "Scenes" folder under "Assets". Then click the play button at the top of the editor, and the program will run.

### Choosing a Function
Open the "Funcs" file in the "Scripts" folder under "Assets" with a code editor. Inside "Funcs", there is a method named "Function", which contains the assignment you will change. Below the commented lines, there is an assignment that looks like "Complex output = (function)". Simply change "(function)" to whatever function you want to graph. Most functions can be found in the "Complex" class, but any arithmetic operator will work too. If the "Complex" class does not contain your desired function, feel free to program it yourself!

### Controls
<strong>Left Click:</strong> Left clicking smoothly brings the input selector to the position of your mouse.

<strong>Arrow Keys:</strong> Clicking any arrow key will move the input selector 0.001 units in its corresponding direction.

<strong>Ctrl:</strong> Holding ctrl while pressing the arrow keys will cause the input selector to smoothly move in the arrow key's corresponding direction.

<strong>Shift:</strong> Holding shift will multiply any movement triggered by the arrow keys by a factor of 5.

<strong>Alt:</strong> Holding alt will lock the input selector to the unit circle.

<strong>X:</strong> While holding x, type in any real co-ordinate you want your input selector to land on, once you are done typing, release x and in input selector will move to the desired real co-ordinate.

<strong>Y:</strong> While holding y, type in any imaginary coordinate (don't add "i" at any point) you want your input selector to land on, once you are done typing, release x and in input, the selector will move to the desired imaginary co-ordinate.

<strong>I:</strong> Clicling i will toggle between controlling the input selector, or the integral start point selector.

<strong>R:</strong> Clicling r will recalculate the heat map, applying any changes made under the "Code" object of the project. Only click r if you make changes while the program is running, any changes made before running the program are automatically applied.
