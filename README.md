# Heat-Transfer
Code for MAE 423 Heat Transfer (Fall 2019-20). 
***

## Code Overview

### P-Set 3
Use Gauss-Seidel Iteration to relax the temperature to reach equilibrium for data points. "PSet-3.p<i></i>y" is the project code. "p-set-3-code.tex" compiles images and the code into a PDF for printing.

<p align="center">
<img alt="PSet-3 Final Image" src="https://raw.githubusercontent.com/sdale28/Heat-Transfer/master/Problem-Set-3/images/pset-3-figure-3-final.png">
</p>

### P-Set 5
Numerically solve a specified grid for temperature across time-steps until reaching steady-state using finite difference approximations. This problem has unsteady heat conduction with an insulated boundary, convective boundary, and two boundaries maintained at constant temperature. 
"PSet-5.p<i></i>y" is the project code. "PSet-5-Print.tex" compiles images and the code into a PDF for printing.

<p align="center">
<img alt="PSet-5 Final Image" src="https://raw.githubusercontent.com/sdale28/Heat-Transfer/master/Problem-Set-5/images/Problem-Set-5-Figure-4.png">
</p>

---
## Software Setup
Setup Python Code: 
- [Virtual Studio Code][3]
- [Python][4] (For Windows: "Windows x86-64 executable installer")
- [Python Extension for VS Code][5]

Setup Github:
- Windows: [Install Git][1]
- Mac: Go to Terminal, enter the command "git --version". If you don't have git, Apple's XCode Extension will automatically prompt to be installed. 
- [Github Desktop][2]
- In VS Code Settings, search for "Post Git Commit", change to "Push".

To Run Python Code:
- Have the desired Python file open in VS Code (e.g. "PSet-3.py").
- Hit the green play button to the right of the tabs.

Recommended VS Code Setting:
- File > Preferences > Keyboard Shortcuts
- Search for "Python: Run Python File in Terminal".
- Set this to "Shift + Enter".
- Now you can just hit "Shift + Enter" to run the Python code.

To Get the Code:
- Open Github Desktop
- File > "Clone Repository..."
- Select URL Tab
- URL of this repo: https://github.com/sdale28/Heat-Transfer.git
- Make sure the local path is where you want it, as you cannot move this folder in the future (if you do, it will no longer sync with Github). 

To Make Sure Code is Up to Date:
- On left menu, click on the three connected points that branch ("Source Control"). 
- Click three dot icon to the right of the checkmark and refresh arrow ("More Actions").
- Select "Pull". 
- (If you want to force reset back to the version that is stored on Github: "Pull (rebase)". Not that this overwrites all changes you've made locally). 

To Post Code to Github: 
- On left menu, click on the three connected points that branch ("Source Control"). 
- Put note in textbox (usually required in Git, example "Added color bar to plot"). 
- Click checkmark (or "Ctrl + Enter" in message text area). 

[1]:https://git-scm.com/download/win
[2]:https://desktop.github.com/
[3]:https://code.visualstudio.com/
[4]:https://www.python.org/downloads/
[5]:https://marketplace.visualstudio.com/items?itemName=ms-python.python
