<div align="center">
  <h1>Research on Nonlinear Reconstruction Method of Particle Velocity in AMS Time-of-Flight Detector</h1>
  
  [![English](https://badgen.net/badge/Language/English/blue?icon=github)](README_EN.md) [![简体中文](https://badgen.net/badge/语言/简体中文/red?icon=github)](README.md)
</div>

---
**Stage Two**
> Expected work content of this stage:
> 1. Implement nonlinear reconstruction
> 2. Add energy fluctuations in Monte Carlo simulation

## Project Introduction

This project aims to study the nonlinear reconstruction method of particle velocity in the Alpha Magnetic Spectrometer (AMS) Time-of-Flight (TOF) detector. AMS is a particle physics experiment operating on the International Space Station, and the TOF detector is used to measure the flight direction and velocity of charged particles. However, due to the ionization energy loss of charged particles in the TOF material, the particles decelerate. If a linear function is used to fit the time-space relationship of the particles, it will lead to a certain bias in the particle velocity reconstruction, which is more pronounced at low energies.

This study introduces the particle energy loss term in the particle velocity fitting process for nonlinear fitting to reduce the bias in velocity reconstruction and improve accuracy. This will help AMS to more accurately identify nuclear isotopes.

## Abstract

The Alpha Magnetic Spectrometer (AMS) is a particle physics experiment operating on the International Space Station. The Time-of-Flight (TOF) detector of AMS is used to measure the flight direction and velocity of charged particles. Due to the ionization energy loss of charged particles in the TOF material, the particles decelerate. If a linear function is used to fit the time-space relationship of the particles, it will lead to a certain bias in the particle velocity reconstruction, which is more pronounced at low energies. This study intends to introduce the particle energy loss term in the particle velocity fitting process for nonlinear fitting to reduce the bias in velocity reconstruction and improve accuracy, contributing to the identification of nuclear isotopes by AMS.

## Related Information
1. Scintillator detector used in AMS-02 TOF: [EJ-200](https://eljentechnology.com/products/plastic-scintillators/ej-200-ej-204-ej-208-ej-212)
2. Properties of some materials can be found in [Atomic and Nuclear Properties of Materials](https://pdg.lbl.gov/2024/AtomicNuclearProperties)
3. [NIST ESTAR Program](https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html) can be used to calculate stopping power, density effect correction parameters, etc., for electrons.

## Usage Instructions (For Beginners)
**Note: This method does not consider robustness, it is only for replicating the author's environment!!**

> As a reference:  
> Author's system: Windows 11 Pro for Workstations (26100.2894)
> Note: After version [de779bbab71c05b85a56a5bb39faab48a7d5aaa9](commit/de779bbab71c05b85a56a5bb39faab48a7d5aaa9), the author switched from VS Code editor to Cursor. You can continue using VS Code editor, or switch to Cursor editor. Since Cursor is a fork version of VS Code, the switch is almost seamless.

1. Go to [this link](https://code.visualstudio.com/Download) to download and install the appropriate version of [Visual Studio Code (VS Code)](https://code.visualstudio.com/).
2. Open VS Code and download the following extensions: [Remote - SSH](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh) and [WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl).
3. Right-click the Windows icon, open Terminal as Administrator, and execute the command `wsl --install` to automatically download [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/). Restart after successful installation. Alternatively, refer to [this link](https://learn.microsoft.com/en-us/windows/wsl/install-manual) for manual download. If any issues occur during automatic download, refer to the [manual download](https://learn.microsoft.com/en-us/windows/wsl/install-manual) steps.
4. (Optional) In `C:\Users\<your username>` (you can directly enter `%USERPROFILE%` in File Explorer and press `Enter`), create a new text document, open it, and type the following content:
   ```
   [experimental]
   autoMemoryReclaim=gradual
   networkingMode=mirrored
   dnsTunneling=true
   firewall=true
   autoProxy=true
   sparseVhd=true
   ```
   Then rename the document to `.wslconfig` (ensure the file extension is `wslconfig` and not `txt`).
5. Right-click the Windows icon, open Terminal, and enter `ubuntu` to access the WSL system.
6. (Optional) [Highly recommended] If you are in China, refer to [this link](https://mirrors.tuna.tsinghua.edu.cn/help/ubuntu/) to change the source.
7. Update the package list: `sudo apt update`
8. (Optional) Upgrade installed packages: `sudo apt upgrade -y`
9. Install necessary packages: `sudo apt install g++ gdb make snap snapd libgif-dev -y`
10. Download [CERN ROOT](https://root.cern.ch/): `sudo snap install root-framework`
11. (Temporarily skip this step) Download [Geant 4](https://geant4.web.cern.ch/) (although the author downloaded Geant 4 and plans to use it, it is not currently used in this project)
12. Execute the command `git clone https://github.com/abwuge/capstoneProject.git` to clone this repository.
13. Enter the capstoneProject folder: `cd capstoneProject`
14. Open the project with VS Code: `code .`

After completing the above steps, you will be able to run or debug the project normally!

If there are any missing steps, feel free to submit Issues or Pull requests!