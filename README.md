# 2nd Semester Project ECE upatras - BJT simulation

A comprehensive Python application that implements the Newton-Raphson and the Bisection to solve the non linear system of a BJT in  *exercise on page 96-101 of Numerical Methods for non-linear Engineering Models, John R. Hauser, Springer*. The project includes a visualization of results and a terminal-based user interface (check the [assets/exercise.pdf](assets/exercise.pdf) for the original exercise description).

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Reports](#reports)
- [LTspice Simulation](#ltspice-simulation)
- [Project Structure](#project-structure)
- [Contributors](#contributors)
- [License](#license)
- [Resources](#resources)

## Features

- **Newton-Raphson Method Implementation** Python code: [src/bjt_sim.py](src/bjt_sim.py)

- **Bisection Method** Python code: [src/bjt_sim.py](src/bjt_sim_bisection.py)

- **Visualization of Results**: Graphical representation of the BJT characteristics and convergence behavior.
![plot](assets/V(R_2).png)


## Requirements

- Python 3.7+
- numpy (for numerical computations)
- scipy (for optimization and root finding) --> **Only** if you want to use the fsolve function for comparison, not required for the Newton-Raphson implementation
- sympy (for symbolic mathematics, optional but useful for deriving equations) --> **Only** if you want to derive the equations symbolically, not required for the numerical implementation
- matplotlib (for visualization)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/konst3/bjt_sim.git
cd bjt_sim
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the application from the project directory:

```bash
cd src
```

- To run the Newton-Raphson method:
```bash
python bjt_sim.py
```

- To run the Bisection method:
```bash
python bjt_sim_bisection.py
```

NOTE: We found the jacobean using the [utils/find_j.py](utils/find_j.py) script, which uses sympy to derive the jacobean matrix symbolically. You can run it to see how we derived the jacobean, but it's not required for the main application. But to run:
```bash
cd ../utils
python find_j.py
```

## LTspice Simulation
We also implemented the same circuit in LTspice to compare the results. The LTspice simulation files can be found in the `ltspice/` directory: [BJT Simulation](ltspice/bjt_sim.asc).
![ltspice](assets/ltspice.png)

## Reports
The requested reports can be found in the `reports/` directory:
- [Project Report](reports/Project_Report.pdf): Detailed documentation of the project, design decisions, and implementation details.

## Project Structure
```
bjt_sim/
├── src/
│   ├── bjt_sim.py
│   ├── bjt_sim_bisection.py
├── utils/
│   ├── find_j.py
```

## Contributors (with Github usernames)
- **Konstantinos Papalamprou** (@konst3) <-- Project Lead (python implementations, reports, and LTspice simulation)
- **Argo Papadopoulou** <-- Bisection method implementation & reports
- **Thomas Tsekas** <-- Bisection method math & reports
- **Adriana Vagena** <-- Newton-Raphson method math & reports
- **Alexandros Voukelatos** <-- LTspice simulation & reports
- **Konstantinos Papadimas** <-- LTspice simulation & reports

## License
This project is released into the public domain under The Unlicense. See the [LICENSE](LICENSE) file for more details.

## Resources
The books in the exercise description were very helpful for understanding the numerical methods and their application to engineering problems. 