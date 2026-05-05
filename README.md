# 2nd Semester Project ECE upatras - BJT simulation

A comprehensive Python application that implements the Newton-Raphson method to simulate the behavior of a Bipolar Junction Transistor (BJT) under various conditions. The project includes a visualization of results and a terminal-based user interface (check the [assets/exercise.pdf](assets/exercise.pdf) for the original exercise description).

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Reports](#reports)
- [Project Structure](#project-structure)
- [Contributors](#contributors)
- [License](#license)
- [Resources](#resources)

## Features

- **Newton-Raphson Method Implementation**

- **Bisection Method**

- **Visualization of Results**: Graphical representation of the BJT characteristics and convergence behavior.


## Requirements

- Python 3.7+
- numpy (for numerical computations)
- scipy (for optimization and root finding) --> Only if you want to use the fsolve function for comparison, not required for the Newton-Raphson implementation
- sympy (for symbolic mathematics, optional but useful for deriving equations) --> Only if you want to derive the equations symbolically, not required for the numerical implementation
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
cd src && python3 main.py
```

NOTE: We found the jacobean using the [utils/find_j.py](utils/find_j.py) script, which uses sympy to derive the jacobean matrix symbolically. You can run it to see how we derived the jacobean, but it's not required for the main application.

### Application Guide (also available in the UI)



## Reports
The requested reports can be found in the `reports/` directory:
- [Project Report](reports/Project_Report.pdf): Detailed documentation of the project, design decisions, and implementation details.

## Project Structure

```
bjt_sim/
├── src/
│   ├── bjt_sim.py
├── utils/
│   ├── find_j.py
```

## Contributors (with Github usernames)

- **Konstantinos Papalamprou** (@konst3)

## License
This project is released into the public domain under The Unlicense. See the [LICENSE](LICENSE) file for more details.

## Resources

