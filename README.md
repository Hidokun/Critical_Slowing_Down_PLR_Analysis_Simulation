# Pupillary Light Reflex (PLR) Simulation

This is a computational neuroscience simulation framework that models the Pupillary Light Reflex (PLR) using a Nonlinear Delay-Differential Equation (DDE). 

The project investigates the dynamics of the pupil as the feedback loop approaches a point of instability (Hopf bifurcation), specifically exploring:
- **Critical Slowing Down:** Exponentially longer recovery times after a light pulse as loop gain increases.
- **Hippus:** Spontaneous pupillary oscillations (1-3 Hz) when the system loses stability.

## Project Structure

- `prl_model.py`: Core mathematical DDE model of the iris sphincter using a Hill function.
- `stimulus.py`: Generates the light pulse experimental protocols.
- `analysis.py`: Data processing, including exponential curve fitting for recovery time and spectral analysis for hippus detection.
- `main.py`: To run parameter sweeps and saves results.
- `figures.py`: Generates  figures from the simulated data.

## Installation

Ensure you have Python installed, then install the required dependencies:

```bash
pip install -r requirements.txt
```

## Usage

1. Run the main simulation to generate the data. This will create an `output/` folder and populate it with CSV files:
```bash
python main.py
```

2. Generate the figures based on the simulated data. The images will be saved in the `output/` folder:
```bash
python figures.py
```
