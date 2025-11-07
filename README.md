# 🔗 Network Synchronization of Lorenz Oscillators

This repository contains **MATLAB** code for simulating and visualizing the synchronization of a network of coupled Lorenz oscillators.
The project demonstrates the application of **graph theory**—specifically weighted Laplacian matrices and a novel weight assignment method—to achieve and analyze network synchronization across **switching topologies**.

-----

## 📖 Project Overview

The core of this project is a set of MATLAB scripts that model and simulate a network of 10 chaotic Lorenz oscillators across **three consecutive time intervals** (using graphs **$G1, G2, G3$**).

  * The **`SyncCouplingAssign`** function implements a method based on **Negative Imbalance Vectors** and **Cycle Basis Vectors** to assign coupling weights that guarantee synchronization.
  * Synchronization error is calculated as the **L2 norm** distance between an oscillator's state and all other oscillators.
  * Simulation data (time and synchronization error) is exported to an Excel file, **`synchronization_data.xlsx`**.

-----

## 📂 File Descriptions

### MATLAB Scripts

| File Name | Description | Key Mechanism |
| :--- | :--- | :--- |
| **`main.m`** | **Main simulation script.** Defines parameters ($\sigma, \rho, \beta$), sets up the three switching connectivity digraphs ($G1, G2, G3$), assigns weights using `SyncCouplingAssign`, runs the full time-multiplexed simulation, and plots/saves the synchronization error. | Switching topology simulation across $G1 \to G2 \to G3$ |
| **`LorenzOscillator.m`** | Defines the uncoupled dynamics ($f(X)$) of a single Lorenz system (state: $[x; y; z]$). | $\sigma=10$, $\rho=28$, $\beta=8/3$ |
| **`SimulateCoupledSystems.m`** | Integrates the network dynamics using `ode45`. It first builds the **weighted Laplacian matrix ($L$)** from the input graph $G$. | Numerical integration via `ode45` |
| **`CoupledDynamics.m`** | Defines the full coupled network dynamics ($\dot{X}$) required by `ode45`. Implements the diffusive coupling $\left( -P L X^T - P X \right)$. | Applies coupling through Projection Matrix **$P$** and Laplacian **$L$** |
| **`SyncCouplingAssign.m`** | **Core weight assignment function.** Calls both `NegativeImbalanceVector` and `CycleBasisVector` sequentially to determine the synchronizing edge weights based on parameter $a$. | Calls `NegativeImbalanceVector` and `CycleBasisVector` |
| **`NegativeImbalanceVector.m`** | Assigns an initial set of positive edge weights to ensure all but one vertex in each Strongly Connected Component (SCC) has a **negative imbalance** with respect to parameter $a$. | Uses **Shortest Path** algorithm to assign weights |
| **`CycleBasisVector.m`** | Computes a final set of positive edge weights based on the **cycle basis** of the SCCs to ensure all edges are active and weighted appropriately for synchronization. | Uses cycles to ensure **non-zero weights** for edges within SCCs |

-----

## 🚀 Getting Started

### Prerequisites

  * **MATLAB** (preferably R2017b or newer, which fully supports the `digraph` object).
  * Required toolboxes are typically included in the base installation (e.g., Graph and Network Algorithms).

### Usage

1.  Place all `.m` files in the same directory.
2.  Open **MATLAB**.
3.  Run the main script from the MATLAB command window:
    ```matlab
    main
    ```

### Outputs

1.  **Multiple Figures:** Plots showing the three directed graphs ($G1, G2, G3$) with their final assigned edge weights.
2.  **Synchronization Error Plot:** A figure showing the time evolution of the synchronization error for all 10 oscillators.
3.  **Data File:** An Excel file named **`synchronization_data.xlsx`** containing the simulation time vector and the error vector for each oscillator ($e_1$ through $e_{10}$).

-----

## 📌 Notes on Synchronization Theory

  * **Coupling Parameter ($a$):** The value `a` is analytically derived in `main.m` to be sufficient for synchronizing the Lorenz system, relating to the one-sided Lipschitz constant of the uncoupled system.
  * **Projection Matrix ($P$):** The coupling is applied only to the $x$-state, as set by `P = diag([1, 0, 0])`.
  * **Generalizability:** This weight assignment framework generalizes beyond Lorenz oscillators to any system whose local dynamics satisfy the required QUAD condition, using the corresponding constant as parameter $a$.
