# ReEntryAI  
**Physics-Grounded Re-entry Simulation with Uncertainty Propagation and AI-Ready Design**

ReEntryAI is a research-grade re-entry trajectory simulation framework that I built to explore how **physics-based modeling, uncertainty propagation, and learning-based guidance** can coexist in a clean, extensible system.

The project starts from first principles (orbital mechanics, atmospheric models, aerodynamics) and is intentionally structured to support **future AI models** for guidance, robustness, and fast surrogate prediction.

This is not a black-box simulator — every modeling decision is explicit, inspectable, and designed with downstream learning in mind.

---

## Why I Built This

Most re-entry simulations fall into one of two extremes:

- **High-fidelity CFD / flight code** that is too slow or opaque for learning
- **Oversimplified toy models** that break down under uncertainty

My goal with ReEntryAI was to build a **middle-ground system**:
- Physically meaningful
- Computationally lightweight
- Explicit about uncertainty
- Structured for ML and control

This mirrors how real aerospace systems are prototyped before high-cost fidelity is introduced.

---

## Core Capabilities

### 1. Physics-Based Re-entry Dynamics
- **2D planar model** for rapid experimentation and intuition building
- **Full 3D spherical Earth model** with:
  - Radial position
  - Latitude / longitude
  - Flight-path angle
  - Heading
- Bank-angle controlled lift for trajectory shaping

The equations of motion are implemented directly, not wrapped from a library, to preserve transparency and control.

---

### 2. US Standard Atmosphere 1976 (USSA-76)
- Accurate up to **86 km altitude**
- Implements:
  - Geometric → geopotential altitude conversion
  - Layer-wise lapse rates
  - Temperature, pressure, and density computation
- Clean separation between **environment modeling** and **vehicle dynamics**

This ensures aerodynamic forces and heating are grounded in accepted aerospace standards rather than ad-hoc approximations.

---

### 3. Interval Arithmetic for Uncertainty Propagation
One of the key design choices in ReEntryAI is **interval-based dynamics**, not just Monte Carlo noise.

I implemented a custom interval arithmetic module that supports:
- Guaranteed enclosures for nonlinear dynamics
- Interval versions of `sin`, `cos`, `exp`, `sqrt`
- Interval Euler integration
- Vectorized “box” operations for state propagation

This allows me to answer questions like:
> “Given bounded uncertainty in density, mass, or initial conditions, what trajectories are *guaranteed* possible?”

This is especially important for safety-critical systems where probabilistic averages are insufficient.

---

### 4. Heating & Control-Relevant Outputs
Rather than flooding a learning system with raw mesh-level data, the simulator is designed to output **control-relevant summaries**, such as:
- Peak heating
- Mean heating
- Fraction of surface exceeding thermal limits
- Spatial non-uniformity metrics

This reflects how real flight systems compress high-dimensional physics into actionable signals.

---

## Project Structure
├── AtmosphereModel.py # USSA-76 atmosphere (float + interval versions)
├── constants.py # Physical constants and base atmospheric layers
├── interval_math.py # Interval arithmetic + inclusion functions
├── math_2d.py # Planar re-entry dynamics
├── math_3d.py # Full 3D spherical Earth dynamics
├── heat.py # Heating and thermal metrics
├── main.py # End-to-end 3D simulation runner
├── main_2d.ipynb # 2D model exploration
└── main_3d_interval.ipynb # Interval-based trajectory analysis


The codebase is modular by intent — each file corresponds to a clear physical or mathematical responsibility.

---

## Design Choices (What I’d Defend in an Interview)

### Why Euler Integration?
I intentionally start with explicit Euler:
- Forces clarity over hidden numerical behavior
- Makes uncertainty growth visible
- Easy to swap for higher-order integrators later

This mirrors how early-stage flight software is often developed.

---

### Why Interval Arithmetic Instead of Just Monte Carlo?
Monte Carlo answers *“what usually happens.”*  
Intervals answer *“what can possibly happen.”*

For re-entry, worst-case envelopes matter. Interval propagation gives:
- Deterministic guarantees
- Early detection of unsafe regions
- A natural bridge to robust control and verification

---

### Why Not Use a Physics Engine?
Because understanding, modifying, and learning from the dynamics matters more here than raw fidelity. Every equation in ReEntryAI is readable, traceable, and explainable.

---

## Future AI Direction

ReEntryAI is explicitly designed to support learning-based methods **without replacing physics**.

Planned extensions include:

### 1. Surrogate Re-entry Models
- Train lightweight neural networks to approximate:
  - Terminal states
  - Peak heating
  - Time-to-deploy
- Millisecond-level inference for onboard or real-time planning
- Physics-generated data, not synthetic labels

---

### 2. Robust Guidance with Uncertainty Awareness
- Policies trained against **interval-bounded dynamics**
- Controllers that learn to respect heating and load constraints
- Bank-angle scheduling as a learned control problem

---

### 3. Hybrid Physics + Learning
- Physics remains the ground truth
- ML models act as:
  - Accelerators
  - Constraint predictors
  - Decision-making layers

This avoids the “black-box autopilot” trap and aligns with how aerospace AI is actually adopted.

---

## What This Project Demonstrates

- Strong fundamentals in **physics-based simulation**
- Comfort working across **math, numerical methods, and software design**
- Thoughtful handling of **uncertainty and safety**
- Ability to design systems that are **AI-ready, not AI-dependent**

If you’re a recruiter or engineer reading this:  
I care deeply about building things that are *correct first*, *clever second*, and *scalable third*.

---

## Status

Active and evolving.  
This project is both a learning platform and a foundation for future research-grade work in aerospace simulation and AI-assisted guidance.
