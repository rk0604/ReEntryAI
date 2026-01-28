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

