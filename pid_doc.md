# PID Control Algorithm Documentation

## Overview

This algorithm implements a **PID-like controller** applied to a time-series of positions (e.g., signal trajectory). It uses past values to smooth the signal and compute control outputs that respond to deviations from the local trend.

---

## Data Structure

* **Input**: `std::map<float, float> positions`
  A mapping of timestamps to position values. Assumes uniform or near-uniform time intervals.

* **Output**:
  A stream of smoothed, controlled outputs for each timestamp.

---

## Parameters and Coefficients

| Parameter         | Symbol | Description                                | Default Value |
| ----------------- | ------ | ------------------------------------------ | ------------- |
| Window Size       | `N`    | Number of past samples used for averaging  | 5             |
| Proportional Gain | `Kp`   | Controls response to current error         | 1.0           |
| Integral Gain     | `Ki`   | Controls response to accumulated error     | 0.1           |
| Smoothing Factor  | `α`    | Blends past output with new control signal | 0.5           |

---

## Control Loop Details

For each timestamped position value:

### 1. **Running Average (Prediction)**

The algorithm calculates the average of the last `N` values:

$$
\text{avg}_t = \frac{1}{N} \sum_{i=t-N+1}^{t} x_i
$$

---

### 2. **Error Calculation**

$$
\text{error}_t = \text{avg}_t - x_t
$$

Where:

* $\text{avg}_t$: running average (prediction)
* $x_t$: actual position at time `t`

---

### 3. **Integral Accumulation**

$$
I_t = I_{t-1} + \text{error}_t
$$

This captures long-term bias or drift in the signal.

---

### 4. **Proportional-Integral Output**

$$
\text{PI}_t = K_p \cdot \text{error}_t + K_i \cdot I_t
$$

---

### 5. **Slope Half-Point**

To smooth changes, the controller computes a midpoint between the current PI output and the last smoothed output:

$$
\text{mid}_t = \frac{\text{PI}_t + \text{output}_{t-1}}{2}
$$

---

### 6. **Smoothed Control Output**

The final output is calculated using exponential smoothing:

$$
\text{output}_t = \alpha \cdot \text{output}_{t-1} + (1 - \alpha) \cdot \text{mid}_t
$$

## Notes

* The derivative term (`Kd`) is **not included** but can be added to improve transient response.
* Initial values (first few samples) are assumed to be stable enough to begin control.
* Use in control systems, signal smoothing, or trajectory estimation.
