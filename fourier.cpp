#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <fstream>
#include <complex>
#include <deque>
#include <numeric>

// Compact Fourier series
struct fourier_component {
    float angular_velocity;
    float amplitude;
    float phase;
    int n_degree;

    float compute_component(const float time) const {
        return amplitude * std::cos(angular_velocity * time + phase);
    }
};

std::map<float, float> compute_trajectory(const float& initial_condition,
                        fourier_component fourier_series[],
                        const int fourier_series_size,
                        const float& time_span,
                        const float& time_step) {
    int num_samples = static_cast<int>(time_span / time_step);
    if (num_samples == 0) return std::map<float, float>();
    std::map<float, float> positions;

    for (int i = 0; i < num_samples; i++) {
        float time = i * time_step;
        float position = initial_condition;
        for (int j = 0; j < fourier_series_size; j++) {
            position += fourier_series[j].compute_component(time);
        }
        positions[time] = position;
    }

    return positions;
}

void report_trajectory(fourier_component fourier_series[], const int& fourier_series_size,
                       const std::map<float, float>& positions, const std::string& filename) {
    // Initialize the output file
    std::ofstream outputFile(filename);
    outputFile << "# Trajectory Report\n\n";

    // Write the Fourier components to the output file
    outputFile << "## Fourier components\n\n";
    outputFile << "| n | Angular Velocity | Amplitude | Phase |\n";
    outputFile << "| --- | --- | --- | --- |\n";

    for (int i = 0; i < fourier_series_size; i++) {
        outputFile << "| " << fourier_series[i].n_degree << " | " << fourier_series[i].angular_velocity << " | " << fourier_series[i].amplitude << " | " << fourier_series[i].phase << " |\n";
    }

    // Write the positions to the output file
    outputFile << "\n## Positions\n\n";
    outputFile << "| Time | Position |\n";
    outputFile << "| --- | --- |\n";

    for (const auto& position : positions) {
        outputFile << "| " << position.first << " | " << position.second << " |\n";
    }

    outputFile.close();
}

void generate_components(fourier_component fourier_series[], const int& fourier_series_size,
                         const float& base_angular_velocity, const float& amplitude_range, const float& phase_range) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<float> amplitude_dist(0.0f, amplitude_range);
    std::uniform_real_distribution<float> phase_dist(0.0f, phase_range);
    std::uniform_real_distribution<float> noise(-0.1f, 0.1f);

    for (int n = 0; n < fourier_series_size; n++) {
        fourier_series[n].n_degree = n;
        fourier_series[n].angular_velocity = base_angular_velocity * n + noise(gen);
        fourier_series[n].amplitude = amplitude_dist(gen);
        fourier_series[n].phase = phase_dist(gen);
    }
}

void discrete_fourier_transform(std::vector<fourier_component>& components,
                                const std::map<float, float>& positions) {
    const int N = positions.size();
    if (N == 0) return;

    // Extract positions into an array
    std::vector<float> samples;
    samples.reserve(N);
    for (const auto& [t, x] : positions) {
        samples.push_back(x);
    }

    // Compute dt assuming uniform sampling
    auto it = positions.begin();
    float t0 = it->first;
    std::advance(it, 1);
    float t1 = it->first;
    float dt = t1 - t0;

    const float T = N * dt; // total duration
    components.clear();
    components.reserve(N);

    for (int n = 0; n < N; ++n) {
        std::complex<float> sum(0.0f, 0.0f);
        for (int k = 0; k < N; ++k) {
            float angle = -2.0f * M_PI * n * k / N;
            sum += std::complex<float>(samples[k], 0.0f) * std::exp(std::complex<float>(0.0f, angle));
        }

        float amplitude = std::abs(sum) / N;
        float phase = std::arg(sum);
        float angular_velocity = 2.0f * M_PI * n / T;

        components.push_back({angular_velocity, amplitude * 2.0f, phase, n});  // multiply by 2 (except DC) if reconstructing real signal
    }
}

std::map<float, float> reconstruct_signal(
    const std::vector<fourier_component>& components,
    const std::map<float, float>& time_stamps){
    std::map<float, float> reconstructed;

    for (const auto& [t, _] : time_stamps) {
        float value = 0.0f;
        for (const auto& comp : components) {
            value += comp.amplitude * 0.5f * std::cos(comp.angular_velocity * t + comp.phase);
        }
        reconstructed[t] = value;
    }

    return reconstructed;
}

void PID_customed(const std::map<float, float>& positions, const std::string file_name) {
    /* PID control logic:
    - Compute running average of a window size (past samples only)
    - Compute the error between the running average and the current position
    - Adjust the proportional term based on the error
    - Take slope half point (between PI value and current value)
    - Smooth value of the PI output
    */

    const int window_size = 7;
    const float Kp = 1.0f;
    const float Ki = 0.01f;
    const float smoothing = 0.5f;

    std::deque<float> window;
    window.resize(window_size);
    float integral = 0.0f;
    float previous_output = 0.0f;
    bool initialized = false;

    std::ofstream outputFile(file_name);
    outputFile << "Time, Current Position, Running Average" << "(" << window_size << "), Error, PI Output, Smoothed Output\n";

    for (const auto& [time, current_position] : positions) {
        // Add to running window
        window.push_back(current_position);
        if (window.size() > window_size) {
            window.pop_front();
        }

        // Compute running average
        float avg = std::accumulate(window.begin(), window.end(), 0.0f) / window.size();

        // Compute error
        float error = avg - current_position;

        // Update integral
        integral += error;

        // PI controller output
        float pi_output = Kp * error + Ki * integral;

        // Midpoint slope between previous and current PI output
        float mid = (pi_output + previous_output) / 2.0f;

        // Smooth output
        float smoothed = initialized ? smoothing * previous_output + (1.0f - smoothing) * mid : mid;
        initialized = true;

        previous_output = smoothed;

        // Write to CSV
        outputFile << time << "," << current_position << "," << avg << "," << error << "," << pi_output << "," << smoothed << "\n";
    }
}

void Kalman_filter(const std::map<float, float>& positions, const std::string file_name) {
    // PID Coefficients
    const int window_size = 5;
    const float Kp = 1.0f;
    const float Ki = 0.1f;
    const float smoothing = 0.5f;

    // Kalman Filter Initialization (1D)
    float x_hat = 0.0f;           // Estimated state
    float P = 1.0f;               // Estimate uncertainty
    const float Q = 1e-5f;        // Process noise covariance
    const float R = 1e-2f;        // Measurement noise covariance

    // Others
    std::deque<float> window;
    float integral = 0.0f;
    float previous_output = 0.0f;
    bool initialized = false;

    std::ofstream outputFile(file_name);
    outputFile << "Time, Current Position, Running Average" << "(" << window_size << "), Error, PI Output, Smoothed Output, Measurement, PID Out, Kalman Estimate, Kalman Gain\n";

    for (const auto& [time, measurement] : positions) {
        // ---- Running Average ----
        window.push_back(measurement);
        if (window.size() > window_size) {
            window.pop_front();
        }

        float average = std::accumulate(window.begin(), window.end(), 0.0f) / window.size();

        // ---- PID Calculation ----
        float error = average - measurement;
        integral += error;
        float pi_output = Kp * error + Ki * integral;

        float mid = (pi_output + previous_output) / 2.0f;
        float smoothed_pid = initialized ? smoothing * previous_output + (1 - smoothing) * mid : mid;
        initialized = true;
        previous_output = smoothed_pid;

        // ---- Kalman Filter ----
        // Prediction
        float x_hat_prior = x_hat;
        float P_prior = P + Q;

        // Measurement update
        float K = P_prior / (P_prior + R);
        x_hat = x_hat_prior + K * (measurement - x_hat_prior);
        P = (1 - K) * P_prior;

        // Write to CSV
        outputFile << time << "," << measurement << "," << average << "," << error << "," << pi_output << "," 
        << smoothed_pid << "," << measurement << "," << smoothed_pid << "," << x_hat << "," << K << "\n";
    }
}

int main() {
    const float pi = 3.14159f;
    const float angular_velocity = 5 * pi;
    const float amplitude = 10;
    const float phase = 2 * pi; // phase range now 0 to pi
    const float time_span = 100;
    const float time_step = 0.1f;
    const int max_degree = 25;

    fourier_component fourier_series[max_degree];
    generate_components(fourier_series, max_degree, angular_velocity, amplitude, phase);
    std::map<float, float> positions = compute_trajectory(0, fourier_series, max_degree, time_span, time_step);
    report_trajectory(fourier_series, max_degree, positions, "trajectory_report.md");

    // ----- Inverting signal to build Fourier model -----
    
    std::vector<fourier_component> discrete_fourier_series;

    discrete_fourier_transform(discrete_fourier_series, positions);
    std::map<float, float> reconstructed_positions = reconstruct_signal(discrete_fourier_series, positions);

    // ----- PID controller -----
    PID_customed(positions, "pid1.csv");
    PID_customed(reconstructed_positions, "pid2.csv");

    // ----- Kalman Filter -----
    Kalman_filter(positions, "kalman1.csv");
    Kalman_filter(reconstructed_positions, "kalman2.csv");

    return 0;
}
