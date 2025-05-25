#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <fstream>
#include <complex>
#include <deque>
#include <numeric>

const float cap_threshold = 10.0f;

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
        positions[time] = abs(position) > cap_threshold ? cap_threshold * (abs(position) / position) : position;
    }

    return positions;
}

void report_trajectory(fourier_component fourier_series[], const int& fourier_series_size,
                       const std::map<float, float>& positions, const std::string& filename) {
    // Component report
    std::ofstream componentFile(filename + "component.csv");
    componentFile << "n,angular_velocity,amplitude,phase\n";

    for (int i = 0; i < fourier_series_size; i++) {
        componentFile << fourier_series[i].n_degree << "," << fourier_series[i].angular_velocity << "," << fourier_series[i].amplitude << "," << fourier_series[i].phase << "\n";
    }

    componentFile.close();

    // Position report
    std::ofstream positionFile(filename + "position.csv");
    positionFile << "time,position\n";

    for (const auto& position : positions) {
        positionFile << position.first << "," << position.second << "\n";
    }

    positionFile.close();
}

void generate_components(fourier_component fourier_series[], const int& fourier_series_size,
                         const float& base_angular_velocity, const float& amplitude_range, const float& phase_range) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<float> amplitude_dist(0.0f, amplitude_range);
    std::uniform_real_distribution<float> phase_dist(0.0f, phase_range);
    std::uniform_real_distribution<float> noise(-0.1f, 0.1f);
    // std::uniform_int_distribution<int> angular_velocity_dist(0, fourier_series_size);
    std::uniform_real_distribution<float> angular_velocity_dist(-base_angular_velocity, base_angular_velocity);

    for (int n = 0; n < fourier_series_size; n++) {
        fourier_series[n].n_degree = n;
        // fourier_series[n].angular_velocity = angular_velocity_dist(gen) * base_angular_velocity;
        fourier_series[n].angular_velocity = angular_velocity_dist(gen);
        fourier_series[n].amplitude = amplitude_dist(gen);
        fourier_series[n].phase = phase_dist(gen);
    }
}

void discrete_fourier_transform(std::vector<fourier_component>& components,
                                const std::map<float, float>& positions) {
    const int N = positions.size();
    if (N == 0) {
        std::cout << "No data points provided." << std::endl;
        return;
    }

    // Extract samples
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

    const float T = N * dt;
    components.clear();
    components.reserve(N);

    for (int n = 0; n < N; ++n) {
        std::complex<float> sum(0.0f, 0.0f);
        for (int k = 0; k < N; ++k) {
            float angle = -2.0f * M_PI * n * k / N;
            sum += std::complex<float>(samples[k], 0.0f) * std::exp(std::complex<float>(0.0f, angle));
        }

        float amplitude = 2 * std::abs(sum) / N;
        float phase = std::arg(sum);
        float angular_velocity = 2.0f * M_PI * n / T;

        float scaled_amplitude = (n == 0) ? amplitude : amplitude * 2.0f;
        components.push_back({angular_velocity, scaled_amplitude, phase, n});
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

    const int window_size = 5;
    const float Kp = 1.0f;
    const float Ki = 0.01f;
    const float smoothing = 0.5f;

    std::deque<float> window;
    window.resize(window_size);
    float integral = 0.0f;
    float previous_output = 0.0f;
    bool initialized = false;

    std::ofstream outputFile(file_name);
    outputFile << "time, running_average" << "(" << window_size << "), position_error, PI_output, smoothed_output, total_output\n";

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
        outputFile << time  << "," << avg << "," << error << "," << pi_output << "," << smoothed << "," << -(pi_output+smoothed) << "\n";
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
    outputFile << "time, running_average" << "(" << window_size << "), position_error, "
               << "PI_output, smoothed_output, PID_out, kalman_estimate, kalman_gain\n";

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
        outputFile << time << "," << average << "," << -error << "," << pi_output << "," 
        << smoothed_pid << "," << smoothed_pid << "," << x_hat << "," << K << "\n";
    }
}

void alternative_Kalman_filter(const std::map<float, float>& positions, const std::string file_name) {
    if (positions.size() < 2) {
        std::cout << "Not enough data for filtering." << std::endl;
        return;
    }

    std::ofstream out(file_name);
    if (!out.is_open()) {
        std::cerr << "Failed to open file for writing: " << file_name << std::endl;
        return;
    }

    // Header for CSV output
    out << "time,predicted_position,filtered_position,position_error,"
        << "kalman_gain_position,kalman_gain_velocity\n";

    // Initial state: [position, velocity]
    float x = 0.0f;       // Estimated position
    float v = 0.0f;       // Estimated velocity

    // Initial estimation covariance matrix
    float p11 = 1.0f, p12 = 0.0f;
    float p21 = 0.0f, p22 = 1.0f;

    // Process noise covariance (Q) and measurement noise (R)
    const float q = 0.01f;
    const float r = 1.0f;

    auto it = positions.begin();
    float prev_time = it->first;
    float z = it->second;

    // Initialize state from first measurement
    x = z;
    v = 0.0f;
    ++it;

    for (; it != positions.end(); ++it) {
        float curr_time = it->first;
        float dt = curr_time - prev_time;
        z = it->second;

        // ---------- Predict ----------
        float x_pred = x + v * dt;
        float v_pred = v;

        float p11_pred = p11 + dt * (p12 + p21 + dt * p22) + q;
        float p12_pred = p12 + dt * p22;
        float p21_pred = p21 + dt * p22;
        float p22_pred = p22 + q;

        // ---------- Update ----------
        float y = z - x_pred; // innovation (error)
        float s = p11_pred + r;
        float k1 = p11_pred / s;
        float k2 = p21_pred / s;

        x = x_pred + k1 * y;
        v = v_pred + k2 * y;

        p11 = (1 - k1) * p11_pred;
        p12 = (1 - k1) * p12_pred;
        p21 = -k2 * p11_pred + p21_pred;
        p22 = -k2 * p12_pred + p22_pred;

        // Output to file
        out << curr_time << "," << x_pred << "," << x << "," << y << "," << k1 << "," << k2 << "\n";

        prev_time = curr_time;
    }

    out.close();
}

void Luenberger_observer(const std::map<float, float>& positions, const std::string file_name) {
    if (positions.size() < 2) return;

    std::ofstream out(file_name);
    out << "time,measured_position,estimated_position\n";

    float x = 0.0f;   // estimated position
    float v = 0.0f;   // estimated velocity

    // Observer gains (tuned manually)
    const float L1 = 0.5f;
    const float L2 = 0.1f;

    auto it = positions.begin();
    float prev_time = it->first;
    float z = it->second;
    x = z;

    ++it;
    for (; it != positions.end(); ++it) {
        float t = it->first;
        float dt = t - prev_time;
        z = it->second;

        // Prediction
        float x_pred = x + v * dt;

        // Correction
        float error = z - x_pred;
        x = x_pred + L1 * error;
        v = v + L2 * error;

        out << t << "," << z << "," << x << "\n";
        prev_time = t;
    }

    out.close();
}

void LQR_observer(const std::map<float, float>& positions, const std::string file_name) {
    if (positions.size() < 2) return;

    std::ofstream out(file_name);
    out << "time,measured_position,estimated_position\n";

    float x = 0.0f;
    float v = 0.0f;

    // Tuned LQR-style gains (derived from solving Riccati in real applications)
    const float K1 = 0.8f;
    const float K2 = 0.2f;

    auto it = positions.begin();
    float prev_time = it->first;
    float z = it->second;
    x = z;

    ++it;
    for (; it != positions.end(); ++it) {
        float t = it->first;
        float dt = t - prev_time;
        z = it->second;

        // Prediction
        float x_pred = x + v * dt;

        // Correction
        float error = z - x_pred;
        x = x_pred + K1 * error;
        v = v + K2 * error;

        out << t << "," << z << "," << x << "\n";
        prev_time = t;
    }

    out.close();
}

void extended_Kalman_filter(const std::map<float, float>& positions, const std::string file_name) {
    if (positions.size() < 2) return;

    std::ofstream out(file_name);
    out << "time,predicted_position,filtered_position,position_error\n";

    float x = 0.0f; // position
    float v = 0.0f; // velocity

    float p11 = 1.0f, p12 = 0.0f;
    float p21 = 0.0f, p22 = 1.0f;

    const float q = 0.01f; // process noise
    const float r = 1.0f;  // measurement noise

    auto it = positions.begin();
    float prev_time = it->first;
    float z = it->second;
    x = z;

    ++it;
    for (; it != positions.end(); ++it) {
        float t = it->first;
        float dt = t - prev_time;
        z = it->second;

        // Jacobian H = [1 0], system F = [[1 dt], [0 1]]

        // Predict
        float x_pred = x + v * dt;
        float v_pred = v;

        float p11_pred = p11 + dt * (p12 + p21 + dt * p22) + q;
        float p12_pred = p12 + dt * p22;
        float p21_pred = p21 + dt * p22;
        float p22_pred = p22 + q;

        // Update
        float y = z - x_pred;
        float s = p11_pred + r;
        float k1 = p11_pred / s;
        float k2 = p21_pred / s;

        x = x_pred + k1 * y;
        v = v_pred + k2 * y;

        p11 = (1 - k1) * p11_pred;
        p12 = (1 - k1) * p12_pred;
        p21 = -k2 * p11_pred + p21_pred;
        p22 = -k2 * p12_pred + p22_pred;

        out << t << "," << x_pred << "," << x << "," << y << "\n";
        prev_time = t;
    }

    out.close();
}

int main() {
    const float pi = 3.14159f;
    const float angular_velocity = 3 * pi;
    const float amplitude = 4;
    const float phase = 3 * pi; // phase range now 0 to pi
    const float time_span = 100;
    const float time_step = 0.5f;
    const int max_degree = 10;

    fourier_component fourier_series[max_degree];
    generate_components(fourier_series, max_degree, angular_velocity, amplitude, phase);
    std::map<float, float> positions = compute_trajectory(0, fourier_series, max_degree, time_span, time_step);
    report_trajectory(fourier_series, max_degree, positions, "trajectory_");

    // ----- Inverting signal to build Fourier model -----
    
    std::vector<fourier_component> discrete_fourier_series;
    discrete_fourier_series.reserve(max_degree);

    discrete_fourier_transform(discrete_fourier_series, positions);
    std::map<float, float> reconstructed_positions = reconstruct_signal(discrete_fourier_series, positions);
    for(int i = 0; i < max_degree; i++){
        fourier_series[i].amplitude = discrete_fourier_series[i].amplitude;
        fourier_series[i].angular_velocity = discrete_fourier_series[i].angular_velocity;
        fourier_series[i].phase = discrete_fourier_series[i].phase;
        fourier_series[i].n_degree = discrete_fourier_series[i].n_degree;
    }
    report_trajectory(fourier_series, max_degree, reconstructed_positions, "inverted_trajectory_");

    // ----- PID controller -----
    PID_customed(positions, "PID1.csv");
    PID_customed(reconstructed_positions, "PID2.csv");

    // ----- Kalman Filter -----
    Kalman_filter(positions, "kalman1.csv");
    Kalman_filter(reconstructed_positions, "kalman2.csv");

    // ----- Alternative Kalman Filter -----
    alternative_Kalman_filter(positions, "alternative_kalman1.csv");
    alternative_Kalman_filter(reconstructed_positions, "alternative_kalman2.csv");

    // ----- Luenberger observer -----
    Luenberger_observer(positions, "luenberger1.csv");
    Luenberger_observer(reconstructed_positions, "luenberger2.csv");

    // ----- LQR observer -----
    LQR_observer(positions, "LQR1.csv");
    LQR_observer(reconstructed_positions, "LQR2.csv");

    // ----- Extended Kalman Filter -----
    extended_Kalman_filter(positions, "extended_kalman1.csv");
    extended_Kalman_filter(reconstructed_positions, "extended_kalman2.csv");

    return 0;
}
