#include "fourier.h"

int main() {
    const float pi = 3.14159f;
    const float angular_velocity = 3 * pi;
    const float amplitude = 4;
    const float phase = 3 * pi; // phase range now 0 to pi
    const float time_span = 100;
    const float time_step = 0.5f;
    const int max_degree = 10;
    const std::string root = "fourier_data/";

    fourier_component fourier_series[max_degree];
    generate_components(fourier_series, max_degree, angular_velocity, amplitude, phase);
    std::map<float, float> positions = compute_trajectory(0, fourier_series, max_degree, time_span, time_step);
    report_trajectory(fourier_series, max_degree, positions, root + "trajectory_");

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
    report_trajectory(fourier_series, max_degree, reconstructed_positions, root + "inverted_trajectory_");

    // ----- PID controller -----
    PID_customed(positions, root + "PID1.csv");
    PID_customed(reconstructed_positions, root + "PID2.csv");

    // ----- Kalman Filter -----
    Kalman_filter(positions, root + "kalman1.csv");
    Kalman_filter(reconstructed_positions, root + "kalman2.csv");

    // ----- Alternative Kalman Filter -----
    alternative_Kalman_filter(positions, root + "alternative_kalman1.csv");
    alternative_Kalman_filter(reconstructed_positions, root + "alternative_kalman2.csv");

    // ----- Luenberger observer -----
    Luenberger_observer(positions, root + "luenberger1.csv");
    Luenberger_observer(reconstructed_positions, root + "luenberger2.csv");

    // ----- LQR observer -----
    LQR_observer(positions, root +  "LQR1.csv");
    LQR_observer(reconstructed_positions, root +  "LQR2.csv");

    // ----- Extended Kalman Filter -----
    extended_Kalman_filter(positions, root +  "extended_kalman1.csv");
    extended_Kalman_filter(reconstructed_positions, root +  "extended_kalman2.csv");

    return 0;
}
