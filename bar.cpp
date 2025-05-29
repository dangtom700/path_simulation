#include <iomanip>

#include "fourier.h"

struct linkage_bar {
    float length;
    float angle_rad;         // relative angle (joint's own angle)
    float angular_velocity;  // rad/s
    float x_pos;
    float y_pos;

    linkage_bar(float length, float angle_rad, float angular_velocity)
        : length(length < 0 ? 0 : length),
          angle_rad(angle_rad),
          angular_velocity(angular_velocity),
          x_pos(0.0f),
          y_pos(0.0f) {}

    void update_angle(float delta_time) {
        angle_rad += angular_velocity * delta_time;
        // Normalize between -PI and PI
        while (angle_rad > M_PI) {
            angle_rad -= 2 * M_PI;
        }
        while (angle_rad < -M_PI) {
            angle_rad += 2 * M_PI;
        }
    }
};

std::vector<linkage_bar> generate_bars(int num_bars, float max_length, float max_angle, float max_velocity) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis_length(1.0f, max_length);  // Avoid zero-length
    std::uniform_real_distribution<float> dis_angle(-max_angle, max_angle);
    std::uniform_real_distribution<float> dis_velocity(-max_velocity, max_velocity);

    std::vector<linkage_bar> bars;
    for (int i = 0; i < num_bars; ++i) {
        bars.emplace_back(dis_length(gen), dis_angle(gen), dis_velocity(gen));
    }
    return bars;
}

float compute_end_effector_distance(float x, float y) {
    return std::sqrt(x * x + y * y);
}

int main() {
    const int num_bars = 3;
    const float max_length = 7.0f;
    const float max_angle = 1.0f * M_PI;
    const float max_velocity = 3.0f;
    const int time_steps = 1000;
    const float delta_time = 0.1f;
    const std::string root = "bar_data/";
    std::map<float, float> positions;

    auto bars = generate_bars(num_bars, max_length, max_angle, max_velocity);

    std::ofstream csv_file(root + "joint_positions.csv");
    if (!csv_file) {
        std::cerr << "Failed to open CSV file." << std::endl;
        return 1;
    }
    std::ofstream angle_file(root + "joint_angles.csv");
    if (!angle_file) {
        std::cerr << "Failed to open angle file." << std::endl;
        return 1;
    }

    // Header
    csv_file << "time,";
    angle_file << "time,";
    for (int i = 0; i < num_bars; ++i) {
        csv_file << "x" << i << ",y" << i << ",";
        angle_file << "angle" << i << ",";
    }
    csv_file << "distance\n";
    angle_file << "distance\n";

    for (int step = 0; step < time_steps; ++step) {
        float time = step * delta_time;

        // 1. Update angles
        for (auto& bar : bars) {
            bar.update_angle(delta_time);
        }

        // 2. Compute absolute angles and positions
        float x = 0.0f;
        float y = 0.0f;

        for (size_t i = 0; i < bars.size(); ++i) {
            x += bars[i].length * std::cos(bars[i].angle_rad);
            y += bars[i].length * std::sin(bars[i].angle_rad);
            bars[i].x_pos = x;
            bars[i].y_pos = y;
        }

        // 3. Write data
        csv_file << std::fixed << std::setprecision(4) << time << ",";
        for (const auto& bar : bars) {
            csv_file << bar.x_pos << "," << bar.y_pos << ",";
        }
        float distance = compute_end_effector_distance(bars.back().x_pos, bars.back().y_pos);
        csv_file << distance << "\n";

        angle_file << std::fixed << std::setprecision(4) << time << ",";
        for (const auto& bar : bars) {
            angle_file << bar.angle_rad << ",";
        }
        angle_file << distance << "\n";

        positions[time] = distance;
    }

    csv_file.close();
    angle_file.close();

    // // ----- Inverting signal to build Fourier model -----
    // const int max_degree = num_bars;

    // std::vector<fourier_component> discrete_fourier_series;
    // discrete_fourier_series.reserve(max_degree);
    // fourier_component fourier_series[max_degree];

    // discrete_fourier_transform(discrete_fourier_series, positions);
    // std::map<float, float> reconstructed_positions = reconstruct_signal(discrete_fourier_series, positions);
    // for(int i = 0; i < max_degree; i++){
    //     fourier_series[i].amplitude = discrete_fourier_series[i].amplitude;
    //     fourier_series[i].angular_velocity = discrete_fourier_series[i].angular_velocity;
    //     fourier_series[i].phase = discrete_fourier_series[i].phase;
    //     fourier_series[i].n_degree = discrete_fourier_series[i].n_degree;
    // }
    // report_trajectory(fourier_series, max_degree, reconstructed_positions, root + "inverted_trajectory_");

    // // ----- PID controller -----
    // PID_customed(positions, root + "PID1.csv");
    // PID_customed(reconstructed_positions, root + "PID2.csv");

    // // ----- Kalman Filter -----
    // Kalman_filter(positions, root + "kalman1.csv");
    // Kalman_filter(reconstructed_positions, root + "kalman2.csv");

    // // ----- Alternative Kalman Filter -----
    // alternative_Kalman_filter(positions, root + "alternative_kalman1.csv");
    // alternative_Kalman_filter(reconstructed_positions, root + "alternative_kalman2.csv");

    // // ----- Luenberger observer -----
    // Luenberger_observer(positions, root + "luenberger1.csv");
    // Luenberger_observer(reconstructed_positions, root + "luenberger2.csv");

    // // ----- LQR observer -----
    // LQR_observer(positions, root +  "LQR1.csv");
    // LQR_observer(reconstructed_positions, root +  "LQR2.csv");

    // // ----- Extended Kalman Filter -----
    // extended_Kalman_filter(positions, root +  "extended_kalman1.csv");
    // extended_Kalman_filter(reconstructed_positions, root +  "extended_kalman2.csv");

    return 0;
}
