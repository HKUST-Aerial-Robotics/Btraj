#include <functional>
#include <arc_utilities/arc_helpers.hpp>
#include <arc_utilities/eigen_helpers.hpp>

#ifndef SHORTCUT_SMOOTHING_HPP
#define SHORTCUT_SMOOTHING_HPP

namespace shortcut_smoothing
{
    EigenHelpers::VectorVector3d ShortcutSmooth(const EigenHelpers::VectorVector3d& input_vector, const size_t first_ind, const size_t second_ind, const double step_size, const std::function<bool(const Eigen::Vector3d&)>& collision_fn)
    {
        const size_t starting_ind = std::min(first_ind, second_ind);
        const size_t ending_ind = std::max(first_ind, second_ind);

        const Eigen::Vector3d& starting_point = input_vector[starting_ind];
        const Eigen::Vector3d& ending_point = input_vector[ending_ind];
        const Eigen::Vector3d delta = ending_point - starting_point;
        const Eigen::Vector3d delta_unit_vec = delta.normalized();

        const double total_dist = delta.norm();

        if (total_dist <= step_size)
        {
            return input_vector;
        }

        bool collision = false;
        for (double dist = step_size; !collision && dist < total_dist; dist += step_size)
        {
            const Eigen::Vector3d point_to_check = starting_point + dist * delta_unit_vec;
            collision = collision_fn(point_to_check);
        }


        if (!collision)
        {
            const size_t num_original_elements = ending_ind - starting_ind - 1;
            const size_t num_new_elements = (size_t)std::ceil(total_dist / step_size);

            EigenHelpers::VectorVector3d output_vector(input_vector.size() - num_original_elements + num_new_elements - 1);

            // Copy in the first unchanged elements of the vector
            std::copy(input_vector.begin(), input_vector.begin() + starting_ind + 1, output_vector.begin());

            // Assign the replaced elements
            for (size_t new_element_ind = 1; new_element_ind < num_new_elements; ++new_element_ind)
            {
                const double dist = (double)new_element_ind * step_size;
                output_vector[starting_ind + new_element_ind] = starting_point + dist * delta_unit_vec;
            }

            // Copy in the last unchanged elements of the vector
            std::copy(input_vector.begin() + ending_ind, input_vector.end(), output_vector.begin() + starting_ind + num_new_elements);

            return output_vector;

        }
        else
        {
            return input_vector;
        }
    }

    template<typename PRNG, typename Configuration, typename ConfigAlloc=std::allocator<Configuration>>
    std::vector<Configuration, ConfigAlloc> ShortcutSmoothPath(const std::vector<Configuration, ConfigAlloc>& path, const uint32_t max_iterations, const uint32_t max_failed_iterations, const double max_shortcut_fraction, const std::function<bool(const Configuration&, const Configuration&)>& edge_validity_check_fn, PRNG& prng)
    {
        std::vector<Configuration, ConfigAlloc> current_path = path;
        uint32_t num_iterations = 0;
        uint32_t failed_iterations = 0;
        while (num_iterations < max_iterations && failed_iterations < max_failed_iterations && current_path.size() > 2)
        {
            // Attempt a shortcut
            const int64_t base_index = std::uniform_int_distribution<size_t>(0, current_path.size() - 1)(prng);
            // Pick an offset fraction
            const double offset_fraction = std::uniform_real_distribution<double>(-max_shortcut_fraction, max_shortcut_fraction)(prng);
            // Compute the offset index
            const int64_t offset_index = base_index + (int64_t)((double)current_path.size() * offset_fraction); // Could be out of bounds
            const int64_t safe_offset_index = arc_helpers::ClampValue(offset_index, (int64_t)0, (int64_t)(current_path.size() - 1));
            // Get start & end indices
            const size_t start_index = (size_t)std::min(base_index, safe_offset_index);
            const size_t end_index = (size_t)std::max(base_index, safe_offset_index);
            if (end_index <= start_index + 1)
            {
                num_iterations++;
                continue;
            }
            // Check if the edge is valid
            const Configuration& start_config = current_path[start_index];
            const Configuration& end_config = current_path[end_index];
            const bool edge_valid = edge_validity_check_fn(start_config, end_config);
            if (edge_valid)
            {
                // Make the shortened path
                std::vector<Configuration, ConfigAlloc> shortened_path;
                shortened_path.insert(shortened_path.end(), current_path.begin(), current_path.begin() + start_index + 1);
                shortened_path.insert(shortened_path.end(), current_path.begin() + end_index, current_path.end());
                current_path = shortened_path;
                num_iterations++;
            }
            else
            {
                num_iterations++;
                failed_iterations++;
            }
        }
        return current_path;
    }

    template<typename Configuration, typename ConfigAlloc=std::allocator<Configuration>>
    std::vector<Configuration, ConfigAlloc> ResamplePath(const std::vector<Configuration, ConfigAlloc>& path, const double resampled_state_distance, const std::function<double(const Configuration&, const Configuration&)>& state_distance_fn, const std::function<Configuration(const Configuration&, const Configuration&, const double)>& state_interpolation_fn)
    {
        if (path.size() <= 1)
        {
            return path;
        }
        std::vector<Configuration, ConfigAlloc> resampled_path;
        // Add the first state
        resampled_path.push_back(path[0]);
        // Loop through the path, adding interpolated states as needed
        for (size_t idx = 1; idx < path.size(); idx++)
        {
            // Get the states from the original path
            const Configuration& previous_state = path[idx - 1];
            const Configuration& current_state = path[idx];
            // We want to add all the intermediate states to our returned path
            const double distance = state_distance_fn(previous_state, current_state);
            const double raw_num_intervals = distance / resampled_state_distance;
            const uint32_t num_segments = (uint32_t)std::ceil(raw_num_intervals);
            // If there's only one segment, we just add the end state of the window
            if (num_segments == 0u)
            {
                // Do nothing because this means distance was exactly 0
            }
            else if (num_segments == 1u)
            {
                // Add a single point for the other end of the segment
                resampled_path.push_back(current_state);
            }
            // If there is more than one segment, interpolate between previous_state and current_state (including the current_state)
            else
            {
                for (uint32_t segment = 1u; segment <= num_segments; segment++)
                {
                    const double interpolation_ratio = (double)segment / (double)num_segments;
                    const Configuration interpolated_state = state_interpolation_fn(previous_state, current_state, interpolation_ratio);
                    resampled_path.push_back(interpolated_state);
                }
            }
        }
        return resampled_path;
    }
}

#endif // SHORTCUT_SMOOTHING_HPP
