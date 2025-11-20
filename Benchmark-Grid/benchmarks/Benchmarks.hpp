#pragma once
#include "json.hpp"
#include "DeoFlops.hpp"

namespace Benchmark
{
    void Memory(nlohmann::json& json_results);
    void SU4(nlohmann::json& json_results);
    void Comms(nlohmann::json& json_results);
    void Latency(nlohmann::json& json_results);
    void P2P(nlohmann::json& json_results);
}
