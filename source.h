#ifndef SOURCE_H
#define SOURCE_H

#include <vector>

class Source {
public:
    std::vector<float> time;
    std::vector<float> rate;
    float interpolated_rate;

    void ReadSource(std::string filename);
    void InterpolateRate(float sim_time);
  Source(void) : interpolated_rate(0.0) {};
};

#endif
