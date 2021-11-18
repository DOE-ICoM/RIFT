#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <exception>
#include "source.h"

using namespace std;

void Source::ReadSource(std::string filename) {
    fstream file; 
    string line;
    int lno = 1;
    
    file.open(filename.c_str());
    if (file.is_open()) {
        while (getline(file, line)) {
            lno++;
            // Remove white space from the beginning of the string.
            /*line.erase(line.begin(), find_if(line.begin(), line.end(), 
              not1(ptr_fun<int, int>(isspace))));
            */
            // If the line of the file begins with '#', skip it.
            if (line[0] == '#') {
                continue;
            }
        
            // Read the comma-delimited source file
            std::string data;
            float value;
            int count = 0;
            stringstream linestream(line);
            while (getline(linestream, data, ',')) {
                stringstream(data) >> value;
                if (count == 0) {
                    time.push_back(value*60*60);
                } else { 
                    rate.push_back(value);
                    //std::cout << value << endl;
                }
                count++;
            }
        }
        file.close();
    } else {
        std::string msg(filename);
        msg += ": error: cannot open";
        throw std::runtime_error(msg);
    }

    // If the source file ends with a nonzero rate, append a zero rate
    if (rate[rate.size()-1] != 0.f) {
        // Estimate the timestep used throughout the source file
        float dt_estimate = 0.f;
        for (int i = 1; i < time.size(); i++) {
            dt_estimate += time[i]-time[i-1];
        }
        dt_estimate /= (float)time.size()-1;
        time.push_back(time[time.size()-1]+dt_estimate);
        rate.push_back(0.f);
    }
}

void Source::InterpolateRate(float sim_time) {
    for (int i = 0; i < time.size()-1; i++) {
        if (sim_time >= time[i] && sim_time < time[i+1]) {
            interpolated_rate = rate[i] + (rate[i+1]-rate[i]) /
			                              (time[i+1]-time[i]) *
			                              (sim_time-time[i]);
        }
    }
}

/* Local variables: */
/* mode: c++ */
/* tab-width: 4 */
/* c-basic-offset: 4 */
/* End: */
