/*

Why aren't we usting an IOT solution? Seems fundamentally like an IOT problem.

just looking at getting a mantid conda setup, got files copied over, good that the SNS system
is mounted so I actually need
that mount to be able to run the workflow. I'm working on adding stories/tasks. And the current
tasks I'm working on getting mantid conda setup locally, to make sure it's repeatable etc, then
next week hopefully will be begin encapsulating it into a galaxy tool. That's when the tasks for
having a nexus file reader will come into play.

add up these 
intersect is instrument focused

message based instead of API to include the custom OS's on the instruments into the
communication interface

limiting flexibility?

*/
#include <iostream>
#include <sstream>

#include "nlohmann/json.hpp"

using json = nlohmann::json;

int main()
{
  // a JSON text
  auto text = R"(
              {
                 "Image": {
                   "Width":  800,
                   "Height": 600,
                   "Title":  "View from 15th Floor",
                   "Thumbnail": {
                   "Url":    "http://www.example.com/image/481989943",
                                                                                                            "Height": 125,
            "Width":  100
            },
            "Animated" : false,
            "IDs": [116, 943, 234, 38793]
            }
            }
            )";

  // fill a stream with JSON text
  std::stringstream ss;
  ss << text;

  // parse and serialize JSON
  json j_complete = json::parse(ss);

  return 0;
}

