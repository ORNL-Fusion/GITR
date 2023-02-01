
# Steps to refactor the code

```ruby 
   gitr.cpp
  
   - task 1: create an input.cpp to parse the gitrInput.cfg
   - task 2: can we add execute_command (strcmp(command)) that will dispatch to the respective class. 
   This will essentially allow us to reduce by half the input file. 
   (see here: https://github.com/lammps/lammps/blob/stable/src/input.cpp)
   - task 3: create an dump.cpp/diagnostic.cpp: this file will contain functions to write output files.
   
Task 1 and 3 are currently being done inside gitr.cpp, it is just a question of moving these functionalities into "separate" files. 

