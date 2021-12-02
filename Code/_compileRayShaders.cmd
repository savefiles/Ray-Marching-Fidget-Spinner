..\Bin\glslangValidator.exe -V ray.vert -o ray.vert.spv
..\Bin\glslangValidator.exe -V ray.frag -o ray.frag.spv
..\Bin\spirv-opt --strip-debug ray.vert.spv -o ray2.vert.spv
..\Bin\spirv-opt --strip-debug ray.frag.spv -o ray2.frag.spv
bin2hex --i ray2.vert.spv --o ray.vert.inc
bin2hex --i ray2.frag.spv --o ray.frag.inc
del ray.vert.spv
del ray.frag.spv
del ray2.vert.spv
del ray2.frag.spv

pause