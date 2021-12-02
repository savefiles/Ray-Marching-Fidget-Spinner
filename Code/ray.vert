/*
Copyright 2019
Original authors: Niko Procopi
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
<http://www.gnu.org/licenses/>.

Special Thanks to Exzap from Team Cemu,
he gave me advice on how to optimize Vulkan
graphics, he is working on a Wii U emulator
that utilizes Vulkan, see more at http://cemu.info
*/

// https://stackoverflow.com/questions/2588875/whats-the-best-way-to-draw-a-fullscreen-quad-in-opengl-3-2

#version 450

layout (location = 0) out vec2 outUV;

void main() 
{	
	// Create the vertices array (have to assign one by one).
	vec2 vertices[3];
	vertices[0] = vec2(-1.0, 3.0);
	vertices[1] = vec2(3.0, -1.0);
	vertices[2] = vec2(-1.0, -1.0);

	gl_Position = vec4(vertices[gl_VertexIndex], 0.0, 1.0);
	outUV = 0.5 * gl_Position.xy + vec2(0.5);
}