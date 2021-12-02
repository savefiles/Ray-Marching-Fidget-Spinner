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

#version 450

layout (location = 0) in vec2 uv;

layout (pixel_center_integer) in vec4 gl_FragCoord;

layout (location = 0) out vec4 outColor;

layout (std140, binding = 0) uniform bufferVals {
    vec3 cameraPos;
	vec3 cameraForward;
	float aspectRatio;
	float time;
};


// - - - - - - - -//
// MATH FUNCTIONS //
// - - - - - - - -//

// params:
// p: arbitrary point in 3D space
// c: the center of our sphere
// r: the radius of our sphere
float sdfSphere(in vec3 p, in vec3 c, float r)
{
	return length(p - c) - r;
}

// params:
// p: point in 3D space
// c: center of the box
// b: dimensions of the box
float sdfBox(in vec3 p, in vec3 c, in vec3 b) 
{
	vec3 q = abs(p - c) - b;
	return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// p: point in 3D space
// c: center of torus
// R: inner radius (radius of inner circle)
// r: outer radius (radius of tube)
float sdfTorus(vec3 p, vec3 c, float R, float r) {
	vec3 rel = p-c;
	vec2 q = vec2(length(rel.xz)-R,rel.y);
  	return length(q)-r;
}

// Same as torus, except with an sc term.
// sc: not sure, controls the angle :shrug:
float sdfCappedTorus(vec3 p, vec3 c, float R, float r, vec2 sc) {
	vec3 rel = p-c;
	rel.x = abs(rel.x);
  	float k = (sc.y*rel.x>sc.x*rel.z) ? dot(rel.xz,sc) : length(rel.xz);
  	return sqrt( dot(rel,rel) + R*R - 2.0*R*k ) - r;
}

float sdfCylinder(vec3 p, vec3 c, float h, float r) {
	vec3 rel = p-c;
	vec2 d = abs(vec2(length(rel.xz),rel.y)) - vec2(r,h);
  	return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// Create a plane.
// n: normal of plane.
// h: height from center.
float sdfPlane(vec3 p, vec3 c, vec3 n, float h) {
	return dot(p-c, n) + h;
}

float sdfPyramid( vec3 p, vec3 c, float h)
{
	float m2 = h*h + 0.25;

	p = p-c;  
	p.xz = abs(p.xz);
	p.xz = (p.z>p.x) ? p.zx : p.xz;
	p.xz -= 0.5;

	vec3 q = vec3( p.z, h*p.y - 0.5*p.x, h*p.x + 0.5*p.y);
	
	float s = max(-q.x,0.0);
	float t = clamp( (q.y-0.5*p.z)/(m2+0.25), 0.0, 1.0 );
		
	float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
	float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);
		
	float d2 = min(q.y,-q.x*m2-q.y*0.5) > 0.0 ? 0.0 : min(a,b);
		
	return sqrt( (d2+q.z*q.z)/m2 ) * sign(max(q.z,-p.y));
}

// Polynomial smooth min function.
// A smoother min function, looks better for
// object blending.
float sminCubic(float a, float b, float k) {
	float h = max(k-abs(a-b), 0.0)/k;
	return min(a, b) - h*h*h*k*(1.0/6.0);
}

vec2 opBlend(vec2 d1, vec2 d2, float k) {
	float d = sminCubic(d1.x, d2.x, k);
	float m = mix(d1.y, d2.y, clamp(d1.x-d, 0.0, 1.0));
	return vec2(d, m);
}

vec2 opUnion(vec2 d1, vec2 d2) {
	return (d1.x < d2.x) ? d1 : d2;
}

// Subtract d2 from d1 (for easier chaining purposes).
vec2 opSubtract(vec2 d1, vec2 d2) {
	return (-d2.x > d1.x) ? vec2(-d2.x, d1.y) : d1;
}
// Rotate the point by the given number of radians around the origin. 
vec3 rotatePointOrigin(vec3 p, float rad) {
	float s = sin(rad);
	float c = cos(rad);
	return vec3(p.x * c - p.z * s, p.y, p.x * s + p.z * c);
}

// OLD SDF FUNCTIONS
// these don't take into account materials.
//float opUnion(float d1, float d2) {
//	return min(d1, d2);
//}
//
//float opIntersect(float d1, float d2) {
//	return max(d1, d2);
//}


// - - - - - - - - - - - - - -//
// COMPLETE MAPPING FUNCTIONS //
// - - - - - - - - - - - - - -//

vec2 mapFidgetSpinner(in vec3 p) {

	p = rotatePointOrigin(p, time);

	// MATERIAL LIST:
	// 1.0 - red
	// 2.0 - black
	vec2 torus1outside = vec2(sdfTorus(p, vec3(0.0, 0.0, 0.77), 0.4, 0.05), 1.0);
	vec2 torus1inside  = vec2(sdfTorus(p, vec3(0.0, 0.0, 0.77), 0.3, 0.05), 2.0);
	vec2 arm1 		   = opUnion(torus1outside, torus1inside);						 	// this is rendered.
	vec2 torus1volume  = vec2(sdfCylinder(p, vec3(0.0, 0.0, 0.77), 0.08, 0.4), 0.0); 	// used for subtract operations

	vec2 torus2outside = vec2(sdfTorus(p, vec3(0.69, 0.0, -0.40), 0.4, 0.05), 1.0);
	vec2 torus2inside  = vec2(sdfTorus(p, vec3(0.69, 0.0, -0.40), 0.3, 0.05), 2.0);
	vec2 arm2 		   = opUnion(torus2outside, torus2inside);
	vec2 torus2volume  = vec2(sdfCylinder(p, vec3(0.69, 0.0, -0.40), 0.08, 0.4), 0.0);	// used for subtract operations

	vec2 torus3outside = vec2(sdfTorus(p, vec3(-0.69, 0.0, -0.40), 0.4, 0.05), 1.0);
	vec2 torus3inside  = vec2(sdfTorus(p, vec3(-0.69, 0.0, -0.40), 0.3, 0.05), 2.0);
	vec2 arm3 		   = opUnion(torus3outside, torus3inside);
	vec2 torus3volume  = vec2(sdfCylinder(p, vec3(-0.67, 0.0, -0.39), 0.08, 0.4), 0.0);	// used for subtract operations

	vec2 middleWideCylinder = vec2(sdfCylinder(p, vec3(0.0, 0.0, 0.0), 0.05, 0.55), 1.0);
	vec2 middleTallCylinder = vec2(sdfCylinder(p, vec3(0.0, 0.0, -0.02), 0.09, 0.30), 1.0);

	// All three are used for subtractions, never rendered.
	vec2 bite1 = vec2(sdfCylinder(p, vec3(0.0, 0.0, -0.8),   0.10, 0.45), 0.0);
	vec2 bite2 = vec2(sdfCylinder(p, vec3(-0.69, 0.0, 0.40), 0.10, 0.45), 0.0);
	vec2 bite3 = vec2(sdfCylinder(p, vec3(0.69, 0.0, 0.40),  0.10, 0.45), 0.0);

	// Combine all of our objects.
	float k = 0.1;
	vec2 arms = opUnion(opUnion(arm1, arm2), arm3);
	vec2 centerObj = opSubtract(opSubtract(opSubtract(middleWideCylinder, bite1), bite2), bite3);
	centerObj = opSubtract(opSubtract(opSubtract(centerObj, torus1volume), torus2volume), torus3volume);
	centerObj = opUnion(centerObj, middleTallCylinder);

	return opUnion(arms, centerObj);
}

vec2 mapMandelbulb(vec3 p) {

	p = rotatePointOrigin(p, time);
	
    vec3 w = p;
    vec3 q = p;

    q.xz = mod( q.xz+1.0, 2.0 ) -1.0;

    float d = sdfBox(p, q, vec3(1.0));
    float s = 1.0;
    for( int m=0; m<7; m++ )
    {
        float h = float(m)/6.0;

        p =  q.yzx - 0.5*sin( 1.5*p.x + 6.0 + p.y*3.0 + float(m)*5.0 + vec3(1.0,0.0,0.0));

        vec3 a = mod( p*s, 2.0 )-1.0;
        s *= 3.0;
        vec3 r = abs(1.0 - 3.0*abs(a));

        float da = max(r.x,r.y);
        float db = max(r.y,r.z);
        float dc = max(r.z,r.x);
        float c = (min(da,min(db,dc))-1.0)/s;
        d = max( c, d );
   }
    
   return vec2(d*0.5, 1.0);
}

vec2 mapTerrain(vec3 p) {
	
	//float y = sin(p.x);
	// Small bumps.
	float y = 3.0 * (sin((p.x + 15.13) / 7.0) * sin((p.z - 22.03) / 11.0));
	// Larger bumps.
	y += 10.0 * sin((p.x - 39.04) / 33.0) * sin((p.z + 49.07) / 37.0);

	vec2 result = vec2(0.0, 3.0);
	result.x = p.y - y;
	//y = max(y, p.y);
	return result;

}

// - - - - - - - - - - - -//
// RAY MARCHING FUNCTIONS //
// - - - - - - - - - - - -//

// sdf function for the scene.
// p: point in 3d space.
// return: vec2.x is the distance, vec2.y is the material id.
vec2 map_the_world(in vec3 p) 
{
	//vec2 sphere1 = vec2(sdfSphere(p, vec3(-0.5, 0.0, -0.5), 1.0), 1.0);
	//vec2 sphere2 = vec2(sdfSphere(p, vec3(0.5, 0.0, 0.5), 1.0), 1.0);
	//return opUnion(sphere1, sphere2);
	//return opBlend(sphere1, sphere2, 1.0);

	//float displacement = (p.x + p.y + p.z) * 0.2;
	//float displacement = sin(5.0 * p.x) * sin(5.0 * p.y) * sin(5.0 * p.z) * 0.25;
	//return opBlend(sphere1, sphere2, 1.0) + vec2(displacement, 0.0);

	//return mapFidgetSpinner(p);
	
	vec2 fidgetSpinner = mapFidgetSpinner(p);
	vec2 planeFloor = vec2(sdfPlane(p, vec3(0.0, -3.0, 0.0), vec3(0.0, 1.0, 0.0), 1.0), 3.0);
	vec2 pyramid = vec2(sdfPyramid(p, vec3(0.0, -4.0, 0.0), 4.0), 3.0);
	return opUnion(opUnion(fidgetSpinner, planeFloor), pyramid);

	//vec2 terrain = mapTerrain(p);
	//return terrain;

	//vec2 sdfPyramid()
	//return mapMandelbulb(p);
}


// calculate the gradient.
vec3 calculate_normal(in vec3 p, in float distToCamera)
{
    float h = distToCamera * 0.000001; // replace by an appropriate value
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*map_the_world( p + k.xyy*h ).x + 
                      k.yyx*map_the_world( p + k.yyx*h ).x + 
                      k.yxy*map_the_world( p + k.yxy*h ).x + 
                      k.xxx*map_the_world( p + k.xxx*h ).x );
}

// Cast a singular ray from the ro in the direction of rd.
vec2 cast_ray(in vec3 ro, in vec3 rd) {
	
	float total_distance_traveled = 0.0;
    const int NUMBER_OF_STEPS = 256;
    const float MINIMUM_HIT_DISTANCE = 0.001;
    const float MAXIMUM_TRACE_DISTANCE = 500.0;
	
	// a material of 0 means no material.
	vec2 result = vec2(-1.0, 0.0);

    for (int i = 0; i < NUMBER_OF_STEPS; ++i)
    {
        // We wrote this function earlier in the tutorial -
        // assume that the sphere is centered at the origin
        // and has unit radius
        vec2 distance_to_closest = map_the_world(ro + total_distance_traveled * rd);

        if (distance_to_closest.x < MINIMUM_HIT_DISTANCE) // hit
        {
			result.x = total_distance_traveled;
			result.y = distance_to_closest.y;
			return result;
        }

        // accumulate the distance traveled thus far
        total_distance_traveled += distance_to_closest.x;

        if (total_distance_traveled > MAXIMUM_TRACE_DISTANCE) // complete miss
        {
			result.x = -1.0;
            return result;
        }
	}
	result.x = -1.0;
	return result;
}


// - - - - - - - - - -//
// LIGHTING FUNCTIONS //
// - - - - - - - - - -//

// Light struct
struct Light 
{
	vec3 pos;
	vec3 diffuse;
	vec3 ambient;
	vec3 specular;
	float shininess;
};

// Basic point light
// http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/#model-transformations
// ambient -> ambient color of light.
// diffuse -> main color of light.
// specular -> specular (shiny point) color of light.
// shine -> shineness coefficient
vec3 phongPointLight(Light light, vec3 L, vec3 litPointPos, float distToCamera) {
	
	vec3 N = calculate_normal(litPointPos, distToCamera);
	vec3 V = normalize(cameraPos - litPointPos);
	vec3 R = normalize(reflect(-L, N));

	float dotLN = dot(L, N);
	float dotRV = dot(R, V);

	if(dotLN < 0.0) {
		// Light is not visible from this surface.
		return light.ambient;
	}

	if(dotRV < 0.0) {
		// Light reflection is in opposite direction as viewer. Show diffuse only.
		return light.ambient + light.diffuse * dotLN;
	}

	return light.ambient + (light.diffuse * dotLN) + (light.specular * pow(dotRV, light.shininess));
}

// Same as phong point, but without any specular.
vec3 phongDirectionalLight(Light light, vec3 L, vec3 N) {
	
	float dotLN = dot(L, N);
	dotLN = clamp(dotLN, 0.0, 1.0);

	return light.ambient + (light.diffuse * dotLN);

}


// Function to create soft shadows at the current point.
// Is essentially a raymarch but from the light instead of camera.
// https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
// ro -> ray origin, which is each point.
// rd -> ray direction, which is direction from point to light.
// minDist -> the initial distance to move by.
// maxDist -> the maximum distance that the ray can go (dist between point and light).
float soft_shadow( in vec3 ro, in vec3 rd, in float minDist, in float maxDist)
{
    const float MINIMUM_HIT_DISTANCE = 0.0005;
	const float k = 4.0;	// controls how hard the shadows are.
    float res = 1.0;
    float ph = 1e20;
    for(float t = minDist; t < maxDist;)
    {
        float h = map_the_world(ro + rd * t).x;
        if(h < MINIMUM_HIT_DISTANCE)
            return 0.0;
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min( res, k*d/max(0.0,t-y) );
        ph = h;
        t += h;
    }
    return res;
}

// https://www.youtube.com/watch?v=XBF-g4y4clg
float fog_value (in float dist) {
	float fogStartDist = 1000.0;
	float fogAmount = 1.0 - exp((8.0-dist)/fogStartDist);
	return fogAmount;
}


// Checker board color pattern from Inigo Quilez.
float checkers(vec2 p)
{
    vec2 w = fwidth(p) + 0.001;
    vec2 i = 2.0*(abs(fract((p-0.5*w)*0.5)-0.5)-abs(fract((p+0.5*w)*0.5)-0.5))/w;
    return 0.5 - 0.5*i.x*i.y;
}

// A material ID to color lookup function for the fidget spinner.
vec3 colorLookupFidgetSpinner(float matID, vec3 current_position) {
	
	if (matID == 1.0) {
		// red
		return vec3(1.0, 0.0, 0.0);
	}
	else if (matID == 2.0) {
		// black
		return vec3(0.0, 1.0, 0.0);
	}
	else if (matID == 3.0) {
		// checkerboard pattern
		return vec3(clamp(checkers(vec2(current_position.x, current_position.z)), 0.45, 0.55));
	}
	else {
		// something went wrong.
		return vec3(1.0, 0.0, 1.0);
	}

}



// - - - - - - - - - - - - //
// MAIN LIGHTING FUNCTIONS //
// - - - - - - - - - - - - //

// This function calculates the point's color, and takes care of lighting too.
vec3 calcPointColorThreeSpot(vec2 result, vec3 current_position, float distToCamera) {

	// We will use the three point lighting setup.
	Light lightArray[3];

	// Key light.
	Light lightKey;
	lightKey.pos = vec3(6.0, 2.0, 6.0);
	lightKey.diffuse = vec3(1.0, 1.0, 1.0);
	//lightKey.ambient = vec3(0.15, 0.15, 0.15);
	lightKey.ambient = vec3(0.05, 0.05, 0.05);
	lightKey.specular = vec3(0.6, 0.6, 0.6);
	lightKey.shininess = 2.0;
	lightArray[0] = lightKey;
	
	// Fill light.
	Light lightFill;
	lightFill.pos = vec3(6.0, 2.0, -6.0);
	lightFill.diffuse = vec3(0.6, 0.6, 0.6);
	//lightFill.ambient = vec3(0.05, 0.05, 0.05);
	lightFill.ambient = vec3(0.0, 0.0, 0.0);
	lightFill.specular = vec3(0.3, 0.3, 0.3);
	lightFill.shininess = 2.0;
	lightArray[1] = lightFill;

	// Back light.
	Light lightBack;
	lightBack.pos = vec3(-6.0, 2.0, 6.0);
	lightBack.diffuse = vec3(0.3, 0.3, 0.3);
	lightBack.ambient = vec3(0.0, 0.0, 0.0);
	lightBack.specular = vec3(0.0, 0.0, 0.0);
	lightBack.shininess = 2.0;
	lightArray[2] = lightBack;

	// Loop through the lights, sum up the light values.
	vec3 colorFromLights = vec3(0.0);
	for(int i = 0; i < 3; ++i) 
	{
		
		// Calculate light information.
		vec3 toLight = lightArray[i].pos - current_position;
		float toLightDist = length(toLight);
		vec3 L = normalize(toLight);	// vector that points in the direction of the light.

		// Calculate color from light.
		vec3 pixelColor = 0.4 * phongPointLight(lightArray[i], L, current_position, distToCamera);
		
		// Calculate shadows.
		vec3 gradient_normal = calculate_normal(current_position, distToCamera);
		vec3 shadowStart = current_position + gradient_normal * 0.1;	// Moving the shadow start out a bit to prevent initial start intersection.
		float shadowVal = soft_shadow(shadowStart, L, 0.1, toLightDist);
		pixelColor = mix( pixelColor*0.01, pixelColor, shadowVal);
		
		// Apply fog.
		float fog = fog_value(distToCamera);
		pixelColor = mix(pixelColor, vec3(0.1, 0.1, 0.1), fog);

		colorFromLights += pixelColor;
	}

	// Lookup the point's surface color.
	vec3 pointSurfaceColor = colorLookupFidgetSpinner(result.y, current_position);
	
	// Return this pixel's color.
	return colorFromLights * pointSurfaceColor;

}

// This function calculates the point's color, and takes care of lighting too.
vec3 calcPointColorSun(vec2 result, vec3 current_position, float distToCamera) {


	// Key light.
	Light lightSun;
	//lightSun.pos = vec3(0.0, 5.0, 0.0);	// It's actually a direction
	lightSun.pos = vec3(-3.0, -3.0, -3.0);
	lightSun.diffuse = vec3(0.8, 0.8, 0.8);
	lightSun.ambient = vec3(0.15, 0.15, 0.15);

	
	// Calculate light information.
	vec3 toLight = -lightSun.pos;
	float toLightDist = 100.0;	// since this is towards the sun, it just has to break a certain threshold.
	vec3 L = normalize(toLight);	// vector that points in the direction of the light.

	// Calculate the normal.
	vec3 N = calculate_normal(current_position, distToCamera);

	// Calculate color from light.
	vec3 lightColor = 0.6 * phongDirectionalLight(lightSun, L, N);
	
	// Calculate shadows.
	vec3 shadowStart = current_position + (N * 0.05);	// Moving the shadow start out a bit to prevent initial start intersection.
	float shadowVal = soft_shadow(shadowStart, L, 0.001, toLightDist);
	lightColor = mix(lightColor*0.1, lightColor, shadowVal);
	
	// Apply fog.
	float fog = fog_value(distToCamera);
	lightColor = mix(lightColor, vec3(0.1, 0.1, 0.1), fog);

	// Lookup the point's surface color.
	vec3 pointSurfaceColor = colorLookupFidgetSpinner(result.y, current_position);
	
	// Return this pixel's color.
	return lightColor * pointSurfaceColor;
	//return lightColor;
}




// - - - - - - - //
// MAIN FUNCTION //
// - - - - - - - //

void main() 
{

	//// RAY MARCHING
	//// Shift UV to be from -1 to 1 on both axes.
	vec2 uv_mod = uv * 2.0 - 1.0;
	uv_mod.x *= aspectRatio;	// Correct for aspect ratio.


	// Calculate the ray direction.
	vec3 cameraRight = normalize(cross(cameraForward, vec3(0.0, 1.0, 0.0)));
	vec3 cameraUp = normalize(cross(cameraRight, cameraForward));
	vec3 rd = normalize((uv_mod.x * cameraRight) + (uv_mod.y * cameraUp) + (cameraForward * 2.0));


	// Get the ray origin.
	vec3 ro = cameraPos;

	// Cast a ray.
	vec2 t = cast_ray(ro, rd);
	if	(t.x == -1.0) {
		// There was no collision, don't render this fragment.	
		discard;
		//outColor = vec4(0.1, 0.5, 0.1, 0.5);
		//return;
	}
	vec3 current_position = (rd * t.x) + ro;


	// ISSUE: modColor.y is not returning correct value for positions under the camera.
	//vec3 modColor = mod(current_position, 1.0);
	//outColor = vec4(0.1, modColor.y, 0.1, 0.0);
	//return; 

	float distToCamera = length(current_position - cameraPos);

	// Color and shade the pixel.
	vec3 pixelColor = calcPointColorThreeSpot(t, current_position, distToCamera);
	//vec3 pixelColor = calcPointColorSun(t, current_position, distToCamera);


	// Gamma correct.
	pixelColor = pow(pixelColor, vec3(0.4545));

	// Output the pixel color.
	outColor = vec4(pixelColor, 1.0);

	
}