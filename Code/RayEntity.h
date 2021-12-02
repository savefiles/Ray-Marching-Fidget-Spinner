#pragma once

#include <glm/glm.hpp>
#include "Mesh.h"
#include "Entity.h"
#include <vector>
#include <memory>
#include "Camera.h"

class RayEntity
{
public:
	
	RayEntity();
	// Since our rect is only used for ray marching, need to pass camera pos.
	RayEntity(Camera* camera, float* totalRunTime);
	~RayEntity();

	void Draw(Demo* d);
	void Update();

private:
	Entity* entity;
};