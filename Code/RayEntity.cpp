#include "RayEntity.h"

RayEntity::RayEntity()
{
}

RayEntity::RayEntity(Camera* camera, float* totalRunTime)
{
	// create an entity for the rectangle
	// we don't need a mesh as the vertex shader handles the
	// vertices by itself.
	entity = new Entity();
	//entity->mesh = new Mesh(vertexList, 4);	// Entity takes care of destruction.
	entity->camera = camera;
	entity->totalRunTime = totalRunTime;
	entity->CreateDescriptorSetRay();
}

RayEntity::~RayEntity()
{
	delete entity;
}

void RayEntity::Draw(Demo* d)
{
	d->ApplyPipelineRay();
	d->DrawEntity(entity);
}

void RayEntity::Update()
{
	entity->Update();
}
