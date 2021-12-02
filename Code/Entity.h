#pragma once
#include "Mesh.h"
#include "Texture.h"
#include "Demo.h"
#include "Camera.h"

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm/glm.hpp>
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>

// forward declaration
class Demo;
class Mesh;
class Texture;

class Entity
{
public:

	// This is only for dynamic font
	// It will be left empty for everything else
	char name[100];

	// for wireframe hitboxes
	glm::vec3 color;

	// Entity data
	// Only one of these will be used
	// per entity, not both at once
	BufferCPU* matrixBufferCPU;
	BufferCPU* spriteBufferCPU;
	BufferCPU* colorBufferCPU;
	BufferCPU* cameraBufferCPU;
	VkDescriptorSet descriptor_set;

	// mesh of the entity
	Mesh* mesh;
	
	// for color and normal
	Texture* texture[2];

	glm::vec3 pos;
	//glm::vec3 rot;
	glm::quat rotQuat;
	glm::vec3 scale;
	glm::mat4 modelMatrix;
	glm::mat4 parentModelMatrix = glm::mat4(1.0);
	float radius;

	// pass in camera position if you are using ray tracing.
	Camera* camera;
	float* totalRunTime;

	glm::mat4 GetModelMatrix();
	glm::vec3 GetWorldPosition();

	Entity(glm::vec3 position = glm::vec3(0));
	~Entity();

	void CreateDescriptorSetColor();
	void CreateDescriptorSetRay();
	void CreateDescriptorSetBasic();
	void CreateDescriptorSetBumpy();
	void CreateDescriptorSet2D();
	void CreateDescriptorSetWire();

	void Draw(VkCommandBuffer cmd, VkPipelineLayout pipeline_layout);
	void Update();

	bool SeperatingAxisTest(Entity* const other);
};

